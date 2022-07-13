use core::fmt;

use crate::{
    gate::{CNotGate, Gate, HadamardGate},
    Instruction, Measurement, PW,
};

pub type BinaryMatrix = Box<[Box<[u64]>]>;

/// Create a binary matrix for an `n` number of qubits.
pub fn binary_matrix(n: usize) -> BinaryMatrix {
    let len = 2 * n + 1;
    let over32 = (n >> 5) + 1;

    (0..len)
        .map(|_| vec![0; over32].into_boxed_slice())
        .collect::<Vec<_>>()
        .into_boxed_slice()
}

/// Quantum state (from [Improved Simulation of Stabilizer Circuits](https://arxiv.org/abs/quant-ph/0406196)
/// by Scott Aaronson and Daniel Gottesman)
pub struct State {
    /// Number of qubits.
    pub n: usize,

    /// floor(n/8)+1
    pub over32: usize,

    /// (2n+1)*n matrix for stabilizer/destabilizer x bits.
    pub x: BinaryMatrix,

    /// (2n+1)*n matrix for stabilizer/destabilizer z bits.
    pub z: BinaryMatrix,

    /// Phase bits (0 for +1, 1 for i, 2 for -1, 3 for -i). Normally either 0 or 2.
    pub r: Box<[i32]>,
}

impl State {
    /// Create a quantum state with `n` number of qubits.
    pub fn new(n: usize) -> Self {
        let len = 2 * n + 1;
        let over32 = (n >> 5) + 1;
        let mut x = binary_matrix(n);
        let mut z = binary_matrix(n);
        let r = vec![0; len].into_boxed_slice();

        for i in 0..len {
            if i < n {
                x[i][i >> 5] = PW[i & 31];
            } else if i < 2 * n {
                let j = i - n;
                z[i][j >> 5] = PW[j & 31];
            }
        }

        Self { n, x, z, r, over32 }
    }

    pub fn run<I>(&mut self, iter: I) -> Measurements<'_, I::IntoIter>
    where
        I: IntoIterator<Item = Instruction>,
    {
        Measurements {
            state: self,
            iter: iter.into_iter(),
        }
    }

    /// Apply the controlled-NOT gate, also known as the controlled-x (CX) gate.
    /// It performs a NOT on the `target` whenever the `control` is in state `|1⟩`.
    pub fn cx(&mut self, target: usize, control: usize) {
        let gate = CNotGate { target, control };
        gate.apply(self);
    }

    /// Apply the Hadamard gate.
    /// Rotates the states `|0⟩` and `|1⟩` to `|+⟩` and `|-⟩`, respectively.
    pub fn h(&mut self, target: usize) {
        let gate = HadamardGate { target };
        gate.apply(self);
    }

    /// Measure the `target` qubit.
    pub fn measure(&mut self, target: usize) -> Measurement {
        let mut is_indeterminate = false;

        let b5 = target >> 5;
        let pw = PW[target & 31];
        let mut p = 0;

        // loop over stabilizer generators
        for a in 0..self.n {
            // if a Zbar does NOT commute with Z_b (the operator being measured), then outcome is random
            if (self.x[a + self.n][b5] & pw) > 0 {
                is_indeterminate = true;
                break;
            }
            p += 1;
        }

        if is_indeterminate {
            // Outcome is indeterminate
            self.rowcopy(p, p + self.n); // Set Xbar_p := Zbar_p
            self.rowset(p + self.n, target + self.n); // Set Zbar_p := Z_b
            self.r[p + self.n] = 2 * (rand::random::<i32>() % 2); // moment of quantum randomness
            for i in 0..2 * self.n {
                // Now update the Xbar's and Zbar's that don't commute with
                if (i != p) && (self.x[i][b5] & pw > 0) {
                    self.rowmult(i, p);
                } // Z_b
            }

            if self.r[p + self.n] > 0 {
                Measurement::random(true)
            } else {
                Measurement::random(false)
            }
        } else {
            let mut m = 0;
            for a in 0..self.n {
                // Before we were checking if stabilizer generators commute
                if self.x[a][b5] & pw > 0 {
                    // with Z_b; now we're checking destabilizer generators
                    break;
                }
                m += 1;
            }
            self.rowcopy(2 * self.n, m + self.n);
            for i in (m + 1)..self.n {
                if self.x[i][b5] & pw > 0 {
                    self.rowmult(2 * self.n, i + self.n);
                }
            }

            if self.r[2 * self.n] > 0 {
                Measurement::fixed(true)
            } else {
                Measurement::fixed(false)
            }
        }
    }

    /// Perform Gaussian elimination and calculate the number of nonzero basis states (in 2^n).
    pub fn nonzero(&mut self) -> usize {
        let mut i = self.n;
        let mut k = i;
        for j in 0..self.n {
            let j5 = j >> 5;
            let pw = PW[j & 31];
            for a in i..2 * self.n {
                // Find a generator containing X in jth column
                if (self.x[a][j5] & pw) > 0 {
                    break;
                }
                k += 1;
            }

            if k < 2 * self.n {
                self.rowswap(i, k);
                self.rowswap(i - self.n, k - self.n);
                for k2 in (i + 1)..2 * self.n {
                    if (self.x[k2][j5] & pw) > 0 {
                        // Gaussian elimination step
                        self.rowmult(k2, i);
                        self.rowmult(i - self.n, k2 - self.n);
                    }
                }

                i += 1;
            }
        }

        let g = i - self.n;

        for j in 0..self.n {
            let j5 = j >> 5;
            let pw = PW[j & 31];
            let mut k = i;
            for a in i..2 * self.n {
                // Find a generator containing Z in jth column
                if (self.z[a][j5] & pw) > 0 {
                    break;
                }
                k += 1;
            }

            if k < 2 * self.n {
                self.rowswap(i, k);
                self.rowswap(i - self.n, k - self.n);
                for k2 in (i + 1)..2 * self.n {
                    if (self.z[k2][j5] & pw) > 0 {
                        self.rowmult(k2, i);
                        self.rowmult(i - self.n, k2 - self.n);
                    }
                }

                i += 1;
            }
        }

        g
    }

    /// Format the current state as a string in bra-ket notation.
    pub fn ket(&mut self) -> String {
        let g = self.nonzero();

        let mut s = String::new();
        self.ket_basis_state(&mut s);

        for t in 0..PW[g] - 1 {
            let t2 = t ^ (t + 1);
            for i in 0..g {
                if t2 & PW[i] > 0 {
                    self.rowmult(2 * self.n, self.n + i);
                }
            }
            self.ket_basis_state(&mut s);
        }

        s
    }

    fn clifford(&mut self, i: usize, k: usize) -> i32 {
        let mut e = 0;

        for j in 0..self.over32 {
            for l in 0..PW.len() {
                let pw = PW[l];
                // X
                if (self.x[k][j] & pw) > 0 && (!(self.z[k][j] & pw)) > 0 {
                    if (self.x[i][j] & pw) > 0 && (self.z[i][j] & pw) > 0 {
                        e += 1; // XY=iZ
                    }

                    if (!(self.x[i][j] & pw) > 0) && (self.z[i][j] & pw) > 0 {
                        e -= 1; // XZ=-iY
                    }
                    if (self.x[k][j] & pw) > 0 && (self.z[k][j] & pw) > 0
                    // Y
                    {
                        if (!(self.x[i][j] & pw) > 0) && (self.z[i][j] & pw) > 0 {
                            e += 1; // YZ=iX
                        }

                        if (self.x[i][j] & pw) > 0 && (!(self.z[i][j] & pw) > 0) {
                            e -= 1; // YX=-iZ
                        }
                    }
                    if (!(self.x[k][j] & pw) > 0) && (self.z[k][j] & pw) > 0
                    // Z
                    {
                        if (self.x[i][j] & pw) > 0 && (!(self.z[i][j] & pw) > 0) {
                            e += 1; // ZX=iY
                        }

                        if (self.x[i][j] & pw) > 0 && (self.z[i][j] & pw) > 0 {
                            e -= 1; // ZY=-iX
                        }
                    }
                }
            }
        }

        e = (e + self.r[i] + self.r[k]) % 4;
        if e >= 0 {
            e
        } else {
            e + 4
        }
    }

    fn ket_basis_state(&self, s: &mut String) {
        let mut e = self.r[2 * self.n];

        for j in 0..self.n {
            let j5 = j >> 5;
            let pw = PW[j & 31];

            // Pauli operator is "Y"
            if (self.x[2 * self.n][j5] & pw) > 0 && (self.z[2 * self.n][j5] & pw) > 0 {
                e = (e + 1) % 4;
            }
        }

        match e {
            0 => s.push_str(" +|"),
            1 => s.push_str("+i|"),
            2 => s.push_str(" -|"),
            3 => s.push_str("-i|"),
            _ => {}
        }

        for j in 0..self.n {
            let j5 = j >> 5;
            let pw = PW[j & 31];

            if (self.x[2 * self.n][j5] & pw) > 0 {
                s.push_str("1")
            } else {
                s.push_str("0")
            }
        }

        s.push_str(">\n");
    }

    fn rowset(&mut self, i: usize, b: usize) {
        for j in 0..self.over32 {
            self.x[i][j] = 0;
            self.z[i][j] = 0;
        }
        self.r[i] = 0;
        if b < self.n {
            let b5 = b >> 5;
            let b31 = b & 31;
            self.x[i][b5] = PW[b31];
        } else {
            let b5 = (b - self.n) >> 5;
            let b31 = (b - self.n) & 31;
            self.z[i][b5] = PW[b31];
        }
    }

    fn rowcopy(&mut self, i: usize, k: usize) {
        for j in 0..self.over32 {
            self.x[i][j] = self.x[k][j];
            self.z[i][j] = self.z[k][j];
        }
        self.r[i] = self.r[k];
    }

    fn rowswap(&mut self, i: usize, k: usize) {
        self.rowcopy(2 * self.n, k);
        self.rowcopy(k, i);
        self.rowcopy(i, 2 * self.n);
    }

    fn rowmult(&mut self, i: usize, k: usize) {
        self.r[i] = self.clifford(i, k);
        for j in 0..self.over32 {
            self.x[i][j] ^= self.x[k][j];
            self.z[i][j] ^= self.z[k][j];
        }
    }
}

impl fmt::Display for State {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..2 * self.n {
            if i == self.n {
                f.write_str("\n")?;
                for _ in 0..self.n + 1 {
                    f.write_str("-")?;
                }
            }
            if self.r[i] == 2 {
                f.write_str("\n-")?;
            } else {
                f.write_str("\n+")?;
            }
            for j in 0..self.n {
                let j5 = j >> 5;
                let pw = PW[j & 31];
                if (!(self.x[i][j5] & pw > 0)) && (!(self.z[i][j5] & pw > 0)) {
                    f.write_str("I")?;
                }
                if (self.x[i][j5] & pw > 0) && (!(self.z[i][j5] & pw > 0)) {
                    f.write_str("X")?;
                }
                if (self.x[i][j5] & pw > 0) && (self.z[i][j5] & pw > 0) {
                    f.write_str("Y")?;
                }
                if (!(self.x[i][j5] & pw) > 0) && (self.z[i][j5] & pw > 0) {
                    f.write_str("Z")?;
                }
            }
        }
        f.write_str("\n")
    }
}

pub struct Measurements<'s, I> {
    state: &'s mut State,
    iter: I,
}

impl<I> Iterator for Measurements<'_, I>
where
    I: Iterator<Item = Instruction>,
{
    type Item = Measurement;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(instruction) = self.iter.next() {
                match instruction {
                    Instruction::Gate(gate) => {
                        gate.apply(self.state);
                    }
                    Instruction::Measure { target } => break Some(self.state.measure(target)),
                }
            } else {
                break None;
            }
        }
    }
}
