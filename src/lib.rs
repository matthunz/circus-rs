pub mod gate;
pub use gate::{CNotGate, Gate, Gates, HadamardGate};

pub enum Instruction {
    Gate(Gates),
    Measure { target: usize },
}

const PW: [u64; 32] = {
    let mut pw = [1; 32];
    let mut i = 1;
    while i < 32 {
        pw[i] = 2 * pw[i - 1];
        i += 1;
    }
    pw
};

pub struct State {
    n: usize,
    x: Box<[Box<[u64]>]>,
    z: Box<[Box<[u64]>]>,
    r: Box<[i32]>,
    over32: usize,
}

impl State {
    /// Create a quantum state with `n` number of qubits.
    pub fn new(n: usize) -> Self {
        let len = 2 * n + 1;
        let over32 = (n >> 5) + 1;
        let mut x = (0..len)
            .map(|_| vec![0; over32].into_boxed_slice())
            .collect::<Vec<_>>()
            .into_boxed_slice();
        let mut z = (0..len)
            .map(|_| vec![0; over32].into_boxed_slice())
            .collect::<Vec<_>>()
            .into_boxed_slice();
        let r = vec![0; len].into_boxed_slice();

        for i in 0..2 * n + 1 {
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

    pub fn cnot(&mut self, target: usize, control: usize) {
        let gate = CNotGate { target, control };
        gate.apply(self);
    }

    pub fn hadamard(&mut self, target: usize) {
        let gate = HadamardGate { target };
        gate.apply(self);
    }

    pub fn measure(&mut self, b: usize) -> u8 {
        let mut ran = false;

        let b5 = b >> 5;
        let pw = PW[b & 31];
        let mut p = 0;
        for a in 0..self.n
        // loop over stabilizer generators
        {
            if (self.x[a + self.n][b5] & pw) > 0 {
                ran = true;
            } // if a Zbar does NOT commute with Z_b (the
            if ran {
                break;
            } // operator being measured), then outcome is random

            p += 1;
        }

        // If outcome is indeterminate
        if ran {
            self.rowcopy(p, p + self.n); // Set Xbar_p := Zbar_p
            self.rowset(p + self.n, b + self.n); // Set Zbar_p := Z_b
            self.r[p + self.n] = 2 * (rand::random::<i32>() % 2); // moment of quantum randomness
            for i in 0..2 * self.n {
                // Now update the Xbar's and Zbar's that don't commute with
                if (i != p) && (self.x[i][b5] & pw > 0) {
                    self.rowmult(i, p);
                } // Z_b
            }

            if self.r[p + self.n] > 0 {
                return 3;
            } else {
                return 2;
            }
        }

        // If outcome is determinate
        if !ran {
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
                return 1;
            } else {
                return 0;
            }
        }

        return 0;
    }

    pub fn print_ket(&mut self) {
        let g = self.gaussian_elimination();
        println!("2^{g} nonzero basis states");

        self.print_basis_state();

        for t in 0..PW[g] - 1 {
            let t2 = t ^ (t + 1);
            for i in 0..g {
                if t2 & PW[i] > 0 {
                    self.rowmult(2 * self.n, self.n + i);
                }
            }
            self.print_basis_state();
        }
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

    fn gaussian_elimination(&mut self) -> usize {
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
                        self.rowmult(k2, i); // Gaussian elimination step
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

    pub fn print(&self) {
        for i in 0..2 * self.n {
            if i == self.n {
                print!("\n");
                for _ in 0..self.n + 1 {
                    print!("-");
                }
            }
            if self.r[i] == 2 {
                print!("\n-");
            } else {
                print!("\n+");
            }
            for j in 0..self.n {
                let j5 = j >> 5;
                let pw = PW[j & 31];
                if (!(self.x[i][j5] & pw > 0)) && (!(self.z[i][j5] & pw > 0)) {
                    print!("I");
                }
                if (self.x[i][j5] & pw > 0) && (!(self.z[i][j5] & pw > 0)) {
                    print!("X");
                }
                if (self.x[i][j5] & pw > 0) && (self.z[i][j5] & pw > 0) {
                    print!("Y");
                }
                if (!(self.x[i][j5] & pw) > 0) && (self.z[i][j5] & pw > 0) {
                    print!("Z");
                }
            }
        }
        print!("\n");
    }

    fn print_basis_state(&self) {
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
            0 => print!(" +|"),
            1 => print!("+i|"),
            2 => print!(" -|"),
            3 => print!("-i|"),
            _ => {}
        }

        for j in 0..self.n {
            let j5 = j >> 5;
            let pw = PW[j & 31];

            if (self.x[2 * self.n][j5] & pw) > 0 {
                print!("1")
            } else {
                print!("0")
            }
        }

        println!(">");
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

pub struct Measurements<'s, I> {
    state: &'s mut State,
    iter: I,
}

impl<I> Iterator for Measurements<'_, I>
where
    I: Iterator<Item = Instruction>,
{
    type Item = u8;

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

pub fn add(left: usize, right: usize) -> usize {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let mut state = State::new(2);
        state.hadamard(0);
        state.cnot(0, 1);
        dbg!(state.measure(1));

        state.print();
        state.gaussian_elimination();
        state.print();
        state.print_ket();
    }
}
