pub struct State {
    n: usize,
    x: Vec<Vec<u64>>,
    z: Vec<Vec<u64>>,
    r: Vec<i32>,
    // TODO maybe i64?
    over32: usize,
    pw: [u64; 32],
}

impl State {
    pub fn new(n: usize) -> Self {
        let len = 2 * n + 1;
        let over32 = (n >> 5) + 1;
        let x = vec![vec![0; over32]; len];
        let z = vec![vec![0; over32]; len];
        let r = vec![0; len];

        let mut pw = [1; 32];
        for i in 1..pw.len() {
            pw[i] = 2 * pw[i - 1];
        }

        Self {
            n,
            x,
            z,
            r,
            over32,
            pw,
        }
    }

    fn clifford(&mut self, i: usize, k: usize) -> i32 {
        let mut e = 0;

        for j in 0..self.over32 {
            for l in 0..self.pw.len() {
                let pw = self.pw[l];
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
            let pw = self.pw[j & 31];
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
            let pw = self.pw[j & 31];
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

    pub fn print_ket(&mut self) {
        let g = self.gaussian_elimination();
        println!("{g} nonzero basis states");

        self.print_basis_state();

        for t in 0..self.pw[g] - 1 {
            let t2 = t ^ (t + 1);
            for i in 0..g {
                if t2 & self.pw[i] > 0 {
                    self.rowmult(2 * self.n, self.n + i);
                }
            }
            self.print_basis_state();
        }

        println!();
    }

    fn print_basis_state(&self) {
        let mut e = self.r[2 * self.n];

        for j in 0..self.n {
            let j5 = j >> 5;
            let pw = self.pw[j & 31];

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
            let pw = self.pw[j & 31];

            if (self.x[2 * self.n][j5] & pw) > 0 {
                print!("1")
            } else {
                print!("0")
            }
        }

        print!(">");
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

pub fn add(left: usize, right: usize) -> usize {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let mut state = State::new(2);
        state.print_ket();
    }
}
