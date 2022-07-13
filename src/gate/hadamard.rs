use super::Gate;
use crate::State;

pub struct HadamardGate {
    pub target: usize,
}

impl Gate for HadamardGate {
    fn apply(&self, state: &mut State) {
        let b5 = self.target >> 5;
        let pw = state.pw[self.target & 31];
        for i in 0..2 * state.n {
            let tmp = state.x[i][b5];
            state.x[i][b5] ^= (state.x[i][b5] ^ state.z[i][b5]) & pw;
            state.z[i][b5] ^= (state.z[i][b5] ^ tmp) & pw;
            if (state.x[i][b5] & pw) > 0 && (state.z[i][b5] & pw) > 0 {
                state.r[i] = (state.r[i] + 2) % 4;
            }
        }
    }
}
