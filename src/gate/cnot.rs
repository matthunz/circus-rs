use super::Gate;
use crate::{State, PW};

pub struct CNotGate {
    pub target: usize,
    pub control: usize,
}

impl Gate for CNotGate {
    fn apply(&self, state: &mut State) {
        let b5 = self.target >> 5;
        let c5 = self.control >> 5;
        let pwb = PW[self.target & 31];
        let pwc = PW[self.control & 31];
        for i in 0..2 * state.n {
            if state.x[i][b5] & pwb > 0 {
                state.x[i][c5] ^= pwc;
            }
            if state.z[i][c5] & pwc > 0 {
                state.z[i][b5] ^= pwb;
            }
            if (state.x[i][b5] & pwb > 0)
                && (state.z[i][c5] & pwc > 0)
                && (state.x[i][c5] & pwc > 0)
                && (state.z[i][b5] & pwb > 0)
            {
                state.r[i] = (state.r[i] + 2) % 4;
            }
            if (state.x[i][b5] & pwb > 0)
                && (state.z[i][c5] & pwc > 0)
                && !(state.x[i][c5] & pwc > 0)
                && !(state.z[i][b5] & pwb > 0)
            {
                state.r[i] = (state.r[i] + 2) % 4;
            }
        }
    }
}
