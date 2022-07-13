use super::Gate;
use crate::{State, PW};

pub struct PhaseGate {
    pub target: usize,
}

impl Gate for PhaseGate {
    fn apply(&self, state: &mut State) {
        let b5 = self.target >> 5;
        let pw = PW[self.target & 31];

        for i in 0..2 * state.n {
            if state.x[i][b5] & pw > 0 && state.z[i][b5] & pw > 0 {
                state.r[i] = (state.r[i] + 2) % 4;
            }
            state.z[i][b5] ^= state.x[i][b5] & pw;
        }
    }
}
