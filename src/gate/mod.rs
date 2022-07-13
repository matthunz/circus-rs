mod cnot;
pub use cnot::CNotGate;

mod hadamard;
pub use hadamard::HadamardGate;

mod phase;
pub use phase::PhaseGate;

use crate::State;

pub trait Gate {
    fn apply(&self, state: &mut State);
}

pub enum Gates {
    CNot(CNotGate),
    Hadamard(HadamardGate),
    Phase(PhaseGate),
}

impl Gate for Gates {
    fn apply(&self, state: &mut State) {
        match self {
            Self::CNot(cx) => cx.apply(state),
            Self::Hadamard(h) => h.apply(state),
            Self::Phase(p) => p.apply(state),
        }
    }
}
