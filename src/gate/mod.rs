mod cnot;
pub use cnot::CNotGate;

mod hadamard;
pub use hadamard::HadamardGate;

use crate::State;

pub trait Gate {
    fn apply(&self, state: &mut State);
}

pub enum Gates {
    Hadamard(HadamardGate),
}

impl Gate for Gates {
    fn apply(&self, state: &mut State) {
        match self {
            Self::Hadamard(h) => h.apply(state),
        }
    }
}
