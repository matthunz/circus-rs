//! ```
//! // Create a bell state, or EPR pair, a superposition of qubits 0 and 1
//! let mut state = State::new(2);
//! state.h(0);
//! state.cx(0, 1);
//! state.measure(1);
//!
//! println!("{state}");
//! println!("{}", state.ket());
//! ```

pub mod gate;
use gate::Gates;

mod measurement;
pub use measurement::Measurement;

pub mod state;
pub use state::State;

pub enum Instruction {
    Gate(Gates),
    Measure { target: usize },
}

// Powers of 2 (PW[i] = 2^i)
const PW: [u64; 32] = {
    let mut pw = [1; 32];
    let mut i = 1;
    while i < 32 {
        pw[i] = 2 * pw[i - 1];
        i += 1;
    }
    pw
};
