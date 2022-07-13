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
        state.nonzero();
        state.print();
        state.print_ket();
    }
}
