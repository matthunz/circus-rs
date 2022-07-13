use circus::State;

fn main() {
    // Create a bell state, a superposition of qubits 0 and 1
    let mut state = State::new(2);
    state.h(0);
    state.cx(1, 0);
    dbg!(state.measure(1));

   println!("{state}");
   println!("{}", state.ket());
}