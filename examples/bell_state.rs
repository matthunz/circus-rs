use circus::State;

fn main() {
    // Create a bell state, or EPR pair, a superposition of qubits 0 and 1
    let mut state = State::new(2);
    state.h(0);
    state.cx(0, 1);
    state.measure(1);

    println!("{state}");
    println!("{}", state.ket());
}
