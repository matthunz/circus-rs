#[derive(Clone, Copy, Debug)]
pub struct Measurement {
    byte: u8,
}

impl Measurement {
    pub const fn new(byte: u8) -> Self {
        Self { byte }
    }

    pub const fn fixed(bit: bool) -> Self {
        Self::new(bit as u8)
    }

    pub const fn random(bit: bool) -> Self {
        Self::new(bit as u8 + 2)
    }

    pub const fn is_zero(self) -> bool {
        self.byte == 0 || self.byte == 2
    }

    pub const fn is_one(self) -> bool {
        self.byte == 1 || self.byte == 3
    }

    pub const fn is_random(self) -> bool {
        self.byte >= 2
    }
}
