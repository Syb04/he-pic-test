pub struct Rng {
    state: u64,
}

impl Rng {
    pub fn new(seed: u64) -> Self {
        let state = if seed == 0 {
            0x9e37_79b9_7f4a_7c15
        } else {
            seed
        };
        Self { state }
    }

    pub fn uniform(&mut self) -> f64 {
        let x = self.next_u64() >> 11;
        (x as f64) * (1.0 / ((1u64 << 53) as f64))
    }

    pub fn normal(&mut self) -> f64 {
        let u1 = (1.0 - self.uniform()).max(1.0e-300);
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }

    pub fn sample_count(&mut self, expected: f64) -> usize {
        if expected <= 0.0 || !expected.is_finite() {
            return 0;
        }
        let whole = expected.floor();
        let frac = expected - whole;
        let mut count = whole as usize;
        if self.uniform() < frac {
            count += 1;
        }
        count
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.state = x;
        x.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_is_in_range() {
        let mut rng = Rng::new(42);
        for _ in 0..1000 {
            let x = rng.uniform();
            assert!((0.0..1.0).contains(&x));
        }
    }
}
