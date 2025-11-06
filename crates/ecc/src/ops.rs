pub const trait Pow<Rhs = Self> {
    type Output;
    #[must_use = "this returns the result of the operation, without modifying the original"]
    fn pow(self, exp: Rhs) -> Self::Output;
}
