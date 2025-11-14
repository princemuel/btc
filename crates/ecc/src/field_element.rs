use core::fmt;
use core::ops::{Add, Div, Mul, Neg, Sub};

use crate::ops::Pow;

/// Field element in prime field GF(P) with compile-time prime `P`.
///
/// Invariant: 0 <= num < P
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct FieldElement<const P: u128>(u128);

impl<const P: u128> FieldElement<P> {
    pub const fn new(value: u128) -> Option<Self> {
        if P == 0 {
            return None;
        }
        let value = ((value % P) + P) % P;
        Some(Self(value))
    }

    /// Construct by reducing modulo P (always succeeds).
    pub const fn new_unchecked(value: u128) -> Self { Self(value % P) }

    pub const fn prime() -> u128 { P }

    pub const fn value(&self) -> u128 { self.0 }
}

impl<const P: u128> FieldElement<P> {
    /// Fast modular exponentiation: base^exp mod modulus.
    /// base and modulus are `u128` here for safe intermediate arithmetic.
    const fn pow_mod(base: u128, mut exp: u128, modulus: Option<u128>) -> u128 {
        let mut result = 1u128;
        let mut base = base;
        let modulus = if let Some(modulus) = modulus { modulus } else { P };

        while exp > 0 {
            if exp & 1 == 1 {
                result = (result * base) % modulus;
            }

            base = (base * base) % modulus;

            exp >>= 1;
        }

        result
    }

    pub const fn recip(&self) -> Self { Self(Self::pow_mod(self.0, P - 2, None)) }
}

impl<const P: u128> fmt::Debug for FieldElement<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { write!(f, "FieldElement_{}({})", P, self.0) }
}

impl<const P: u128> fmt::Display for FieldElement<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { write!(f, "{}", self.0) }
}

impl<const P: u128> Add for FieldElement<P> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output { Self((self.value() + rhs.value()) % P) }
}

impl<const P: u128> Sub for FieldElement<P> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        // compute (a + P - b) % P to avoid negative numbers
        Self((self.value() + P - rhs.value()) % P)
    }
}

impl<const P: u128> Mul for FieldElement<P> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output { Self((self.value() * rhs.value()) % P) }
}

impl<const P: u128> Mul<u128> for FieldElement<P> {
    type Output = Self;

    fn mul(self, rhs: u128) -> Self::Output { Self((self.value() * rhs) % P) }
}

impl<const P: u128> Mul<FieldElement<P>> for u128 {
    type Output = FieldElement<P>;

    fn mul(self, rhs: FieldElement<P>) -> Self::Output { rhs * self }
}

impl<const P: u128> Div for FieldElement<P> {
    type Output = Self;

    // use Fermat's little theorem: a / b = a * b ^ (p - 2)
    // (n ** (p - 1)) % p == 1
    // this means: 1 / n == pow(n, p - 2, p)
    fn div(self, rhs: Self) -> Self::Output { self.mul(rhs.recip()) }
}

impl<const P: u128> Neg for FieldElement<P> {
    type Output = Self;

    fn neg(self) -> Self {
        if self.value() == 0 {
            self
        } else {
            Self(P - self.value())
        }
    }
}

/// Exponentiation with a unsigned exponent.
impl<const P: u128> Pow<u128> for FieldElement<P> {
    type Output = Self;

    fn pow(self, exp: u128) -> Self::Output {
        let exp = exp % (P - 1);
        Self(Self::pow_mod(self.value(), exp, Some(P)))
    }
}

/// Exponentiation with signed exponent.
impl<const P: u128> Pow<i128> for FieldElement<P> {
    type Output = Self;

    fn pow(self, exp: i128) -> Self {
        let prime = P - 1;
        let exp = if exp < 0 {
            prime - ((-exp) as u128) % prime
        } else {
            (exp as u128) % prime
        };
        Self(Self::pow_mod(self.value(), exp, Some(P)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    type F31 = FieldElement<31>;
    type F12 = FieldElement<12>;
    type F13 = FieldElement<13>;
    type F19 = FieldElement<19>;
    // type F7 = FieldElement<7>;

    #[test]
    fn test_ne() {
        let a = F31::new(2).unwrap();
        let b = F31::new(2).unwrap();
        let c = F31::new(15).unwrap();

        assert_eq!(a, b);
        assert_ne!(a, c);
        assert!(a == b);
    }

    #[test]
    fn test_add() {
        let a = F31::new(2).unwrap();
        let b = F31::new(15).unwrap();
        assert_eq!(a + b, F31::new(17).unwrap());

        let a = F31::new(17).unwrap();
        let b = F31::new(21).unwrap();
        assert_eq!(a + b, F31::new(7).unwrap());

        let a = F12::new(3).unwrap();
        let b = F12::new(47).unwrap();
        assert_eq!(a + b, F12::new(2,).unwrap());
    }

    #[test]
    fn test_sub() {
        let a = F31::new(29).unwrap();
        let b = F31::new(4).unwrap();
        assert_eq!(a - b, F31::new(25).unwrap());

        let a = F31::new(15).unwrap();
        let b = F31::new(30).unwrap();
        assert_eq!(a - b, F31::new(16).unwrap());

        let a = F12::new(3).unwrap();
        let b = F12::new(16).unwrap();
        assert_eq!(a - b, F12::new(11,).unwrap());
    }

    #[test]
    fn test_mul() {
        let a = F31::new(24).unwrap();
        let b = F31::new(19).unwrap();
        assert_eq!(a * b, F31::new(22).unwrap());

        let a = F13::new(3).unwrap();
        assert_eq!((a * 4).value(), 12);
        assert_eq!((4 * a).value(), 12);
    }

    #[test]
    fn test_pow() {
        let a = F31::new(17).unwrap();
        assert_eq!(a.pow(3i128), F31::new(15).unwrap());

        let a = F31::new(5).unwrap();
        let b = F31::new(18).unwrap();
        assert_eq!(a.pow(5i128) * b, F31::new(16).unwrap());
    }

    #[test]
    fn test_div() {
        let a = F31::new(3).unwrap();
        let b = F31::new(24).unwrap();
        assert_eq!(a / b, F31::new(4).unwrap());

        let a = F31::new(17).unwrap();
        assert_eq!(a.pow(-3i128), F31::new(29).unwrap());

        let a = F31::new(4).unwrap();
        let b = F31::new(11).unwrap();
        assert_eq!(a.pow(-4i128) * b, F31::new(13).unwrap());
    }

    #[test]
    fn test_pow_u128() {
        let a = F13::new(3).unwrap();
        assert_eq!(a.pow(3u128).value(), 1);
        assert_eq!(a.pow(0u128).value(), 1);
    }

    #[test]
    fn test_pow_i128() {
        let a = F19::new(7).unwrap();
        let b = a.pow(-1i128); // multiplicative inverse
        assert_eq!((a * b).value(), 1);

        let c = F13::new(3).unwrap();
        assert_eq!(c.pow(-3i128).value(), c.pow(3i128).recip().value());
    }
}
