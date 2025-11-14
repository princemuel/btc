use core::fmt;
use core::ops::{Add, Div, Mul, Neg, Sub};
use std::hash::Hash;

use crate::ops::Pow;

/// A field
///
/// <https://en.wikipedia.org/wiki/Field_(mathematics)>
pub trait Field:
    Neg<Output = Self>
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Eq
    + Copy
    + fmt::Debug
{
    const CHARACTERISTIC: u128;
    const ZERO: Self;
    const ONE: Self;

    /// Multiplicative inverse
    fn inverse(self) -> Self;

    /// Z-mod structure
    fn int_mul(self, a: i128) -> Self;
    fn from_int(a: i128) -> Self { Self::ONE.int_mul(a) }

    /// Iterate over all elements in this field
    ///
    /// The iterator finishes only for finite fields.
    type ElementsIter: Iterator<Item = Self>;
    fn elements() -> Self::ElementsIter;
}

/// Field element in the prime field GF(P) = ℤ/Pℤ, where P is a prime.
///
/// This represents an element in a finite field of prime order P. The field
/// operations (addition, multiplication, etc.) are performed modulo P.
///
/// # Type Parameters
/// * `P` - The prime modulus defining the field. Must be prime for correct
///   field semantics.
///
/// # Invariants
/// * The inner value is always in the range [0, P)
/// * P should be prime (not enforced at compile time)
///
/// # Representation
/// Uses `u128` internally, supporting primes up to 2^128 - 1, though practical
/// usage may be limited by performance considerations for very large primes.
///
/// # Examples
/// ```
/// // Field element in GF(17)
/// let elem = PrimeField::<17>::new(15);
/// ```
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PrimeField<const P: u128>(i128);
impl<const P: u128> PrimeField<P> {
    /// Creates a new field element from a signed integer.
    ///
    /// Automatically reduces the value into the canonical range [0, P) using
    /// modular arithmetic. Negative values wrap around correctly: -1 becomes
    /// P-1, -2 becomes P-2, etc.
    ///
    /// # Errors
    ///
    /// Returns `Err` if:
    /// * `P == 0` - Field modulus cannot be zero (division by zero)
    /// * `P == 1` - Field modulus must be at least 2 for a valid field
    /// * `P > i128::MAX` - Modulus too large to fit in i128
    ///
    /// # Examples
    ///
    /// ```
    /// let a = PrimeField::<17>::new(-3)?; // -3 ≡ 14 (mod 17)
    /// let b = PrimeField::<17>::new(20)?; // 20 ≡ 3 (mod 17)
    /// assert_eq!(a.value(), 14);
    /// assert_eq!(b.value(), 3);
    /// ```
    pub const fn new(value: i128) -> Result<Self, &'static str> {
        if P == 0 {
            return Err("Field modulus cannot be zero");
        }

        if P == 1 {
            return Err("Field modulus must be at least 2");
        }

        if P > i128::MAX as u128 {
            return Err("Modulus too large for i128");
        }

        let prime: i128 = P as i128;

        let reduced = value.rem_euclid(prime);
        debug_assert!(reduced >= 0 && reduced < prime);
        Ok(Self(reduced))
    }

    /// Returns the prime modulus P defining this field.
    pub const fn prime() -> u128 { P }

    /// Returns the canonical representative of this field element in [0, P).
    ///
    /// Note: The internal representation may be any i128 value before
    /// reduction, but this always returns a value in the valid range.
    pub const fn value(&self) -> u128 {
        debug_assert!(self.0 >= 0 && self.0 < P as i128);
        self.0 as u128
    }

    /// Reduces the representation into the range [0, P)
    fn reduce(self) -> Self {
        let p: i128 = P.try_into().expect("Modulus too large for i128");
        let reduced = self.0.rem_euclid(p);
        debug_assert!(reduced >= 0 && reduced < p);
        Self(reduced)
    }
}

impl<const P: u128> PrimeField<P> {
    /// Computes the multiplicative inverse using the Extended Euclidean
    /// Algorithm.
    ///
    /// Panics if the element is zero (which has no inverse).
    pub fn inverse(&self) -> Self {
        assert!(self.0 != 0, "Cannot invert zero element");

        let p = P as i128;
        let (gcd, x, _) = extended_gcd(self.0, p);

        assert_eq!(gcd, 1, "Element has no inverse (P might not be prime)");

        // x might be negative, reduce to [0, P)
        Self(x).reduce()
    }
}

impl<const P: u128> fmt::Debug for PrimeField<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { write!(f, "PrimeField_{}({})", P, self.0) }
}

impl<const P: u128> fmt::Display for PrimeField<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { write!(f, "{}", self.0) }
}

impl<const P: u128> From<PrimeField<P>> for u128 {
    fn from(value: PrimeField<P>) -> Self { value.value() }
}

impl<const P: u128> From<i128> for PrimeField<P> {
    fn from(value: i128) -> Self { Self(value.rem_euclid(P as i128)) }
}

impl<const P: u128> Add for PrimeField<P> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let p = P as i128;
        let sum = self.0 + rhs.0;
        // Since both operands are in [0, P), sum is in [0, 2P)
        // So we only need one conditional subtraction
        Self(if sum >= p { sum - p } else { sum })
    }
}

impl<const P: u128> Sub for PrimeField<P> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        // compute (a + P - b) % P to avoid negative numbers
        Self(self.0 - rhs.0).reduce()
    }
}

impl<const P: u128> Mul for PrimeField<P> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output { Self(self.0 * rhs.0).reduce() }
}

impl<const P: u128> Mul<i128> for PrimeField<P> {
    type Output = Self;

    fn mul(self, rhs: i128) -> Self::Output { Self(self.0 * rhs).reduce() }
}

impl<const P: u128> Mul<u128> for PrimeField<P> {
    type Output = Self;

    fn mul(self, rhs: u128) -> Self::Output { self * (rhs as i128) }
}

impl<const P: u128> Mul<PrimeField<P>> for u128 {
    type Output = PrimeField<P>;

    fn mul(self, rhs: PrimeField<P>) -> Self::Output { rhs * self }
}

impl<const P: u128> Div for PrimeField<P> {
    type Output = Self;

    // use Fermat's little theorem: a / b = a * b ^ (p - 2)
    // (n ** (p - 1)) % p == 1
    // this means: 1 / n == pow(n, p - 2, p)
    fn div(self, rhs: Self) -> Self::Output { self.mul(rhs.inverse()) }
}

impl<const P: u128> Neg for PrimeField<P> {
    type Output = Self;

    fn neg(self) -> Self {
        if self.0 == 0 {
            self // -0 = 0
        } else {
            Self((P as i128) - self.0)
        }
    }
}

/// Exponentiation with a unsigned exponent.
impl<const P: u128> Pow<u128> for PrimeField<P> {
    type Output = Self;

    fn pow(self, exp: u128) -> Self::Output {
        if exp == 0 {
            return Self::new(1).expect("P must be >= 2");
        }

        let mut base = self;
        let mut exponent = exp;
        let mut result = Self::new(1).expect("P must be >= 2");

        // Exponentiation by squaring: O(log exp)
        while exponent > 0 {
            if exponent & 1 == 1 {
                result = result * base;
            }
            base = base * base;
            exponent >>= 1;
        }

        result
    }
}

/// Exponentiation with signed exponent.
impl<const P: u128> Pow<i128> for PrimeField<P> {
    type Output = Self;

    fn pow(self, exp: i128) -> Self {
        if exp >= 0 {
            // Positive exponent: use normal exponentiation
            self.pow(exp as u128)
        } else {
            // Negative exponent: (a^-n) = (a^-1)^n
            self.inverse().pow((-exp) as u128)
        }
    }
}

/// Extended Euclidean Algorithm
/// Returns (gcd, x, y) where ax + by = gcd(a, b)
fn extended_gcd(a: i128, b: i128) -> (i128, i128, i128) {
    if a == 0 {
        return (b, 0, 1);
    }

    let (gcd, x1, y1) = extended_gcd(b % a, a);
    let x = y1 - (b / a) * x1;
    let y = x1;

    (gcd, x, y)
}

#[cfg(test)]
mod tests {
    use super::*;

    type F31 = PrimeField<31>;
    type F12 = PrimeField<12>;
    type F13 = PrimeField<13>;
    type F19 = PrimeField<19>;
    // type F7 = PrimeField<7>;

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
        assert_eq!((a * 4i128).value(), 12);
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
        assert_eq!(c.pow(-3i128).value(), c.pow(3i128).inverse().value());
    }
}
