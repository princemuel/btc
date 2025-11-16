use core::fmt;
use core::marker::PhantomData;
use core::ops::{Add, Div, Mul, Neg, Sub};
use std::hash::Hash;

use crate::helpers::{extended_euclidean_algorithm, mulmod, powmod};
use crate::macros::*;
use crate::ops::Pow;

/// A field supporting arithmetic operations and trait-based exponentiation.
///
/// <https://en.wikipedia.org/wiki/Field_(mathematics)>
pub trait Field:
    Neg<Output = Self>
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Pow<u128, Output = Self>
    + Pow<i128, Output = Self>
    + From<i128>
    + From<u128>
    + Eq
    + Copy
    + fmt::Debug
{
    const CHARACTERISTIC: u128;
    const ZERO: Self;
    const ONE: Self;

    /// Computes the multiplicative inverse using the Extended Euclidean
    /// Algorithm.
    ///
    /// Panics if the element is zero (which has no inverse).
    fn inverse(&self) -> Self;

    /// Scalar multiplication by an integer
    fn int_mul(&self, a: i128) -> Self;

    /// Converts an integer to a field element
    fn from_int(a: i128) -> Self { Self::from(a) }

    /// Iterate over all elements in this field
    ///
    /// The iterator finishes only for finite fields.
    type Iter: Iterator<Item = Self>;
    fn iter() -> Self::Iter;
}

/// Represents possible errors that can occur when creating a `PrimeField`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FieldError {
    ModulusIsLessThanTwo,
    ModulusIsTooLarge,
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
    /// let a = PrimeField::<17>::new(-3).unwrap(); // -3 ≡ 14 (mod 17)
    /// let b = PrimeField::<17>::new(20).unwrap(); // 20 ≡ 3 (mod 17)
    /// assert_eq!(a.value(), 14);
    /// assert_eq!(b.value(), 3);
    /// ```
    pub const fn new(value: i128) -> Result<Self, FieldError> {
        if P < 2 {
            return Err(FieldError::ModulusIsLessThanTwo);
        }

        if P > i128::MAX as u128 {
            return Err(FieldError::ModulusIsTooLarge);
        }

        let prime: i128 = P as i128;

        let reduced = value.rem_euclid(prime);
        debug_assert!(reduced >= 0 && reduced < prime, "wasn't true");
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
    const fn reduce(self) -> Self { Self(self.0.rem_euclid(P as i128)) }
}

impl<const P: u128> Field for PrimeField<P> {
    type Iter = PrimeFieldElements<P>;

    const CHARACTERISTIC: u128 = P;
    const ONE: Self = Self(1);
    const ZERO: Self = Self(0);

    fn inverse(&self) -> Self {
        assert!(self.0 != 0, "Cannot invert zero element");

        let p = P as i128;
        let (gcd, x, _) = extended_euclidean_algorithm(self.0, p);

        assert_eq!(gcd, 1, "Element has no inverse (P might not be prime)");

        Self(x).reduce()
    }

    fn int_mul(&self, a: i128) -> Self { *self * Self::from(a) }

    fn iter() -> Self::Iter {
        PrimeFieldElements {
            ptr:     0,
            _marker: PhantomData,
        }
    }
}

/// Iterator over all elements in a prime field GF(P)
pub struct PrimeFieldElements<const P: u128> {
    ptr:     u128,
    _marker: PhantomData<PrimeField<P>>,
}

impl<const P: u128> Iterator for PrimeFieldElements<P> {
    type Item = PrimeField<P>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.ptr < P {
            let value = PrimeField::from(self.ptr);
            self.ptr += 1;
            Some(value)
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = (P - self.ptr) as usize;
        (remaining, Some(remaining))
    }
}

impl<const P: u128> ExactSizeIterator for PrimeFieldElements<P> {
    fn len(&self) -> usize { (P - self.ptr) as usize }
}
impl<const P: u128> fmt::Debug for PrimeField<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { write!(f, "PrimeField_{P}({})", self.0) }
}

impl<const P: u128> fmt::Display for PrimeField<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { write!(f, "{}", self.0) }
}

impl<const P: u128> From<PrimeField<P>> for i128 {
    fn from(value: PrimeField<P>) -> Self { value.0 }
}

impl<const P: u128> From<PrimeField<P>> for u128 {
    fn from(value: PrimeField<P>) -> Self { value.0 as u128 }
}

impl_from_int_for_primefield_signed!(i8, i16, i32, i64);
impl_from_int_for_primefield_unsigned!(u8, u16, u32, u64);
impl_tryfrom_primefield_for_small_ints!(i8, i16, i32, i64, u8, u16, u32, u64);

impl<const P: u128> From<i128> for PrimeField<P> {
    fn from(value: i128) -> Self { Self(value.rem_euclid(P as i128)) }
}

impl<const P: u128> From<u128> for PrimeField<P> {
    fn from(value: u128) -> Self { Self((value % P) as i128) }
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
        let p = P as i128;
        let diff = self.0 - rhs.0;
        Self(if diff < 0 { diff + p } else { diff })
    }
}

impl<const P: u128> Mul for PrimeField<P> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        // let p = P as i128;
        // // Use checked multiplication to detect overflow
        // match self.0.checked_mul(rhs.0) {
        //     Some(product) => Self(product.rem_euclid(p)),
        //     None => {
        //         // Overflow occurred - use modular multiplication
        //         Self(mulmod(self.0, rhs.0, p))
        //     },
        // }

        let a = self.0 as u128;
        let b = rhs.0 as u128;

        Self(mulmod(a, b, P) as i128)
    }
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
impl<const P: u128> Pow for PrimeField<P> {
    type Output = Self;

    fn pow(self, exp: u128) -> Self::Output {
        if exp == 0 {
            return Self::ONE;
        }
        Self(powmod(self.0, exp, P as i128))
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

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use super::*;

    type F31 = PrimeField<31>;
    type F12 = PrimeField<12>;
    type F13 = PrimeField<13>;
    type F19 = PrimeField<19>;
    type F60 = PrimeField<60>;
    // type F7 = PrimeField<7>;

    #[test]
    fn test_point_negation() {
        let a = F31::new(2).unwrap();
        let b = F31::new(2).unwrap();
        let c = F31::new(15).unwrap();

        assert_eq!(a, b);
        assert_ne!(a, c);
        assert!(a == b);
    }

    #[test]
    fn test_point_addition() {
        let a = F31::new(2).unwrap();
        let b = F31::new(15).unwrap();
        assert_eq!(a + b, F31::new(17).unwrap());

        let a = F31::new(17).unwrap();
        let b = F31::new(21).unwrap();
        assert_eq!(a + b, F31::new(7).unwrap());

        let a = F12::new(3).unwrap();
        let b = F12::new(47).unwrap();
        assert_eq!(a + b, F12::new(2).unwrap());

        let a = F60::new(12).unwrap();
        let b = F60::new(843).unwrap();
        assert_eq!(a + b, F60::new(15).unwrap());

        let a = F60::new(23).unwrap();
        let b = F60::new(97).unwrap();
        assert_eq!(a + b, F60::new(0).unwrap());
    }

    #[test]
    fn test_point_subtraction() {
        let a = F31::new(29).unwrap();
        let b = F31::new(4).unwrap();
        assert_eq!(a - b, F31::new(25).unwrap());

        let a = F31::new(15).unwrap();
        let b = F31::new(30).unwrap();
        assert_eq!(a - b, F31::new(16).unwrap());

        let a = F12::new(3).unwrap();
        let b = F12::new(16).unwrap();
        assert_eq!(a - b, F12::new(11).unwrap());
    }

    #[test]
    fn test_point_mul() {
        let a = F31::new(24).unwrap();
        let b = F31::new(19).unwrap();
        assert_eq!(a * b, F31::new(22).unwrap());

        let a = F13::new(3).unwrap();
        assert_eq!((a * 4i128).value(), 12);
        assert_eq!((4 * a).value(), 12);
    }

    #[test]
    fn test_point_pow() {
        let a = F31::new(17).unwrap();
        assert_eq!(a.pow(3i128), F31::new(15).unwrap());

        let a = F31::new(5).unwrap();
        let b = F31::new(18).unwrap();
        assert_eq!(a.pow(5i128) * b, F31::new(16).unwrap());

        let a = F13::new(3).unwrap();
        assert_eq!(a.pow(3u128).value(), 1);
        assert_eq!(a.pow(0u128).value(), 1);

        let a = F19::new(7).unwrap();
        let b = a.pow(-1i128); // multiplicative inverse
        assert_eq!((a * b).value(), 1);

        let a = F13::new(3).unwrap();
        assert_eq!(a.pow(-3i128).value(), a.pow(3i128).inverse().value());
    }

    #[test]
    fn test_point_division() {
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
    fn test_additive_group_full_cycles() {
        use std::collections::HashSet;

        fn test<const P: u128>() {
            let expected = (0..P as i128).map(Into::into).collect();

            for item in 1..P {
                let item = PrimeField::<P>::from(item as i128);

                let mut x = PrimeField::ZERO;
                let mut actual = HashSet::with_capacity(P as usize);

                for _ in 0..P {
                    actual.insert(x);
                    x = x + item;
                }

                assert_eq!(actual, expected);
            }
        }

        test::<5>();
        test::<7>();
        test::<11>();
        test::<13>();
        test::<17>();
        test::<19>();
        test::<23>();
        test::<71>();
        test::<101>();
    }

    #[test]
    fn test_multiplicative_group_orders() {
        fn test<const P: u128>() {
            for a in PrimeField::<P>::iter() {
                if a == PrimeField::ZERO {
                    continue;
                }

                // compute order of a in multiplication
                let mut x = PrimeField::ONE;

                for k in 1..=P - 1 {
                    x = x * a;

                    if x == PrimeField::ONE {
                        assert!((P - 1).is_multiple_of(k), "multiplicative order must divide P-1");
                        break;
                    }

                    // must return to 1 within P-1 steps
                    assert!(k < P - 1, "invalid cycle detected");
                }
            }
        }

        test::<5>();
        test::<7>();
        test::<11>();
        test::<13>();
        test::<17>();
        test::<19>();
        test::<23>();
        test::<71>();
        test::<101>();
    }

    #[test]
    fn test_field_axioms() {
        fn test<const P: u128>() {
            let elems: HashSet<_> = PrimeField::<P>::iter().collect();

            for a in &elems {
                for b in &elems {
                    for c in &elems {
                        // additive associativity
                        assert_eq!((*a + *b) + *c, *a + (*b + *c));

                        // multiplicative associativity
                        assert_eq!((*a * *b) * *c, *a * (*b * *c));

                        // distributivity
                        assert_eq!(*a * (*b + *c), (*a * *b) + (*a * *c));
                    }
                }
            }

            for a in &elems {
                // additive identity
                assert_eq!(*a + PrimeField::ZERO, *a);

                // multiplicative identity
                assert_eq!(*a * PrimeField::ONE, *a);

                // additive inverse exists
                let neg = -*a;
                assert_eq!(*a + neg, PrimeField::ZERO);

                if *a != PrimeField::ZERO {
                    let inv = a.inverse();
                    assert_eq!(*a * inv, PrimeField::ONE);
                }
            }
        }

        test::<5>();
        test::<7>();
        test::<11>();
        test::<13>();
        test::<17>();
        test::<19>();
        test::<23>();
        test::<71>();
        test::<101>();
    }

    #[test]
    fn test_find_generator() {
        fn multiplicative_order<const P: u128>(a: PrimeField<P>) -> u128 {
            let mut x = PrimeField::ONE;
            for k in 1..=P - 1 {
                x = x * a;
                if x == PrimeField::ONE {
                    return k;
                }
            }
            unreachable!()
        }

        fn test<const P: u128>() {
            // search for the first generator
            let gnrt = PrimeField::iter()
            .skip(1) // skip 0
            .find(|&x| multiplicative_order::<P>(x) == P-1)
            .expect("GF(P)* must have at least one generator");

            // verify generator actually hits everything
            let mut set = HashSet::new();
            let mut x = PrimeField::ONE;

            for _ in 0..P - 1 {
                set.insert(x);
                x = x * gnrt;
            }

            assert_eq!(set.len() as u128, P - 1);
        }

        test::<5>();
        test::<7>();
        test::<11>();
        test::<13>();
        test::<17>();
        test::<19>();
        test::<23>();
        test::<71>();
        test::<101>();
    }

    #[test]
    fn test_additive_and_multiplicative_inverses() {
        fn test<const P: u128>() {
            for x in -7i128..=7 {
                let x = PrimeField::<P>::from(x);
                // Additive inverse
                assert_eq!(x + (-x), PrimeField::ZERO, "x + (-x) = 0 for x={:?}", x);
                assert_eq!((-x) + x, PrimeField::ZERO, "(-x) + x = 0 for x={:?}", x);
                assert_eq!(x - x, PrimeField::ZERO, "x - x = 0 for x={:?}", x);

                // Multiplicative inverse
                if x != PrimeField::ZERO {
                    let inv = x.inverse();

                    assert_eq!(x * inv, PrimeField::ONE, "x * x⁻¹ = 1 for x={:?}", x);
                    assert_eq!(inv * x, PrimeField::ONE, "x⁻¹ * x = 1 for x={:?}", x);
                    assert_eq!(
                        (x.0 * inv.0).rem_euclid(P as i128),
                        1,
                        "integer product matches 1 mod P for x={x:?}",
                    );
                    assert_eq!(x / x, PrimeField::ONE, "x / x = 1 for x={:?}", x);
                }
            }
        }

        test::<5>();
        test::<7>();
        test::<11>();
        test::<13>();
        test::<17>();
        test::<19>();
        test::<23>();
        test::<71>();
        test::<101>();
    }

    #[test]
    fn test_large_prime_field() {
        const P: u128 = 2_u128.pow(127) - 1; // Largest prime fitting into i128 (Mersenne prime M127)
        type F = PrimeField<P>;

        let x = F::from((P - 1) as i128);
        let y = x.inverse();
        assert_eq!(x * y, F::ONE);
    }

    #[test]
    fn integer_mul() {
        type F = PrimeField<23>;
        for x in 0..23u128 {
            let x = F::from(x);
            for n in -7..7 {
                assert_eq!(x.int_mul(n), F::from(n * x.0));
            }
        }
    }

    #[test]
    fn from_integer() {
        type F = PrimeField<23>;

        assert_eq!(F::from(0u128), F::ZERO);
        assert_eq!(F::from(1u128), F::ONE);
        assert_eq!(F::from(23u128), F::ZERO);
        assert_eq!(F::from(24u128), F::ONE);

        for x in -100..100 {
            assert_eq!(F::from_int(x), F::from(x));
        }
    }

    #[test]
    fn test_iterator_yields_correct_elements() {
        type F = PrimeField<13>;

        let elems: Vec<_> = F::iter().collect();
        assert_eq!(elems.len(), 13);

        for (i, &field) in elems.iter().enumerate().take(13) {
            assert_eq!(field, F::from(i as i128));
        }
    }

    #[test]
    fn test_iterator_is_finite_and_exhaustive() {
        type F = PrimeField<17>;
        let elems: Vec<_> = F::iter().collect();
        assert_eq!(elems.len(), 17);

        let mut set = HashSet::with_capacity(elems.len());
        for e in elems {
            assert!(set.insert(e));
        }
    }

    #[test]
    fn test_iterator_additive_identity_behavior() {
        type F = PrimeField<11>;

        for x in F::iter() {
            assert_eq!(x + F::ZERO, x);
            assert_eq!(F::ZERO + x, x);
        }
    }

    #[test]
    fn test_iterator_multiplicative_identity_behavior() {
        type F = PrimeField<11>;

        for x in F::iter() {
            assert_eq!(x * F::ONE, x);
            assert_eq!(F::ONE * x, x);
        }
    }

    #[test]
    fn test_iterator_finds_all_inverses() {
        type F = PrimeField<17>;

        for x in F::iter() {
            if x != F::ZERO {
                let inv = x.inverse();
                assert_eq!(x * inv, F::ONE);
            }
        }
    }

    #[test]
    fn test_iterator_fermats_little_theorem() {
        type F = PrimeField<19>;

        for a in F::iter() {
            if a != F::ZERO {
                assert_eq!(a.pow((19 - 1) as u128), F::ONE);
            }
        }
    }

    #[test]
    fn test_iterator_additive_closure() {
        type F = PrimeField<23>;
        let elems: HashSet<_> = F::iter().collect();

        for a in &elems {
            for b in &elems {
                let sum = *a + *b;
                assert!(sum.0 < 23, "sum must also be a field element");
                assert!(
                    elems.contains(&sum),
                    "Addition is not closed: {a:?} + {b:?} = {sum:?} not in field",
                );
            }
        }
    }

    #[test]
    fn test_iterator_multiplicative_closure() {
        type F = PrimeField<23>;
        let elems: HashSet<_> = F::iter().collect();

        for a in &elems {
            for b in &elems {
                let product = *a * *b;
                assert!(product.0 < 23, "product must also be a field element");
                assert!(
                    elems.contains(&product),
                    "Multiplication is not closed: {a:?} + {b:?} = {product:?} not in field",
                );
            }
        }
    }

    #[test]
    fn test_iterator_produces_correct_scalar_multiplication() {
        type F = PrimeField<7>;

        for x in F::iter() {
            for k in -10..10 {
                assert_eq!(x.int_mul(k), x * F::from(k));
            }
        }
    }

    #[test]
    fn test_iterator_zero_positioned_first() {
        type F = PrimeField<101>;
        let mut iter = F::iter();
        assert_eq!(iter.next(), Some(F::ZERO));
    }

    #[test]
    fn test_iterator_last_element() {
        type F = PrimeField<101>;
        let elems: Vec<_> = F::iter().collect();
        assert_eq!(elems.last(), Some(&F::from(100)));
    }

    #[test]
    fn test_iterator_size_hint() {
        type F = PrimeField<31>;
        let mut it = F::iter();

        let (low, high) = it.size_hint();
        assert_eq!(low, 31);
        assert_eq!(high, Some(31));

        it.next();
        let (low, high) = it.size_hint();
        assert_eq!(low, 30);
        assert_eq!(high, Some(30));
    }

    #[test]
    fn test_iterator_full_cycle() {
        type F = PrimeField<29>;

        let elems: Vec<_> = F::iter().collect();
        assert_eq!(elems.len(), 29);

        let mut seen = [false; 29];

        for e in elems {
            seen[e.0 as usize] = true;
        }

        assert!(seen.iter().all(|&x| x));
    }

    #[test]
    fn test_discrete_log_correctness_for_small_prime() {
        type F = PrimeField<11>;
        let g = F::from(2u128);

        for h in F::iter() {
            let mut x = F::ONE;

            // brute force discrete log
            let mut found = None;
            for k in 0..11 {
                if x == h {
                    found = Some(k);
                    break;
                }
                x = x * g;
            }

            if let Some(k) = found {
                assert_eq!(g.pow(k as i128), h);
            }
        }
    }

    #[test]
    fn test_iterator_random_sampling() {
        use rand::Rng;
        type F = PrimeField<97>;

        let elems: Vec<_> = F::iter().collect();

        let mut rng = rand::rng();
        for _ in 0..500 {
            let idx = rng.random_range(0..97) as u128;
            assert_eq!(elems[idx as usize], F::from(idx));
        }
    }
}
