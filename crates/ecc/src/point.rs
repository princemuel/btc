use core::fmt;
use core::ops::Add;

use crate::field::Field;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Point<F, const A: i128, const B: i128>
where
    F: Field,
{
    Point(F, F),
    Infinity,
}

impl<F, const A: i128, const B: i128> Point<F, A, B>
where
    F: Field,
{
    pub fn new(x: F, y: F) -> Result<Self, &'static str> {
        // make sure that the elliptic curve equation is satisfied
        // y**2 == x**3 + a*x + b
        let y_sqr = y.pow(2u128);
        let x_cbd = x.pow(3u128);
        let ax = F::from(A) * x;
        let b = F::from(B);

        if y_sqr != (x_cbd + ax + b) {
            return Err("Invalid parameters to Point::new(). Point is not on the curve");
        }

        Ok(Self::Point(x, y))
    }

    pub const fn x(&self) -> Option<&F> {
        match self {
            Self::Point(x, _) => Some(x),
            Self::Infinity => None,
        }
    }

    pub const fn y(&self) -> Option<&F> {
        match self {
            Self::Point(_, y) => Some(y),
            Self::Infinity => None,
        }
    }

    pub const fn a(&self) -> i128 { A }

    pub const fn b(&self) -> i128 { B }

    pub const fn infinity() -> Self { Self::Infinity }

    pub const fn is_infinity(&self) -> bool { matches!(self, Self::Infinity) }
}

impl<F, const A: i128, const B: i128> fmt::Display for Point<F, A, B>
where
    F: Field + fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Infinity => write!(f, "Point(infinity)"),
            Self::Point(x, y) => {
                let pre = f.precision().unwrap_or(2);
                write!(f, "Point({x:.pre$},{y:.pre$})_{A:.pre$}_{B:.pre$}")
            },
        }
    }
}

impl<F, const A: i128, const B: i128> Add for Point<F, A, B>
where
    F: Field,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (Self::Infinity, Self::Infinity) => Self::infinity(),
            // self is the point at infinity, return rhs
            (Self::Infinity, Self::Point(..)) => rhs,
            // rhs is the point at infinity, return self
            (Self::Point(..), Self::Infinity) => self,

            // self.x ≠ rhs.x; (x, y) == (x1, y1) + (x2, y2)
            // s = (y2 - y1) / (x2 - x1)
            // x3 = s ^ 2 - x1 - x2
            // y3 = s * (x1 - x3) - y1
            (Self::Point(x1, y1), Self::Point(x2, y2)) if x1 != x2 => {
                let s = (y2 - y1) / (x2 - x1);
                let x3 = s.pow(2u128) - x1 - x2;
                let y3 = s * (x1 - x3) - y1;
                Self::Point(x3, y3)
            },

            // self.x == rhs.x, self.y != rhs.y; the vertical line results in a point at infinity
            (Self::Point(x1, y1), Self::Point(x2, y2)) if x1 == x2 && y1 != y2 => Self::infinity(),

            // if we are tangent to the vertical line,
            (Self::Point(x1, y1), Self::Point(..)) if y1 == F::ZERO => Self::infinity(),

            // self == rhs; (x, y) = (x1, y1) + (x1, y1)
            // s = (3 * x1 ^ 2 + a) / (2 * y1)
            // x3 = s ^ 2 - 2 * x1
            // y3 = s * (x1 - x3) - y1
            (Self::Point(x1, y1), Self::Point(x2, y2)) if x1 == x2 && y1 == y2 => {
                let s = (F::from(3u128) * x1.pow(2u128) + F::from(A)) / (F::from(2u128) * y1);
                let x3 = s.pow(2u128) - F::from(2u128) * x1;
                let y3 = s * (x1 - x3) - y1;
                Self::Point(x3, y3)
            },
            _ => unreachable!(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::PrimeField;
    use crate::ops::Pow as _;

    // Helper type alias for testing with regular integers (treating them as a "field")
    // For these tests, we'll use a large prime field to simulate integer arithmetic
    type TestField = PrimeField<223>; // Large enough to avoid most wraparound issues

    #[test]
    fn test_create_valid_point() {
        type P = Point<TestField, 5, 7>;
        let x = TestField::new(-1).unwrap();
        let y = TestField::new(-1).unwrap();
        assert!(P::new(x, y).is_ok());
    }

    #[test]
    fn test_create_invalid_point() {
        type P = Point<TestField, 5, 7>;
        let x = TestField::new(-1).unwrap();
        let y = TestField::new(-2).unwrap();
        assert!(P::new(x, y).is_err());
    }

    #[test]
    fn test_points_with_same_cords_are_equal() {
        type P = Point<TestField, 5, 7>;
        let a = P::new(TestField::new(-1).unwrap(), TestField::new(-1).unwrap()).unwrap();
        let b = P::new(TestField::new(-1).unwrap(), TestField::new(-1).unwrap()).unwrap();

        assert_eq!(a, b);
    }

    #[test]
    fn test_points_with_different_cords_are_different() {
        type P = Point<TestField, 5, 7>;
        let a = P::new(TestField::new(-1).unwrap(), TestField::new(-1).unwrap()).unwrap();
        let b = P::new(TestField::new(18).unwrap(), TestField::new(77).unwrap()).unwrap();

        assert_ne!(a, b);
    }

    #[test]
    fn test_sum_infinity_and_another_point_returns_the_point() {
        type P = Point<TestField, 5, 7>;
        let inf = P::infinity();
        let point = P::new(TestField::new(-1).unwrap(), TestField::new(-1).unwrap()).unwrap();

        assert_eq!(point + inf, point);
        assert_eq!(inf + point, point);
    }

    #[test]
    fn test_two_different_points_from_the_same_curve_can_be_added() {
        type P = Point<TestField, 5, 7>;
        let a = P::new(TestField::new(2).unwrap(), TestField::new(5).unwrap()).unwrap();
        let b = P::new(TestField::new(-1).unwrap(), TestField::new(-1).unwrap()).unwrap();

        let expected = P::new(TestField::new(3).unwrap(), TestField::new(-7).unwrap()).unwrap();
        assert_eq!(a + b, expected);
        assert_eq!(b + a, expected);
    }

    #[test]
    fn test_adding_opposite_points_results_infinity() {
        type P = Point<TestField, 5, 7>;
        let a = P::new(TestField::new(-1).unwrap(), TestField::new(-1).unwrap()).unwrap();
        let b = P::new(TestField::new(-1).unwrap(), TestField::new(1).unwrap()).unwrap();

        let expected = P::infinity();
        assert_eq!(a + b, expected);
    }

    #[test]
    fn test_finite_field_points() {
        // Test with actual finite field arithmetic
        type F223 = PrimeField<223>;
        type P = Point<F223, 0, 7>;

        // Valid points on y² = x³ + 7 over F_223
        let valid_points = [(192, 105), (17, 56), (1, 193)];
        for (x, y) in valid_points {
            let point = P::new(F223::new(x).unwrap(), F223::new(y).unwrap());
            assert!(point.is_ok(), "Point ({}, {}) should be on the curve", x, y);
        }

        // Invalid points
        let invalid_points = [(200, 119), (42, 99)];
        for (x, y) in invalid_points {
            let point = P::new(F223::new(x).unwrap(), F223::new(y).unwrap());
            assert!(point.is_err(), "Point ({}, {}) should not be on the curve", x, y);
        }
    }

    #[test]
    fn test_point_addition_in_finite_field() {
        // Test point addition with finite field arithmetic
        type F223 = PrimeField<223>;
        type P = Point<F223, 0, 7>;

        let p1 = P::new(F223::new(192).unwrap(), F223::new(105).unwrap()).unwrap();
        let p2 = P::new(F223::new(17).unwrap(), F223::new(56).unwrap()).unwrap();

        // Adding two points should produce another point on the curve
        let p3 = p1 + p2;
        assert!(!p3.is_infinity());

        // Verify the result is on the curve by extracting coordinates
        if let Point::Point(x, y) = p3 {
            let y_sqr = y.pow(2u128);
            let x_cbd = x.pow(3u128);
            let rhs = x_cbd + F223::from(7u128);
            assert_eq!(y_sqr, rhs);
        }
    }

    #[test]
    fn test_point_doubling() {
        type F223 = PrimeField<223>;
        type P = Point<F223, 0, 7>;

        let p = P::new(F223::new(192).unwrap(), F223::new(105).unwrap()).unwrap();

        // Doubling a point (adding it to itself)
        let doubled = p + p;

        if let Point::Point(x, y) = doubled {
            // Verify it's still on the curve
            let y_sqr = y.pow(2u128);
            let x_cbd = x.pow(3u128);
            let rhs = x_cbd + F223::from(7u128);
            assert_eq!(y_sqr, rhs);
        }
    }
}
