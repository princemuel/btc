use core::fmt;
use core::ops::Add;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Point<const A: i128, const B: i128> {
    Point(i128, i128),
    Infinity,
}

impl<const A: i128, const B: i128> Point<A, B> {
    pub fn new(x: i128, y: i128) -> Result<Self, &'static str> {
        // if y.pow(2) != x.pow(3) + (a * x) + b
        if y.pow(2) != x.pow(3) + A * x + B {
            return Err("Invalid parameters to Point::<5, 7>::new(). Point is not on the curve");
        }

        Ok(Self::Point(x, y))
    }

    pub const fn x(&self) -> Option<i128> {
        match self {
            Self::Point(x, _) => Some(*x),
            Self::Infinity => None,
        }
    }

    pub const fn y(&self) -> Option<i128> {
        match self {
            Self::Point(_, y) => Some(*y),
            Self::Infinity => None,
        }
    }

    pub const fn a(&self) -> i128 { A }

    pub const fn b(&self) -> i128 { B }
}

impl<const A: i128, const B: i128> Point<A, B> {
    pub const fn infinity() -> Self { Self::Infinity }

    pub const fn is_infinity(&self) -> bool { matches!(self, Self::Infinity) }
}

impl<const A: i128, const B: i128> Add for Point<A, B> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (Self::Infinity, Self::Infinity) => Self::infinity(),
            // Case 0.0: self is the point at infinity, return rhs
            (Self::Infinity, Self::Point(..)) => rhs,
            // Case 0.1: rhs is the point at infinity, return self
            (Self::Point(..), Self::Infinity) => self,
            // Case 1: self.x == rhs.x, self.y != rhs.y
            // The vertical line results in a point at infinity
            (Self::Point(x1, y1), Self::Point(x2, y2)) if x1 == x2 && y1 != y2 => Self::infinity(),
            // Case 2: self.x ≠ rhs.x
            // Formula (x, y) == (x1, y1) + (x2, y2)
            // s = (y2 - y1) / (x2 - x1)
            // x = s ** 2 - x1 - x2
            // y = s * (x1 - x) - y1
            (Self::Point(x1, y1), Self::Point(x2, y2)) if x1 != x2 => {
                let s = (y2 - y1) / (x2 - x1);
                let x = s.pow(2) - (x1 + x2);
                let y = s * (x1 - x) - y1;
                Self::Point(x, y)
            },
            // if we are tangent to the vertical line,
            (Self::Point(x1, y1), Self::Point(x2, y2)) if x1 == x2 && (y1 == 0 || y2 == 0) => {
                Self::infinity()
            },
            // Case 3: self == rhs
            // Formula (x, y) = (x1, y1) + (x1, y1)
            // s = (3 * x1 ** 2 + a) / (2 * y1)
            // x = s ** 2 - 2 * x1
            // y = s * (x1 - x) - y1
            (Self::Point(x1, y1), Self::Point(x2, y2)) if x1 == x2 && y1 == y2 => {
                let s = (3 * x1.pow(2) + self.a()) / (2 * y1);
                let x = s.pow(2) - 2 * x1;
                let y = s * (x1 - x) - y1;
                Self::Point(x, y)
            },
            _ => unreachable!(),
        }
    }
}

impl<const A: i128, const B: i128> fmt::Display for Point<A, B> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.x().is_none() || self.y().is_none() {
            return write!(f, "Point(infinity)");
        }
        let p = f.precision().unwrap_or(3);
        write!(
            f,
            "Point({:.*},{:.*})_{:.*}_{:.*}",
            p,
            self.x().unwrap_or_default(),
            p,
            self.y().unwrap_or_default(),
            p,
            self.a(),
            p,
            self.b()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_valid_point() {
        assert!(Point::<5, 7>::new(-1, -1).is_ok());
    }

    #[test]
    fn test_create_invalid_point() {
        assert!(Point::<5, 7>::new(-1, -2).is_err());
    }

    #[test]
    fn test_points_with_same_cords_are_equal() {
        let a = Point::<5, 7>::new(-1, -1).unwrap();
        let b = Point::<5, 7>::new(-1, -1).unwrap();

        assert_eq!(a, b);
    }

    #[test]
    fn test_points_with_different_cords_are_different() {
        let a = Point::<5, 7>::new(-1, -1).unwrap();
        let b = Point::<5, 7>::new(18, 77).unwrap();

        assert_ne!(a, b);
    }

    #[test]
    fn test_sum_inifity_and_another_point_returns_the_point() {
        let inf = Point::<5, 7>::infinity();
        let point = Point::<5, 7>::new(-1, -1).unwrap();

        assert_eq!(point + inf, point);
        assert_eq!(inf + point, point);
    }

    #[test]
    fn test_two_different_points_from_the_same_curve_can_be_added() {
        let a = Point::<5, 7>::new(2, 5).unwrap();
        let b = Point::<5, 7>::new(-1, -1).unwrap();

        let expected = Point::<5, 7>::new(3, -7).unwrap();
        assert_eq!(a + b, expected);
        assert_eq!(b + a, expected);
    }

    #[test]
    fn test_adding_the_same_point_results_infinity() {
        let a = Point::<5, 7>::new(-1, -1).unwrap();
        let b = Point::<5, 7>::new(-1, 1).unwrap();

        let expected = Point::<5, 7>::infinity();
        assert_eq!(a + b, expected);
    }
}
