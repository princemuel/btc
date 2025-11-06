use core::fmt;
use core::ops::Add;

#[derive(Default, Debug, Clone, Copy)]
pub struct Point(Option<i128>, Option<i128>, i128, i128);

impl Point {
    pub fn new(x: Option<i128>, y: Option<i128>, a: i128, b: i128) -> Result<Self, String> {
        match (x, y) {
            (Some(x), Some(y)) => {
                if y.pow(2) != x.pow(3) + (a * x) + b {
                    return Err(format!("({x}, {y}) is not on the curve"));
                }
                Ok(Self(Some(x), Some(y), a, b))
            },
            (None, None) => Ok(Self(None, None, a, b)),
            _ => Err("Invalid parameters to Point::new()".to_owned()),
        }
    }

    pub const fn x(&self) -> Option<i128> { self.0 }

    pub const fn y(&self) -> Option<i128> { self.1 }

    pub const fn a(&self) -> i128 { self.2 }

    pub const fn b(&self) -> i128 { self.3 }
}

impl Point {
    pub const fn infinity(&self) -> Self { Self(None, None, self.a(), self.b()) }
}

impl PartialEq for Point {
    fn eq(&self, rhs: &Self) -> bool {
        self.x() == rhs.x() && self.y() == rhs.y() && self.a() == rhs.a() && self.b() == rhs.b()
    }
}

impl Add for Point {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        if self.a() != rhs.a() || self.b() != rhs.b() {
            panic!("Points {self}, {rhs} are not on the same curve");
        }

        // Case 0.0: self is the point at infinity, return rhs
        if self.x().is_none() {
            return rhs;
        }

        // Case 0.1: rhs is the point at infinity, return self
        if rhs.x().is_none() {
            return self;
        }

        // We can unwrap thegit se as we know they're not none
        let [x1, y1, x2, y2] = [
            self.x().unwrap(),
            self.y().unwrap(),
            rhs.x().unwrap(),
            rhs.y().unwrap(),
        ];

        // Case 1: self.x == rhs.x, self.y != rhs.y. The vertical
        // The vertical line results in a point at infinity
        if x1 == x2 && y1 != y2 {
            return self.infinity();
        }

        // Case 2: self.x ≠ rhs.x
        // Formula (x, y) == (x1, y1) + (x2, y2)
        // s = (y2 - y1) / (x2 - x1)
        // x = s ** 2 - x1 - x2
        // y = s * (x1 - x) - y1
        if x1 != x2 {
            let s = (y2 - y1) / (x2 - x1);
            let x = s.pow(2) - x1 - x2;
            let y = s * (x1 - x) - y1;
            return Self(Some(x), Some(y), self.a(), self.b());
        }

        // Case 3: self == rhs
        // Formula (x, y) = (x1, y1) + (x1, y1)
        // s = (3 * x1 ** 2 + a) / (2 * y1)
        // x = s ** 2 - 2 * x1
        // y = s * (x1 - x) - y1
        if self == rhs {
            // if we are tangent to the vertical line,
            if y1 == 0 {
                return self.infinity();
            }

            let s = (3 * x1.pow(2) + self.a()) / (2 * y1);
            let x = s.pow(2) - 2 * x1;
            let y = s * (x1 - x) - y1;
            return Self(Some(x), Some(y), self.a(), self.b());
        }

        unreachable!();
    }
}

impl fmt::Display for Point {
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
    fn test_ne() {
        let a = Point::new(Some(3), Some(-7), 5, 7).unwrap();
        let b = Point::new(Some(18), Some(77), 5, 7).unwrap();
        assert!(a != b);
        assert!(!(a != a));
    }

    #[test]
    fn test_on_curve() {
        let _ = Point::new(Some(3), Some(-7), 5, 7).unwrap();
        let _ = Point::new(Some(18), Some(77), 5, 7).unwrap();
    }

    #[test]
    #[should_panic]
    fn test_not_on_curve() { let _ = Point::new(Some(-2), Some(4), 5, 7).unwrap(); }

    #[test]
    fn test_add0() {
        let a = Point::new(None, None, 5, 7).unwrap();
        let b = Point::new(Some(2), Some(5), 5, 7).unwrap();
        let c = Point::new(Some(2), Some(-5), 5, 7).unwrap();

        assert_eq!(a + b, b);
        assert_eq!(b + a, b);
        assert_eq!(b + c, a);
    }

    #[test]
    fn test_add1() {
        let a = Point::new(Some(3), Some(7), 5, 7).unwrap();
        let b = Point::new(Some(-1), Some(-1), 5, 7).unwrap();
        assert_eq!(a + b, Point(Some(2), Some(-5), 5, 7));
    }

    #[test]
    fn test_add2() {
        let a = Point::new(Some(-1), Some(1), 5, 7).unwrap();
        assert_eq!(a + a, Point(Some(18), Some(-77), 5, 7));
    }
}
