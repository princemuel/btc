#![feature(const_trait_impl)]
#![feature(const_ops)]
#![feature(stmt_expr_attributes)]

mod field;
pub(crate) mod helpers;
mod macros;
mod ops;
mod point;

pub use field::{Field, PrimeField};
pub use point::Point;

#[cfg(test)]
mod tests {
    use crate::{Point, PrimeField};

    #[test]
    fn test_should_be_on_curve() {
        // tests the following points whether they are on the curve or not on curve
        // y^2=x^3-7 over F_223: (192,105) (17,56) (200,119) (1,193) (42,99).
        // the ones that aren't should fail
        type F223 = PrimeField<223>;
        type P = Point<F223, 0, 7>;

        let values = [(192, 105), (17, 56), (1, 193)];
        for (x, y) in values {
            let fx = F223::new(x).unwrap();
            let fy = F223::new(y).unwrap();

            let point = P::new(fx, fy);
            assert!(point.is_ok());
        }

        let values = [(200, 119), (42, 99)];
        for (x, y) in values {
            let fx = F223::new(x).unwrap();
            let fy = F223::new(y).unwrap();

            let point = P::new(fx, fy);
            assert!(point.is_err());
        }
    }

    #[test]
    fn test_point_addition_over_finite_fields() {}
}
