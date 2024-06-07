use num::{traits::Euclid, BigInt};

use crate::{extended_euclidean::extended_euclidean, CurvePoint};

pub struct Curve {
    pub generator: CurvePoint,
    pub field_modulus: BigInt,
    pub parameter_a: BigInt,
}

impl Curve {
    pub fn new(
        generator: CurvePoint,
        field_modulus: impl Into<BigInt>,
        parameter_a: impl Into<BigInt>,
    ) -> Self {
        Self {
            generator,
            field_modulus: field_modulus.into(),
            parameter_a: parameter_a.into(),
        }
    }

    /// Add two points on the elliptic curve.
    ///
    /// Formulas taken from https://en.wikipedia.org/wiki/Elliptic_curve_point_multiplication.
    pub fn add(&self, p: CurvePoint, q: CurvePoint) -> CurvePoint {
        match (&p, &q) {
            (CurvePoint::PointAtInfinity, CurvePoint::PointAtInfinity) => {
                CurvePoint::PointAtInfinity
            }
            (CurvePoint::PointAtInfinity, CurvePoint::Point { .. }) => q,
            (CurvePoint::Point { .. }, CurvePoint::PointAtInfinity) => p,
            (CurvePoint::Point { x: x_p, y: y_p }, CurvePoint::Point { x: x_q, y: y_q }) => {
                println!("Computing {p:?} + {q:?}");

                if p == q {
                    let lambda = (3 * x_p.pow(2) + &self.parameter_a)
                        * (Self::modular_multiplicative_inverse(
                            2 * y_p,
                            self.field_modulus.clone(),
                        ));
                    let lambda = Euclid::rem_euclid(&lambda, &self.field_modulus);

                    let x_r = lambda.pow(2) - 2 * x_p;
                    let x_r = Euclid::rem_euclid(&x_r, &self.field_modulus);

                    let y_r = lambda * (-&x_r + x_p) - y_p;
                    let y_r = Euclid::rem_euclid(&y_r, &self.field_modulus);

                    CurvePoint::Point { x: x_r, y: y_r }
                } else if x_p == x_q {
                    // If the x-coordinates match, there will be no intersection with a third point,
                    // so we return the point at infinity.
                    CurvePoint::PointAtInfinity
                } else {
                    let lambda = (y_q - y_p)
                        * Self::modular_multiplicative_inverse(
                            x_q - x_p,
                            self.field_modulus.clone(),
                        );
                    let lambda = Euclid::rem_euclid(&lambda, &self.field_modulus);

                    let x_r = lambda.pow(2) - x_p - x_q;
                    let x_r = Euclid::rem_euclid(&x_r, &self.field_modulus);

                    let y_r = lambda * (x_p - &x_r) - y_p;
                    let y_r = Euclid::rem_euclid(&y_r, &self.field_modulus);

                    CurvePoint::Point { x: x_r, y: y_r }
                }
            }
        }
    }

    /// Computes the modular multiplicative inverse `i` of `a mod b`, such that `a * i mod b = 1`.
    pub fn modular_multiplicative_inverse(a: BigInt, b: BigInt) -> BigInt {
        extended_euclidean(a, b).bezout_coefficient_a
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn can_generate_all_points() {
        let expected_points = [
            CurvePoint::new(4, 10),
            CurvePoint::new(7, 7),
            CurvePoint::new(1, 9),
            CurvePoint::new(0, 6),
            CurvePoint::new(8, 8),
            CurvePoint::new(2, 0),
            CurvePoint::new(8, 3),
            CurvePoint::new(0, 5),
            CurvePoint::new(1, 2),
            CurvePoint::new(7, 4),
            CurvePoint::new(4, 1),
            CurvePoint::PointAtInfinity,
            CurvePoint::new(4, 10),
        ]
        .to_vec();

        // Generator Point for curve y^2 = x^3 + 3 mod 11.
        let curve = Curve::new(CurvePoint::new(4, 10), 11, 0);

        let mut computed_points = vec![curve.generator.clone()];
        let mut new_point = curve.generator.clone();
        loop {
            new_point = curve.add(curve.generator.clone(), new_point);
            computed_points.push(new_point.clone());
            if new_point == curve.generator {
                break;
            }
        }

        assert_eq!(expected_points, computed_points);
    }
}
