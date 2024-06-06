use num::{bigint::Sign, BigInt};

use crate::{extended_euclidean::extended_euclidean, CurvePoint};

pub struct Curve {
    generator: CurvePoint,
    field_modulus: i64,
    parameter_a: i64,
}

impl Curve {
    pub fn new(generator: CurvePoint, field_modulus: i64, parameter_a: i64) -> Self {
        Self {
            generator,
            field_modulus,
            parameter_a,
        }
    }

    /// Add two points on the elliptic curve.
    ///
    /// Formulas taken from https://en.wikipedia.org/wiki/Elliptic_curve_point_multiplication.
    pub fn add(&self, p: &CurvePoint, q: &CurvePoint) -> CurvePoint {
        match (p, q) {
            (CurvePoint::PointAtInfinity, CurvePoint::PointAtInfinity) => {
                CurvePoint::PointAtInfinity
            }
            (CurvePoint::PointAtInfinity, CurvePoint::Point { .. }) => *q,
            (CurvePoint::Point { .. }, CurvePoint::PointAtInfinity) => *p,
            (CurvePoint::Point { x: x_p, y: y_p }, CurvePoint::Point { x: x_q, y: y_q }) => {
                println!("Computing {p:?} + {q:?}");
                if p == q {
                    let lambda = (3 * x_p.pow(2) + self.parameter_a)
                        * (Self::modular_multiplicative_inverse(2 * y_p, self.field_modulus));
                    let lambda = Self::modulo(lambda, self.field_modulus);

                    let x_r = lambda.pow(2) - 2 * x_p;
                    let x_r = Self::modulo(x_r, self.field_modulus);

                    let y_r = lambda * (-x_r + x_p) - y_p;
                    let y_r = Self::modulo(y_r, self.field_modulus);

                    CurvePoint::Point { x: x_r, y: y_r }
                } else if x_p == x_q {
                    // If the x-coordinates match, there will be no intersection with a third point,
                    // so we return the point at infinity.
                    CurvePoint::PointAtInfinity
                } else {
                    let lambda = (y_q - y_p)
                        * Self::modular_multiplicative_inverse(x_q - x_p, self.field_modulus);
                    let lambda = Self::modulo(lambda, self.field_modulus);

                    let x_r = lambda.pow(2) - x_p - x_q;
                    let x_r = Self::modulo(x_r, self.field_modulus);

                    let y_r = lambda * (x_p - x_r) - y_p;
                    let y_r = Self::modulo(y_r, self.field_modulus);

                    CurvePoint::Point { x: x_r, y: y_r }
                }
            }
        }
    }

    /// Computes the modular multiplicative inverse `i` of `a mod b`, such that `a * i mod b = 1`.
    pub fn modular_multiplicative_inverse(a: i64, b: i64) -> i64 {
        let (sign, result) = extended_euclidean(BigInt::from(a), BigInt::from(b))
            .bezout_coefficient_a
            .to_u64_digits();

        let result = result.first().unwrap();
        if let Sign::Minus = sign {
            *result as i64 * -1
        } else {
            *result as i64
        }
    }

    /// Implements truncated division modulo.
    ///
    /// In the standard Rust modulo operator: -a % b = -a, but this function returns a result between 0 and `modulus`.
    pub fn modulo(a: i64, modulus: i64) -> i64 {
        if a >= 0 {
            a % modulus
        } else {
            let mut remainder = a;
            loop {
                let quotient = remainder / modulus;
                remainder = remainder - (quotient * modulus);

                if quotient == 0 {
                    if remainder == 0 {
                        return 0;
                    } else {
                        return modulus + remainder;
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn modulo() {
        let tests = [
            (12, 11, 1),
            (11, 11, 0),
            (0, 11, 0),
            (-11, 11, 0),
            (-10, 11, 1),
            (-15, 11, 7),
            (-37, 11, 7),
            (-111, 11, 10),
        ];

        for (a, modulus, expected_result) in tests {
            let actual_result = Curve::modulo(a, modulus);
            assert_eq!(
                expected_result, actual_result,
                "expected {a} mod {modulus} to be {expected_result} but was {actual_result}"
            );
        }
    }

    #[test]
    fn can_generate_all_points() {
        let expected_points = [
            CurvePoint::Point { x: 4, y: 10 },
            CurvePoint::Point { x: 7, y: 7 },
            CurvePoint::Point { x: 1, y: 9 },
            CurvePoint::Point { x: 0, y: 6 },
            CurvePoint::Point { x: 8, y: 8 },
            CurvePoint::Point { x: 2, y: 0 },
            CurvePoint::Point { x: 8, y: 3 },
            CurvePoint::Point { x: 0, y: 5 },
            CurvePoint::Point { x: 1, y: 2 },
            CurvePoint::Point { x: 7, y: 4 },
            CurvePoint::Point { x: 4, y: 1 },
            CurvePoint::PointAtInfinity,
            CurvePoint::Point { x: 4, y: 10 },
        ]
        .to_vec();

        // Generator Point for curve y^2 = x^3 + 3 mod 11.
        let curve = Curve::new(CurvePoint::new(4, 10), 11, 0);

        let mut computed_points = vec![curve.generator];
        let mut new_point = curve.generator;
        loop {
            new_point = curve.add(&curve.generator, &new_point);
            computed_points.push(new_point);
            if new_point == curve.generator {
                break;
            }
        }

        assert_eq!(expected_points, computed_points);
    }
}
