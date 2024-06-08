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

    /// Multiplies `scalar` with `p` in logarithmic time.
    pub fn multiply(&self, mut scalar: BigInt, p: CurvePoint) -> CurvePoint {
        if scalar == BigInt::ZERO {
            return CurvePoint::PointAtInfinity;
        }

        // The number of doublings we need to efficiently compute the multiplication.
        let doublings = (scalar.bits() - 1) as usize;

        // Create the doubling cache.
        // The ith entry in the cache is the result of 2^i * p.
        // +1 capacity to account for the 0th entry we add manually.
        let mut double_cache = Vec::with_capacity(doublings + 1);
        double_cache.push(p);

        for i in 1..=doublings {
            double_cache.push(self.add(double_cache[i - 1].clone(), double_cache[i - 1].clone()));
        }

        let mut result = CurvePoint::PointAtInfinity;
        let two = BigInt::from(2);
        while scalar != BigInt::ZERO {
            // SAFETY: Subtracting 1 is fine since we never enter the loop
            // if scalar is zero, and hence bits() always returns > 0.
            let next_smaller_power_of_two = scalar.bits() - 1;
            scalar -= two.pow(next_smaller_power_of_two as u32);
            result = self.add(
                result,
                double_cache[next_smaller_power_of_two as usize].clone(),
            );
        }

        result
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

    fn create_test_curve() -> Curve {
        // Generator Point for curve y^2 = x^3 + 3 mod 11.
        Curve::new(CurvePoint::new(4, 10), 11, 0)
    }

    fn create_bn128_curve() -> Curve {
        // Definition from https://eips.ethereum.org/EIPS/eip-197.
        let modulus = BigInt::parse_bytes(
            b"21888242871839275222246405745257275088696311157297823662689037894645226208583",
            10,
        )
        .expect("should be a valid base 10 number");
        Curve::new(CurvePoint::new(1, 2), modulus, 0)
    }

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

        let curve = create_test_curve();

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

    #[test]
    fn scalar_point_multiplication() {
        let curve = create_test_curve();
        // Equal to 2^8 + 2^4 + 2^2 + 2^0 to test doubling implementation.
        let scalar = 256 + 16 + 4 + 1;

        let multiplication = curve.multiply(scalar.into(), curve.generator.clone());

        let mut addition_result = CurvePoint::PointAtInfinity;
        for _ in 0..scalar {
            addition_result = curve.add(addition_result, curve.generator.clone());
        }

        assert_eq!(addition_result, multiplication);
    }

    #[test]
    fn scalar_point_multiplication_bn128() {
        // Computed with py_ecc.
        let scalar = 300_000_000;
        let expected_result = CurvePoint::new(
            BigInt::parse_bytes(
                b"12600240597266143967986535800884193324885833839429757878922176041119260815197",
                10,
            )
            .unwrap(),
            BigInt::parse_bytes(
                b"21411986724719982918952311537408507205322239197649094947485347628796002057456",
                10,
            )
            .unwrap(),
        );

        let bn128 = create_bn128_curve();
        let actual_result = bn128.multiply(scalar.into(), bn128.generator.clone());

        assert_eq!(expected_result, actual_result);
    }

    #[test]
    fn multiplication_is_associative() {
        let curve = create_bn128_curve();
        let left_first = curve.add(
            curve.add(
                curve.multiply(5.into(), curve.generator.clone()),
                curve.multiply(15.into(), curve.generator.clone()),
            ),
            curve.multiply(7.into(), curve.generator.clone()),
        );
        let right_first = curve.add(
            curve.multiply(5.into(), curve.generator.clone()),
            curve.add(
                curve.multiply(15.into(), curve.generator.clone()),
                curve.multiply(7.into(), curve.generator.clone()),
            ),
        );
        assert_eq!(left_first, right_first);
    }

    /// Test that 9/7 * G + 5/7 *G = 14/7 * G = 2 * G.
    /// This works because G * x = G * (x + curve_order), which implies (x + y) * G mod curve_order = x * G + y * G.
    #[test]
    fn encoding_rational_numbers() {
        let curve = create_bn128_curve();
        let curve_order_bn128 = BigInt::parse_bytes(
            b"21888242871839275222246405745257275088548364400416034343698204186575808495617",
            10,
        )
        .unwrap();

        // Encodes the rational number 9 * 1/7.
        let nine_over_seven = BigInt::from(9)
            * Curve::modular_multiplicative_inverse(7.into(), curve_order_bn128.clone());
        let nine_over_seven = Euclid::rem_euclid(&nine_over_seven, &curve_order_bn128);

        // Encodes the rational number 5 * 1/7.
        let one_half = BigInt::from(5)
            * Curve::modular_multiplicative_inverse(7.into(), curve_order_bn128.clone());
        let five_over_seven = Euclid::rem_euclid(&one_half, &curve_order_bn128);

        let whole_multiplication = curve.multiply(2.into(), curve.generator.clone());
        let rational_multiplication = curve.add(
            curve.multiply(nine_over_seven, curve.generator.clone()),
            curve.multiply(five_over_seven, curve.generator.clone()),
        );
        assert_eq!(whole_multiplication, rational_multiplication);
    }
}
