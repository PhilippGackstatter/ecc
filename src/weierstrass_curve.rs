use crate::CurvePoint;
use num::BigInt;

/// Parameter definitions for Weierstrass elliptic curves.
pub trait WeierstrassCurve {
    /// Returns the generator point of the curve.
    fn generator() -> CurvePoint<Self>
    where
        Self: Sized;
    /// Returns the parameter `a` of the curve.
    fn a() -> BigInt;
    /// Returns the field modulus of the curve.
    fn field_modulus() -> BigInt;
}

#[cfg(test)]
mod tests {
    use std::ops::Add;

    use k256::elliptic_curve::{self, sec1::ToEncodedPoint};
    use num::{bigint::Sign, traits::Euclid};
    use once_cell::sync::Lazy;

    use crate::{
        curves::{Bn128, Secp256k1},
        mod_mul_inverse, Point,
    };

    use super::*;

    static GENERATOR: Lazy<CurvePoint<TestCurve>> =
        Lazy::new(|| CurvePoint::new(BigInt::from(4), BigInt::from(10)));
    static FIELD_MODULUS: Lazy<BigInt> = Lazy::new(|| BigInt::from(11));
    /// A test curve for initial testing with a small modulus.
    #[derive(Debug, PartialEq, Eq)]
    struct TestCurve;
    impl WeierstrassCurve for TestCurve {
        fn generator() -> CurvePoint<Self>
        where
            Self: Sized,
        {
            GENERATOR.clone()
        }

        fn a() -> BigInt {
            BigInt::ZERO
        }

        fn field_modulus() -> BigInt {
            FIELD_MODULUS.clone()
        }
    }

    #[test]
    fn can_generate_all_test_curve_points() {
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
            CurvePoint::point_at_infinity(),
            CurvePoint::new(4, 10),
        ]
        .to_vec();

        let mut computed_points = vec![TestCurve::generator().clone()];
        let mut new_point = TestCurve::generator().clone();
        loop {
            new_point = TestCurve::generator() + &new_point;
            computed_points.push(new_point.clone());
            if new_point == TestCurve::generator() {
                break;
            }
        }

        assert_eq!(expected_points, computed_points);
    }

    #[test]
    fn scalar_point_multiplication() {
        // Equal to 2^8 + 2^4 + 2^2 + 2^0 to test doubling implementation.
        let scalar = 256 + 16 + 4 + 1;

        let multiplication = TestCurve::generator() * &scalar.into();

        let mut addition_result = CurvePoint::point_at_infinity();
        for _ in 0..scalar {
            addition_result = addition_result + &TestCurve::generator();
        }

        assert_eq!(addition_result, multiplication);
    }

    #[test]
    fn scalar_point_multiplication_bn128() {
        // Expected result computed with py_ecc.
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

        let actual_result = Bn128::generator() * &scalar.into();

        assert_eq!(expected_result, actual_result);
    }

    #[test]
    fn negate_points() {
        let random_point = Bn128::generator() * &5000.into();
        let negated = random_point.negate();

        assert_eq!(random_point + &negated, Point::PointAtInfinity.into());

        assert_eq!(
            CurvePoint::<Bn128>::from(Point::PointAtInfinity)
                + &(CurvePoint::from(Point::PointAtInfinity).negate()),
            Point::PointAtInfinity.into()
        );
    }

    #[test]
    fn multiplication_is_associative() {
        let left_first = Add::add(
            Bn128::generator() * &5.into() + &(Bn128::generator() * &15.into()),
            &(Bn128::generator() * &7.into()),
        );

        let right_first = Add::add(
            Bn128::generator() * &5.into(),
            &(Bn128::generator() * &15.into() + &(Bn128::generator() * &7.into())),
        );

        assert_eq!(left_first, right_first);
    }

    /// Test that 9/7 * G + 5/7 *G = 14/7 * G = 2 * G.
    /// This works because G * x = G * (x + curve_order), which implies (x + y) * G mod curve_order = x * G + y * G.
    #[test]
    fn encoding_rational_numbers() {
        let curve_order_bn128 = BigInt::parse_bytes(
            b"21888242871839275222246405745257275088548364400416034343698204186575808495617",
            10,
        )
        .unwrap();

        // Encodes the rational number 9 * 1/7.
        let nine_over_seven =
            BigInt::from(9) * mod_mul_inverse(7.into(), curve_order_bn128.clone());
        let nine_over_seven = Euclid::rem_euclid(&nine_over_seven, &curve_order_bn128);

        // Encodes the rational number 5 * 1/7.
        let one_half = BigInt::from(5) * mod_mul_inverse(7.into(), curve_order_bn128.clone());
        let five_over_seven = Euclid::rem_euclid(&one_half, &curve_order_bn128);

        let whole_multiplication = Bn128::generator() * &2.into();
        let rational_multiplication =
            Bn128::generator() * &nine_over_seven + &(Bn128::generator() * &five_over_seven);

        assert_eq!(whole_multiplication, rational_multiplication);
    }

    /// Test that the public key computation run by k256 and this libary are equivalent.
    #[test]
    fn elliptic_curve_public_key_calculation() {
        let secret: [u8; 32] = rand::random();

        let secret_int = BigInt::from_bytes_be(Sign::Plus, &secret);
        let secret_key = k256::SecretKey::from_slice(&secret).unwrap();
        let public_key = secret_key.public_key();

        let curve_pk = Secp256k1::generator() * &secret_int;

        let Point::Point { x, y } = &curve_pk.point() else {
            panic!("expected regular point");
        };

        let encoded_point = public_key.to_encoded_point(false);
        let public_key_x = BigInt::from_bytes_be(Sign::Plus, encoded_point.x().unwrap());
        let public_key_y = BigInt::from_bytes_be(Sign::Plus, encoded_point.y().unwrap());

        assert_eq!(x, &public_key_x);
        assert_eq!(y, &public_key_y);
    }

    /// Test that ECDH run by k256 and this libary are equivalent.
    #[test]
    fn elliptic_curve_diffie_hellman() {
        let sk1: [u8; 32] = rand::random();
        let sk2: [u8; 32] = rand::random();

        let sk1_int = BigInt::from_bytes_be(num::bigint::Sign::Plus, &sk1);
        let sk2_int = BigInt::from_bytes_be(num::bigint::Sign::Plus, &sk2);

        // Run ECDH with the k256 library.
        let sk1 = k256::SecretKey::from_slice(&sk1).unwrap();
        let sk2 = k256::SecretKey::from_slice(&sk2).unwrap();

        let pk1 = sk1.public_key();
        let pk2 = sk2.public_key();

        let shared_secret1 =
            elliptic_curve::ecdh::diffie_hellman(sk1.to_nonzero_scalar(), pk2.as_affine());
        let shared_secret2 =
            elliptic_curve::ecdh::diffie_hellman(sk2.to_nonzero_scalar(), pk1.as_affine());

        assert_eq!(
            shared_secret1.raw_secret_bytes(),
            shared_secret2.raw_secret_bytes()
        );

        // Run ECDH with this library.
        let curve_pk1 = Secp256k1::generator() * &sk1_int;
        let curve_pk2 = Secp256k1::generator() * &sk2_int;

        let curve_shared_secret1 = curve_pk1 * &sk2_int;
        let curve_shared_secret2 = curve_pk2 * &sk1_int;

        assert_eq!(curve_shared_secret1, curve_shared_secret2);

        // Assert that ECDH was executed the same by both libraries.
        let shared_secret_x = shared_secret1.raw_secret_bytes().as_slice();
        let curve_shared_secret_x = curve_shared_secret1
            .as_coordinates()
            .unwrap()
            .0
            .to_bytes_be()
            .1;

        assert_eq!(shared_secret_x, curve_shared_secret_x);
    }
}
