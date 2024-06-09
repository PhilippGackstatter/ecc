use num::BigInt;

use crate::{CurvePoint, WeierstrassCurve};
use once_cell::sync::Lazy;

static GENERATOR: Lazy<CurvePoint<Secp256k1>> = Lazy::new(|| {
    let x = BigInt::parse_bytes(
        b"79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798",
        16,
    )
    .unwrap();
    let y = BigInt::parse_bytes(
        b"483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8",
        16,
    )
    .unwrap();

    CurvePoint::new(x, y)
});

static FIELD_MODULUS: Lazy<BigInt> = Lazy::new(|| {
    BigInt::parse_bytes(
        b"fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f",
        16,
    )
    .unwrap()
});

/// Curve secp256k1 as defined in <http://www.secg.org/sec2-v2.pdf>.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Secp256k1;

impl WeierstrassCurve for Secp256k1 {
    fn generator() -> CurvePoint<Self> {
        GENERATOR.clone()
    }

    fn a() -> BigInt {
        BigInt::ZERO
    }

    fn field_modulus() -> BigInt {
        FIELD_MODULUS.clone()
    }
}
