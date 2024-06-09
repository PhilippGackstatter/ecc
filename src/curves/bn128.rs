use num::BigInt;

use crate::{CurvePoint, WeierstrassCurve};
use once_cell::sync::Lazy;

static GENERATOR: Lazy<CurvePoint<Bn128>> =
    Lazy::new(|| CurvePoint::new(BigInt::from(1), BigInt::from(2)));

static FIELD_MODULUS: Lazy<BigInt> = Lazy::new(|| {
    BigInt::parse_bytes(
        b"21888242871839275222246405745257275088696311157297823662689037894645226208583",
        10,
    )
    .unwrap()
});

/// Curve `bn128` as defined in https://eips.ethereum.org/EIPS/eip-197.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bn128;

impl WeierstrassCurve for Bn128 {
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
