use num::{traits::Euclid, BigInt};

use crate::extended_euclidean;

/// Computes the modular multiplicative inverse `i` of `a mod b`, such that `a * i mod b = 1`.
///
/// ## Example
///
/// For example, the inverse of `2 mod 11` is `6` because `2 * 6 mod 11 = 1`.
///
/// ```
/// # use num::BigInt;
/// # use ecc::mod_mul_inverse;
/// assert_eq!(mod_mul_inverse(2.into(), 11.into()), BigInt::from(6));
/// ```
pub fn mod_mul_inverse(a: BigInt, b: BigInt) -> BigInt {
    Euclid::rem_euclid(&extended_euclidean(a, b.clone()).bezout_coefficient_a, &b)
}
