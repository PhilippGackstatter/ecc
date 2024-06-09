use std::{
    marker::PhantomData,
    ops::{Add, Mul},
};

use num::{traits::Euclid, BigInt};

use crate::{mod_mul_inverse, WeierstrassCurve};

#[derive(Debug, PartialEq, Eq)]
pub struct CurvePoint<C: WeierstrassCurve> {
    point: Point,
    phantom: PhantomData<C>,
}

#[derive(Debug, PartialEq, Eq, Clone)]
pub enum Point {
    PointAtInfinity,
    Point { x: BigInt, y: BigInt },
}

impl<C: WeierstrassCurve> From<Point> for CurvePoint<C> {
    fn from(point: Point) -> Self {
        Self {
            point,
            phantom: PhantomData,
        }
    }
}

impl<C: WeierstrassCurve> Clone for CurvePoint<C> {
    fn clone(&self) -> Self {
        self.point.clone().into()
    }
}

impl<C: WeierstrassCurve> CurvePoint<C> {
    /// Creates a new point on the curve with the given coordinates.
    ///
    /// At present, it does not check whether the point is actually on the curve.
    pub fn new(x: impl Into<BigInt>, y: impl Into<BigInt>) -> Self {
        Point::Point {
            x: x.into(),
            y: y.into(),
        }
        .into()
    }

    /// Returns the underlying point's x, y coordinates in that order
    /// or `None` if it's the point at infinity.
    pub fn as_coordinates(&self) -> Option<(&BigInt, &BigInt)> {
        match &self.point {
            Point::PointAtInfinity => None,
            Point::Point { x, y } => Some((x, y)),
        }
    }

    /// Returns the underlying [`Point`].
    pub fn point(&self) -> &Point {
        &self.point
    }

    /// Creates the `CurvePoint` representing the point at infinity, i.e. the identity element.
    pub fn point_at_infinity() -> Self {
        Self::from(Point::PointAtInfinity)
    }

    /// Multiplies `scalar` with `p` in logarithmic time.
    fn multiply(&self, scalar: &BigInt) -> CurvePoint<C> {
        if scalar == &BigInt::ZERO {
            return Point::PointAtInfinity.into();
        }

        // The number of doublings we need to efficiently compute the multiplication.
        let doublings = (scalar.bits() - 1) as usize;

        // Create the doubling cache.
        // The ith entry in the cache is the result of 2^i * p.
        // +1 capacity to account for the 0th entry we add manually.
        let mut double_cache: Vec<CurvePoint<C>> = Vec::with_capacity(doublings + 1);
        double_cache.push(self.clone());

        for i in 1..=doublings {
            double_cache.push(&double_cache[i - 1] + &double_cache[i - 1]);
        }

        let mut result = Point::PointAtInfinity.into();
        let mut scalar = scalar.clone();
        let two = BigInt::from(2);

        while scalar != BigInt::ZERO {
            // SAFETY: Subtracting 1 is fine since we never enter the loop
            // if scalar is zero, and hence bits() always returns > 0.
            let next_smaller_power_of_two = scalar.bits() - 1;
            scalar -= two.pow(next_smaller_power_of_two as u32);
            result = &result + &double_cache[next_smaller_power_of_two as usize];
        }

        result
    }

    /// Adds `q` to `self` on the elliptic curve.
    ///
    /// Formulas taken from https://en.wikipedia.org/wiki/Elliptic_curve_point_multiplication.
    fn add(&self, q: &CurvePoint<C>) -> CurvePoint<C> {
        let p = &self.point;
        let q = &q.point;

        match (p, q) {
            (Point::PointAtInfinity, Point::PointAtInfinity) => Point::PointAtInfinity.into(),
            (Point::PointAtInfinity, Point::Point { .. }) => q.clone().into(),
            (Point::Point { .. }, Point::PointAtInfinity) => p.clone().into(),
            (Point::Point { x: x_p, y: y_p }, Point::Point { x: x_q, y: y_q }) => {
                if p == q {
                    let lambda =
                        (3 * x_p.pow(2) + C::a()) * (mod_mul_inverse(2 * y_p, C::field_modulus()));
                    let lambda = Euclid::rem_euclid(&lambda, &C::field_modulus());

                    let x_r = lambda.pow(2) - 2 * x_p;
                    let x_r = Euclid::rem_euclid(&x_r, &C::field_modulus());

                    let y_r = lambda * (-&x_r + x_p) - y_p;
                    let y_r = Euclid::rem_euclid(&y_r, &C::field_modulus());

                    Point::Point { x: x_r, y: y_r }.into()
                } else if x_p == x_q {
                    // If the x-coordinates match, there will be no intersection with a third point,
                    // so we return the point at infinity.
                    Point::PointAtInfinity.into()
                } else {
                    let lambda = (y_q - y_p) * mod_mul_inverse(x_q - x_p, C::field_modulus());
                    let lambda = Euclid::rem_euclid(&lambda, &C::field_modulus());

                    let x_r = lambda.pow(2) - x_p - x_q;
                    let x_r = Euclid::rem_euclid(&x_r, &C::field_modulus());

                    let y_r = lambda * (x_p - &x_r) - y_p;
                    let y_r = Euclid::rem_euclid(&y_r, &C::field_modulus());

                    Point::Point { x: x_r, y: y_r }.into()
                }
            }
        }
    }

    /// Returns the inverse `inv` of `self` such that `self` + `inv` equals the [`Point::PointAtInfinity`].
    pub fn negate(&self) -> CurvePoint<C> {
        let Point::Point { x, y } = &self.point else {
            // The inverse of the point at infinity is itself.
            return Point::PointAtInfinity.into();
        };

        // The inverse of a point has the same x-coordinate,
        // and the new y value satisfies old_y * new_y mod field_modulus = 1, i.e. the modular multiplicate inverse.
        Point::Point {
            x: x.clone(),
            y: mod_mul_inverse(y.clone(), C::field_modulus()),
        }
        .into()
    }
}

impl<C: WeierstrassCurve> Add<&CurvePoint<C>> for &CurvePoint<C> {
    type Output = CurvePoint<C>;

    fn add(self, q: &CurvePoint<C>) -> Self::Output {
        CurvePoint::add(self, q)
    }
}

// Additional implementation for convenience.
impl<C: WeierstrassCurve> Add<&CurvePoint<C>> for CurvePoint<C> {
    type Output = CurvePoint<C>;

    fn add(self, q: &CurvePoint<C>) -> Self::Output {
        CurvePoint::add(&self, q)
    }
}

impl<C: WeierstrassCurve> Mul<&BigInt> for &CurvePoint<C> {
    type Output = CurvePoint<C>;

    fn mul(self, scalar: &BigInt) -> Self::Output {
        CurvePoint::multiply(self, scalar)
    }
}

// Additional implementation for convenience.
impl<C: WeierstrassCurve> Mul<&BigInt> for CurvePoint<C> {
    type Output = CurvePoint<C>;

    fn mul(self, scalar: &BigInt) -> Self::Output {
        CurvePoint::multiply(&self, scalar)
    }
}
