use num::BigInt;

#[derive(Debug, PartialEq, Eq, Clone)]
pub enum CurvePoint {
    PointAtInfinity,
    Point { x: BigInt, y: BigInt },
}

impl CurvePoint {
    pub fn new(x: impl Into<BigInt>, y: impl Into<BigInt>) -> Self {
        Self::Point {
            x: x.into(),
            y: y.into(),
        }
    }

    pub fn as_coordinates(&self) -> Option<(&BigInt, &BigInt)> {
        match self {
            CurvePoint::PointAtInfinity => None,
            CurvePoint::Point { x, y } => Some((x, y)),
        }
    }
}
