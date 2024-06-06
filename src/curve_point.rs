#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub enum CurvePoint {
    PointAtInfinity,
    Point { x: i64, y: i64 },
}

impl CurvePoint {
    pub fn new(x: i64, y: i64) -> Self {
        Self::Point { x, y }
    }
}
