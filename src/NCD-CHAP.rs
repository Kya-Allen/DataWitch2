use crate::linear;
use std::vec::Vec;
use linear::Matrix;

pub struct Point2d {
    x: f64,
    y: f64
}

pub struct Point3d {
    x: f64,
    y: f64,
    z: f64
}

pub struct IndexedPoint {
    index: u32,
    point: Point2d
}

pub struct ConvexSet {
    data: Vec<Point2d>,
    boundary: Vec<Point2d>,
}

impl ConvexSet {
    pub fn graham_scan() {

    }
}