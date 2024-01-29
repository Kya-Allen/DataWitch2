use crate::zipsort;
use crate::linear;
use crate::rng;
use crate::estimators;
use std::vec::Vec;
use linear::Matrix;

pub struct Point {
    position: Vec<f64>,
    id: usize,
    clusters: Vec<Vec<f64>>
}

pub struct Distance {
    point: Point,
    distance: f64
}

impl Distance {
    pub fn new(point: Point, distance: f64) -> Self {
        Self {
            point,
            distance,
        }
    }
}

pub fn k_nearest_neighbors(k: usize, data: Vec<Point>, query: Vec<f64>) -> Vec<Distance> {
    let mut nearest_neighbors: Vec<Distance> = Vec::with_capacity(k);
    let mut distances: Vec<Distance> = Vec::with_capacity(data.len());
    for point in data {
        let num = linear::euclidean_distance(&query, &point.position);
        let distance: Distance = Distance::new(point, num);
        nearest_neighbors.push(distance)
    }
    distances.sort_by(|a, b| b.distance.partial_cmp(&a.distance).unwrap());
    for index in 0..k {
        nearest_neighbors.push(distances.swap_remove(distances.len()-1))
    }
    return nearest_neighbors
}

pub fn k_means(k: usize, data: Vec<Point>) -> Vec<Point> {
    return data
}