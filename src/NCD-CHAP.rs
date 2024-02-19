use crate::linear;
use std::vec::Vec;
use linear::Matrix;


pub enum Point<T> {
    Point2d {x: T, y: T},
    Point3d {x: T, y: T, z: T}
}

pub struct DataPoint<T> {
    index: usize,
    point: Point<T>
}

pub struct ConvexSet<T> {
    data: Vec<Point<T>>,
    size: usize,
    boundary: Vec<Point<T>>
}

impl<T> ConvexSet<T> {
    /// Constructor 
    pub fn new(data: Vec<Point<T>>) -> Self {
        Self {
            data,
            size: data.len(),
            boundary: Vec::with_capacity(data.len() / 3)
        }
    }

    pub fn compute_boundary(&self, subset_size: usize) -> &self::boundary {
        let num_partitions: usize = size / subset_size; 
        let mut sub_boundaries: Vec<Vec<Point<T>>> = Vec::with_capacity(partitions);
        let starting_point = starting_point();
        let subsets: Vec<Vec<Point<T>>> = partition_dataset(&self, num_partitions);
        let mut index: usize = 0;
        sub_boundaries = multi_graham_scan(&num_partitions, index, &subsets, sub_boundaries)
    }

    fn partition_dataset() -> Vec<Vec<Point>> {
        // arbitrarily partition the set of points into K = n/m subsets, with at most m points each.
    }

    pub fn single_graham_scan(self) {
        let starting_point: Point<T> = self.starting_point();
        // let sorted_set: Vec<Point<T>> = sort set in increasing order of the angle points moake with P
    }

    /// Picks a point in the dataset guaranteed to be on the convex hull
    fn starting_point(&self) -> Point<T> {
        let mut starting_point: Point<T> = self.data[0];
        for (index, point) in self.data.enumerate() {
            if (point.x < starting_point.x) { starting_point = point.copy() }
        }
        return starting_point
    }
}

/// recursive function computes graham scan for each subset of points
pub fn multi_graham_scan(num_partitions: &usize, index: usize, subsets: &Vec<Vec<Point<T>>>, sub_boudaries: Vec<Vec<Point<T>>>) -> Vec<Vec<Point<T>>> {
    //do graham scan
    if (index+1 < num_partitions) { 
        index = index + 1;
        return multi_graham_scan(num_partitions, index, sub_boundaries)
    }
    return sub_boundaries
}
