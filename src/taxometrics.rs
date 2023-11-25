use crate::zipsort;
use crate::linear;
use crate::rng;
use crate::estimators;
use std::vec::Vec;
use linear::Matrix;

fn ccfi(empirical_curve: &[f64], categorical_simulation: &[f64], dimensional_simulation: &[f64]) -> f64 {
    let dimensional_fit = estimators::rmse(empirical_curve, dimensional_simulation); 
    let categorical_fit = estimators::rmse(empirical_curve, categorical_simulation);
    let index: f64 = dimensional_fit / (categorical_fit + dimensional_fit);
    return index
}

pub trait TaxometricFit {
    fn fit(&self, partition_variable: Vec<f64>, num_partitions: usize) -> Vec<f64>;
}

pub struct Mambac {
    data_matrix: linear::Matrix,
    curve: Vec<f64>,
    ccfi: f64,
}

impl Mambac {
    fn fit(&mut self, partition_variable: Vec<f64>, comparison_variable: Vec<f64>, num_partitions: usize) -> Vec<f64> {
        // defined for two variables at a time
        // where the number of partitions is the number of divisions, not the number of groups
        let partition_space: f64 = linear::vector_norm(&partition_variable, -1) - linear::minimum(&partition_variable);
        let (sorted_partition, sorted_comparison): (Vec<f64>, Vec<f64>) = zipsort::vecsort(&partition_variable[..], &comparison_variable[..]);
        let partition_interval: f64 = partition_space / (num_partitions + 1) as f64;
        let mut results = Vec::with_capacity(num_partitions);
        let mut partition_pointer: f64 = 0.0;
        for partition in 0..num_partitions {
            partition_pointer += partition_interval;
            let mut index_pointer: usize = 0;
            while sorted_partition[index_pointer] <= partition_pointer {
                index_pointer += index_pointer;
            }
            let result = estimators::mean_difference(&sorted_comparison[0..index_pointer], &sorted_comparison[index_pointer..sorted_comparison.len()]);
            results.push(result)
        }
        self.curve = results;
        return self.curve
    }
}

pub struct MaxEig {
    data_matrix: linear::Matrix,
    curve: Vec<f64>,
    ccfi: f64
}

impl MaxEig {
    fn covariance_matrix(&self) {
        // diagonals all zero
    }

    /// returns key and complement, sorted by key
    fn sort_data(&self, indicator_index_p: usize) -> (Vec<f64>, Vec<Vec<f64>>) {
        // switch key to front of data matrix

        //return zipsort::vecsort(indicator, &complement);
        let data_matrix_transpose: Matrix = linear::get_transpose(&self.data_matrix);
        let indicator: &[f64] = data_matrix_transpose.get_col(indicator_index_p);
        let mut complement: Vec<&[f64]> = Vec::with_capacity(data_matrix_transpose.data.len() - data_matrix_transpose.num_rows);
        for index in 0..data_matrix_transpose.num_columns {
            if index != indicator_index_p {
                complement.push(data_matrix_transpose.get_col(index))
            }
        }; 
        let (x, y) = zipsort::matsort(indicator, complement);
        return  (x, y)
    }

    fn fit(&mut self, partition_variable: usize, num_partitions: u32) -> Vec<f64> {
        // based on finding the maximum eigenvalue
            // perfect use case for linalg::power_iteration()
        let (indicator, complement) = self.sort_data(partition_variable);
        let partition_space: f64 = indicator.iter().max() - indicator.iter().min();
        let partition_interval: f64 = partition_space / num_partitions + 1;
        let mut results = Vec::with_capacity(num_partitions);
        // for each parition subspace of the indicator variable, do:
        let mut partition_pointer: f64 = 0;
        for partition in 0..num_partitions {
            partition_pointer = partition_pointer + partition_interval
        }
            // generate covariance matrix &self::covariance_matrix()
            // find maximum eigenvalue linalg::power_iteration() -> rayleigh_quotient()
            // associate with mean value of partition
    }
}

