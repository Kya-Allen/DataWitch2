
use crate::zipsort;
use crate::linear;
use crate::estimators;
use std::vec::Vec;
use linear::Matrix;

// Dynamic Programming solution to the L1 Convex Clustering problem. Computes a full clusterpath, each cluster solution corresponding to one of the penalty parameters input as an argument to the function. 
// Based on the conference paper "Dynamic Visualization for L1 Fusion Convex Clustering in Near-Linear Time" by Binguan Zhang, Jie Chen, and Yoshikazu Terada
// published in Proceedings of the Thirty-Seventh Conference on Uncertainty in Artificial Intelligence (UAI 2021), PMLR 161:515–524
// 
// assume matrix columns are observations. 

pub fn c_paint(design_matrix: Matrix, path_length: usize) -> Vec<Matrix> {
    let n: usize = design_matrix.num_columns;
    let p: usize = design_matrix.num_rows;
    let design_transpose: Matrix = linear::get_transpose(&design_matrix); 
    let penalty_max: f64 =  1.0; // max_penalty(&design_transpose);
    let mut penalty_path: Vec<f64> = vec![0.0; path_length];
    penalty_path = penalty_path.iter().enumerate().map(|(i, x)| return (penalty_max * ((i+1) as f64 /path_length as f64))).collect();
    println!("penalty path: {:?}", penalty_path);
    // iterate over each dimension
    // for each dimension, we will recieve k vectors of partial centroids. k = path_length
    let mut cluster_solutions: Vec<Matrix> = vec![Matrix::new(Vec::with_capacity(n * p), n); path_length];
    for index in 0..p {
        // sort the values in a given dimension, and then save the ordering
        let mut sorted_dimension = design_transpose.get_col(index).to_vec();
        let mut ordering: Vec<f64> = (0..n).enumerate().map(|(i, _)|  i as f64).collect();
        (sorted_dimension, ordering) = zipsort::vecsort(&sorted_dimension, &ordering);


        // iterate over each penalty in the path
        // for each penalty, we recieve a vector of partial centroids
        let mut dimension_path: Vec<Vec<f64>> = Vec::with_capacity(path_length);
        for penalty in penalty_path.iter() {
            // use the dynamic programming solver
            // recieve a vector of partial centroids
            // let boop: f64 = 0.05;
            let partial_centroids: Vec<f64> = partial_clusterpath(&sorted_dimension, penalty, &n);
            println!("{:?}", partial_centroids);
            dimension_path.push(partial_centroids); 
        }
        // dimension_path: Vec<Vec<partial_centroids>>
        for path_step in 0..dimension_path.len() {
            // for each 'partial_centroids' in dimension path
            let temp_order: Vec<f64>;
            (dimension_path[path_step], temp_order) = zipsort::vecsort(&ordering, &dimension_path[path_step]);
        }
        for i in 0..cluster_solutions.len() {
            // start constructing cluster solution matrices where each column is a dimension of the data. Will need to be transposed
            // dimension_path[i].clone().into_iter().enumerate().map(|(i, x)| {cluster_solutions[i].data[i] = x; 0});
            cluster_solutions[i].data.append(&mut dimension_path[i])
        }

    }
    for matrix_index in 0..cluster_solutions.len() {
        cluster_solutions[matrix_index] = linear::get_transpose(&cluster_solutions[matrix_index]);
    }
    return cluster_solutions
}

/// Finds that maximum penalty value with a non-trivial solution, to determine where the path starts
pub fn max_penalty(transposed_matrix: &Matrix) -> f64 { 
    let nobs: usize = transposed_matrix.num_rows;
    let mut partial_maximums: Vec<f64> = Vec::with_capacity(transposed_matrix.num_columns);
    for dimension in 0..transposed_matrix.num_columns {
        let dimension: Vec<f64> = transposed_matrix.get_col(dimension).to_vec();
        let mean: f64 = estimators::mean(&dimension);
        let sum: f64 = dimension.iter().sum();
        let candidates: Vec<f64> = (0..nobs-1).enumerate().map(|(index, _)| {
            let partial_mean: f64 = dimension[0..index].iter().sum();
            (mean - (1 as f64/ (index+1) as f64) * partial_mean) / ((nobs - (index+1)) as f64) 
        }).collect();
        let partial_max: f64 = linear::vector_norm(&candidates, -1); // -1 get inf norm, which = max()
        println!("dimension candidates: {:?}", candidates);
        println!("dimension max: {:}", partial_max);
        partial_maximums.push(partial_max);
    };
    println!("interdimensional candidates: {:?}", partial_maximums);
    return linear::vector_norm(&partial_maximums[..], -1)
}

/// algorithm to find the partial centroids (The values of the centroids in a single dimension) given some penalty
pub fn partial_clusterpath(sorted_dimension: &[f64], penalty: &f64, n: &usize) -> Vec<f64> {
    let penalty_adjustment = |index: usize, penalty: &f64, n: &usize| -> f64 {
        (index+1) as f64 * (*n - (index+1)) as f64 * penalty
    };

    // PC = Partial Centroid
    let mut potential_pc_list: Vec<f64> = Vec::with_capacity(sorted_dimension.len());
    let adjusted_penalties: Vec<f64> = (0..n-1).enumerate().map(|(i, _)| penalty_adjustment(i, penalty, n)).collect();
    //println!("adjusted penalties: {:?}", adjusted_penalties);

    potential_pc_list = optimize_potential_pc(sorted_dimension, &adjusted_penalties, potential_pc_list);
    // solve a_n such that h_n(a_n) = 0
    let mut partial_centroids: Vec<f64> = vec![0.0; sorted_dimension.len()];
    partial_centroids[n-1] = last_partial_centroid(sorted_dimension, &adjusted_penalties, n);
    // find remaining a_i
    for k in (0..n-1).rev() {
        // a_k = max(a_k+1, U_k)
        partial_centroids[k] = partial_centroids[k+1].min(potential_pc_list[k])
    }
    return partial_centroids
}

pub fn first_potential_pc(first_value: &f64, penalty: &f64) -> f64 {
    return first_value + penalty
}

/// finds U_k as the value of beta that makes the derivative of the piecewise function that gives the centroid, equal to a penalty
pub fn optimize_potential_pc(sorted_dimension: &[f64], penalties: &[f64], mut potential_pc_list: Vec<f64>) -> Vec<f64> {
    for index in 0..potential_pc_list.capacity()-1 {
        let mut beta: f64;
        if index == 0 {
            potential_pc_list.push(first_potential_pc(&sorted_dimension[0], &penalties[0]));
        }
        else {
            // check right side of knot for the b > U_index-1 case
            beta = penalties[index] - penalties[index - 1] + sorted_dimension[index];
            if beta >= potential_pc_list[index-1] {
                potential_pc_list.push(beta);
            }
            else {
                // our U_index must be on the line to the left of the knot
                beta = penalties[index] + (sorted_dimension[index-1] + sorted_dimension[index]) / 2.0;
                potential_pc_list.push(beta);
            };
        }
    }
    return potential_pc_list
}

pub fn last_partial_centroid(sorted_dimension: &[f64], penalties: &[f64], n: &usize) -> f64 {
    let mut beta: f64 = sorted_dimension[n-1] - penalties[n-1-1]; 
    return beta
}

