use crate::zipsort;
use crate::linear;
use crate::rng;
use crate::estimators;
use std::vec::Vec;
use linear::Matrix;

/// treat columns of design_matrix as observations. Each column consists of every dimension for a single observation
pub fn clusterpath_l2(design_matrix: Matrix, weight: f64, mut penalty: f64) -> (Vec<Matrix>, Vec<f64>) {
    let mut centroids: Matrix = design_matrix.clone();
    let n: usize = centroids.num_columns;
    let mut clusters: Vec<usize> = (0..n).enumerate().map(|(x, _)| x).collect();

    let mut solutions: Vec<Matrix> = Vec::with_capacity(10);
    let mut penalties: Vec<f64> = Vec::with_capacity(10);

    while clusters.len() > 1 {
        // (centroids, clusters) = solve_l2()

        solutions.push(centroids.clone());
        penalties.push(penalty);
        penalty = penalty * 1.5;
        clusters = (0..n).enumerate().map(|(x, _)| x).collect();
    }
    return (solutions, penalties)
}

fn solve_l2(centroids: &mut Matrix, clusters: Vec<usize>, design_matrix: &Matrix, weight: &f64, penalty: &f64) -> () {
    let tolerance: f64 = 0.01;
    let step_size: f64 = 0.5;
    let mut gradients: Matrix = Matrix::new(Vec::with_capacity(clusters.len() * centroids.num_rows), clusters.len());
    for cluster_index in 0..clusters.len() {
        let mut cluster_gradient = subgradient_l2(&centroids, &clusters, design_matrix, &cluster_index, penalty);
        gradients.data.append(&mut cluster_gradient);
    } 
    let mut iter: usize = 1;
    while gradients.frobenius().powf(2.0) > tolerance {
        let step: f64 = step_size / (iter as f64);
        let direction: Vec<f64> = gradients.data.iter().map(|x| x*step).collect();
        centroids.data = linear::vector_subtract(&centroids.data, &direction);

        let mut temporary_grad: Vec<f64> = Vec::with_capacity(clusters.len());
        for cluster_index in 0..clusters.len() {
            let mut cluster_gradient = subgradient_l2(&centroids, &clusters, design_matrix, &cluster_index, penalty);
            temporary_grad.append(&mut cluster_gradient);
        } 
        gradients.data = temporary_grad.clone();
        iter += 1;
    }
    //finish this section
}

fn subgradient_l2(centroids: &Matrix, clusters: &[usize], design_matrix: &Matrix, cluster_index: &usize, penalty: &f64) -> Vec<f64> {
    // composed of gradient for the fission term + subgradient for the fusion 
    let centroid = centroids.get_col(cluster_index.clone());
    let fission_term = fission_gradient(centroid, design_matrix, clusters);
    let fusion_term = fusion_subgradient(centroid, centroids, clusters, cluster_index, penalty, design_matrix);
    let subgradient: Vec<f64> = fusion_term.iter().enumerate().map(|(i, x)| x + fission_term[i]).collect();

    return subgradient
}

fn fusion_subgradient(centroid: &[f64], centroids: &Matrix, clusters: &[usize], cluster_index: &usize, penalty: &f64, design_matrix: &Matrix) -> Vec<f64> {
    let cluster_cardinality: usize = clusters.len();
    let penalty_scalar: f64 = penalty / (*cluster_index as f64);
    let mut penalty_term: Vec<f64> = vec![0.0; cluster_cardinality]; 

    for cluster_jindex in 0..cluster_cardinality {
        if cluster_jindex == *cluster_index { continue };
        let numerator = linear::vector_subtract(centroid, centroids.get_col(cluster_jindex));
        let denomenator: f64 = linear::l2_norm(&numerator);
        let weight: f64 = (0..cluster_cardinality).enumerate().map(|(i, _)| {
            exponential_weight(design_matrix.get_col(i), design_matrix.get_col(cluster_jindex))
        }).collect::<Vec<f64>>().iter().sum();

        let divided: Vec<f64> = numerator.iter().map(|x| x/denomenator).collect();
        let weighted: Vec<f64> = divided.iter().map(|x| x*weight).collect();
        linear::vector_add_inplace(&mut penalty_term, &weighted); 
    }
    let result: Vec<f64> = penalty_term.iter().map(|x| x*penalty_scalar).collect();
    return result
}

/// The fission term of the fused lasso model. Responsible for the minimization of the centroids, from their respective data points
fn fission_gradient(centroid: &[f64], design_matrix: &Matrix, clusters: &[usize]) -> Vec<f64> {
    let design_mean = clusterpath_mean(design_matrix, clusters);
    let result = linear::vector_subtract(centroid, &design_mean);
    return result
}

/// special mean function, taking the sum of original data points, divided by cluster-set cardinality
fn clusterpath_mean(design_matrix: &Matrix, clusters: &[usize]) -> Vec<f64> {
    let cluster_cardinality: usize = clusters.len();
    let mut mean: Vec<f64> = Vec::with_capacity(cluster_cardinality);
    for index in 0..cluster_cardinality {
        linear::vector_add_inplace(&mut mean, design_matrix.get_col(index))
    }
    for element in mean.iter_mut() {
        *element = *element / (cluster_cardinality as f64);
    }
    return mean
}

// weights that decrease exponentially with distance, inducing greater densisty importance in the convex clustering solution
fn exponential_weight(datapoint_1: &[f64], datapoint_2: &[f64]) -> f64 {
    let gamma = 1.0;
    let diff: Vec<f64> = linear::vector_subtract(datapoint_1, datapoint_2);
    let squared_norm: f64 = linear::l2_norm(&diff).powf(2.0);
    return (-gamma * squared_norm).exp();
}

