#![allow(warnings)]
mod linear; 
mod rng;
//mod taxometrics;
mod L2_clusterpath;
mod zipsort;
mod estimators;
mod convex_clustering;
mod approximate_nearest_neighbors;
use std::vec;

fn main() {
    let my_matrix = linear::Matrix::new(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], 3);
    let my_matrix2 = linear::Matrix{num_columns: 3, num_rows: 3, data: vec![2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]};
    let my_matrix_3 = linear::Matrix::new(vec![6.0, 4.0, 5.0, 5.0], 2);
    let my_matrix_3 = linear::Matrix::new(vec![0.0, 1.0, 1.0, 1.0], 2);
    //let my_matrix3 = fake::naive_matmul(&my_matrix, &my_matrix2);

    let x1 = vec![5.0, 2.0, 12.0];
    let x2 = vec![4.0, 5.0, 6.0];
    let x3 = vec![6.0, 2.0, 9.0];
    let m1: Vec<&[f64]> = vec![&x2[..], &x3[..]];
    //let my_dot = fake::vector_dot_product(&x1[0..x1.len()], &x2[0..x2.len()]);

    //let vec_mat_product = fake::matrix_vector_product(&my_matrix2, &x2[0..x2.len()]);
    let l1_norm = linear::vector_norm(&x2[0..x2.len()], 1);
    let l2_norm = linear::vector_norm(&x2[0..x2.len()], 2);
    //let powit = fake::power_iteration(&my_matrix_3);
    //let powit_value = fake::rayleigh_quotient(&my_matrix_3, &powit[..]);
    //println!("power iteration eigenvector {:?}", powit);
    //println!("power iteration eigenvalue {:?}", powit_value);
    //let ok = fake::linf_norm(&x2[0..x2.len()]);
    //println!("inf norm {:?}", ok);
    //let boop = rng::uniform(5, 0.0, 10.0, 3.0537);
    //let boop2 = rng::exprand(5, 0.2, 3.0537);
    //let (beep, boop) = zipsort::matsort(&x1[..], m1);
    //println!("result 1 {:?}", beep);
    //println!("result 2 {:?}", boop);
    let dataset: linear::Matrix = linear::Matrix::new(vec![30.0, 40.0, 20.0, 25.0, 30.0, 60.0, 2.0, 2.0, 4.0, 1.0, 1.0, 1.0], 6);
    let path_length: usize = 6;
    //let penalty_result = convex_clustering::max_penalty(&linear::get_transpose(&dataset));
    //println!("max penality {:}", penalty_result);
    let result: Vec<linear::Matrix> = convex_clustering::c_paint(dataset, path_length);
    for (i, matrix) in result.iter().enumerate() {
        println!("{:}th Solution", i);
        println!("{:?}", matrix.display());
    }
    let sorted_dim = vec![1.0, 2.0, 4.0, 4.0, 6.0, 6.0];
    let new_pen = 0.05;
    let en: usize = 6;
    let more_results = convex_clustering::partial_clusterpath(&sorted_dim, &new_pen, &en);
    println!("partial test {:?}", more_results);
    let fun_data: linear::Matrix = linear::Matrix::new(vec![5.0, 3.0, -2.0, -3.0], 4);
    let even_more = convex_clustering::c_paint(fun_data, 6);
    
}
