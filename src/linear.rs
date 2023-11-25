use std::vec;

#[derive(Clone)]
pub struct Matrix {
    pub num_columns: usize,
    pub num_rows: usize,
    pub data: Vec<f64>,
}

pub fn make_matrix(data: Vec<f64>, num_columns: usize) -> Matrix {
    return Matrix {
        num_columns,
        num_rows: data.len() / num_columns,
        data,
    }
}

pub fn matadd(matrix1: Matrix, matrix2: Matrix) -> Matrix {
    let checker = ((matrix1.num_columns == matrix2.num_columns) && (matrix1.num_rows == matrix2.num_rows));
    assert!(checker, "ERROR: matrices must be of equal shape");

    let mut new_vec: Vec<f64> = vec![0.0; matrix1.data.len()];
    new_vec = new_vec.into_iter().enumerate().map(|(i, mut x)| {x = (matrix1.data[i] + matrix2.data[i]); return x}).collect();
    //for i in 0..matrix1.data.len() {
    //    new_vec[i] = matrix1.data[i] + matrix2.data[i];
    //}
    return make_matrix(new_vec, matrix1.num_columns)

}

/// Elementwise Subtraction of matrices
pub fn matsub(matrix1: Matrix, matrix2: Matrix) -> Matrix {
    let checker = ((matrix1.num_columns == matrix2.num_columns) && (matrix1.num_rows == matrix2.num_rows));
    assert!(checker, "ERROR: matrices must be of equal shape");
    
    let mut new_vec: Vec<f64> = vec![0.0; matrix1.data.len()];
    new_vec = new_vec.into_iter().enumerate().map(|(i, mut x)| {x = (matrix1.data[i] - matrix2.data[i]); return x}).collect();
    return make_matrix(new_vec, matrix1.num_columns)

}

pub fn get_transpose(matrix: &Matrix) -> Matrix {
    let mut new_flat_matrix: Vec<f64> = Vec::with_capacity(matrix.data.len());
    for i in 0..matrix.col_size() {
        for j in 0..matrix.row_size() {
            new_flat_matrix.push(matrix.index(j, i));
        };
    };
    //return Matrix {data: new_flat_matrix, num_columns: matrix.num_rows, num_rows: matrix.num_columns}
    return Matrix::new(new_flat_matrix, matrix.num_rows)
}

pub fn elementwise_matmul(matrix1: &Matrix, matrix2: &Matrix) -> Matrix {
    let checker = ((matrix1.num_columns == matrix2.num_columns) && (matrix1.num_rows == matrix2.num_rows));
    assert!(checker, "ERROR: matrices must be of equal shape");

    let mut new_vec: Vec<f64> = vec![0.0; matrix1.data.len()];
    new_vec = new_vec.into_iter().enumerate().map(|(i, mut x)| {x = matrix1.data[i] * matrix2.data[i]; return x}).collect();
    return Matrix::new(new_vec, matrix1.num_columns)
}

pub fn vector_dot_product(vector_1: &[f64], vector_2: &[f64]) -> f64 {
    let mut result: f64 = 0.0;
    for index in 0..vector_1.len() {
        let elementwise_product = vector_1[index] * vector_2[index];
        result = result + elementwise_product;
    };
    return result
}

/// Proper Matrix Multiplication
pub fn naive_matmul(matrix_1: &Matrix, matrix_2: &Matrix) -> Matrix {
    assert!(matrix_1.num_columns == matrix_2.num_rows, "ERROR: The num_columns of the first matrix must equal num_rows of the second matrix");
    let mut result_flat_mat: Vec<f64> = Vec::with_capacity(matrix_1.data.len());
    let new_matrix_1: Matrix = get_transpose(matrix_1);
    let matrix_2_col_size: usize = matrix_2.col_size();
    let new_mat_1_col_size: usize = new_matrix_1.col_size();
    for (right_column_index, right_vector) in matrix_2.data.chunks_exact(matrix_2_col_size).enumerate() {
        for (left_column_index, left_vector) in new_matrix_1.data.chunks_exact(new_mat_1_col_size).enumerate() {
            result_flat_mat.push(vector_dot_product(left_vector, right_vector));
        };
    };
    return Matrix {data: result_flat_mat, num_columns: matrix_2.num_columns, num_rows: matrix_1.num_rows}
}

pub fn matrix_vector_product(matrix_1: &Matrix, vector_2: &[f64]) -> Vec<f64> {
    assert!(vector_2.len() == matrix_1.col_size(), "Error: The dimension of the vector must equal the dimension of the matrix columns");
    let new_matrix: Matrix = get_transpose(matrix_1);
    let mut result_flat_mat: Vec<f64> = Vec::with_capacity(matrix_1.col_size());
    for (left_column_index, left_vector) in new_matrix.data.chunks_exact(new_matrix.col_size()).enumerate() {
        result_flat_mat.push(vector_dot_product(left_vector, vector_2));
    };
    return result_flat_mat
}

pub fn vector_norm(vector_1: &[f64], order: isize) -> f64 {
    match order {
        1 => return l1_norm(vector_1),
        2 => return l2_norm(vector_1),
        -1 => return linf_norm(vector_1),
        _ => panic!("The order of the norm is out of range. only l1, l2, and l-inf are supported")
    };
    return 0.0
}

pub fn linf_norm(vector_1: &[f64]) -> f64 {
    let mut value: f64 = vector_1[0].clone();
    for (index, element) in vector_1.iter().enumerate() {
        value = value.max(*element);
    }
    return value
}

pub fn minimum(vector_1: &[f64]) -> f64 {
    let mut value: f64 = vector_1[0].clone();
    for (index, element) in vector_1.iter().enumerate() {
        value = value.min(*element);
    }
    return value
}

pub fn l1_norm(vector_1: &[f64]) -> f64 {
    let mut value: f64 = 0.0;
    for (index, element) in vector_1.iter().enumerate() {
        value = value + element.abs()
    };
    return value
}

pub fn l2_norm(vector_1: &[f64]) -> f64 {
    let mut value: f64 = 0.0;
    for (index, element) in vector_1.iter().enumerate() {
        value = value + (element.abs()).powf(2.0)
    };
    return value.sqrt()
}

/// Iterative method converges on the vector with the ||largest|| eigenvalue of a matrix 
pub fn power_iteration(matrix: &Matrix) -> Vec<f64> {
    let max_iters: usize = 20;
    let tolerance: f64 = 0.05;
    let mut current_vector: Vec<f64> = vec![1.0, 1.0];
    let mut previous_vector: Vec<f64> = current_vector.clone();
    for iteration in 0..max_iters {
        // matrix-by-vector dot product
        current_vector = matrix_vector_product(matrix, &current_vector[0..current_vector.len()]);
        // calculate the norm of above result
        let current_norm = vector_norm(&current_vector[0..current_vector.len()], -1);
        // normalize the vector by taking the product / norm
        current_vector = current_vector.into_iter().map(|x| x / current_norm).collect();

        if euclidean_distance(&current_vector[..], &previous_vector[..]) <= tolerance {
            return current_vector
        }
        previous_vector = current_vector.clone();
    }
    return current_vector
}

/// gets the associated eigenvalue for an eigenvector
pub fn rayleigh_quotient(matrix: &Matrix, eigenvector: &[f64]) -> f64 {
    let product: Vec<f64> = matrix_vector_product(matrix, eigenvector);
    let top_dot: f64 = vector_dot_product(eigenvector, &product[..]);
    let bottom_dot: f64 = vector_dot_product(eigenvector, eigenvector);
    return top_dot / bottom_dot; 
}

pub fn euclidean_distance(vector_1: &[f64], vector_2: &[f64]) -> f64 {
    let difference_vector: Vec<f64> = vector_subtract(vector_1, vector_2);
    return l2_norm(&difference_vector[..]);  
}

pub fn vector_subtract(vector_1: &[f64], vector_2: &[f64]) -> Vec<f64> {
    let mut new_vec: Vec<f64> = Vec::with_capacity(vector_1.len());
    for (index, element) in vector_1.into_iter().enumerate() {
        new_vec.push(element - vector_2[index]);
    }
    return new_vec
}

pub fn vector_add_inplace(vector_1: &mut Vec<f64>, vector_2: &[f64]) {
    for (index, element) in vector_1.into_iter().enumerate() {
        *element = *element + vector_2[index];
    }
}

impl Matrix {
    /// Constructor
    pub fn new(input_data: Vec<f64>, num_columns: usize) -> Self {
        if input_data.len() > 0 {
            Self {
                num_columns,
                num_rows: input_data.len() / num_columns,
                data: input_data,
            }
        }
        else {
            Self {
                num_columns,
                num_rows: input_data.capacity() / num_columns,
                data: input_data,
            }
        }
    }
    /// Sum of all Matrix Elements
    pub fn sum(&self) -> f64 {
        return self.data.iter().map(|x| *x).sum::<f64>()
    }

    /// length of Columns
    pub fn col_size(&self) -> usize {
        return self.data.len() / self.num_columns
    }

    /// length of Rows
    pub fn row_size(&self) -> usize {
        return self.data.len() / self.num_rows
    }
    
    /// query Matrix by index (column, row)
    pub fn index(&self, i: usize, j: usize) -> f64 {
        assert!(i <= self.num_columns, "ERROR: index i out of range, columns");
        assert!(j <= self.num_rows, "ERROR: index j out of range, rows");
        let mut pointer: usize = 0;
        let col_size = self.col_size();

        pointer = pointer + (col_size * i);
        pointer = pointer + j;
        
        return self.data[pointer]
    }
    
    pub fn get_col(&self, column: usize) -> &[f64] {
        assert!(column <= self.num_columns, "ERROR: column index out of range");
        let mut pointer: usize = 0;
        let col_size: usize = self.col_size();
        pointer = pointer + (col_size * column);

        return &self.data[pointer..(pointer + col_size)]
    }

    pub fn display(&self) -> Vec<&[f64]> {
        let mut outer_vec: Vec<&[f64]> = Vec::with_capacity(self.num_columns);
        let mut pointer: usize = 0;
        let col_size = self.col_size();

        for column in 0..self.num_columns {
           outer_vec.push(&self.data[pointer..(pointer + col_size)]); 
           pointer = pointer + col_size
        }

        return outer_vec
    }

    pub fn frobenius(&self) -> f64 {
        let mut collection: Vec<f64> = self.data.iter().map(|x| x.abs()).collect();
        collection = collection.into_iter().map(|mut x| x.powf(2.0)).collect();
        let sum: f64 = collection.iter().sum();
        return sum.sqrt()
    }
}