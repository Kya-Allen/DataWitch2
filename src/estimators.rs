

///Gaussian Maximum Likelihood Estimator, f64
pub fn mean(data: &[f64]) -> f64 {
    let sum: f64 = data.iter().sum();
    let mean: f64 = sum / data.len() as f64;
    return mean
}

pub fn mean_difference(sample_1: &[f64], sample_2: &[f64]) -> f64 {
    let mut difference: f64 = mean(sample_1) - mean(sample_2);
    difference = difference.abs();
    return difference
}

/// Mean Squared Error
pub fn mse(data1: &[f64], data2: &[f64]) -> f64 {
   let mut squared_errors: Vec<f64> = Vec::with_capacity(data1.len());
   let paired_data = data1.iter().zip(data2.iter());

   for (index, (x, y)) in paired_data.enumerate() {
        let error: f64 = x - y;
        let squared_error: f64 = error.powf(2.0);
        squared_errors.push(squared_error);
   }

   let mean_squared_error = mean(&squared_errors[..]);
   return mean_squared_error;
}

/// Root Mean Squared Error
pub fn rmse(data1: &[f64], data2: &[f64]) -> f64 {
    let mean_squared_error: f64 = mse(data1, data2);
    return mean_squared_error.sqrt();
}