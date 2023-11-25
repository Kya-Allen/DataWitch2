use std::vec;

pub fn uniform(n: usize, min: f64, max: f64, seed: f64) -> Vec<f64> {
    let modulus: f64 = (2.0_f64).powf(48.0);
    let uniform_variates: Vec<f64> = linear_congruential_generator(n, seed);
    let standard_uniform_variates: Vec<f64> = uniform_variates.into_iter().map(|x| x/modulus).collect();
    if (min == 0.0) && (max == 1.0) {
        return standard_uniform_variates
    }
    let scaled_variates: Vec<f64> = standard_uniform_variates.into_iter().map(|x| min + (x * (max - min))).collect();
    return scaled_variates
}

pub fn exprand(n: usize, lambda: f64, seed: f64) -> Vec<f64> {
    let pre_transform: Vec<f64> = uniform(n, 0.0, 1.0, seed);
    let result: Vec<f64> = pre_transform.into_iter().map(|u| exponential_dist_inverse_cdf(&u, &lambda)).collect();
    return result
}

pub fn linear_congruential_generator(n: usize, seed: f64) -> Vec<f64> {
    let mut result: Vec<f64> = Vec::with_capacity(n);
    //for i in 0..n {result[i] = linear_congruential_equation(result[i-1])}
    let mut previous_value = seed;
    for index in 0..n {
        previous_value = linear_congruential_equation(&previous_value);
        //previous_value = scaled_value;
        result.push(previous_value)
    }
    return result;
}

pub fn linear_congruential_equation(previous_value: &f64) -> f64 {
    // parameters borrowed from Newlib
    let modulus: f64 = (2.0_f64).powf(48.0);
    let multiplier: f64 = 25214903917.0; 
    let increment: f64 = 11.0;
    let result = ((multiplier * previous_value) + increment).rem_euclid(modulus);
    return result
}

pub fn exponential_dist_inverse_cdf(value: &f64, lambda: &f64) -> f64{
    let left = -(1.0/lambda); let right = (1.0-value).ln();
    return left * right;
}

//fn Linear_feedback_shift_register() {
//    
//}