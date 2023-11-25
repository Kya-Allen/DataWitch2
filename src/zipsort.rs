use std::vec::Vec;

/// Datastructure holding an observation of 2 variables
struct Observation {
    key_variable: f64,
    complement_variables: f64,
}

/// Datastructure holding an observation of more than 2 variables
struct MatObservation {
    key_variable: f64,
    complement_variables: Vec<f64>,
}



pub fn vecsort(key: &[f64], complement: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let length: usize = key.len();
    assert!(key.len() == complement.len(), "ERROR: key and complement must have the same length");
    let paired_data = key.iter().zip(complement.iter());
    let mut sorted: Vec<Observation> = Vec::with_capacity(length);
    for (index, (x, y)) in paired_data.enumerate() {
        sorted.push(Observation {key_variable: *x, complement_variables: *y})
    }
    sorted.sort_by(|a, b| a.key_variable.partial_cmp(&b.key_variable).unwrap());
    let mut new_key: Vec<f64> = Vec::with_capacity(length);
    let mut new_complement: Vec<f64> = Vec::with_capacity(length);
    for obs in sorted {
        new_key.push(obs.key_variable);
        new_complement.push(obs.complement_variables)
    }
    return (new_key, new_complement)
}

pub fn matsort(key: &[f64], complement: Vec<&[f64]>) ->(Vec<f64>, Vec<Vec<f64>>) {
    let length: usize = key.len();
    assert!(key.len() == complement[0].len(), "ERROR: key and complement must have the same length");
    let mut sorted: Vec<MatObservation> = Vec::with_capacity(length);
    for index in 0..length {
        let current_key: f64 = key[index];
        let mut current_complements: Vec<f64> = Vec::with_capacity(complement.len());
        for jindex in 0..complement.len() {
            current_complements.push(complement[jindex][index]);
        }
        sorted.push(MatObservation {key_variable: current_key, complement_variables: current_complements})

    }
    sorted.sort_by(|a, b| a.key_variable.partial_cmp(&b.key_variable).unwrap());
    let mut new_key: Vec<f64> = Vec::with_capacity(length);
    let mut new_complement: Vec<Vec<f64>> = Vec::with_capacity(complement.len());
    
    new_key = sorted.iter().map(|x| x.key_variable.clone()).collect();
    for index in 0..complement.len() {
        let mut current_complement: Vec<f64> = Vec::with_capacity(length);
        for jindex in 0..length {
            current_complement.push(sorted[jindex].complement_variables[index].clone());
        };
        new_complement.push(current_complement)
    }
    return (new_key, new_complement)
}