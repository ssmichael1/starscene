//! The poission distribution module in the rand_dist crate
//! does not allow for easily changin the rate constant parameter
//! lambda, so we have to implement our own poisson distribution
//! here.
//!

use rand::Rng;

/// Log-gamma function
///
/// The algorithm used here is from the python numpy library:
///  The algorithm comes from SPECFUN by Shanjie Zhang and Jianming Jin
///  and their book
///  "Computation of Special Functions", 1996, John Wiley & Sons, Inc.
fn random_loggam(x: f64) -> f64 {
    let a = [
        8.333_333_333_333_333e-2,
        -2.777_777_777_777_778e-3,
        7.936_507_936_507_937e-4,
        -5.952_380_952_380_952e-4,
        8.417_508_417_508_418e-4,
        -1.917_526_917_526_918e-3,
        6.410_256_410_256_41e-3,
        -2.955_065_359_477_124e-2,
        1.796_443_723_688_307e-1,
        -1.39243221690590e00,
    ];
    let n: u64;

    if x == 1.0 || x == 2.0 {
        return 0.0;
    } else if x < 7.0 {
        n = 7 - x as u64;
    } else {
        n = 0;
    }
    let mut x0 = x + n as f64;
    let x2 = 1.0 / (x0 * x0);
    let lg2pi = (2.0 * std::f64::consts::PI).ln();
    let mut gl0 = a[9];
    for i in (0..9).rev() {
        gl0 = gl0 * x2 + a[i];
    }
    let mut gl = gl0 / x0 + 0.5 * lg2pi + (x0 - 0.5) * x0.ln() - x0;
    if x < 7.0 {
        for _ in 1..(n + 1) {
            gl -= (x0 - 1.0).ln();
            x0 -= 1.0;
        }
    }
    gl
}

/// Generate a random number from a poisson distribution
/// with a given rate constant lambda.
///
/// # Arguments
/// * `rng` - A mutable reference to a random number generator
/// * `lambda` - The rate constant of the poisson distribution
///
/// # Returns
/// * A random number from the poisson distribution
///
/// # Notes:
/// * The algorithm is copied from the python numpy library:
/// * The algorithm uses
///   "The transformed rejection method for generating Poisson random variables"
///    by W. Hoermann, Insurance: Mathematics and Economics 12, 39-45 (1993)
///    for lambda >= 10.0
/// * For lambda < 10.0, the algorithm uses the Knuth algorithm
///
///
/// # Example
/// ```
/// use rand::thread_rng;
/// use starscene::poisson::poissrand;
/// p = poissrand(&mut thread_rng(), 1.0);
/// ```
///
pub fn poissrand(rng: &mut impl Rng, lambda: f64) -> u64 {
    // Use the algorithm used by the python numpy library
    // to generate a poisson random number
    if lambda >= 10.0 {
        /*
         * The transformed rejection method for generating Poisson random variables
         * W. Hoermann
         * Insurance: Mathematics and Economics 12, 39-45 (1993)
         */
        let slam = lambda.sqrt();
        let loglam = lambda.ln();
        let b = 0.931 + 2.53 * slam;
        let a = -0.059 + 0.02483 * b;
        let invalpha = 1.1239 + 1.1328 / (b - 3.4);
        let vr = 0.9277 - 3.6224 / (b - 2.0);
        loop {
            let u = rng.gen::<f64>() - 0.5;
            let v = rng.gen::<f64>();
            let us = 0.5 - u.abs();
            let k = ((2.0 * a / us + b) * u + lambda + 0.43).floor() as i64;
            if (us > 0.07) && (v <= vr) {
                return k as u64;
            }
            if k < 0 || (us < 0.013 && v > us) {
                continue;
            }
            if (v.ln() + invalpha.ln() - (a / (us * us) + b).ln())
                <= (-lambda + k as f64 * loglam - random_loggam(k as f64 + 1.0))
            {
                return k as u64;
            }
        }
    } else if lambda > 0.0 {
        let enlam = (-lambda).exp();
        let mut x = 0;
        let mut prod = 1.0;
        loop {
            prod *= rng.gen::<f64>();
            if prod > enlam {
                x += 1;
            } else {
                return x;
            }
        }
    } else {
        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::thread_rng;

    #[test]
    fn test_poisson() {
        let mut rng = thread_rng();
        let lambda = 20.0;
        let mut sum = 0;
        let mut sumsq = 0;
        let nsample = 10000;
        for _ in 0..nsample {
            let p = poissrand(&mut rng, lambda);
            sum += p;
            sumsq += p * p;
        }
        let mean = sum as f64 / nsample as f64;
        let var = sumsq as f64 / nsample as f64 - mean * mean;
        println!("mean = {}", mean);
        println!("var = {}", var);
        assert!((mean - lambda).abs() < 0.2);
        assert!((var - lambda).abs() < 0.5);
    }
}
