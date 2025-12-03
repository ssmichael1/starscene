/// Approximation of the error function using the Hastings approximation.
/// Accurate to about 7 decimal places.
///
/// # Arguments
/// * `x` - The input value
///
/// # Returns
/// The approximate value of erf(x)
///
///
#[inline]
pub fn erf(x: f64) -> f64 {
    // Approximation with max abs error ~1.2e-7
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let xx = x.abs();
    let xx_sq = xx * xx;

    if xx <= 0.5 {
        return small_range_erf(x);
    }

    let t = 1.0 / (1.0 + 0.5 * xx);

    // Polynomial in t, evaluated with Horner's rule
    let exponent = if cfg!(target_feature = "fma") {
        hastings_poly_fma(t, xx_sq)
    } else {
        hastings_poly_scalar(t, xx_sq)
    };

    let tau = t * exponent.exp();

    sign * (1.0 - tau)
}

#[inline(always)]
fn hastings_poly_scalar(t: f64, xx_sq: f64) -> f64 {
    // Classic Horner evaluation (portable path)
    let mut poly = 0.17087277f64;
    poly = poly * t + -0.82215223;
    poly = poly * t + 1.48851587;
    poly = poly * t + -1.13520398;
    poly = poly * t + 0.27886807;
    poly = poly * t + -0.18628806;
    poly = poly * t + 0.09678418;
    poly = poly * t + 0.37409196;
    poly = poly * t + 1.00002368;
    poly * t + (-xx_sq - 1.26551223)
}

#[inline(always)]
fn hastings_poly_fma(t: f64, xx_sq: f64) -> f64 {
    // Uses fused multiply-adds when the target enables FMA
    let mut poly = 0.17087277f64;
    poly = poly.mul_add(t, -0.82215223);
    poly = poly.mul_add(t, 1.48851587);
    poly = poly.mul_add(t, -1.13520398);
    poly = poly.mul_add(t, 0.27886807);
    poly = poly.mul_add(t, -0.18628806);
    poly = poly.mul_add(t, 0.09678418);
    poly = poly.mul_add(t, 0.37409196);
    poly = poly.mul_add(t, 1.00002368);
    poly.mul_add(t, -xx_sq - 1.26551223)
}

#[inline]
fn small_range_erf(x: f64) -> f64 {
    if cfg!(target_feature = "fma") {
        small_range_erf_fma(x)
    } else {
        small_range_erf_scalar(x)
    }
}

#[inline(always)]
fn small_range_erf_scalar(x: f64) -> f64 {
    // Maclaurin series truncated at x^11, scaled by 2/sqrt(pi)
    let x2 = x * x;
    let mut poly = -0.0008556239969770373f64; // coefficient for x^10 term before multiplying by x
    poly = poly * x2 + 0.005223977625442187;
    poly = poly * x2 + -0.026866170645131252;
    poly = poly * x2 + 0.11283791670955126;
    poly = poly * x2 + -0.37612638903183754;
    poly = poly * x2 + 1.1283791670955126;
    x * poly
}

#[inline(always)]
fn small_range_erf_fma(x: f64) -> f64 {
    let x2 = x * x;
    let mut poly = -0.0008556239969770373f64;
    poly = poly.mul_add(x2, 0.005223977625442187);
    poly = poly.mul_add(x2, -0.026866170645131252);
    poly = poly.mul_add(x2, 0.11283791670955126);
    poly = poly.mul_add(x2, -0.37612638903183754);
    poly = poly.mul_add(x2, 1.1283791670955126);
    x * poly
}

/// Inverse error function erf^{-1}(x) for x in (-1, 1).
/// Uses Winitzki initial guess + Newton-Raphson refinement.
///
/// Panics if |x| >= 1.0 (since erf^{-1} diverges).
pub fn erf_inv(x: f64) -> f64 {
    assert!(x > -1.0 && x < 1.0, "erf_inv domain is (-1, 1)");

    if x == 0.0 {
        return 0.0;
    }

    // ---- Winitzki approximation for initial guess ----
    // erf^{-1}(x) ≈ sign(x) * sqrt( sqrt(t^2 - ln/a) - t )
    // with a ≈ 0.147
    let a = 0.147;
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let xx = x * x;

    let ln = (1.0 - xx).ln(); // ln(1 - x^2), negative
    let t = 2.0 / (std::f64::consts::PI * a) + ln / 2.0;
    let inside = (t * t - ln / a).sqrt() - t;
    let mut y = sign * inside.sqrt();

    // ---- Newton refinement ----
    // y_{n+1} = y_n - (erf(y_n) - x) / (2/sqrt(pi) * exp(-y_n^2))
    let two_over_sqrt_pi = 2.0 / std::f64::consts::PI.sqrt();

    // usually 2 iterations is enough; 3 is very safe
    for _ in 0..3 {
        let err = erf(y) - x;
        let deriv = two_over_sqrt_pi * (-y * y).exp();
        y -= err / deriv;
    }

    y
}
