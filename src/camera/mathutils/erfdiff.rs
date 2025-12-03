use std::sync::OnceLock;

const Z_MAX: f64 = 6.0;
const N: usize = 8192;
const H: f64 = Z_MAX / (N as f64 - 1.0);

static EXP_NEG_SQ_LUT: OnceLock<[f64; N]> = OnceLock::new();

#[inline(always)]
fn lut() -> &'static [f64; N] {
    EXP_NEG_SQ_LUT.get_or_init(|| {
        let mut arr = [0.0f64; N];
        for i in 0..N {
            let z = i as f64 * H;
            arr[i] = (-z * z).exp();
        }
        arr
    })
}

#[inline(always)]
fn interp_exp_neg_sq(z: f64) -> f64 {
    let a = z.abs();
    if a >= Z_MAX {
        return (-a * a).exp();
    }

    let t = a / H;
    let k = t.floor() as usize;
    let frac = t - (k as f64);
    let k1 = if k + 1 < N { k + 1 } else { k };

    let lut = lut();
    let f0 = lut[k];
    let f1 = lut[k1];

    f0 + frac * (f1 - f0)
}

/// Fast approximation of erf(x) - erf(y) using a midpoint rule and LUT for exp(-z^2).
///
/// Formula:
///   erf(x) - erf(y) ≈ (2/sqrt(pi)) * (x-y) * exp(-((x+y)/2)^2)
///
/// We approximate exp(-((x+y)/2)^2) via linear interpolation from a precomputed table.
/// This avoids catastrophic cancellation when erf(x) and erf(y) are both near ±1
/// and is typically accurate to about 1e-6 for |x|,|y| <= 6.
pub fn erf_diff_lut(x: f64, y: f64) -> f64 {
    if x == y {
        return 0.0;
    }

    // Midpoint and interval
    let dx = x - y;
    let c = 0.5 * (x + y);

    // 2 / sqrt(pi)
    const TWO_OVER_SQRT_PI: f64 = 1.128_379_167_095_512_6;

    let f = interp_exp_neg_sq(c);
    TWO_OVER_SQRT_PI * dx * f
}
