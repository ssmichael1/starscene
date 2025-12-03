use rand::Rng;
use std::sync::OnceLock;

/// Precomputed ln(k!) table up to KMAX, Stirling beyond.
pub struct LogFactTable {
    logf: Vec<f64>,
    kmax: usize,
}

impl LogFactTable {
    pub fn new(kmax: usize) -> Self {
        let mut logf = vec![0.0; kmax + 1];
        let mut acc = 0.0;
        for k in 1..=kmax {
            acc += (k as f64).ln();
            logf[k] = acc;
        }
        Self { logf, kmax }
    }

    #[inline]
    pub fn log_factorial(&self, k: u32) -> f64 {
        let ku = k as usize;
        if ku <= self.kmax {
            unsafe { *self.logf.get_unchecked(ku) }
        } else {
            stirling_log_factorial(k as f64)
        }
    }
}

#[inline]
fn stirling_log_factorial(x: f64) -> f64 {
    // Stirling series to O(1/x^3)
    x * x.ln() - x + 0.5 * (2.0 * std::f64::consts::PI * x).ln() + 1.0 / (12.0 * x)
        - 1.0 / (360.0 * x * x * x)
}

// ---- Global OnceLock table ----
static LOGFACT: OnceLock<LogFactTable> = OnceLock::new();

#[inline]
fn logfact() -> &'static LogFactTable {
    // tune KMAX as you like; 65536 ~ 32 KB of f64s
    LOGFACT.get_or_init(|| LogFactTable::new(65_536))
}

/// 3-regime Poisson sampler: exact for small/mid, normal approx for huge.
#[inline]
pub fn poisson_u32<R: Rng + ?Sized>(rng: &mut R, lambda: f64) -> u32 {
    if lambda <= 0.0 {
        return 0;
    }

    if lambda < 12.0 {
        knuth(rng, lambda)
    } else if lambda < 1.0e5 {
        ptrs(rng, lambda, logfact())
    } else {
        normal_approx(rng, lambda)
    }
}

// --- Small λ: Knuth product method (exact) ---
#[inline]
fn knuth<R: Rng + ?Sized>(rng: &mut R, lambda: f64) -> u32 {
    let l = (-lambda).exp();
    let mut k = 0_u32;
    let mut p = 1.0_f64;
    loop {
        k += 1;
        p *= rng.random::<f64>();
        if p <= l {
            return k - 1;
        }
    }
}

// --- Mid λ: Hörmann PTRS transformed rejection (exact) ---
#[inline]
fn ptrs<R: Rng + ?Sized>(rng: &mut R, lambda: f64, lf: &LogFactTable) -> u32 {
    let slam = lambda.sqrt();
    let loglam = lambda.ln();
    let b = 0.931 + 2.53 * slam;
    let a = -0.059 + 0.02483 * b;
    let inv_alpha = 1.1239 + 1.1328 / (b - 3.4);
    let v_r = 0.9277 - 3.6224 / (b - 2.0);

    loop {
        let u = rng.random::<f64>() - 0.5;
        let v = rng.random::<f64>();
        let us = 0.5 - u.abs();
        let k = ((2.0 * a / us + b) * u + lambda + 0.43).floor();
        if k < 0.0 {
            continue;
        }
        let ku = k as u32;

        // quick accept/reject
        if us >= 0.07 && v <= v_r {
            return ku;
        }
        if us < 0.013 && v > us {
            continue;
        }

        // log acceptance test
        let log_v = v.ln();
        let lhs = log_v + inv_alpha.ln() - (a / (us * us + a)).ln();
        let rhs = -lambda + k * loglam - lf.log_factorial(ku);

        if lhs <= rhs {
            return ku;
        }
    }
}

// --- Huge λ: Normal approximation (very fast, accurate when λ large) ---
#[inline]
fn normal_approx<R: Rng + ?Sized>(rng: &mut R, lambda: f64) -> u32 {
    // Box–Muller for standard normal Z
    let u1 = rng.random::<f64>().max(1e-300);
    let u2 = rng.random::<f64>();
    let z = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();

    // Continuity correction: Poisson ≈ N(λ, λ)
    let k = (lambda + lambda.sqrt() * z + 0.5).floor();
    if k <= 0.0 {
        0
    } else {
        k as u32
    }
}

// Sample a big lambda field in parallel.
// pub fn poisson_field(lambdas: &[f64], seed: u64) -> Vec<u32> {
//     // ensure LOGFACT initialized on main thread (optional; avoids one-time init on worker)
//     let _ = logfact();

//     let mut out = vec![0u32; lambdas.len()];
//     out.par_chunks_mut(4096)
//         .zip(lambdas.par_chunks(4096))
//         .enumerate()
//         .for_each(|(cid, (o, l))| {
//             let mut rng = SmallRng::seed_from_u64(
//                 (seed ^ (cid as u64).wrapping_mul(0x9E3779B97F4A7C15)).into()
//             );
//             for (oi, &lam) in o.iter_mut().zip(l.iter()) {
//                 *oi = poisson_u32(&mut rng, lam);
//             }
//         });

//     out
// }

// fn main() {
//     let n = 4_000_000;

//     // Example λ field spanning 0 .. >1e6
//     let lambdas: Vec<f64> = (0..n)
//         .map(|i| {
//             let t = i as f64 / n as f64;
//             if t < 0.3 { 5.0 * (20.0 * t).sin().abs() }
//             else if t < 0.8 { 500.0 * (10.0 * t).cos().abs() + 20.0 }
//             else { 1.0e6 * (t - 0.8).powi(2) + 1.0e5 }
//         })
//         .collect();

//     let samples = poisson_field(&lambdas, 12345);
//     println!("first 10 samples: {:?}", &samples[..10]);
// }
