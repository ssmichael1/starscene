use crate::Star;
use pyo3::prelude::*;

#[pyclass(name = "Star")]
#[derive(Clone)]
pub struct PyStar {
    inner: Star,
}

#[pymethods]
impl PyStar {
    #[new]
    pub fn new(j2000_vec: (f64, f64, f64), v_mag: f64) -> Self {
        let (x, y, z) = j2000_vec;
        Self {
            inner: Star {
                j2000_vec: [x, y, z],
                v_mag,
            },
        }
    }

    #[getter]
    pub fn j2000_vec(&self) -> (f64, f64, f64) {
        let [x, y, z] = self.inner.j2000_vec;
        (x, y, z)
    }

    #[getter]
    pub fn v_mag(&self) -> f64 {
        self.inner.v_mag
    }

    #[getter]
    pub fn ra_dec_rad(&self) -> (f64, f64) {
        let [x, y, z] = self.inner.j2000_vec;
        let ra_rad = y.atan2(x).rem_euclid(2.0 * std::f64::consts::PI);
        let dec_rad = z.asin();
        (ra_rad, dec_rad)
    }

    #[getter]
    pub fn ra_rad(&self) -> f64 {
        let [x, y, _z] = self.inner.j2000_vec;
        y.atan2(x).rem_euclid(2.0 * std::f64::consts::PI)
    }

    #[getter]
    pub fn dec_rad(&self) -> f64 {
        self.inner.j2000_vec[2].asin()
    }
}

impl From<Star> for PyStar {
    fn from(value: Star) -> Self {
        Self { inner: value }
    }
}

impl From<PyStar> for Star {
    fn from(value: PyStar) -> Self {
        value.inner
    }
}
