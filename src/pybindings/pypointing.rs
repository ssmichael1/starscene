use nalgebra::Quaternion as RawQuaternion;
use numpy::{
    ndarray::{arr1, arr2},
    PyArray1, PyArray2,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::camera::pointing::FixedInertialPointing;
use crate::camera::pointing::PointingProvider;

#[pyclass(name = "FixedInertialPointing")]
pub struct PyFixedInertialPointing {
    pointing: FixedInertialPointing,
}

impl PyFixedInertialPointing {
    pub(crate) fn as_provider(&self) -> &FixedInertialPointing {
        &self.pointing
    }
}

#[pymethods]
impl PyFixedInertialPointing {
    #[new]
    pub fn new(ra_rad: f64, dec_rad: f64, roll_rad: f64) -> Self {
        Self {
            pointing: FixedInertialPointing::from_ra_dec_roll_rad(ra_rad, dec_rad, roll_rad),
        }
    }

    #[staticmethod]
    pub fn from_ra_dec_roll_rad(ra_rad: f64, dec_rad: f64, roll_rad: f64) -> Self {
        Self {
            pointing: FixedInertialPointing::from_ra_dec_roll_rad(ra_rad, dec_rad, roll_rad),
        }
    }

    #[staticmethod]
    pub fn from_ra_dec_roll_deg(ra_deg: f64, dec_deg: f64, roll_deg: f64) -> Self {
        let ra_rad = ra_deg.to_radians();
        let dec_rad = dec_deg.to_radians();
        let roll_rad = roll_deg.to_radians();
        Self {
            pointing: FixedInertialPointing::from_ra_dec_roll_rad(ra_rad, dec_rad, roll_rad),
        }
    }

    #[staticmethod]
    pub fn from_quaternion(w: f64, x: f64, y: f64, z: f64) -> PyResult<Self> {
        let raw = RawQuaternion::new(w, x, y, z);
        let unit = nalgebra::Unit::try_new(raw, 1e-12)
            .ok_or_else(|| PyErr::new::<PyValueError, _>("Quaternion has near-zero norm"))?;
        Ok(Self {
            pointing: FixedInertialPointing::new(unit),
        })
    }

    pub fn right_ascension_rad(&self) -> f64 {
        self.pointing.right_ascension_rad()
    }

    pub fn declination_rad(&self) -> f64 {
        self.pointing.declination_rad()
    }

    pub fn roll_rad(&self) -> f64 {
        self.pointing.roll_rad()
    }

    pub fn boresight_j2000<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        let v = self.pointing.boresight_j2000();
        let arr = arr1(&[v.x, v.y, v.z]);
        PyArray1::from_owned_array(py, arr)
    }

    pub fn camera_from_j2000<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let rot = self.pointing.camera_from_j2000().to_rotation_matrix();
        let m = rot.matrix();
        let data = arr2(&[
            [m[(0, 0)], m[(0, 1)], m[(0, 2)]],
            [m[(1, 0)], m[(1, 1)], m[(1, 2)]],
            [m[(2, 0)], m[(2, 1)], m[(2, 2)]],
        ]);
        PyArray2::from_owned_array(py, data)
    }
}
