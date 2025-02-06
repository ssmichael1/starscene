use crate::StarScene;

use pyo3::prelude::*;
use pyo3::types::PyDateTime;
use pyo3::types::PyDict;

#[pyclass(name = "StarScene")]
#[derive(Debug)]
pub struct PyStarScene(StarScene);

#[pymethods]
impl PyStarScene {
    #[new]
    #[pyo3(signature=(**kwargs))]
    fn new(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let mut ss = StarScene::default();
        if let Some(kw) = kwargs {
            if let Some(v) = kw.get_item("pixel_pitch")? {
                ss.pixel_pitch = v.extract()?;
            }
            if let Some(v) = kw.get_item("focal_length")? {
                ss.focal_length = v.extract()?;
            }
            if let Some(v) = kw.get_item("pedestal")? {
                ss.pedestal = v.extract()?;
            }
            if let Some(v) = kw.get_item("bit_depth")? {
                ss.bit_depth = v.extract()?;
            }
            if let Some(v) = kw.get_item("gain")? {
                ss.gain = v.extract()?;
            }
            if let Some(v) = kw.get_item("exposure_time")? {
                ss.exposure_time = v.extract()?;
            }
            if let Some(v) = kw.get_item("dark_current")? {
                ss.dark_current = v.extract()?;
            }
            if let Some(v) = kw.get_item("read_noise")? {
                ss.read_noise = v.extract()?;
            }
            if let Some(v) = kw.get_item("quantum_efficiency")? {
                ss.solar_weighted_qe = v.extract()?;
            }
            if let Some(v) = kw.get_item("diameter")? {
                ss.diameter = v.extract()?;
            }
            if let Some(v) = kw.get_item("obscuration")? {
                ss.obscuration = v.extract()?;
            }
            if let Some(v) = kw.get_item("zodiacal")? {
                ss.zodiacal_light = v.extract()?;
            }
            if let Some(v) = kw.get_item("psf_1sigma")? {
                ss.psf_1sigma = v.extract()?;
            }
            if let Some(v) = kw.get_item("time")? {
                let tm: Py<PyDateTime> = v.extract().unwrap();
                pyo3::Python::with_gil(|py| {
                    let ts: f64 = tm
                        .call_method(py, "timestamp", (), None)
                        .unwrap()
                        .extract::<f64>(py)
                        .unwrap();
                    ss.time = satkit::Instant::from_unixtime(ts);
                });
            }
        }
        Ok(PyStarScene(ss))
    }

    #[getter]
    fn pixel_pitch(&self) -> f64 {
        self.0.pixel_pitch
    }

    #[getter]
    fn focal_length(&self) -> f64 {
        self.0.focal_length
    }

    #[getter]
    fn pedestal(&self) -> f64 {
        self.0.pedestal
    }

    #[getter]
    fn bit_depth(&self) -> u8 {
        self.0.bit_depth
    }

    #[getter]
    fn gain(&self) -> f64 {
        self.0.gain
    }

    #[getter]
    fn exposure_time(&self) -> f64 {
        self.0.exposure_time
    }

    #[getter]
    fn dark_current(&self) -> f64 {
        self.0.dark_current
    }

    #[getter]
    fn read_noise(&self) -> f64 {
        self.0.read_noise
    }

    #[getter]
    fn quantum_efficiency(&self) -> f64 {
        self.0.solar_weighted_qe
    }

    #[getter]
    fn diameter(&self) -> f64 {
        self.0.diameter
    }

    #[getter]
    fn obscuration(&self) -> f64 {
        self.0.obscuration
    }

    #[getter]
    fn zodiacal(&self) -> f64 {
        self.0.zodiacal_light
    }

    #[getter]
    fn psf_1sigma(&self) -> f64 {
        self.0.psf_1sigma
    }

    #[getter]
    fn time(&self) -> PyResult<Py<PyDateTime>> {
        let ts = self.0.time.as_unixtime();
        pyo3::Python::with_gil(|py| {
            let dt = PyDateTime::from_timestamp(py, ts, None)?.unbind();
            Ok(dt)
        })
    }

    #[setter]
    fn set_pixel_pitch(&mut self, value: f64) {
        self.0.pixel_pitch = value;
    }

    #[setter]
    fn set_focal_length(&mut self, value: f64) {
        self.0.focal_length = value;
    }

    #[setter]
    fn set_pedestal(&mut self, value: f64) {
        self.0.pedestal = value;
    }

    #[setter]
    fn set_bit_depth(&mut self, value: u8) {
        self.0.bit_depth = value;
    }

    #[setter]
    fn set_gain(&mut self, value: f64) {
        self.0.gain = value;
    }

    #[setter]
    fn set_exposure_time(&mut self, value: f64) {
        self.0.exposure_time = value;
    }

    #[setter]
    fn set_dark_current(&mut self, value: f64) {
        self.0.dark_current = value;
    }

    #[setter]
    fn set_read_noise(&mut self, value: f64) {
        self.0.read_noise = value;
    }

    #[setter]
    fn set_quantum_efficiency(&mut self, value: f64) {
        self.0.solar_weighted_qe = value;
    }

    #[setter]
    fn set_diameter(&mut self, value: f64) {
        self.0.diameter = value;
    }

    #[setter]
    fn set_obscuration(&mut self, value: f64) {
        self.0.obscuration = value;
    }

    #[setter]
    fn set_zodiacal(&mut self, value: f64) {
        self.0.zodiacal_light = value;
    }

    #[setter]
    fn set_psf_1sigma(&mut self, value: f64) {
        self.0.psf_1sigma = value;
    }
}
