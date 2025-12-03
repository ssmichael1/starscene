use numpy::{ndarray::Array2, PyArray2};
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};
use pyo3::Bound;

use super::pypointing::PyFixedInertialPointing;

use crate::StarCamera;

#[pyclass(name = "StarCamera")]
pub struct PyStarCamera {
    camera: StarCamera,
}

impl PyStarCamera {
    fn extract_pixel_format(value: &Bound<'_, PyAny>) -> PyResult<[usize; 2]> {
        if let Ok((width, height)) = value.extract::<(usize, usize)>() {
            return Ok([width, height]);
        }

        if let Ok(arr) = value.extract::<[usize; 2]>() {
            return Ok(arr);
        }

        Err(PyValueError::new_err(
            "pixel_format must be a length-2 sequence of integers",
        ))
    }

    fn apply_kwargs(camera: &mut StarCamera, kwargs: &Bound<'_, PyDict>) -> PyResult<()> {
        for (key, value) in kwargs.iter() {
            let key_str: &str = key.extract()?;
            match key_str {
                "pixel_width" => camera.pixel_format[0] = value.extract()?,
                "pixel_height" => camera.pixel_format[1] = value.extract()?,
                "pixel_format" => camera.pixel_format = Self::extract_pixel_format(&value)?,
                "focal_length_mm" => camera.focal_length_mm = value.extract()?,
                "pixel_pitch_um" => camera.pixel_pitch_um = value.extract()?,
                "ensquared" => camera.ensquared = value.extract()?,
                "read_noise" => camera.read_noise = value.extract()?,
                "dark_current_e_per_s" => camera.dark_current_e_per_s = value.extract()?,
                "solar_weighted_qe" => camera.solar_weighted_qe = value.extract()?,
                "aperture_diameter_mm" => camera.aperture_diameter_mm = value.extract()?,
                "obscuration_diameter_mm" => camera.obscuration_diameter_mm = value.extract()?,
                "transmission" => camera.transmission = value.extract()?,
                "bit_depth" => camera.bit_depth = value.extract()?,
                "dn_per_electron" => camera.dn_per_electron = value.extract()?,
                "full_well_electrons" | "full_well" => {
                    camera.full_well_electrons = value.extract()?;
                }
                "pedestal_electrons" | "pedestal" => {
                    camera.pedestal_electrons = value.extract()?;
                }
                "integration_time" | "integration_time_seconds" => {
                    camera.integration_time_seconds = value.extract()?;
                }
                _ => {
                    return Err(PyValueError::new_err(format!(
                        "Unknown keyword argument '{}'.",
                        key_str
                    )))
                }
            }
        }
        Ok(())
    }
}

#[pymethods]
impl PyStarCamera {
    #[new]
    #[pyo3(signature = (**kwargs))]
    pub fn new(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let mut camera = StarCamera::default();
        if let Some(kwargs) = kwargs {
            Self::apply_kwargs(&mut camera, kwargs)?;
        }

        Ok(Self { camera })
    }

    #[getter]
    pub fn pixel_format(&self) -> (usize, usize) {
        (self.camera.pixel_format[0], self.camera.pixel_format[1])
    }

    #[setter]
    pub fn set_pixel_format(&mut self, value: &Bound<'_, PyAny>) -> PyResult<()> {
        self.camera.pixel_format = Self::extract_pixel_format(value)?;
        Ok(())
    }

    #[getter]
    pub fn pixel_width(&self) -> usize {
        self.camera.pixel_format[0]
    }

    #[setter]
    pub fn set_pixel_width(&mut self, width: usize) -> PyResult<()> {
        self.camera.pixel_format[0] = width;
        Ok(())
    }

    #[getter]
    pub fn pixel_height(&self) -> usize {
        self.camera.pixel_format[1]
    }

    #[setter]
    pub fn set_pixel_height(&mut self, height: usize) -> PyResult<()> {
        self.camera.pixel_format[1] = height;
        Ok(())
    }

    #[getter]
    pub fn focal_length_mm(&self) -> f64 {
        self.camera.focal_length_mm as f64
    }

    #[setter]
    pub fn set_focal_length_mm(&mut self, focal_length_mm: f64) -> PyResult<()> {
        self.camera.focal_length_mm = focal_length_mm as f32;
        Ok(())
    }

    #[getter]
    pub fn pixel_pitch_um(&self) -> f64 {
        self.camera.pixel_pitch_um as f64
    }

    #[setter]
    pub fn set_pixel_pitch_um(&mut self, pixel_pitch_um: f64) -> PyResult<()> {
        self.camera.pixel_pitch_um = pixel_pitch_um as f32;
        Ok(())
    }

    #[getter]
    pub fn ensquared(&self) -> f64 {
        self.camera.ensquared as f64
    }

    #[setter]
    pub fn set_ensquared(&mut self, ensquared: f64) -> PyResult<()> {
        self.camera.ensquared = ensquared as f32;
        Ok(())
    }

    #[getter]
    pub fn read_noise(&self) -> f64 {
        self.camera.read_noise as f64
    }

    #[setter]
    pub fn set_read_noise(&mut self, read_noise: f64) -> PyResult<()> {
        if read_noise < 0.0 {
            return Err(PyValueError::new_err("read_noise cannot be negative"));
        }
        self.camera.read_noise = read_noise as f32;
        Ok(())
    }

    #[getter]
    pub fn dark_current_e_per_s(&self) -> f64 {
        self.camera.dark_current_e_per_s as f64
    }

    #[setter]
    pub fn set_dark_current_e_per_s(&mut self, value: f64) -> PyResult<()> {
        if value < 0.0 {
            return Err(PyValueError::new_err(
                "dark_current_e_per_s cannot be negative",
            ));
        }
        self.camera.dark_current_e_per_s = value as f32;
        Ok(())
    }

    #[getter]
    pub fn solar_weighted_qe(&self) -> f64 {
        self.camera.solar_weighted_qe as f64
    }

    #[setter]
    pub fn set_solar_weighted_qe(&mut self, qe: f64) -> PyResult<()> {
        self.camera.solar_weighted_qe = qe as f32;
        Ok(())
    }

    #[getter]
    pub fn aperture_diameter_mm(&self) -> f64 {
        self.camera.aperture_diameter_mm as f64
    }

    #[setter]
    pub fn set_aperture_diameter_mm(&mut self, value: f64) -> PyResult<()> {
        if value <= 0.0 {
            return Err(PyValueError::new_err(
                "aperture_diameter_mm must be positive",
            ));
        }
        self.camera.aperture_diameter_mm = value as f32;
        Ok(())
    }

    #[getter]
    pub fn obscuration_diameter_mm(&self) -> f64 {
        self.camera.obscuration_diameter_mm as f64
    }

    #[setter]
    pub fn set_obscuration_diameter_mm(&mut self, value: f64) -> PyResult<()> {
        if value < 0.0 {
            return Err(PyValueError::new_err(
                "obscuration_diameter_mm cannot be negative",
            ));
        }
        self.camera.obscuration_diameter_mm = value as f32;
        Ok(())
    }

    #[getter]
    pub fn transmission(&self) -> f64 {
        self.camera.transmission as f64
    }

    #[setter]
    pub fn set_transmission(&mut self, transmission: f64) -> PyResult<()> {
        self.camera.transmission = transmission as f32;
        Ok(())
    }

    #[getter]
    pub fn bit_depth(&self) -> u8 {
        self.camera.bit_depth
    }

    #[setter]
    pub fn set_bit_depth(&mut self, bits: u8) -> PyResult<()> {
        self.camera.bit_depth = bits;
        Ok(())
    }

    #[getter]
    pub fn dn_per_electron(&self) -> f64 {
        self.camera.dn_per_electron as f64
    }

    #[setter]
    pub fn set_dn_per_electron(&mut self, value: f64) -> PyResult<()> {
        self.camera.dn_per_electron = value as f32;
        Ok(())
    }

    #[getter]
    pub fn full_well_electrons(&self) -> f64 {
        self.camera.full_well_electrons as f64
    }

    #[setter]
    pub fn set_full_well_electrons(&mut self, value: f64) -> PyResult<()> {
        if value <= 0.0 {
            return Err(PyValueError::new_err(
                "full_well_electrons must be positive",
            ));
        }
        self.camera.full_well_electrons = value as f32;
        Ok(())
    }

    #[getter]
    pub fn pedestal_electrons(&self) -> f64 {
        self.camera.pedestal_electrons as f64
    }

    #[setter]
    pub fn set_pedestal_electrons(&mut self, value: f64) -> PyResult<()> {
        if value < 0.0 {
            return Err(PyValueError::new_err(
                "pedestal_electrons cannot be negative",
            ));
        }
        self.camera.pedestal_electrons = value as f32;
        Ok(())
    }

    #[getter]
    pub fn integration_time_seconds(&self) -> f64 {
        self.camera.integration_time_seconds as f64
    }

    #[setter]
    pub fn set_integration_time_seconds(&mut self, value: f64) -> PyResult<()> {
        if value <= 0.0 {
            return Err(PyValueError::new_err(
                "integration_time_seconds must be positive",
            ));
        }
        self.camera.integration_time_seconds = value as f32;
        Ok(())
    }

    #[getter]
    pub fn max_dn(&self) -> u32 {
        self.camera.max_dn()
    }

    #[getter]
    pub fn ifov_rad(&self) -> f64 {
        self.camera.ifov_rad() as f64
    }

    pub fn fov_rad(&self) -> (f64, f64) {
        let (width_rad, height_rad) = self.camera.fov_rad();
        (width_rad as f64, height_rad as f64)
    }

    pub fn fov_deg(&self) -> (f64, f64) {
        let (w, h) = self.camera.fov_deg();
        (w as f64, h as f64)
    }

    pub fn aperture_area_m2(&self) -> f64 {
        self.camera.aperture_area_m2() as f64
    }

    pub fn load_star_catalog_from_file(&mut self, path: &str) -> PyResult<()> {
        self.camera
            .load_star_catalog_from_file(path)
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }

    #[getter]
    pub fn has_star_catalog(&self) -> bool {
        self.camera.has_star_catalog()
    }

    #[getter]
    pub fn star_count(&self) -> usize {
        self.camera.star_count()
    }

    pub fn frame<'py>(
        &self,
        py: Python<'py>,
        pointing: &PyFixedInertialPointing,
    ) -> PyResult<Bound<'py, PyArray2<u16>>> {
        let image = self
            .camera
            .frame(pointing.as_provider())
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

        let width = image.width() as usize;
        let height = image.height() as usize;
        let mut data = Vec::with_capacity(width * height);
        for pixel in image.pixels() {
            data.push(pixel[0]);
        }

        let array = Array2::from_shape_vec((height, width), data)
            .map_err(|_| PyRuntimeError::new_err("Failed to reshape rendered image buffer"))?;
        Ok(PyArray2::from_owned_array(py, array))
    }

    pub fn render_star_image<'py>(
        &self,
        py: Python<'py>,
        pointing: &PyFixedInertialPointing,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let image = self
            .camera
            .render_star_image(pointing.as_provider())
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

        let width = image.width() as usize;
        let height = image.height() as usize;
        let mut data = Vec::with_capacity(width * height);
        for pixel in image.pixels() {
            data.push(pixel[0]);
        }

        let array = Array2::from_shape_vec((height, width), data)
            .map_err(|_| PyRuntimeError::new_err("Failed to reshape rendered image buffer"))?;
        Ok(PyArray2::from_owned_array(py, array))
    }

    fn __str__(&self) -> PyResult<String> {
        Ok(format!(
            "StarCamera:\n  resolution: {}x{} px\n  focal_length_mm: {:.2}\n  pixel_pitch_um: {:.2}\n  ensquared: {:.3}\n  integration_time_s: {:.3}\n  read_noise_e: {:.2}\n  dark_current_e_per_s: {:.3}\n  pedestal_electrons: {:.2}\n  full_well_electrons: {:.0}\n  star_count: {}",
            self.camera.pixel_format[0],
            self.camera.pixel_format[1],
            self.camera.focal_length_mm,
            self.camera.pixel_pitch_um,
            self.camera.ensquared,
            self.camera.integration_time_seconds,
            self.camera.read_noise,
            self.camera.dark_current_e_per_s,
            self.camera.pedestal_electrons,
            self.camera.full_well_electrons,
            self.camera.star_count()
        ))
    }
}
