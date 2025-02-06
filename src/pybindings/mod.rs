use pyo3::prelude::*;
use pyo3::{wrap_pyfunction, wrap_pymodule};

mod pystarscene;

#[pymodule]
pub fn starscene(py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    let satkit = py.import("satkit")?;
    let instant = satkit.getattr("time")?;
    m.add_class::<satkit::pybindings::PyInstant>()?;

    m.add_class::<pystarscene::PyStarScene>()?;
    Ok(())
}
