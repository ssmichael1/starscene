mod pycamera;
mod pycatalog;
mod pypointing;
mod pystar;

use pyo3::prelude::*;

#[pymodule]
pub fn starscene(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<pystar::PyStar>()?;
    m.add_class::<pycatalog::PyStarKdTree>()?;
    m.add_class::<pypointing::PyFixedInertialPointing>()?;
    m.add_class::<pycamera::PyStarCamera>()?;
    m.add_function(wrap_pyfunction!(pycatalog::load_tycho_star_list, m)?)?;
    Ok(())
}
