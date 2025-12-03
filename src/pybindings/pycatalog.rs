use std::fs::File;
use std::io::Read;

use pyo3::exceptions::{PyIOError, PyRuntimeError};
use pyo3::prelude::*;

use crate::catalogs::{tycho2, Star, StarKdTree};

use super::pystar::PyStar;

#[pyclass(name = "StarKdTree")]
pub struct PyStarKdTree {
    tree: StarKdTree,
}

#[pymethods]
impl PyStarKdTree {
    #[staticmethod]
    pub fn from_rkyv(path: &str) -> PyResult<Self> {
        let mut file = File::open(path).map_err(|e| PyErr::new::<PyIOError, _>(e.to_string()))?;
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer)
            .map_err(|e| PyErr::new::<PyIOError, _>(e.to_string()))?;

        let tree = rkyv::from_bytes::<StarKdTree, rkyv::rancor::Error>(&buffer).map_err(|e| {
            PyErr::new::<PyRuntimeError, _>(format!("Failed to deserialize KD-tree: {e}"))
        })?;

        Ok(Self { tree })
    }

    #[staticmethod]
    pub fn from_tycho_catalog(
        main_path: Option<&str>,
        suppl1_path: Option<&str>,
    ) -> PyResult<Self> {
        let main = tycho2::load_tycho2_main_catalog(main_path.map(|s| s.to_string()))
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(e.to_string()))?;
        let suppl1 = tycho2::load_tycho2_suppl1_catalog(suppl1_path.map(|s| s.to_string()))
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(e.to_string()))?;

        let stars: Vec<Star> = main
            .into_iter()
            .chain(suppl1.into_iter())
            .map(|s| {
                let ra_rad = s.mean_ra_deg.to_radians();
                let dec_rad = s.mean_dec_deg.to_radians();
                let x = dec_rad.cos() * ra_rad.cos();
                let y = dec_rad.cos() * ra_rad.sin();
                let z = dec_rad.sin();
                Star {
                    j2000_vec: [x, y, z],
                    v_mag: s.johnson_v_mag().unwrap_or(99.9) as f64,
                }
            })
            .collect();

        Ok(Self {
            tree: StarKdTree::build(stars),
        })
    }

    pub fn nearest_neighbor(&self, point: (f64, f64, f64)) -> Option<(PyStar, f64)> {
        let query = [point.0, point.1, point.2];
        self.tree
            .nearest_neighbor(query)
            .map(|(idx, dist)| (self.tree.items[idx].clone().into(), dist))
    }

    pub fn radius_search(&self, point: (f64, f64, f64), radius: f64) -> Vec<PyStar> {
        let query = [point.0, point.1, point.2];
        self.tree
            .radius_search(query, radius)
            .into_iter()
            .map(|idx| self.tree.items[idx].clone().into())
            .collect()
    }

    pub fn len(&self) -> usize {
        self.tree.items.len()
    }
}

#[pyfunction]
pub fn load_tycho_star_list(
    main_path: Option<&str>,
    suppl1_path: Option<&str>,
) -> PyResult<Vec<PyStar>> {
    let main = tycho2::load_tycho2_main_catalog(main_path.map(|s| s.to_string()))
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(e.to_string()))?;
    let suppl1 = tycho2::load_tycho2_suppl1_catalog(suppl1_path.map(|s| s.to_string()))
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(e.to_string()))?;

    let stars = main
        .into_iter()
        .chain(suppl1.into_iter())
        .map(|s| {
            let ra_rad = s.mean_ra_deg.to_radians();
            let dec_rad = s.mean_dec_deg.to_radians();
            let x = dec_rad.cos() * ra_rad.cos();
            let y = dec_rad.cos() * ra_rad.sin();
            let z = dec_rad.sin();
            PyStar::from(Star {
                j2000_vec: [x, y, z],
                v_mag: s.johnson_v_mag().unwrap_or(99.9) as f64,
            })
        })
        .collect();

    Ok(stars)
}
