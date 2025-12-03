//! Module to load stars from the Gaia catalog.
//!

use anyhow::{Context, Result};
use rkyv::{Archive, Deserialize, Serialize};
use std::path::Path;

#[derive(Archive, Serialize, Deserialize, Debug, Clone)]
pub struct GaiaStar {
    pub source_id: u64,
    pub ra_deg: f64,
    pub dec_deg: f64,
    pub phot_g_mean_mag: f64,
    pub phot_bp_mean_mag: Option<f64>,
    pub phot_rp_mean_mag: Option<f64>,
    pub parallax_mas: Option<f64>,
    pub pmra_mas_yr: Option<f64>,
    pub pmdec_mas_yr: Option<f64>,
}

impl GaiaStar {
    pub fn j2000_vec(&self) -> [f64; 3] {
        let ra_rad = self.ra_deg.to_radians();
        let dec_rad = self.dec_deg.to_radians();
        let x = dec_rad.cos() * ra_rad.cos();
        let y = dec_rad.cos() * ra_rad.sin();
        let z = dec_rad.sin();
        [x, y, z]
    }
}

pub fn load_gaia_csv<T: AsRef<Path>>(path: T) -> Result<Vec<GaiaStar>> {
    let mut rdr = csv::Reader::from_path(path.as_ref())
        .with_context(|| format!("Failed to open Gaia CSV file: {}", path.as_ref().display()))?;

    let mut stars = Vec::new();
    for result in rdr.records() {
        let record = result.context("Failed to read a record from Gaia CSV")?;
        let source_id: u64 = record
            .get(0)
            .ok_or_else(|| anyhow::anyhow!("Missing source_id field"))?
            .parse()
            .context("Failed to parse source_id")?;
        let ra_deg: f64 = record
            .get(1)
            .ok_or_else(|| anyhow::anyhow!("Missing ra_deg field"))?
            .parse()
            .context("Failed to parse ra_deg")?;
        let dec_deg: f64 = record
            .get(2)
            .ok_or_else(|| anyhow::anyhow!("Missing dec_deg field"))?
            .parse()
            .context("Failed to parse dec_deg")?;
        let phot_g_mean_mag: f64 = record
            .get(3)
            .ok_or_else(|| anyhow::anyhow!("Missing phot_g_mean_mag field"))?
            .parse()
            .context("Failed to parse phot_g_mean_mag")?;

        let phot_bp_mean_mag = record.get(4).and_then(|s| s.parse().ok());
        let phot_rp_mean_mag = record.get(5).and_then(|s| s.parse().ok());
        let parallax_mas = record.get(6).and_then(|s| s.parse().ok());
        let pmra_mas_yr = record.get(7).and_then(|s| s.parse().ok());
        let pmdec_mas_yr = record.get(8).and_then(|s| s.parse().ok());

        stars.push(GaiaStar {
            source_id,
            ra_deg,
            dec_deg,
            phot_g_mean_mag,
            phot_bp_mean_mag,
            phot_rp_mean_mag,
            parallax_mas,
            pmra_mas_yr,
            pmdec_mas_yr,
        });
    }

    Ok(stars)
}

pub fn create_catalog(gaia_stars: &[GaiaStar]) -> Result<()> {
    // Placeholder for catalog creation logic

    let stars: Vec<crate::Star> = gaia_stars
        .iter()
        .map(|gs| crate::Star {
            j2000_vec: gs.j2000_vec(),
            v_mag: gs.phot_g_mean_mag,
        })
        .collect();
    let kd_tree = crate::catalogs::StarKdTree::build(stars);
    let bytes = rkyv::to_bytes::<rkyv::rancor::Error>(&kd_tree).unwrap();
    println!(
        "Serialized KD-tree with {} stars into {} bytes",
        kd_tree.items.len(),
        bytes.len()
    );

    // Write to file
    use std::fs::File;
    use std::io::Write;
    let mut file = File::create("data/gaia_catalog.rkyv").expect("creating output file");
    file.write_all(&bytes).expect("writing rkyv data to file");

    Ok(())
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    #[ignore]
    fn test_create_catalog() {
        let gaia_stars = load_gaia_csv("data/gaia_bright_stars.csv").expect("loading Gaia CSV");
        assert!(!gaia_stars.is_empty());
        println!("loaded {} Gaia stars", gaia_stars.len());

        create_catalog(&gaia_stars).expect("creating Gaia catalog");
    }
}
