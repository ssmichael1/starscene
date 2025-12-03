//! Module to load tycho2 star catalog from raw data files
//!
//! Tycho 2 catalog format at
//! http://tdc-www.harvard.edu/catalogs/tycho2.format.html
//!
//! Files are at:
//! https://archive.eso.org/ASTROM/TYC-2/data/
//!
//! Files:
//!
//! catalog.dat: main catalog file (~ 2.5 Million stars)
//! suppl_1.dat: supplementary file 1 (Stars with high-quality data missing from main catalog)
//! suppl_2.dat: supplementary file 2 (Stars with astrometry only; poor photometry)
//!
//!
//!

use anyhow::{bail, Context};
use rkyv::{Archive, Deserialize, Serialize};

#[derive(Archive, Serialize, Deserialize, Debug, Clone)]
pub struct Tycho2Star {
    pub id: (u16, u16, u8),
    pub mean_ra_deg: f64,
    pub mean_dec_deg: f64,
    pub pm_ra_cosdec_mas_yr: Option<f32>,
    pub pm_dec_mas_yr: Option<f32>,
    pub mean_ra_cosdec_err_mas: Option<f32>,
    pub mean_dec_err_mas: Option<f32>,
    pub pm_ra_cosdec_err_mas_yr: Option<f32>,
    pub pm_dec_err_mas_yr: Option<f32>,
    pub bt_mag: Option<f32>,
    pub vt_mag: Option<f32>,
}

impl Tycho2Star {
    pub fn catalog_id(&self) -> String {
        format!("TYC {}-{}-{}", self.id.0, self.id.1, self.id.2)
    }

    pub fn johnson_v_mag(&self) -> Option<f32> {
        if let (Some(bt), Some(vt)) = (self.bt_mag, self.vt_mag) {
            // Approximate conversion from Tycho VT to Johnson V
            return Some(vt - 0.090 * (bt - vt));
        }
        None
    }
}

pub mod loader {}
fn parse_ids(field: &str) -> anyhow::Result<(u16, u16, u8)> {
    let mut parts = field.split_whitespace();
    let id1 = parts
        .next()
        .context("missing TYC1")?
        .parse::<u16>()
        .context("parsing TYC1")?;
    let id2 = parts
        .next()
        .context("missing TYC2")?
        .parse::<u16>()
        .context("parsing TYC2")?;
    let id3 = parts
        .next()
        .context("missing TYC3")?
        .parse::<u8>()
        .context("parsing TYC3")?;
    Ok((id1, id2, id3))
}

fn parse_required_f64(field: Option<&str>, label: &str) -> anyhow::Result<f64> {
    let value = field.context(format!("missing {}", label))?.trim();
    Ok(value.parse::<f64>().context(format!("parsing {}", label))?)
}

fn parse_optional_f32(field: Option<&str>, label: &str) -> anyhow::Result<Option<f32>> {
    if let Some(raw) = field {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            Ok(None)
        } else {
            Ok(Some(
                trimmed
                    .parse::<f32>()
                    .context(format!("parsing {}", label))?,
            ))
        }
    } else {
        Ok(None)
    }
}

fn parse_tycho_line(line: &str) -> anyhow::Result<Option<Tycho2Star>> {
    let fields: Vec<&str> = line.split('|').collect();
    if fields.len() < 20 {
        bail!(
            "unexpected field count {} in Tycho2 main line",
            fields.len()
        );
    }

    let (id1, id2, id3) = parse_ids(fields[0])?;
    let pflag = fields.get(1).map(|f| f.trim()).unwrap_or("");
    if pflag == "X" {
        return Ok(None);
    }

    let mean_ra_deg = parse_required_f64(fields.get(2).copied(), "mean RA deg")?;
    let mean_dec_deg = parse_required_f64(fields.get(3).copied(), "mean Dec deg")?;
    let pm_ra_cosdec_mas_yr = parse_optional_f32(fields.get(4).copied(), "pm RA*cos(dec)")?;
    let pm_dec_mas_yr = parse_optional_f32(fields.get(5).copied(), "pm Dec")?;
    let mean_ra_cosdec_err_mas =
        parse_optional_f32(fields.get(6).copied(), "mean RA*cos(dec) err")?;
    let mean_dec_err_mas = parse_optional_f32(fields.get(7).copied(), "mean Dec err")?;
    let pm_ra_cosdec_err_mas_yr = parse_optional_f32(fields.get(8).copied(), "pm RA*cos(dec) err")?;
    let pm_dec_err_mas_yr = parse_optional_f32(fields.get(9).copied(), "pm Dec err")?;
    let bt_mag = parse_optional_f32(fields.get(17).copied(), "BT magnitude")?;
    let vt_mag = parse_optional_f32(fields.get(19).copied(), "VT magnitude")?;

    if let Some(vt) = vt_mag {
        if vt < 0.0_f32 {
            println!(
                "Warning: star {} has negative VT mag {}",
                format!("TYC {}-{}-{}", id1, id2, id3),
                vt
            );
        }
    }

    Ok(Some(Tycho2Star {
        id: (id1, id2, id3),
        mean_ra_deg,
        mean_dec_deg,
        pm_ra_cosdec_mas_yr,
        pm_dec_mas_yr,
        mean_ra_cosdec_err_mas,
        mean_dec_err_mas,
        pm_ra_cosdec_err_mas_yr,
        pm_dec_err_mas_yr,
        bt_mag,
        vt_mag,
    }))
}

fn parse_tycho_suppl_1_line(line: &str) -> anyhow::Result<Option<Tycho2Star>> {
    let fields: Vec<&str> = line.split('|').collect();
    if fields.len() < 15 {
        bail!(
            "unexpected field count {} in Tycho2 suppl1 line",
            fields.len()
        );
    }

    let (id1, id2, id3) = parse_ids(fields[0])?;
    let pflag = fields.get(1).map(|f| f.trim()).unwrap_or("");
    if pflag == "X" {
        return Ok(None);
    }

    let mean_ra_deg = parse_required_f64(fields.get(2).copied(), "mean RA deg")?;
    let mean_dec_deg = parse_required_f64(fields.get(3).copied(), "mean Dec deg")?;
    let pm_ra_cosdec_mas_yr = parse_optional_f32(fields.get(4).copied(), "pm RA*cos(dec)")?;
    let pm_dec_mas_yr = parse_optional_f32(fields.get(5).copied(), "pm Dec")?;
    let bt_mag = parse_optional_f32(fields.get(11).copied(), "BT magnitude")?;
    let vt_mag = parse_optional_f32(fields.get(13).copied(), "VT magnitude")?;

    if let Some(vt) = vt_mag {
        if vt < 0.0_f32 {
            println!(
                "Warning: star {} has negative VT mag {}",
                format!("TYC {}-{}-{}", id1, id2, id3),
                vt
            );
        }
    }

    Ok(Some(Tycho2Star {
        id: (id1, id2, id3),
        mean_ra_deg,
        mean_dec_deg,
        pm_ra_cosdec_mas_yr,
        pm_dec_mas_yr,
        mean_ra_cosdec_err_mas: None,
        mean_dec_err_mas: None,
        pm_ra_cosdec_err_mas_yr: None,
        pm_dec_err_mas_yr: None,
        bt_mag,
        vt_mag,
    }))
}

pub fn load_tycho2_suppl1_catalog(filename: Option<String>) -> anyhow::Result<Vec<Tycho2Star>> {
    use anyhow::Context;
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    let filename = filename.unwrap_or_else(|| "data/tycho2_suppl_1.dat".to_string());

    let file = File::open(filename).context("opening tycho2 suppl_1.dat")?;
    let reader = BufReader::new(file);

    let stars = reader
        .lines()
        .enumerate()
        .into_iter()
        .map(|(line_num, line_res)| {
            let line = line_res.context(format!("reading line {}", line_num + 1))?;
            parse_tycho_suppl_1_line(&line).context(format!("parsing line {}", line_num + 1))
        })
        .collect::<anyhow::Result<Vec<Option<Tycho2Star>>>>()?;

    // Remove None entries (stars with pflag 'X')
    let stars: Vec<Tycho2Star> = stars.into_iter().filter_map(|s| s).collect();

    Ok(stars)
}

pub fn load_tycho2_main_catalog(filename: Option<String>) -> anyhow::Result<Vec<Tycho2Star>> {
    use anyhow::Context;
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    let filename = filename.unwrap_or_else(|| "data/tycho2_catalog.dat".to_string());

    let file = File::open(filename).context("opening tycho2 catalog.dat")?;
    let reader = BufReader::new(file);

    let stars = reader
        .lines()
        .enumerate()
        .into_iter()
        .map(|(line_num, line_res)| {
            let line = line_res.context(format!("reading line {}", line_num + 1))?;
            parse_tycho_line(&line).context(format!("parsing line {}", line_num + 1))
        })
        .collect::<anyhow::Result<Vec<Option<Tycho2Star>>>>()?;

    // Remove None entries (stars with pflag 'X')
    let stars: Vec<Tycho2Star> = stars.into_iter().filter_map(|s| s).collect();

    Ok(stars)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_load_tycho2_catalog() {
        let stars = load_tycho2_main_catalog(None).expect("loading tycho2 catalog");
        assert!(stars.len() > 2_000_000, "expected over 2 million stars");
        let first_star = &stars[0];
        println!("First star = {:?}", first_star);
    }

    #[test]
    fn test_load_suppl_1_catalog() {
        let stars = load_tycho2_suppl1_catalog(None).expect("loading tycho2 suppl1 catalog");
        let first_star = &stars[0];
        println!("First suppl1 star = {:?}", first_star);
    }

    #[test]
    #[ignore]

    fn create_catalog() {
        let mut all_stars = load_tycho2_main_catalog(None).expect("loading tycho2 main catalog");
        let suppl1_stars = load_tycho2_suppl1_catalog(None).expect("loading tycho2 suppl1 catalog");
        println!("Loaded {} suppl1 stars", suppl1_stars.len());
        all_stars.extend(suppl1_stars);
        println!("Total stars loaded: {}", all_stars.len());

        // Map to generic star struct
        let stars: Vec<crate::Star> = all_stars
            .into_iter()
            .map(|s| crate::Star {
                j2000_vec: {
                    let ra_rad = s.mean_ra_deg.to_radians();
                    let dec_rad = s.mean_dec_deg.to_radians();
                    let x = dec_rad.cos() * ra_rad.cos();
                    let y = dec_rad.cos() * ra_rad.sin();
                    let z = dec_rad.sin();
                    [x, y, z]
                },
                //v_mag: s.vt_mag.unwrap_or(99.9) as f64,
                v_mag: s.johnson_v_mag().unwrap_or(99.9) as f64,
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
        let mut file = File::create("data/tycho2_catalog.rkyv").expect("creating output file");
        file.write_all(&bytes).expect("writing rkyv data to file");
    }

    #[test]
    #[ignore]
    /// Load Tycho2 catalog and serialize to rkyv format
    fn load_all() {
        let main_stars = load_tycho2_main_catalog(None).expect("loading tycho2 main catalog");
        println!("Loaded {} main catalog stars", main_stars.len());

        let suppl1_stars = load_tycho2_suppl1_catalog(None).expect("loading tycho2 suppl1 catalog");
        println!("Loaded {} suppl1 catalog stars", suppl1_stars.len());

        let stars = [main_stars, suppl1_stars].concat();

        // Map to generic star struct
        let stars: Vec<crate::Star> = stars
            .into_iter()
            .map(|s| crate::Star {
                j2000_vec: {
                    let ra_rad = s.mean_ra_deg.to_radians();
                    let dec_rad = s.mean_dec_deg.to_radians();
                    let x = dec_rad.cos() * ra_rad.cos();
                    let y = dec_rad.cos() * ra_rad.sin();
                    let z = dec_rad.sin();
                    [x, y, z]
                },
                v_mag: s.johnson_v_mag().unwrap_or(99.9) as f64,
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
        let mut file = File::create("data/tycho2_catalog.rkyv").expect("creating output file");
        file.write_all(&bytes).expect("writing rkyv data to file");
    }
}
