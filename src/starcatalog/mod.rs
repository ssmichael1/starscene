use crate::SSResult;

use std::io::Read;

#[derive(Debug, Clone)]
pub struct Star {
    pub ra: f64,
    pub dec: f64,
    pub ra_pm: f32,
    pub dec_pm: f32,
    pub bt_mag: f32,
    pub vt_mag: f32,
}

impl Star {
    /// Calculate the apparent visual magnitude of the star.
    ///
    /// # Returns
    /// Apparent visual magnitude of the star.
    ///
    /// # Notes:
    /// See: https://heasarc.gsfc.nasa.gov/w3browse/all/tycho2.html
    pub fn mv(&self) -> f32 {
        self.vt_mag - 0.09 * (self.bt_mag - self.vt_mag)
    }
}

#[derive(Debug)]
pub struct StarCatalog {
    pub depth: usize,
    pub nindex: usize,
    pub starhash: Vec<Vec<Star>>,
}

impl StarCatalog {
    /// Query the star catalog for stars within a given radius of a given RA and Dec.
    ///
    ///
    /// # Arguments
    ///
    /// * `ra` - The right ascension in radians.
    /// * `dec` - The declination in radians.
    /// * `radius` - The radius in radians.
    ///
    /// # Returns
    ///
    /// A vector of stars within the given radius.
    ///
    pub fn query(&self, ra: f64, dec: f64, radius: f64) -> Vec<Star> {
        let v = cdshealpix::nested::get(self.depth as u8);
        let res = v.cone_coverage_approx(ra, dec, radius);
        let mut stars = Vec::new();
        res.flat_iter().for_each(|i| {
            stars.extend(self.starhash[i as usize].iter().cloned());
        });
        stars
    }

    /// Load a star catalog from a file.
    ///
    /// # Arguments
    ///
    /// * `filename` - The name of the file to load.
    ///
    /// # Returns
    ///
    /// A StarCatalog object.
    ///
    pub fn load(filename: &str) -> SSResult<StarCatalog> {
        let file = std::fs::File::open(filename)?;
        let mut reader = std::io::BufReader::new(file);
        let mut buffer = [0; 4];
        let mut buffer8 = [0; 8];
        reader.read_exact(&mut buffer)?;
        let depth = u32::from_le_bytes(buffer) as usize;
        if depth > 12 {
            return Err("nside too large".into());
        }
        reader.read_exact(&mut buffer)?;
        let nindex = u32::from_le_bytes(buffer) as usize;
        if nindex > 1000000 {
            return Err("index too large".into());
        }
        Ok(StarCatalog {
            depth,
            nindex,
            starhash: {
                let mut v = Vec::with_capacity(nindex);
                for _ in 0..nindex {
                    reader.read_exact(&mut buffer)?;
                    let len = u32::from_le_bytes(buffer);
                    let mut stars = Vec::with_capacity(len as usize);
                    for _ in 0..len {
                        reader.read_exact(&mut buffer8)?;
                        let ra = f64::from_le_bytes(buffer8);
                        reader.read_exact(&mut buffer8)?;
                        let dec = f64::from_le_bytes(buffer8);
                        reader.read_exact(&mut buffer)?;
                        let ra_pm = f32::from_le_bytes(buffer);
                        reader.read_exact(&mut buffer)?;
                        let dec_pm = f32::from_le_bytes(buffer);
                        reader.read_exact(&mut buffer)?;
                        let bt_mag = f32::from_le_bytes(buffer);
                        reader.read_exact(&mut buffer)?;
                        let vt_mag = f32::from_le_bytes(buffer);
                        stars.push(Star {
                            ra,
                            dec,
                            ra_pm,
                            dec_pm,
                            bt_mag,
                            vt_mag,
                        });
                    }
                    v.push(stars);
                }
                v
            },
        })
    }
}

pub fn the_tycho2_catalog() -> &'static StarCatalog {
    static INSTANCE: once_cell::sync::OnceCell<StarCatalog> = once_cell::sync::OnceCell::new();
    INSTANCE.get_or_init(|| StarCatalog::load("tycho2_hash.dat").unwrap())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_star_catalog() {
        let catalog = StarCatalog::load("tycho2_hash.dat").unwrap();
        let s = catalog.query(
            30.0_f64.to_radians(),
            20.0_f64.to_radians(),
            0.1_f64.to_radians(),
        );
        println!("s = {:?}", s);
    }

    #[test]
    fn test_healpix() {
        let nside = 5;
        let v = cdshealpix::nested::get(nside);
        let res = v.cone_coverage_approx(
            90.0_f64.to_radians(),
            20.0_f64.to_radians(),
            1.0_f64.to_radians(),
        );
        println!("res = {:?}", res);
    }
}
