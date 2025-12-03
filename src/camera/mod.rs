pub mod pointing;
use pointing::PointingProvider;

use rand::rngs::SmallRng;
use rand::SeedableRng;
use rand_distr::Distribution;
use rayon::prelude::*;
use rayon::slice::ParallelSlice;
use rayon::slice::ParallelSliceMut;

use crate::Image;
use crate::StarKdTree;

mod mathutils;
use mathutils::{erf_diff_lut, erf_inv};
pub struct StarCamera {
    pub pixel_format: [usize; 2],
    pub focal_length_mm: f32,
    pub pixel_pitch_um: f32,
    pub ensquared: f32,
    pub read_noise: f32,
    pub dark_current_e_per_s: f32,
    pub solar_weighted_qe: f32,
    pub aperture_diameter_mm: f32,
    pub obscuration_diameter_mm: f32,
    pub transmission: f32,
    pub bit_depth: u8,
    pub dn_per_electron: f32,
    pub pedestal_electrons: f32,
    pub full_well_electrons: f32,
    pub integration_time_seconds: f32,
    star_catalog: Option<crate::StarKdTree>,
}

impl Default for StarCamera {
    fn default() -> Self {
        StarCamera {
            pixel_format: [2048, 2048],
            focal_length_mm: 50.0,
            pixel_pitch_um: 3.76,
            ensquared: 0.7,
            read_noise: 5.0,
            dark_current_e_per_s: 0.1,
            solar_weighted_qe: 0.9,
            aperture_diameter_mm: 25.0,
            obscuration_diameter_mm: 0.0,
            transmission: 0.9,
            bit_depth: 12,
            dn_per_electron: 1.0,
            pedestal_electrons: 0.0,
            full_well_electrons: 1.0e8,
            integration_time_seconds: 1.0,
            star_catalog: None,
        }
    }
}

/// SplitMix64 hash function for seeding RNGs
fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E3779B97F4A7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

impl StarCamera {
    pub fn load_star_catalog_from_file(&mut self, rkyv_file: &str) -> anyhow::Result<()> {
        use std::fs::File;
        use std::io::Read;

        let mut file = File::open(rkyv_file)?;
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer)?;

        let stars = rkyv::from_bytes::<StarKdTree, rkyv::rancor::Error>(&buffer)
            .map_err(|e| anyhow::anyhow!("Failed to deserialize star catalog: {}", e))?;

        self.star_catalog = Some(stars);

        Ok(())
    }

    /// Returns true when a catalog has been loaded
    pub fn has_star_catalog(&self) -> bool {
        self.star_catalog.is_some()
    }

    /// Returns the number of stars in the loaded catalog
    pub fn star_count(&self) -> usize {
        self.star_catalog
            .as_ref()
            .map(|catalog| catalog.items.len())
            .unwrap_or(0)
    }

    /// Instantaneous field of view per pixel in radians
    ///
    ///
    /// # Returns
    /// Instantaneous field of view per pixel in radians
    ///
    pub fn ifov_rad(&self) -> f32 {
        let pixel_pitch_mm = self.pixel_pitch_um * 1e-3;
        pixel_pitch_mm / self.focal_length_mm
    }

    /// Maximum digital number (DN) for the camera
    /// given its bit depth
    pub fn max_dn(&self) -> u32 {
        (1 << self.bit_depth) as u32 - 1
    }

    /// Field of view in radians
    ///
    /// # Returns
    ///
    /// Field of view in radians (x, y)
    pub fn fov_rad(&self) -> (f32, f32) {
        (
            (self.pixel_pitch_um * 1.0e-3 * self.pixel_format[0] as f32
                / self.focal_length_mm
                / 2.0)
                .atan()
                * 2.0,
            (self.pixel_pitch_um * 1.0e-3 * self.pixel_format[1] as f32
                / self.focal_length_mm
                / 2.0)
                .atan()
                * 2.0,
        )
    }

    /// Field of view in degrees
    ///
    /// # Returns
    ///
    /// Field of view in degrees (x, y)
    pub fn fov_deg(&self) -> (f32, f32) {
        let (fov_x_rad, fov_y_rad) = self.fov_rad();
        (fov_x_rad.to_degrees(), fov_y_rad.to_degrees())
    }

    /// Aperture area in m^2
    ///
    /// # Returns
    /// Aperture area in m^2
    pub fn aperture_area_m2(&self) -> f32 {
        let r_outer = self.aperture_diameter_mm * 0.5 * 1e-3;
        let r_inner = self.obscuration_diameter_mm * 0.5 * 1e-3;
        std::f32::consts::PI * (r_outer * r_outer - r_inner * r_inner)
    }

    /// Total detected photon flux for a mv=0 star in photons per second
    fn mv0_flux_photons_per_s(&self) -> f32 {
        // Total number of solar photons
        const SOLAR_PHOTONS_M2_S: f32 = 6.8e21;
        const MV0_MAG: f32 = -26.74;

        // Zero magnitude star flu7x at Earth in photons/m^2/s
        SOLAR_PHOTONS_M2_S
            * 10.0_f32.powf(MV0_MAG * 0.4)
            * self.aperture_area_m2()
            * self.transmission
            * self.solar_weighted_qe
    }

    fn ensq_to_1sigma_radius_pixels(&self) -> f32 {
        // Approximate conversion from ensquared energy to 1-sigma radius in pixels
        // Assumes a Gaussian PSF
        1.0 / (8.0_f32.sqrt() * erf_inv(self.ensquared as f64) as f32)
    }

    pub fn frame(&self, pointing: &impl PointingProvider) -> anyhow::Result<Image<u16>> {
        let mut expected_electrons =
            crate::Image::<f64>::new(self.pixel_format[0] as u32, self.pixel_format[1] as u32);

        // Add dark current
        //expected_electrons.fill((self.dark_current_e_per_s * self.integration_time_seconds) as f64);

        self.add_stars_to_image(&mut expected_electrons, pointing);

        let maxdn = self.max_dn() as f32;

        let seed = 0x12345678u64;
        let mut out = vec![0_u16; expected_electrons.len()];

        out.par_chunks_exact_mut(self.pixel_format[0])
            .zip(
                expected_electrons
                    .as_raw()
                    .par_chunks_exact(self.pixel_format[0]),
            )
            .enumerate()
            .for_each(|(cid, (o, l))| {
                let row_seed = splitmix64(seed ^ (cid as u64));
                let mut rng = SmallRng::seed_from_u64(row_seed);
                let rdist = rand_distr::Normal::new(0.0, self.read_noise as f32).unwrap();

                for ix in 0..o.len() {
                    let pval = mathutils::poisson_u32(&mut rng, l[ix]);
                    let elec = (pval as f32 + rdist.sample(&mut rng)) as f32;
                    //let elec = l[0] as f32;
                    o[ix] = ((elec + self.pedestal_electrons) * self.dn_per_electron)
                        .min(maxdn)
                        .max(0_f32) as u16;
                }
            });

        Ok(Image::<u16>::from_raw(
            self.pixel_format[0] as u32,
            self.pixel_format[1] as u32,
            out,
        )
        .unwrap())
    }

    // Add expected star photon counts to image buffer
    fn add_stars_to_image(&self, buf: &mut crate::Image<f64>, pointing: &impl PointingProvider) {
        let psf_sigma_pixels = self.ensq_to_1sigma_radius_pixels();
        let mv0_signal = self.mv0_flux_photons_per_s() * self.integration_time_seconds;

        // Implementation goes here
        let qcamfromj2000 = pointing.camera_from_j2000();
        let boresight = pointing.boresight_j2000();
        let fov_rad = self.fov_rad();
        let search_radius = (fov_rad.0.max(fov_rad.1) as f64 * 2.0_f64.sqrt()) + 0.01; // Add small margin
        let stars = match &self.star_catalog {
            Some(catalog) => catalog,
            None => {
                tracing::warn!("No star catalog loaded");
                return;
            }
        };

        // Get stars that may be in the field of view
        let nearby_star_indices =
            stars.radius_search([boresight[0], boresight[1], boresight[2]], search_radius);

        // PSF range over which to sum star signal
        let base_psf_range = (psf_sigma_pixels * 3.0).ceil() as isize;
        let max_psf_range = (base_psf_range * 4).max(base_psf_range + 2);

        let normfac = 1.0 / (2.0_f64.sqrt() * psf_sigma_pixels as f64);

        let width = self.pixel_format[0] as isize;
        let height = self.pixel_format[1] as isize;

        let reference_flux_per_s = self.mv0_flux_photons_per_s() as f64;

        for &star_idx in nearby_star_indices.iter() {
            let star = &stars.items[star_idx];

            let star_vec_cam = qcamfromj2000.transform_vector(&crate::Vector3::new(
                star.j2000_vec[0],
                star.j2000_vec[1],
                star.j2000_vec[2],
            ));

            if star_vec_cam.z <= 0.0 {
                continue;
            }

            let x_angle_tan = star_vec_cam.x / star_vec_cam.z;
            let y_angle_tan = star_vec_cam.y / star_vec_cam.z;

            let xpix = -x_angle_tan * self.focal_length_mm as f64
                / (self.pixel_pitch_um as f64 * 1.0e-3)
                + (width as f64) * 0.5;
            let ypix = -y_angle_tan * self.focal_length_mm as f64
                / (self.pixel_pitch_um as f64 * 1.0e-3)
                + (height as f64) * 0.5;

            if !(0.0..(width as f64)).contains(&xpix) || !(0.0..(height as f64)).contains(&ypix) {
                continue;
            }

            let star_signal = mv0_signal as f64 * 10f64.powf(-0.4 * star.v_mag);

            let brightness_reference = reference_flux_per_s.max(f64::EPSILON);
            let brightness_multiplier = (star_signal / brightness_reference).max(1.0);
            let extra_range = brightness_multiplier.log2().max(0.0).ceil() as isize;
            let psf_range = (base_psf_range + extra_range).min(max_psf_range);

            let xpix_int = xpix.round() as isize;
            let ypix_int = ypix.round() as isize;

            for dy in -psf_range..=psf_range {
                for dx in -psf_range..=psf_range {
                    let x0 = (xpix - (xpix_int + dx) as f64) * normfac;
                    let x1 = (xpix - (xpix_int + dx + 1) as f64) * normfac;
                    let y0 = (ypix - (ypix_int + dy) as f64) * normfac;
                    let y1 = (ypix - (ypix_int + dy + 1) as f64) * normfac;

                    //let erf_x = erf(x1) - erf(x0);
                    //let erf_y = erf(y1) - erf(y0);
                    let erf_x = erf_diff_lut(x1, x0);
                    let erf_y = erf_diff_lut(y1, y0);
                    let psf_fraction = 0.25 * erf_x * erf_y;

                    let img_x = xpix_int + dx;
                    let img_y = ypix_int + dy;
                    if img_x >= 0 && img_x < width && img_y >= 0 && img_y < height {
                        buf[(img_x as u32, img_y as u32)][0] += star_signal * psf_fraction;
                    }
                }
            }
        }
    }

    /// Render a simulated star image for the supplied pointing and integration time
    pub fn render_star_image(
        &self,
        pointing: &impl PointingProvider,
    ) -> anyhow::Result<crate::Image<f64>> {
        if !self.has_star_catalog() {
            anyhow::bail!("Star catalog has not been loaded");
        }

        let mut image =
            crate::Image::<f64>::new(self.pixel_format[0] as u32, self.pixel_format[1] as u32);
        self.add_stars_to_image(&mut image, pointing);
        Ok(image)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[ignore]
    fn test_load_catalog() {
        let mut camera = StarCamera::default();
        camera
            .load_star_catalog_from_file("data/tycho2_catalog.rkyv")
            .expect("loading star catalog");
        assert!(
            camera.star_catalog.is_some(),
            "star catalog should be loaded"
        );
        println!(
            "Loaded star catalog with {} stars",
            camera.star_catalog.as_ref().unwrap().items.len()
        );
    }

    #[test]
    fn test_betelgeuse() {
        let ra_deg = 88.792937_f64;
        let dec_deg = 7.407064_f64;
        let ra_rad = ra_deg.to_radians();
        let dec_rad = dec_deg.to_radians();
        let pointing = crate::camera::pointing::FixedInertialPointing::from_ra_dec_roll_rad(
            ra_rad, dec_rad, 0.0,
        );
        let boresight = pointing.boresight_j2000();
        println!("Boresight J2000 vector: {:?}", boresight);

        let mut camera = StarCamera::default();
        camera
            .load_star_catalog_from_file("data/tycho2_catalog.rkyv")
            .expect("loading star catalog");

        let cat = camera.star_catalog.as_ref().unwrap();
        let nearby_stars = cat.radius_search(
            [boresight[0], boresight[1], boresight[2]],
            5.0 * std::f64::consts::PI / 180.0 * 2.0_f64.sqrt(),
        );
        // Sort stars by magnitude (smallest magnitude first)
        let mut nearby_stars_sorted: Vec<_> =
            nearby_stars.iter().map(|&idx| &cat.items[idx]).collect();
        nearby_stars_sorted.sort_by(|a, b| a.v_mag.partial_cmp(&b.v_mag).unwrap());
        assert!(
            nearby_stars_sorted.len() > 0,
            "should find at least one nearby star"
        );
        // Print the 10 brightest stars: right ascensin (deg), declination (deg), V mag
        for star in nearby_stars_sorted.iter().take(10) {
            let ra = star.j2000_vec[1].atan2(star.j2000_vec[0]).to_degrees();
            let dec = star.j2000_vec[2].asin().to_degrees();
            println!(
                "Star: RA: {:.6} deg, Dec: {:.6} deg, V mag: {:.2}",
                ra, dec, star.v_mag
            );
        }
    }

    #[test]
    #[ignore]
    fn test_render_image() {
        let mut camera = StarCamera::default();
        camera
            .load_star_catalog_from_file("data/tycho2_catalog.rkyv")
            .expect("loading star catalog");
        camera.ensquared = 0.05;
        let pointing =
            crate::camera::pointing::FixedInertialPointing::from_ra_dec_roll_rad(0.0, 0.0, 0.0);
        let _image = camera
            .render_star_image(&pointing)
            .expect("rendering star image");
    }
}
