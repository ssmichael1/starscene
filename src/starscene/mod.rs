use crate::CamState;
use crate::CameraFrame;
use crate::SSResult;

use satkit as sk;

use rand::thread_rng;
use rand_distr::{Distribution, Normal};

use std::f64::consts::PI;

/// The apparent magnitude of the sun
const SOLAR_MV: f64 = -26.74;

// Convert arcseconds to radians
const ARCSEC: f64 = std::f64::consts::PI / 648000.0;

#[derive(Debug)]
pub struct StarScene {
    /// Pixel format of the camera  Rows, then columns
    pub pixel_format: [usize; 2],
    /// Bit depth of the camera
    pub bit_depth: u8,
    /// Well depth of a pixel, electrons
    pub well_depth: f64,
    /// Read noise, electrons
    pub read_noise: f64,
    /// Dark current, electrons per second per pixel
    pub dark_current: f64,
    /// Gain, electrons per ADU
    pub gain: f64,
    /// Exposure time, seconds
    pub exposure_time: f64,
    /// Pedestal value, electrons
    pub pedestal: f64,
    /// Focal length, meters
    pub focal_length: f64,
    /// Pixel pitch, meters
    pub pixel_pitch: f64,
    /// Aperture diameter, meters
    pub diameter: f64,
    /// Aperture obscuration, meters
    pub obscuration: f64,
    /// Solar-weighted quantum efficiency
    pub solar_weighted_qe: f64,
    /// 1-sigma of point-spread function, pixels
    pub psf_1sigma: f64,
    /// zodiacal light, mV / arcsec^2
    pub zodiacal_light: f64,
    /// Time of center of integration
    pub time: sk::Instant,

    pub state: Box<dyn CamState>,
}

impl Default for StarScene {
    fn default() -> Self {
        StarScene {
            pixel_format: [512, 512],
            bit_depth: 12,
            well_depth: 10000.0,
            read_noise: 1.0,
            dark_current: 0.0,
            gain: 1.0,
            exposure_time: 10.0,
            pedestal: 50.0,
            focal_length: 0.05,
            pixel_pitch: 10.0e-6,
            diameter: 0.05,
            obscuration: 0.0,
            solar_weighted_qe: 1.0,
            psf_1sigma: 2.0,
            zodiacal_light: 22.0,
            time: sk::Instant::J2000,
            state: Box::new(crate::CamStateStaticITRF::from_azelroll(
                sk::ITRFCoord::from_geodetic_deg(0.0, 0.0, 0.0),
                0.0,
                std::f64::consts::PI / 4.0,
                0.0,
            )),
        }
    }
}

/// Cumulative distribution function for the standard normal distribution
/// # Arguments
/// * `x` - The value to evaluate the CDF at
/// # Returns
/// * The value of the CDF at `x`
///
/// # Notes
/// Approximation from:
/// https://m-hikari.com/ams/ams-2014/ams-85-88-2014/epureAMS85-88-2014.pdf
///
fn normcdf(x: f64) -> f64 {
    // m = -sqrt(pi/8)
    const M: f64 = -0.6266570686577501;
    0.5 + 0.5 * (1.0 - (M * x.powi(2)).exp()).sqrt() * f64::signum(x)
}

/// Inverse of the standard normal cumulative distribution function
/// # Arguments
/// * `a` - The value to evaluate the inverse CDF at
///
/// # Returns
/// * The value of the inverse CDF at `a`
///
/// # Notes
/// Approximation from:
/// https://m-hikari.com/ams/ams-2014/ams-85-88-2014/epureAMS85-88-2014.pdf
fn normcdfinv(a: f64) -> f64 {
    // m = -sqrt(pi/8)
    const M: f64 = -0.6266570686577501;
    (1.0 / M * (1.0 - 4.0 * (a - 0.5).powi(2)).ln()).sqrt()
}

impl StarScene {
    /// Get the field of view of the camera in radians.
    ///
    /// # Returns
    /// * A 2-element array with the field of view in radians.
    ///
    pub fn fov(&self) -> [f64; 2] {
        self.pixel_format
            .map(|x| (x as f64 * self.pixel_pitch / self.focal_length / 2.0).atan() * 2.0)
    }

    /// Get the solar flux in photons per second per square meter
    /// over the entire solar spectrum
    pub fn solar_flux() -> f64 {
        #[allow(non_upper_case_globals)]
        const c: f64 = 299792458.0;
        #[allow(non_upper_case_globals)]
        const planck: f64 = 6.62607015e-34;

        static INSTANCE: once_cell::sync::OnceCell<f64> = once_cell::sync::OnceCell::new();
        *INSTANCE.get_or_init(|| {
            crate::am0::AM0
                .iter()
                .skip(1)
                .zip(crate::am0::AM0.iter())
                .fold(0.0, |acc, (&cur, &prev)| {
                    // delta wavelength in nanometers
                    let dnm = cur[0] - prev[0];
                    // Wavelength in meters
                    let lambda = (cur[0] + prev[0]) / 2.0 * 1.0e-9;
                    // Flux in W/m^2/nm
                    let flux = (cur[1] + prev[1]) / 2.0;
                    // Photon energy in Joules
                    let energy = planck * c / lambda;
                    // Photon flux in photons/m^2/s/nm
                    let photon_flux = flux / energy;
                    acc + photon_flux * dnm
                })
        })
    }

    /// Get the field of view of the camera in degrees.
    /// # Returns
    /// * A 2-element array with the field of view in degrees.
    ///
    pub fn fov_deg(&self) -> [f64; 2] {
        self.fov().map(|x| x.to_degrees())
    }

    /// Get the instantaneous field of view of the camera in radians.
    ///
    /// # Returns
    /// * The instantaneous field of view in radians.
    ///
    pub fn ifov(&self) -> f64 {
        (self.pixel_pitch / self.focal_length).atan()
    }

    /// Map a right ascension and declination to pixel space
    ///
    /// # Arguments
    /// * `ra` - The right ascension in radians.
    /// * `dec` - The declination in radians.
    ///
    /// # Returns
    /// * The pixel coordinates of the given right ascension and declination,
    ///   or None if the point is not in the field of view.
    ///
    pub fn radec_to_pixel(
        &self,
        ra: f64,
        dec: f64,
        barycentric_velocity: Option<crate::Vector3>,
    ) -> Option<[f64; 2]> {
        let (_pos, q) = self.state.camstate(&self.time);
        let mut vec = crate::Vector3::new(ra.cos() * dec.cos(), ra.sin() * dec.cos(), dec.sin());

        // If barycentric velocity is passed in, account for stellar aberration
        if let Some(bv) = barycentric_velocity {
            let vscaled = vec + bv / satkit::consts::C;
            vec = vscaled / vscaled.norm();
        }

        let vcam = q * vec;
        let xangle = vcam.x.atan2(vcam.z);
        let yangle = vcam.y.atan2(vcam.z);
        let xpix =
            xangle.tan() * self.focal_length / self.pixel_pitch + self.pixel_format[1] as f64 / 2.0;
        let ypix =
            yangle.tan() * self.focal_length / self.pixel_pitch + self.pixel_format[0] as f64 / 2.0;
        if xpix < 0.0
            || xpix >= self.pixel_format[1] as f64
            || ypix < 0.0
            || ypix >= self.pixel_format[0] as f64
        {
            None
        } else {
            Some([xpix, ypix])
        }
    }

    /// Return the camera boresight unit vector in inertial space.
    ///
    /// # Returns
    /// * The boresight vector in inertial space.
    ///
    pub fn boresight_vec(&self) -> crate::Vector3 {
        let (_pos, q) = self.state.camstate(&self.time);
        // This vector points out along the optical axis of the camera
        let vec = crate::Vector3::new(0.0, 0.0, -1.0);
        q.conjugate() * vec
    }

    pub fn boresight_radec(&self) -> (f64, f64) {
        let vec = self.boresight_vec();
        let ra = vec.y.atan2(vec.x);
        let dec = vec.z.atan2((vec.x.powi(2) + vec.y.powi(2)).sqrt());
        (ra, dec)
    }

    /// Get the f-number of the camera.
    ///
    /// # Returns
    /// * The f-number of the camera.
    ///
    pub fn fnum(&self) -> f64 {
        self.focal_length / self.diameter
    }

    pub fn frame(&self) -> SSResult<CameraFrame<u16>> {
        let mv0_signal = StarScene::solar_flux()
            * 10.0_f64.powf(0.4 * SOLAR_MV)
            * self.solar_weighted_qe
            * self.exposure_time
            * (self.diameter.powi(2) - self.obscuration.powi(2))
            * PI
            / 4.0;

        let (ra, dec) = self.boresight_radec();
        let fov = self.fov();
        let radius = (fov[0].max(fov[1]) / 2.0) * 2.0f64.sqrt();
        let stars = crate::the_tycho2_catalog().query(ra, dec, radius);

        let readnoise = Normal::new(0.0, self.read_noise)?;

        let mut starframe = CameraFrame::<f64>::new(self.pixel_format[0], self.pixel_format[1]);

        // Add stars to frame
        let (_psun, vsun) =
            satkit::jplephem::geocentric_state(satkit::SolarSystem::Sun, &self.time)?;
        let barycentric_velocity = -vsun;

        stars.iter().for_each(|star| {
            let (ra, dec) = (star.ra.to_radians(), star.dec.to_radians());

            // Convert ra, dec to pixel coordinates, accounting for stellar aberration
            // and handle pixel if it is in the field of view
            if let Some(pix) = self.radec_to_pixel(ra, dec, Some(barycentric_velocity)) {
                let cenrow = pix[1];
                let cencol = pix[0];
                let irow = cenrow.floor() as i32;
                let icol = cencol.floor() as i32;

                // Scale the signal by the magnitude of the star
                let total_signal = mv0_signal * 10.0_f64.powf(-0.4 * star.mv() as f64);

                // Capture 99.5% of the total signal, or all but 5 e-, whichever is greater
                let frac = (1.0 - 5.0 / total_signal).max(0.995);

                // How many pixels out to go to capture the signal
                // This is conservative, but it's better to overestimate the PSF than underestimate it
                let radius = (self.psf_1sigma * normcdfinv(frac)).ceil() as i32;

                for row in (irow - radius)..=(irow + radius) {
                    for col in (icol - radius)..=(icol + radius) {
                        if row >= 0
                            && row < starframe.rows() as i32
                            && col >= 0
                            && col < starframe.cols() as i32
                        {
                            // Integrate the gaussian PDF over the pixel by taking
                            // the difference of the CDF at the edges of the pixel
                            // Since the PSF is gaussian, x and y are separable (yay!)
                            starframe[(row as usize, col as usize)] += total_signal
                                * (normcdf((row as f64 + 0.5 - cenrow) / self.psf_1sigma)
                                    - normcdf((row as f64 - 0.5 - cenrow) / self.psf_1sigma))
                                * (normcdf((col as f64 + 0.5 - cencol) / self.psf_1sigma)
                                    - normcdf((col as f64 - 0.5 - cencol) / self.psf_1sigma));
                        }
                    }
                }
            } // end of iterating over star
        });

        // Create the frame
        let mut frame = CameraFrame::<u16>::new(self.pixel_format[0], self.pixel_format[1]);

        // Random number generator
        let mut rng = thread_rng();

        // Expected dark current
        let idark = self.dark_current * self.exposure_time;
        // expected zodiacal light
        let zodiacal_light =
            mv0_signal * 10.0_f64.powf(-0.4 * self.zodiacal_light) * (self.ifov() / ARCSEC).powi(2);

        // Total background
        let background = idark + zodiacal_light;

        frame
            .data
            .iter_mut()
            .zip(starframe.data.iter())
            .for_each(|(pix, &ph)| {
                let mut em: f64 = 0.0;
                em += readnoise.sample(&mut rng);
                em += self.pedestal;

                *pix = ((em + crate::poissrand(&mut rng, ph + background) as f64)
                    .clamp(0.0, self.well_depth)
                    / self.gain)
                    .clamp(0.0, 65535.0) as u16;
            });

        Ok(frame)
    }
}

impl std::fmt::Display for StarScene {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "StarScene")?;
        writeln!(f, "  pixel_format: {:?}", self.pixel_format)?;
        writeln!(
            f,
            "  field of view: {:0.2} X {:0.2} deg",
            self.fov_deg()[0],
            self.fov_deg()[1]
        )?;
        writeln!(f, "  bit_depth: {}", self.bit_depth)?;
        writeln!(f, "  read_noise: {} e-", self.read_noise)?;
        writeln!(f, "  dark_current: {} e- / sec / pixel", self.dark_current)?;
        writeln!(f, "  gain: {} e- / ADU", self.gain)?;
        writeln!(f, "  exposure_time: {} s", self.exposure_time)?;
        writeln!(f, "  pedestal: {} e-", self.pedestal)?;
        writeln!(f, "  focal_length: {} m", self.focal_length)?;
        writeln!(f, "  pixel_pitch: {} µm", self.pixel_pitch * 1.0e6)?;
        writeln!(f, "  diameter: {} m", self.diameter)?;
        writeln!(f, "  obscuration: {} m", self.obscuration)?;
        writeln!(f, "  solar_weighted_qe: {}", self.solar_weighted_qe)?;
        writeln!(f, "  psf_1sigma: {}", self.psf_1sigma)?;
        writeln!(f, "  zodiacal_light: {} mV/arcsec^2", self.zodiacal_light)?;
        writeln!(f, "  time: {}", self.time)?;
        writeln!(f, "  state: {:?}", self.state)
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::Quaternion;
    use crate::Vector3;

    #[test]
    fn test_star() {
        // For star beatelgeuse
        let ra = 88.7929_f64.to_radians();
        let dec = 7.4071_f64.to_radians();
        let v = Vector3::new(ra.cos() * dec.cos(), ra.sin() * dec.cos(), dec.sin());
        let q = Quaternion::rotation_between(&v, &-Vector3::z_axis()).unwrap();

        let scene = StarScene {
            gain: 0.01,
            pixel_format: [2048, 2048],
            exposure_time: 0.1,
            psf_1sigma: 2.8,
            diameter: 0.0254,
            focal_length: 0.059,
            pixel_pitch: 3.46e-6,
            bit_depth: 16,
            solar_weighted_qe: 0.5,
            state: Box::new(crate::CamStateStaticGCRF {
                qgcrf2cam: q,
                ..crate::CamStateStaticGCRF::default()
            }),
            ..StarScene::default()
        };
        let frame = scene.frame().unwrap();
        println!("camera  = {}", scene);
        frame.save("test_frame.png").unwrap();
    }

    #[test]
    fn test_frame() {
        let mut scene = StarScene {
            gain: 10.0,
            exposure_time: 1.0,
            psf_1sigma: 0.8,
            ..StarScene::default()
        };
        let time = sk::Instant::J2000 + sk::Duration::from_hours(16.0);
        scene.time = time;

        let frame = scene.frame().unwrap();
        frame.save("test_frame.png").unwrap();
    }

    #[test]
    fn test_solar_flux() {
        let flux = StarScene::solar_flux();
        println!("flux = {}", flux);
    }

    #[test]
    fn test_normcdf() {
        let cdf = normcdf(-1.0);
        println!("cdf = {}", cdf);
    }

    #[test]
    fn test_normcdfinv() {
        let a = normcdf(1.0);
        println!("a = {}", a);
        let x = normcdfinv(0.99);
        println!("x = {}", x);
    }

    #[test]
    fn test_scene() {
        let scene = StarScene::default();
        let time = sk::Instant::now();
        println!("scene = {:?}", scene);
        let (pos, q) = scene.state.camstate(&time);
        println!("pos = {:?}", pos);
        println!("q = {:?}", q);
    }
}
