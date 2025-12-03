use crate::Quaternion;

pub trait PointingProvider {
    fn camera_from_j2000(&self) -> crate::Quaternion;

    fn boresight_j2000(&self) -> crate::Vector3 {
        let q = self.camera_from_j2000();
        q.inverse().transform_vector(&crate::Vector3::z_axis())
    }

    fn right_ascension_rad(&self) -> f64 {
        let boresight = self.boresight_j2000();
        boresight
            .y
            .atan2(boresight.x)
            .rem_euclid(2.0 * std::f64::consts::PI)
    }

    fn declination_rad(&self) -> f64 {
        let boresight = self.boresight_j2000();
        boresight.z.asin()
    }

    fn roll_rad(&self) -> f64 {
        let qxout2zout =
            Quaternion::from_axis_angle(&crate::Vector3::z_axis(), -std::f64::consts::FRAC_PI_2)
                * Quaternion::from_axis_angle(
                    &crate::Vector3::y_axis(),
                    -std::f64::consts::FRAC_PI_2,
                );

        let q = qxout2zout.conjugate() * self.camera_from_j2000();

        let euler = q.conjugate().euler_angles();
        euler.0
    }
}

pub struct FixedInertialPointing {
    fixed_cam_from_j2000: Quaternion,
}

impl PointingProvider for FixedInertialPointing {
    fn camera_from_j2000(&self) -> Quaternion {
        self.fixed_cam_from_j2000
    }
}

impl FixedInertialPointing {
    pub fn new(orientation: Quaternion) -> Self {
        Self {
            fixed_cam_from_j2000: orientation,
        }
    }

    pub fn from_ra_dec_roll_rad(ra_rad: f64, dec_rad: f64, roll_rad: f64) -> Self {
        // Constant quaternion that rotates from frame with xhat out of boresight to zhat out of boresight
        // (the camera frame)
        let qxout2zout =
            Quaternion::from_axis_angle(&crate::Vector3::z_axis(), -std::f64::consts::FRAC_PI_2)
                * Quaternion::from_axis_angle(
                    &crate::Vector3::y_axis(),
                    -std::f64::consts::FRAC_PI_2,
                );

        let q_ra = Quaternion::from_axis_angle(&crate::Vector3::z_axis(), -ra_rad);
        let q_dec = Quaternion::from_axis_angle(&crate::Vector3::y_axis(), dec_rad);
        let q_roll = Quaternion::from_axis_angle(&crate::Vector3::x_axis(), -roll_rad);
        let orientation = qxout2zout * q_roll * q_dec * q_ra;
        Self {
            fixed_cam_from_j2000: orientation,
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    #[test]
    fn test_fixed_inertial_pointing() {
        let ra = 10.0_f64.to_radians();
        let dec = 20.0_f64.to_radians();
        let roll = 50.0_f64.to_radians();
        let pointing = FixedInertialPointing::from_ra_dec_roll_rad(ra, dec, roll);
        let q = pointing.camera_from_j2000();

        assert!((pointing.right_ascension_rad() - ra).abs() < 1e-10);
        assert!((pointing.declination_rad() - dec) < 1e-10);
        assert!((pointing.roll_rad() - roll).abs() < 1e-10);

        let v = crate::Vector3::new(ra.cos() * dec.cos(), ra.sin() * dec.cos(), dec.sin());
        let v_rotated = q.transform_vector(&v);
        let angle = v_rotated.angle(&crate::Vector3::z_axis());
        assert!(angle < 1e-10);

        // Check roll via signed angle between zero-roll and roll-adjusted camera x-axis
        let boresight = q.inverse().transform_vector(&crate::Vector3::z_axis());
        let cam_x = q.inverse().transform_vector(&crate::Vector3::x_axis());

        let ref_pointing = FixedInertialPointing::from_ra_dec_roll_rad(ra, dec, 0.0);
        let ref_x = ref_pointing
            .camera_from_j2000()
            .inverse()
            .transform_vector(&crate::Vector3::x_axis());

        let project = |v: crate::Vector3| {
            let parallel = boresight * v.dot(&boresight);
            (v - parallel).normalize()
        };

        let cam_x_proj = project(cam_x);
        let ref_x_proj = project(ref_x);

        let sin_angle = boresight.dot(&ref_x_proj.cross(&cam_x_proj));
        let cos_angle = ref_x_proj.dot(&cam_x_proj);
        let roll_angle = sin_angle.atan2(cos_angle);
        assert!((roll_angle - roll).abs() < 1e-10);

        assert!((pointing.right_ascension_rad() - ra).abs() < 1e-10);
        assert!((pointing.declination_rad() - dec).abs() < 1e-10);
        assert!((pointing.roll_rad() - roll).abs() < 1e-10);
    }
}
