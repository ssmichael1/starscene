mod siderealtrack;
mod staticgcrf;
mod staticitrf;

use crate::Quaternion;
use crate::Vector3;

pub use siderealtrack::CamStateSiderealTrack;
pub use staticgcrf::CamStateStaticGCRF;
pub use staticitrf::CamStateStaticITRF;

/// Camera Coordinate Frame Definition:
/// * X axis goes to the right along camera columns
/// * Y axis goes up along the camera rows
/// * Z axis goes into the camera along the optical axis (-Zhat is out)
///
pub trait CamState: Send + Sync {
    /// Return camera state
    /// GCRF position of camera in meters,
    /// quaternion from GCRF to camera frame
    fn camstate(&self, time: &crate::Instant) -> (Vector3, Quaternion);
}

impl std::fmt::Debug for dyn CamState {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "CamState")
    }
}

#[cfg(test)]
mod tests {
    use nalgebra as na;

    #[test]
    fn test_euler() {
        let azimuth = (30.0_f64).to_radians();
        let elevation = (45.0_f64).to_radians();
        let roll: f64 = std::f64::consts::PI / 2.0;

        let ned = na::vector![
            f64::cos(azimuth) * f64::cos(elevation),
            f64::sin(azimuth) * f64::cos(elevation),
            -1.0 * f64::sin(elevation)
        ];
        let qcart2cam = na::UnitQuaternion::from_axis_angle(
            &na::Vector3::<f64>::y_axis(),
            std::f64::consts::PI / 2.0,
        );

        let q = na::UnitQuaternion::from_euler_angles(roll, elevation, azimuth);
        println!("q = {}", q);
        println!("enu = {}", ned);
        println!("q * enu = {}", qcart2cam * q.conjugate() * ned);
    }
}
