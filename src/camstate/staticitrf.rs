use super::*;

pub struct CamStateStaticITRF {
    pub pos: satkit::ITRFCoord,
    pub qned2cam: Quaternion,
}

impl CamStateStaticITRF {
    pub fn new(pos: satkit::ITRFCoord, qned2cam: Quaternion) -> Self {
        Self { pos, qned2cam }
    }

    pub fn from_azelroll(pos: satkit::ITRFCoord, azimuth: f64, elevation: f64, roll: f64) -> Self {
        let qcart2cam = Quaternion::from_axis_angle(&Vector3::y_axis(), std::f64::consts::PI / 2.0);
        let qned2cam =
            qcart2cam * Quaternion::from_euler_angles(roll, elevation, azimuth).conjugate();
        Self::new(pos, qned2cam)
    }
}

impl CamState for CamStateStaticITRF {
    fn camstate(&self, time: &crate::Instant) -> (Vector3, Quaternion) {
        let qgcrf2itrf = satkit::frametransform::qgcrf2itrf(time);
        (
            qgcrf2itrf.conjugate() * self.pos.itrf,
            self.qned2cam * self.pos.q_ned2itrf().conjugate() * qgcrf2itrf,
        )
    }
}
