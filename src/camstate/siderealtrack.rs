use super::*;

pub struct CamStateSiderealTrack {
    pub pos: satkit::ITRFCoord,
    pub qgcrf2cam: Quaternion,
}

impl CamStateSiderealTrack {
    pub fn new(pos: satkit::ITRFCoord, qgcrf2cam: Quaternion) -> Self {
        Self { pos, qgcrf2cam }
    }

    pub fn from_azelroll(pos: satkit::ITRFCoord, azimuth: f64, elevation: f64, roll: f64) -> Self {
        let qcart2cam = Quaternion::from_axis_angle(&Vector3::y_axis(), std::f64::consts::PI / 2.0);
        let qned2cam =
            qcart2cam * Quaternion::from_euler_angles(roll, elevation, azimuth).conjugate();
        Self::new(pos, qned2cam)
    }
}

impl CamState for CamStateSiderealTrack {
    fn camstate(&self, time: &crate::Instant) -> (Vector3, Quaternion) {
        let qgcrf2itrf = satkit::frametransform::qgcrf2itrf(time);
        (qgcrf2itrf.conjugate() * self.pos.itrf, self.qgcrf2cam)
    }
}
