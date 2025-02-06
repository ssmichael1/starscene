use super::*;

use crate::Quaternion;
use crate::Vector3;

pub struct CamStateStaticGCRF {
    pub pos: Vector3,
    pub qgcrf2cam: Quaternion,
}

impl Default for CamStateStaticGCRF {
    fn default() -> Self {
        Self {
            pos: Vector3::zeros(),
            qgcrf2cam: Quaternion::identity(),
        }
    }
}

impl CamState for CamStateStaticGCRF {
    fn camstate(&self, _time: &crate::Instant) -> (Vector3, Quaternion) {
        (self.pos, self.qgcrf2cam)
    }
}
