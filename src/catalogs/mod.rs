pub mod gaia;
pub mod kdtree;
pub mod tycho2;

use rkyv::{Archive, Deserialize, Serialize};

/// A simple representation of a star in a catalog
/// with proper motion already accounted for in
/// the RA and Dec provided
#[derive(Archive, Serialize, Deserialize, Debug, Clone)]
pub struct Star {
    pub j2000_vec: [f64; 3],
    pub v_mag: f64,
}
impl kdtree::KdPoint<3> for Star {
    fn point(&self) -> [f64; 3] {
        self.j2000_vec
    }
}

pub type StarKdTree = kdtree::KdTree<Star, 3>;
pub type StarKdNode = kdtree::KdNode<3>;
