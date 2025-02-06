pub mod am0;
mod cameraframe;
mod camstate;
pub mod consts;
mod poisson;
mod starcatalog;
mod starscene;

// Stars
pub use starcatalog::the_tycho2_catalog;
pub use starcatalog::Star;
pub use starcatalog::StarCatalog;

// Math Types
pub type Quaternion = nalgebra::Unit<nalgebra::Quaternion<f64>>;
pub type Vector3 = nalgebra::Vector3<f64>;

// Random
pub use poisson::poissrand;

// Time type
pub type Instant = satkit::Instant;

// Camera State Trait
pub use camstate::CamState;
pub use camstate::CamStateSiderealTrack;
pub use camstate::CamStateStaticGCRF;
pub use camstate::CamStateStaticITRF;

pub use cameraframe::CameraFrame;

// The Scene
pub use starscene::StarScene;

// Error type
pub type SSResult<T> = Result<T, Box<dyn std::error::Error + Send + Sync>>;

// Python bindings?
#[cfg(feature = "pybindings")]
mod pybindings;
