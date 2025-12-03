mod camera;
pub mod catalogs;

pub use camera::StarCamera;
pub use catalogs::Star;
pub use catalogs::StarKdTree;

#[cfg(feature = "python")]
mod pybindings;

/// Common types used in the library
pub type Quaternion = nalgebra::Unit<nalgebra::Quaternion<f64>>;
pub type Vector3 = nalgebra::Vector3<f64>;
pub type Matrix3 = nalgebra::Matrix3<f64>;

pub type Image<T> = image::ImageBuffer<image::Luma<T>, Vec<T>>;
