// Trait alias
pub trait Pixel: num_traits::Zero + Clone + Copy {}
impl<T: num_traits::Zero + Clone + Copy> Pixel for T {}

/// A frame of camera data
///
#[derive(Debug, Clone)]
pub struct CameraFrame<T: Pixel> {
    pub(crate) data: Vec<T>,
    rows: usize,
    cols: usize,
}

impl<T: Pixel> CameraFrame<T> {
    pub fn new(rows: usize, cols: usize) -> Self {
        Self {
            data: vec![T::zero(); rows * cols],
            rows,
            cols,
        }
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Return an iterator over the pixels at the given row
    pub fn row(&self, row: usize) -> impl Iterator<Item = T> + '_ {
        self.data[row * self.cols..(row + 1) * self.cols]
            .iter()
            .copied()
    }

    /// Return an iterator over mutable references to the pixels at the given row
    pub fn row_mut(&mut self, row: usize) -> impl Iterator<Item = &mut T> + '_ {
        let ptr = self.data.as_mut_ptr();
        let cols = self.cols;
        (0..self.cols).map(move |col| unsafe { &mut *ptr.add(row * cols + col) })
    }

    pub fn col(&self, col: usize) -> impl Iterator<Item = T> + '_ {
        (0..self.rows).map(move |row| self.data[row * self.cols + col])
    }

    pub fn col_mut(&mut self, col: usize) -> impl Iterator<Item = &mut T> + '_ {
        let ptr = self.data.as_mut_ptr();
        let cols = self.cols;
        (0..self.rows).map(move |row| unsafe { &mut *ptr.add(row * cols + col) })
    }

    pub fn get(&self, row: usize, col: usize) -> T {
        self.data[row * self.cols + col]
    }

    pub fn set(&mut self, row: usize, col: usize, val: T) {
        self.data[row * self.cols + col] = val;
    }
}

impl CameraFrame<u16> {
    /// Save frame to a PNG file
    ///
    /// # Arguments
    /// * `filename` - The name of the file to save the frame to
    ///
    /// # Returns
    /// * A Result containing the success or failure of the save operation
    pub fn save(&self, filename: &str) -> std::io::Result<()> {
        let mut encoder = png::Encoder::new(
            std::fs::File::create(filename)?,
            self.cols as u32,
            self.rows as u32,
        );
        encoder.set_color(png::ColorType::Grayscale);
        encoder.set_depth(png::BitDepth::Sixteen);
        let mut writer = encoder.write_header()?;
        writer.write_image_data(
            &self
                .data
                .iter()
                .flat_map(|x| x.to_be_bytes())
                .collect::<Vec<u8>>(),
        )?;
        Ok(())
    }
}

/// Implement indexing for CameraFrame
impl<T: Pixel> std::ops::Index<(usize, usize)> for CameraFrame<T> {
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        &self.data[index.0 * self.cols + index.1]
    }
}

/// Implement mutable indexing for CameraFrame
impl<T: Pixel> std::ops::IndexMut<(usize, usize)> for CameraFrame<T> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        &mut self.data[index.0 * self.cols + index.1]
    }
}
