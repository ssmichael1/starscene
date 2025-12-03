use std::cmp::Ordering;

use rkyv::{Archive, Deserialize, Serialize};

/// Trait for types that can be located in N-dimensional Cartesian space.
pub trait KdPoint<const N: usize> {
    fn point(&self) -> [f64; N];
}

/// Axis-aligned KD-Tree node storing a single item index and optional children.
#[derive(Archive, Serialize, Deserialize, Debug, Clone)]
pub struct KdNode<const N: usize> {
    pub axis: u8,
    pub point_index: usize,
    pub left: Option<usize>,
    pub right: Option<usize>,
    pub bbox_min: [f64; N],
    pub bbox_max: [f64; N],
}

#[derive(Archive, Serialize, Deserialize, Debug, Clone)]
pub struct KdTree<T, const N: usize>
where
    T: KdPoint<N>,
{
    pub items: Vec<T>,
    pub nodes: Vec<KdNode<N>>,
    pub root: Option<usize>,
}

impl<T, const N: usize> KdTree<T, N>
where
    T: KdPoint<N>,
{
    pub fn build(items: Vec<T>) -> Self {
        if items.is_empty() {
            return Self {
                items,
                nodes: Vec::new(),
                root: None,
            };
        }

        assert!(N > 0, "KDTree dimension must be greater than zero");

        let mut indices: Vec<usize> = (0..items.len()).collect();
        let mut builder = TreeBuilder::new(&items);
        let root = builder.build_recursive(&mut indices, 0);
        let nodes = builder.nodes;

        Self { items, nodes, root }
    }

    /// Find the closest item to the provided query point.
    /// Returns the item index and Euclidean distance.
    pub fn nearest_neighbor(&self, query_point: [f64; N]) -> Option<(usize, f64)> {
        let root = self.root?;
        let mut best = NearestMatch {
            index: None,
            dist2: f64::INFINITY,
        };

        self.nearest_recursive(root, &query_point, &mut best);
        best.index.map(|idx| (idx, best.dist2.sqrt()))
    }

    fn nearest_recursive(&self, node_idx: usize, query_point: &[f64; N], best: &mut NearestMatch) {
        let node = &self.nodes[node_idx];
        let point = self.items[node.point_index].point();
        let dist2 = squared_distance(&point, query_point);
        if dist2 < best.dist2 {
            best.dist2 = dist2;
            best.index = Some(node.point_index);
        }

        let axis = node.axis as usize;
        let diff = query_point[axis] - point[axis];
        let (near, far) = if diff <= 0.0 {
            (node.left, node.right)
        } else {
            (node.right, node.left)
        };

        if let Some(child) = near {
            self.nearest_recursive(child, query_point, best);
        }

        if diff * diff < best.dist2 {
            if let Some(child) = far {
                self.nearest_recursive(child, query_point, best);
            }
        }
    }

    /// Collect indices of all items lying within `radius` of the query point.
    pub fn radius_search(&self, query_point: [f64; N], radius: f64) -> Vec<usize> {
        if radius < 0.0 {
            return Vec::new();
        }

        let root = match self.root {
            Some(idx) => idx,
            None => return Vec::new(),
        };

        let radius2 = radius * radius;
        let mut hits = Vec::new();
        self.radius_recursive(root, &query_point, radius2, &mut hits);
        hits
    }

    fn radius_recursive(
        &self,
        node_idx: usize,
        query_point: &[f64; N],
        radius2: f64,
        hits: &mut Vec<usize>,
    ) {
        let node = &self.nodes[node_idx];
        let point = self.items[node.point_index].point();
        let dist2 = squared_distance(&point, query_point);
        if dist2 <= radius2 {
            hits.push(node.point_index);
        }

        let axis = node.axis as usize;
        let diff = query_point[axis] - point[axis];

        let (near, far) = if diff <= 0.0 {
            (node.left, node.right)
        } else {
            (node.right, node.left)
        };

        if let Some(child) = near {
            self.radius_recursive(child, query_point, radius2, hits);
        }

        if diff * diff <= radius2 {
            if let Some(child) = far {
                self.radius_recursive(child, query_point, radius2, hits);
            }
        }
    }
}

struct NearestMatch {
    index: Option<usize>,
    dist2: f64,
}

fn squared_distance<const N: usize>(a: &[f64; N], b: &[f64; N]) -> f64 {
    let mut sum = 0.0;
    for axis in 0..N {
        let diff = a[axis] - b[axis];
        sum += diff * diff;
    }
    sum
}

struct TreeBuilder<'a, T, const N: usize>
where
    T: KdPoint<N>,
{
    items: &'a [T],
    nodes: Vec<KdNode<N>>,
}

impl<'a, T, const N: usize> TreeBuilder<'a, T, N>
where
    T: KdPoint<N>,
{
    fn new(items: &'a [T]) -> Self {
        Self {
            items,
            nodes: Vec::with_capacity(items.len()),
        }
    }

    fn build_recursive(&mut self, indices: &mut [usize], depth: usize) -> Option<usize> {
        if indices.is_empty() {
            return None;
        }

        let axis = (depth % N) as u8;
        indices.sort_unstable_by(|a, b| {
            let lhs = self.items[*a].point()[axis as usize];
            let rhs = self.items[*b].point()[axis as usize];
            lhs.partial_cmp(&rhs).unwrap_or(Ordering::Equal)
        });

        let (bbox_min, bbox_max) = self.bounds(indices);

        let median = indices.len() / 2;
        let (left_slice, rest) = indices.split_at_mut(median);
        let (median_idx_slice, right_slice) =
            rest.split_first_mut().expect("non-empty slice after split");
        let point_index = *median_idx_slice;

        let node_index = self.nodes.len();
        self.nodes.push(KdNode {
            axis,
            point_index,
            left: None,
            right: None,
            bbox_min,
            bbox_max,
        });

        let left = self.build_recursive(left_slice, depth + 1);
        let right = self.build_recursive(right_slice, depth + 1);
        self.nodes[node_index].left = left;
        self.nodes[node_index].right = right;
        Some(node_index)
    }

    fn bounds(&self, indices: &[usize]) -> ([f64; N], [f64; N]) {
        let mut min = [f64::INFINITY; N];
        let mut max = [f64::NEG_INFINITY; N];

        for &idx in indices {
            let point = self.items[idx].point();
            for axis in 0..N {
                min[axis] = min[axis].min(point[axis]);
                max[axis] = max[axis].max(point[axis]);
            }
        }

        (min, max)
    }
}
