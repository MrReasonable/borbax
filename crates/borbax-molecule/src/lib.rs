//! Molecule graphs, canonical form, geometry and shape signatures (spec §8).
//!
//! This crate carries the conceptual weight of Borbax: a molecule becomes a
//! graph, the graph becomes a canonical form, the canonical form becomes a 3D
//! embedding, and the embedding becomes a shape signature that
//! [`borbax_universe`]'s chemistry can score against another one. It is split
//! into small files because that is where nearly all iteration happens.
//!
//! **Everything here is per-species.** Spec §8.6: canonicalisation, embedding,
//! signature construction, folding and cavity extraction are pure functions of
//! the species and are computed once, at intern time. No code reachable from a
//! simulation step may call any of them.

pub mod canonical;
pub mod geodesic;
pub mod graph;
pub mod layout;

pub use canonical::{CanonForm, CanonMol, Capped, SEARCH_LEAF_CAP, SearchStats, canonicalise};
pub use geodesic::{
    GeoError, Geodesic, Mat3, N_ROTATIONS, Rotation, Vec3, apply_mat, is_identity,
    rotation_matrices, vertex_count,
};
pub use graph::{BondError, MAX_ATOMS, Mol12, N_ORDERS};
pub use layout::{Embedding, ITERATIONS, embed, raw_stress, stress};
