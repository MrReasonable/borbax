//! The 3D scene: one molecule from the loaded universe, drawn as shaded
//! spheres.
//!
//! Every file under here is in the seam's **strictest** tier by default — no
//! `egui`, no `eframe`, no `epaint`, and no `wgpu` without an explicit entry in
//! `xtask`'s `VIEWER_GPU_FILES`. That is what keeps the scene's arithmetic
//! testable with no window and no GPU adapter.

pub mod camera;
