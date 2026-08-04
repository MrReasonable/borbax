//! The binary.
//!
//! **One statement, and that is the point.** Everything the window is made of
//! lives in [`borbax_ui::app::viewer_app`], where a test can reach it — a
//! binary crate root cannot be imported, so wiring kept here is wiring nothing
//! can check. That is not hypothetical: an app missing its camera opened an
//! empty window on 2026-08-03 with all 509 tests green.

/// Run the window, and report what it exited with.
///
/// **`impl Termination` rather than the engine's own exit type**, so this file
/// still names no engine at all. Discarding `run()`'s return — which the first
/// version did — makes the process report success however it ended, and that
/// was the one decision left in a file whose whole doc is that it holds none.
fn main() -> impl std::process::Termination {
    borbax_ui::app::viewer_app().run()
}
