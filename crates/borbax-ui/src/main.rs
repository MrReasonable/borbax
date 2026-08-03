//! The binary.
//!
//! **One statement, and that is the point.** Everything the window is made of
//! lives in [`borbax_ui::app::viewer_app`], where a test can reach it — a
//! binary crate root cannot be imported, so wiring kept here is wiring nothing
//! can check. That is not hypothetical: an app missing its camera opened an
//! empty window on 2026-08-03 with all 509 tests green.

fn main() {
    borbax_ui::app::viewer_app().run();
}
