//! The Borbax viewer — the window a universe is looked at through.
//!
//! # The rule this crate exists under
//!
//! **The viewer computes no physics.** It calls into the `borbax-*` crates and
//! displays what comes back; it must never reimplement a calculation "just for
//! display", because the moment it does, the picture and the simulation can
//! disagree and the picture is the thing being trusted.
//!
//! The discriminator, stated so a reviewer can check it rather than believe it:
//! **every number on screen must be traceable to a call into a `borbax-*`
//! crate.** Today that is exactly three — [`borbax_universe::Universe::generate`],
//! `Universe::table().len()` and `u8::from(universe.physics)` — and a grep of
//! this crate for arithmetic should find only layout.
//!
//! # The seam, which is per-file and greppable
//!
//! Three files, and the rule is a property of each file's imports rather than a
//! convention anyone has to remember:
//!
//! | File | May import | Holds |
//! |---|---|---|
//! | [`state`] | no UI type at all | every decision, and every test that can reach one |
//! | [`panel`] | `egui` | the drawing, and no `format!` |
//! | `main.rs` | `eframe` | the window, and nothing else |
//!
//! `egui` is a pure CPU layout library and cannot open a window, which is what
//! lets `egui_kittest` drive [`panel::draw`] headlessly in CI. `eframe` is the
//! shell that can, so it is confined to the one file no test reaches.
//!
//! Formatting lives in [`state::ViewerState::status_line`] rather than in
//! [`panel`] for a reason with a failure attached: a test that asserts on the
//! state while the painted string is built inline in the panel stays green over
//! a window it has stopped describing.
//!
//! # What §13.1 does and does not ask of this crate
//!
//! The bans in §13.1 exist to protect simulation *results*, and this crate
//! produces none — it is a pure consumer. So `f32`, `HashMap`, wall-clock time
//! and unordered iteration are all permissible **inside its own rendering**.
//! What is not permissible is feeding any of them back into simulation state,
//! or letting a display choice reach a `borbax-*` call.
#![expect(
    clippy::multiple_crate_versions,
    reason = "`eframe`'s dependency tree carries 31 duplicated transitive \
              dependencies — `windows-sys` at four versions, the nine \
              `windows_*` target shims and `windows-targets` at two each, \
              `objc2` plus five `objc2-*` crates at two each, `hashbrown` and \
              `redox_syscall` at three, `rustix`/`linux-raw-sys`/`calloop`/\
              `smithay-client-toolkit` straddling a major bump, and `thiserror` \
              1.0 alongside our own 2.0. This workspace chooses none of them \
              and can unify none of them. Scoped to this crate rather than \
              allowed at workspace level so the lint keeps working on the four \
              physics crates, where a duplicated runtime dependency is the \
              expensive case the dependency doctrine is about — allowing it \
              workspace-wide would trade an enforced invariant across those for \
              a suppression in the one crate that cannot affect results. \
              \
              `#[expect]` rather than `#[allow]` so the suppression cannot rot: \
              it fails the build the day the tree unifies. That is not \
              hypothetical here — the count and the placement were BOTH first \
              recorded wrongly (~34, and `required in both crate roots`), and \
              it was an unfulfilled expectation on `main.rs` that caught it. \
              Measured on this tree: the lint is evaluated once per *package*, \
              so this single attribute covers the library and the binary, and a \
              second one in `main.rs` fails the gate as unfulfilled. The count \
              is re-measured on any `eframe` bump, never copied forward"
)]

pub mod panel;
pub mod state;

/// The window's title, and the application's permanent name.
///
/// **Pinned by a test, and the first value was a review decision rather than a
/// test outcome** — `the_window_title_names_the_program_and_not_its_subject`
/// has no power against whatever string is written here first, only against the
/// second. So the rule it enforces is stated rather than implied: **the title
/// may name the program; it may not name the program's subject.** An invented
/// proper noun asserts nothing. Every descriptive alternative — "Chemistry
/// Simulator", "Molecular Explorer", "Atomic Sandbox" — names a real-world
/// subject matter, which is the claim G6 (§5) forbids.
///
/// It is also deliberately **constant**, and an earlier draft that interpolated
/// the seed would have defeated the guard being built: a title that varies
/// cannot be pinned by an equality assertion.
pub const WINDOW_TITLE: &str = "Borbax";

/// The universe the window opens on.
///
/// **An arbitrary constant, and deliberately a boring one.** Any seed would do
/// — that is the point of the whole exercise — so this one is chosen to carry
/// no meaning at all. A seed picked because its universe "looks good" would be
/// a thumb on the scale, and the first thing anyone does with this program is
/// change it.
///
/// It exists because the window should open on a universe rather than on the
/// words "type a seed". Step 1b replaces the number with a name.
pub const OPENING_SEED: u64 = 1;
