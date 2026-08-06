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
//! crate.** The table is the claim; there is deliberately no count in front of
//! it, because a count has gone stale twice now — "exactly three", then
//! "exactly four" while the table had five rows — each time in the commit that
//! widened the table. A reviewer also pointed out the number was measuring the
//! wrong thing: the rows are *provenances*, and three of them are the same
//! field (the seed) arriving three ways.
//!
//! What is on screen, and where each value comes from:
//!
//! | On screen | Comes from |
//! |---|---|
//! | the element count | `universe.table.len()` — `table` is a public **field**, not a method |
//! | the physics version | `u8::from(universe.physics)` |
//! | the seed, when typed | the user |
//! | the seed, after "surprise me" | `borbax_rng::Stream::next_range` |
//! | the seed, after a name is typed | [`borbax_universe::seed_from_phrase`] |
//!
//! and the universe behind the first two comes from
//! [`borbax_universe::Universe::generate`], which has **exactly one call site**
//! in this crate — checked by `cargo xtask`, not merely asserted here.
//!
//! An early version of this list also named
//! `Universe::table()` (a method that does not exist), omitted the seed, and
//! claimed a grep for arithmetic finds only layout. **The arithmetic is not
//! zero**: `state.rs` has a `wrapping_mul`/`wrapping_add` pair (the clock
//! reading, an entropy source where no value is wrong) and a `saturating_add`
//! (the regeneration counter). Neither reaches the screen. The rule that
//! actually holds is narrower and worth stating properly — **no arithmetic on a
//! value that came out of a `borbax-*` call**, because that is the viewer
//! recomputing physics. There is none.
//!
//! # The seam, which is per-file and checked
//!
//! Three files, and the rule is a property of each file's imports rather than a
//! convention anyone has to remember:
//!
//! | File | May name in code | Holds |
//! |---|---|---|
//! | [`state`], [`molecule`], [`orbit`] | no engine or UI type at all | every decision, the chemistry it asks for, the camera's arithmetic — and every test that can reach them with no app at all |
//! | [`panel`] | `egui` (via `bevy_egui`) | the drawing, and no `format!` |
//! | [`app`], [`scene`] | `bevy` and `egui` | the window, the entities and the frame loop — and no `format!` either |
//! | `main.rs` | nothing | one statement |
//!
//! **`main.rs` is on that table with an empty licence, and that is the point.**
//! It held the engine licence until a review removed it and measured the gate
//! still green. A licence nobody uses is permission for the next person to put a
//! decision back into the one file no test can import, which is how an empty
//! window shipped with 509 tests passing.
//!
//! **"In code" is load-bearing and a plain `grep` will not tell you this.** The
//! prose necessarily quotes what it forbids — this paragraph names `bevy`, and
//! [`panel`]'s module doc says "there is no `format!` in this file" — so a raw
//! `grep -rl bevy src/` reports **four** files against the two that name it in
//! code, and a `grep format! panel.rs` reports one, all correctly and all
//! uselessly. (The count is deliberately stated: an earlier version said three,
//! which was the number before this file's own table grew a row.) `cargo xtask` strips
//! whole-line comments before matching, which is the check that means anything.
//!
//! **`bevy` is matched as a whole identifier and `egui` as a substring**, and
//! the asymmetry is deliberate. `bevy_egui` is what the drawing tier is
//! *allowed* — it is where `egui` comes from now — so a substring `bevy` ban
//! would fire on the one import the tier exists to permit. A substring `egui`
//! ban is what catches `bevy_egui` in a file with no drawing licence.
//!
//! `egui` is a pure CPU layout library and cannot open a window, which is what
//! lets `egui_kittest` drive [`panel::draw`] headlessly in CI — measured to
//! survive the move to Bevy completely, because neither it nor [`panel`] ever
//! named `eframe`. Bevy is the shell that opens a window and owns the frame
//! loop, so it is confined to the one file no test reaches.
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
    reason = "the engine's dependency tree emits these — `windows-sys` at \
              four versions, the eight `windows_*` target shims and \
              `windows-targets` at two each, `objc2` plus five `objc2-*` crates \
              at two each, `hashbrown` and `redox_syscall` at three, \
              `rustix`/`linux-raw-sys`/`calloop`/`smithay-client-toolkit` each \
              straddling a major bump, `thiserror` 1.0 alongside our own 2.0, \
              and `bitflags`/`block2`/`calloop-wayland-source`/`foldhash`/\
              `jni-sys`/`rustc-hash`/`syn`/`thiserror-impl`. This workspace \
              chooses none of them and can unify none of them. Scoped to this \
              crate rather than allowed at workspace level so the lint keeps \
              working on the four physics crates, where a duplicated runtime \
              dependency is the expensive case the dependency doctrine is about \
              — allowing it workspace-wide would trade an enforced invariant \
              across those for a suppression in the one crate that cannot \
              affect results. \
              \
              `#[expect]` rather than `#[allow]` so the suppression cannot rot: \
              it fails the build the day the tree unifies. That has already paid \
              twice here. The count was first recorded as ~34, and the \
              placement as `required in both crate roots` — an unfulfilled \
              expectation on `main.rs` refuted the second. The explanation then \
              written for THAT was also wrong: `once per package`. It is not. \
              This lint is a per-crate pass, and `lib` and `bin` are separate \
              crates; the binary escapes because clippy resolves the local \
              package by matching the crate name against `cargo metadata`'s \
              package names, and this binary is called `borbax`, which matches \
              no package. Verified by renaming `[[bin]] name` to `borbax-ui`, \
              at which point the bin lints too and emits the same 31. So if the \
              binary is ever renamed to match the package, or a second bin or \
              example is added whose crate name does, this attribute stops \
              covering it and 31 errors arrive under `-D warnings`. The count \
              is re-measured on any engine bump, never copied forward. The \
              count and the crate list below were measured against `eframe` and \
              are stale for Bevy; `#[expect]` is what makes that a build \
              failure to be re-measured rather than a comment to be believed"
)]

pub mod app;
pub mod molecule;
pub mod orbit;
pub mod palette;
pub mod panel;
pub mod scene;
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
/// words "type a seed".
///
/// **Step 1b did not replace it, and an earlier version of this sentence said
/// it would.** The name box opens *empty*, hinted rather than pre-filled, and
/// this constant is unchanged. Pre-filling a name would mean choosing the
/// universe that name produces — a thumb on the scale in the one place the
/// whole program is about the seed being arbitrary — while an empty box with
/// no hint teaches nothing about what the box is for. The hint does both jobs
/// and neither costs a universe.
pub const OPENING_SEED: u64 = 1;

/// The size the window **requests** at startup, in logical points (§14.0).
///
/// Requested, not guaranteed: the engine hands this to the platform, which may
/// scale it or clamp it to the monitor. So the frame test that reads this
/// constant is asserting about the geometry the program *asks for*, which is the
/// thing this repository controls — not about a physical window it cannot
/// observe headlessly.
///
/// **Shared with the tests deliberately, and that is the whole reason it is a
/// constant.** `every_cell_is_reachable_at_the_size_the_window_actually_opens`
/// is the only test that says anything about the shipped program rather than
/// about a harness someone sized generously — and a retyped `800.0` in the test
/// is a number somebody guessed, which stops tracking the window the first time
/// the window changes.
///
/// It matters more at Step 2 than it did at Step 1. The widest measured table is
/// **74 columns**, whose content extent is ~3531 px and does not shrink to fit:
/// `egui`'s `Grid` sizes columns to content, so there is no reflow to rescue a
/// window that is too narrow. Most universes therefore do *not* fit, which is
/// why the grid scrolls horizontally rather than assuming width.
/// **`u16`, so there is one source of truth and no cast.** The window wants
/// `u32` and the test harness wants `f32`, and `as` conversions are denied at
/// the workspace (§13.1). Both are reachable from `u16` through a *lossless*
/// `From`, so this constant serves both without a second constant that could
/// drift from it and without an `#[expect]` sitting over a narrowing cast.
pub const WINDOW_SIZE: [u16; 2] = [1200, 700];
