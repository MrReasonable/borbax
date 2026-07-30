//! Measurement harnesses for the two open questions that gate V0.
//!
//! This crate is not part of the simulation. It exists because two design
//! premises are stated in the spec but have never been measured, and neither
//! is settleable by argument:
//!
//! - **G2 — does signature space have locality at all?** §8.4's evolvability
//!   argument assumes a one-atom edit moves a molecule's signature a small
//!   distance. An earlier measurement put 32% of one-atom mutants *further*
//!   from their parent than an unrelated molecule sits, which would falsify
//!   the premise rather than reveal a bug.
//! - **Exit criterion 9 — can the acceptance instrument return a negative?**
//!   An instrument that cannot fail is not an instrument.
//!
//! Every harness here ships with a positive and a negative control. A
//! measurement that cannot separate a descriptor known to have locality from
//! one known to have none is measuring itself, and its number must not be
//! reported.
//!
//! [`geodesic`] was the exception to "not part of the simulation": it is
//! Task 7's content, written here first because the measurement needed a
//! rotation table before Task 7 was due. Task 7 has since lifted it into
//! `borbax-molecule`, and the re-export below is deliberately not a second
//! copy — if the table G2 measures ever diverged from the table the binding
//! kernel searches, the measurement would be invalidated without anything
//! failing.

pub mod embed;
pub mod g2;
/// Sample directions and the rotation table, re-exported from where they live.
///
/// See [`borbax_molecule::geodesic`]. Kept as `crate::geodesic` so the
/// harnesses read the same either side of the lift.
pub use borbax_molecule::geodesic;
pub mod molecule;
pub mod openended;
pub mod rng;
pub mod signature;
