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
//! [`geodesic`] is the exception to "not part of the simulation": it is
//! Task 7's content, written here first because the measurement needs it, and
//! it lifts into `borbax-molecule` unchanged whatever G2 returns.

pub mod geodesic;
