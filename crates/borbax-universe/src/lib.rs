//! Layer 0 — the generated chemistry of one universe (spec §7).
//!
//! A 64-bit seed becomes a complete periodic table. Nothing here is read from
//! a data file and nothing is modelled on real chemistry (§5, G1–G6).

pub mod element;
pub mod naming;
/// **`pub(crate)`, not `pub`, and that is §8.6 enforced by the compiler rather
/// than by discipline.** Everything here runs at intern time and its results are
/// already materialised on [`element::Element`]; a step loop has no reason to
/// call any of it. `contacts_upto` costs 11.5 ns, so a redundant call from a hot
/// path would never show up in a profile — visibility is the only guard that
/// costs nothing.
#[expect(
    clippy::redundant_pub_crate,
    reason = "the items are `pub(crate)` *and* the module is, deliberately: if the \
              module is ever widened to `pub`, its contents must not silently become \
              public with it. The redundancy is the belt to the module's braces"
)]
pub(crate) mod packing;
