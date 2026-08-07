# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

Borbax simulates chemical evolution — abiogenesis up to protocells. Molecules are
shapes; binding is a geometric complementarity test; catalysis, membranes and
heredity are all downstream consequences of that single mechanism.

**The chemistry is currently invented from a 64-bit seed rather than modelled on
the real world, and that is a starting point rather than a principle.** It was
adopted because another model refused to help build a chemistry simulation at
all; the design it produced is worth keeping on its own merits, and none of it
depends on the chemistry being *unlike* the real world. Issue #26 is the live
direction: rebuild on approximated real physics with the seed perturbing the
constants, "our universe" being the values left unperturbed. §5's G5 was
withdrawn on 2026-08-06 for standing in the way of exactly that.

**V0 Tasks 1-10 of 21 are complete**, and **Step 1 of the viewer** on top of
them — the workspace and fiction-guarantee gate, `borbax-units` with the
portable-transcendental chokepoint, `borbax-rng`, element generation by fusion,
bond energies with the assembled `Universe`, molecule graphs with canonical
labelling, the geodesic directions with the 60-rotation permutation table, 3D
embedding by stress majorization, the shape signature, and the binding kernel.
Everything else is spec and plan. Read them before proposing anything:

| Document | What it is |
|---|---|
| `docs/superpowers/specs/2026-07-26-borbax-prd.md` | The PRD. Numbered sections (§n) are referenced everywhere — cite them in code comments and commit messages. V0 exit criteria are §23. |
| `docs/superpowers/plans/2026-07-26-borbax-v0.md` | V0 plan, Tasks 1–10. Holds the Global Constraints and File Structure that both plan files share. |
| `docs/superpowers/plans/2026-07-26-borbax-v0-chemistry.md` | V0 plan, Tasks 11–21 (folding onward). Complete. |

Work is plan-*ordered*, not plan-transcribed: execute tasks in order, ticking the
`- [ ]` checkboxes. Commit after every task, conventional-commit style, with the
trailer `MrReasonable <4990954+MrReasonable@users.noreply.github.com>`.

**Start every task with a three-amigos session, and never drop the plan's code
in.** Three roles, told to argue rather than defer:

| Role | Asks |
|---|---|
| **product** | What does the spec actually require here? What is in and out of scope? Which routed requirements from earlier reviews bind this task? |
| **developer** | What does the *current* API make possible? Where does the plan conflict with what actually shipped? What is the smallest honest design? |
| **QA** | What are the test cases? For each, what discriminates a real implementation from a plausible one? |

They must agree the **test list before any implementation**. Then deliver
outside-in: the acceptance-level failing test first, verify it fails for the
stated reason, work inward.

**Why the plan's code blocks are history and not a starting point.** A step
written at Task 1 is specified against an API that later tasks change. Measured
on Task 9: the plan's Step 3 code did not compile against Task 8's shipped
`Embedding`/`Mol12` at all — 18 errors — and after mechanical repair its
arithmetic carried **six** §13.1/§13.4 determinism violations into a crate with
zero. Its Step 1 test block called `embed(&Mol12, ..)`, three tasks after `embed`
started taking a `&CanonMol`, and named two methods its own routed requirements
forbid.

The subtle part, because it caught me: the rule below — *justify a plan deviation
with a verdict a tool can issue* — is still right, and obeying it is what
produced the behaviour this section now forbids. I dropped the plan's code in
**to measure it**, then reverted and wrote tests first. The measurement was
honest and the deviation was correct; the mistake was treating the plan's code as
a candidate at all. Once you accept that plans go stale, quantifying the
staleness is ceremony. Read `embed`'s signature instead.

The plan step keeps the two jobs it is good for: the record of intent, and the
carrier of routed requirements and discriminators from earlier tasks' reviews.
Those are the parts to read closely — Task 9's six routed requirements were all
load-bearing and all correct.

`superpowers:subagent-driven-development` and `superpowers:executing-plans` still
apply for sequencing and review checkpoints; the three-amigos session replaces
their "read the plan and follow the steps" opening.

## How work reaches `main`

**Worktree → commit → `/review-pr` → push → PR. Never push to `main`.**

`/review-pr` runs **before** the push, not after it, and not "if the change
looks risky". Always. Its Step 0 gate must be green before any specialist is
dispatched — a red suite means fix the suite, not review it. All six:
`fmt`, `clippy -D warnings`, `test`, **`test --release`**, `cargo doc`,
`xtask`. The release leg is the one an earlier version of this sentence left
out, which is precisely the leg whose absence let two defects ship.

**Why it is unconditional, with the number attached.** Task 2 was reviewed by
two specialists before it was committed, their findings were applied, and the
*fixes* went to `main` having been read by nobody. A five-specialist review of
that commit then found nine more defects — and **the two worst were in the
fixes**, not the original code:

- `impl Add for Mass` documented a bound that was false by 1024×; two values
  the constructor accepted summed to −2 under release wrapping.
- `canonical_cmp` documented "NaN sorts last", which was true on aarch64 and
  false on x86-64.

Both shipped **with tests that passed while the defect stood**. That is the
recorded pattern — 6 of 7 proposed fixes rejected or materially amended in an
earlier round — happening one level up. Finding a defect and repairing it are
different activities, and only the first one had a process.

So when reviewing work that has already been reviewed, **point the reviewers at
the fixes** and route each to a specialist *other* than the one who proposed
it. That framing is what surfaced both of the above.

Branch protection enforces the rest server-side: PR required, all CI checks
green, **one approving review**, all conversations resolved, linear history, no
force-push, `enforce_admins` on. Bypassing is a deliberate, visible act:
`gh api -X DELETE repos/MrReasonable/borbax/branches/main/protection/enforce_admins`.

**"Rebase and merge" is the only merge button**, from 2026-08-02 — squash and
merge-commits are both disabled at the repository. So a PR's commits land on
`main` individually and each one has to stand on its own. Task 10 landed thirteen
that way, one per review round.

**The approving review comes from CodeRabbit**, configured in `.coderabbit.yaml`
with `request_changes_workflow: true`. That makes its verdict a real gate rather
than a comment: it posts `CHANGES_REQUESTED` while findings are open, which sets
`mergeStateStatus` to `BLOCKED`, and flips to `APPROVED` once they are resolved.
Two consequences worth knowing before they surprise you:

- **`dismiss_stale_reviews` is on, so every push kills the approval.** The last
  commit before merge needs its own verdict.
- **`auto_pause_after_reviewed_commits: 2` means CodeRabbit stops reviewing
  partway through a long PR**, and a green check can therefore be stale. Check
  *which commit* it last reviewed, not just the check colour, and ask with
  `@coderabbitai review` when the branch is ready.

**Zero open threads is not zero open comments.** CodeRabbit folds nitpicks into
the review *body*, where they never become review threads — a thread query
returns clean while a finding sits unaddressed. Sweep both: the thread list, and
every review body for `Nitpick comments (n)` / `Actionable comments posted: n`.

Next action if starting fresh: **viewer Step 4**, sticks between bonded atoms and colour from generated properties. See
`docs/superpowers/plans/2026-08-01-borbax-viewer.md`. Task 10 merged as PR #19
on 2026-08-01, bringing V0 Tasks 1-10 in; viewer Step 1 followed it and created
`crates/borbax-ui`; **Step 1b landed a name instead of a number** — typing
`emily` gives that name's universe, with the seed it produced written into the
box beside it so it can be copied down.

Step 1b also closed §13.1's hasher ban, which existed only as a comment:
`DefaultHasher`, `RandomState` and `SipHasher` passed all six gate legs before
it. The ban now has four enforcement points, and each exists because the
previous one could go silent — `clippy.toml`'s `disallowed-types`, `xtask`'s
textual `BANNED_TYPES`, per-type `#[expect]` liveness anchors (a wrong module
prefix leaves clippy resolving nothing at exit 0), and a check that the lint is
`deny` in the manifest (an `#[expect]` sets the level for its own item, so it
cannot see the level being dropped).

**The crate is `borbax-ui` and the binary is `borbax`.** The viewer plan
originally said `borbax-viewer`; three reviewed documents say `borbax-ui`,
including Task 19 Step 3, which is written against `borbax-render` and
`borbax-ui` as a *pair* of backends over one geometry. The plan was amended, not
the spec.

**The viewer is V1 work pulled forward, and that is a deliberate deviation from
the spec's sequencing rather than an oversight.** §14.5 puts static headless SVG
in V0 (Task 19) and the interactive viewer in V1. Ian's instruction on
2026-08-01 was to bring the visuals forward, because a project whose point is
showing a child how chemistry works had ten more tasks to run with nothing
to look at. Task 19's SVG goldens are **not** replaced — they are how §13.6
proves the geometry is identical across platforms, which a screenshot cannot
do. The plan file says all of this; read it before re-sequencing anything back.

**The viewer's seam is per-file and `cargo xtask` checks it — but a plain grep
will lie to you.** `state.rs`, `molecule.rs` and `orbit.rs` name no engine or UI
type at all; `panel.rs` names `egui` and contains no `format!`; `app.rs` and
`scene.rs` may name `bevy` and `egui` and may build no strings; and `main.rs` is
one statement with no licence at all — **in code**. The prose necessarily quotes
what it forbids, so a raw `grep -rl bevy src/` returns four files against the two
that name it in code, and `grep format! panel.rs` returns one, all correctly and
all uselessly; the guard strips comments **and the contents of string
literals**, and refuses a file that ends inside either rather than scanning a
half-lexed one. The stronger strip is not fussiness: the weaker `code_only`
read an `#![expect(..., reason = "…eframe…")]` as an import and failed the seam
on correct code, which `xtask/src/main.rs` records at the fixture for it. That is
what lets `egui_kittest` drive the real drawing code headlessly in CI — `egui`
lays out on the CPU and cannot open a window. It also means the **computed
status string** comes from `state.rs`, so a test asserting on it is asserting on
what is shown; a `format!` migrating into `panel.rs` leaves those tests green
over a window they have stopped describing. (Fixed labels like `ui.label("seed")`
are painted in `panel.rs` and are not that string — the rule is about strings
built from state, which is what a test can be wrong about.)

**Step 3's §13.1 question is CLOSED, and it was closed in the direction the rule
required.** The camera's trig routes through `det_math`; no exception was
granted and no path exemption was added. The cost of compliance turned out to be
one import and the spelling `det_math::cos(t)`, against a policy change to a
hard invariant — and the exception would not even have worked, since an
`#[expect]` silences clippy and leaves `xtask`'s textual scan firing.
`check_libm_has_one_home` landed with it, closing a hole where a bare `libm::`
call passed all six gate legs.

The argument that the viewer produces no results and so buys nothing by
complying is **still not pre-approved** for anything else. It was not the reason
this went the way it did: the reason is that the guard refuses the alternative,
which is a verdict a tool issues rather than a judgement.

**Task 5b is scaffolding, not a V0 task.** It respecifies `abundance` as a
fusion/fission process and is fully written up, but it gates nothing, the shipped
one-line profile already satisfies §7.2 in 2000 of 2000 universes (`f_w` median
2.44, p5 2.10), and its own steps need a beaker that does not exist until Task
19. Sequence it after the beaker runs, when its three free parameters can be
located by §7.2's battery instead of chosen by eye.

## Commands

These exist once Task 1 lands; none of them work before that.

```bash
proto install                                                   # toolchain from .prototools (Rust only)
cargo xtask setup                                               # prek + git hooks (once per clone)
cargo test --locked --workspace                                 # full suite
cargo test --locked --workspace --release                       # again, in the profile that mints goldens
cargo test --locked -p borbax-molecule signature                # one crate, one filter
cargo test --locked -p borbax-units --doc                       # doc tests (compile_fail unit-mixing tests live here)
cargo xtask                                                     # fiction-guarantee checks (spec §5) — CI gate
cargo run --locked -p borbax-experiments --release --bin g2     # G2 locality measurement (docs/experiments/)
cargo run --locked -p borbax-ui --release                       # the viewer — a window, a seed, a universe
cargo fmt --all --check
cargo clippy --locked --workspace --all-targets -- -D warnings
RUSTDOCFLAGS=-Dwarnings cargo doc --locked --workspace --no-deps --document-private-items
cargo run --locked -p borbax-cli --release -- goldens --emit    # golden state hashes (cross-platform matrix)
UPDATE_GOLDENS=1 cargo test --locked -p borbax-render           # regenerate SVG goldens after a deliberate render change
```

**`--locked` everywhere is not decoration.** From Task 2 the workspace has a
runtime dependency in a result-affecting path (`libm`, a caret range), so a
stale `Cargo.lock` lets the three CI legs resolve different versions — and the
golden matrix would report that as a portability failure in the simulation,
which is the one diagnosis it must never give wrongly. The `cargo xtask` alias
carries it too.

**One command in that list cannot take it, and one used to be missing it.**
`cargo fmt` rejects `--locked` outright (it is `cargo-fmt`, not a cargo build
command) — yet it shells `cargo metadata`, which will happily rewrite
`Cargo.lock`. That is the standing exception, and the reason the CI `fmt` job
carries the same gap at `.github/workflows/ci.yml:37`. The one that *was*
missing it is worse and is now fixed: `UPDATE_GOLDENS=1 cargo test -p
borbax-render` is the only command in the block whose **output is a golden
artefact**, and it was the one permitted to silently re-resolve the lockfile
first. Every result-affecting CI leg — `xtask`, `clippy`, `test`,
`test --release`, `goldens --emit` — was already correct.

**`cargo doc` is part of the gate, not a convenience.** `[workspace.lints.rustdoc]`
denies `broken_intra_doc_links`, and *nothing else enforces it* — clippy does
not run rustdoc lints, measured with a planted broken link that passed
`clippy -D warnings` cleanly. Same shape as the `[lints] workspace = true`
trap: present, believed, inert. It matters here because nearly every type
documents its contract with intra-doc links to its siblings, and a rename rots
them silently. It runs as its own CI job and in the pre-push hook.

**Run the suite in release as well as debug.** `debug_assert!` and
`overflow-checks` both key off `debug_assertions`, so a guard can be present in
`cargo test` and absent from `goldens --emit`. Task 2's review found that
defect twice, in adjacent functions. CI runs both legs.

`cargo xtask` is an alias for `cargo run -p xtask --`, from the checked-in
`.cargo/config.toml`. CI and the git hooks deliberately keep the long form:
flags belonging to `cargo run` — the hook's `-q` — must sit before `run`, and
an alias cannot splice them in.

`borbax-cli` subcommands: `beaker`, `sweep`, `band`, `battery`, `render`, `goldens`.

Transcendentals route through `crates/borbax-units/src/det_math.rs`, backed by
the `libm` crate — rust-lang's pure-Rust MUSL port. It lives in `borbax-units` and nowhere else because `borbax-rng` and
`borbax-universe` both call it and both sit below `borbax-molecule` in the crate
order; anywhere higher and the one file guaranteeing portability is unreachable
from two crates that need it.

**`libm` does dispatch on architecture, and "it doesn't" is the wrong reason to
trust it.** Measured on 0.2.16 — the version `Cargo.lock` currently resolves,
not a pin: `Cargo.toml` carries the caret range `0.2`. Two things hold it, and
naming only the first invites someone to reopen this thread. The **lockfile**,
enforced by `--locked`, keeps the three CI legs on one version. **Renovate TIER
B** is what makes a bump deliberate rather than merely visible: the cargo
manager's `rangeStrategy: auto` resolves to `update-lockfile` for a bare `0.2`,
so a libm 0.2.17 lands as a lockfile-only PR matching `matchDepTypes:
["dependencies", "workspace.dependencies"]` — `automerge: false`, labelled
`deps-runtime` and `determinism-review`, with `lockFileMaintenance` separately
disabled for this exact hazard. That chain already delivers what the dependency
doctrine asks, so `=0.2.16` would add nothing but a resolution failure the first
time something wants `^0.2.17`. `arch` is a *default* feature, and it
routes `sqrt`, `fma`, `rint`, `ceil` and `floor` to hardware — with `fma` on
x86-64 doing **runtime CPU feature detection** between FMA3, FMA4 and soft.

**"The only FMA instructions in our release binaries are inside `libm::cbrt`"
was true and stopped being true at viewer Step 4, on aarch64 specifically.**
`palette.rs`'s `sector.rem_euclid(2.0)` is a `frem` by a power of two. Measured
by compiling `palette.rs` verbatim — it imports nothing, so it compiles
standalone — for the three triples the CI matrix actually runs
(`.github/workflows/ci.yml`'s `macos-latest`/`ubuntu-latest`/`windows-latest`
are `aarch64-apple-darwin`, `x86_64-unknown-linux-gnu` and
`x86_64-pc-windows-msvc`), at this workspace's actual release profile
(`lto = "thin"`, `codegen-units = 1`):

| leg | `rem_euclid(2.0)` | `rem_euclid(360.0)` |
|---|---|---|
| `aarch64-apple-darwin` | inlined, **fused**: `frintz` / `fmsub d1, d1, d2, d0` | `bl _fmod` |
| `x86_64-unknown-linux-gnu` | `callq *fmod@GOTPCREL` | `callq *fmod@GOTPCREL` |
| `x86_64-pc-windows-msvc` | `callq fmod` | `callq fmod` |

**An earlier version of this paragraph said x86-64 also inlines the power-of-two
case, "matching every CI leg", and named `x86_64-apple-darwin` as its evidence.**
`x86_64-apple-darwin` is not a CI leg — GitHub's `macos-latest` runners are
`aarch64` — and its default `target-cpu` (`penryn`, SSE4.1) is not the two real
x86-64 legs' default (the generic `x86-64` baseline, SSE2 only), which is
exactly the difference that makes `roundsd` unavailable and pushes LLVM to a
libcall. Measuring the wrong triple produced the wrong conclusion; this table
is against the three triples that actually run.

So aarch64 is the only leg with an FMA instruction here — a real change from
"the only FMA is inside `libm::cbrt`" — and it is also the only leg that inlines
the power-of-two `rem_euclid` at all. Both x86-64 legs take a genuine `fmod`
libcall for **both** calls, same as aarch64's `360.0` call, which is a `frem`
by a non-power-of-two and is a libcall on every leg without exception.

**Why this is still safe, and the argument is provable rather than merely
measured.** `fmod`'s mathematical result — `x - n·y` for the integer `n`
nearest zero — is *exactly* representable in `f64` for any finite non-zero
divisor: no rounding step exists for a correct implementation to differ on.
So any two implementations that are each individually correct — glibc's,
`msvcrt`'s, and LLVM's inlined bit-fused expansion on aarch64 — must produce
the identical bit pattern, by definition of "correct", not by coincidence
measured over a finite sample. `rem_euclid`'s sign-fixup wrapper around the raw
remainder (`if r < 0 { r + y.abs() } else { r }`) is the same source
expression on every leg, whatever instructions each compiler chooses to
implement it with — LLVM's own range analysis can and does eliminate the
fixup entirely at one of the two call sites here, since `hue` is provably
non-negative before it reaches `rem_euclid(360.0)`, and neither x86-64 leg's
disassembly needs to match aarch64's instruction-for-instruction for the
argument to hold. What has to be true, and is, is that each leg's fixup
correctly implements the same source-level conditional over an
already-identical raw remainder.

What *is* measured, as supporting evidence that the inline expansion itself has
no hidden bug (a sign-of-zero mismatch would be exactly this class of thing):
hashing `aarch64-apple-darwin`'s fused expansion against `x86_64-apple-darwin`'s
inlined-but-unfused SSE4.1 expansion (`mulsd`/`roundsd`/`subsd`/`andpd`/`orpd`)
over 20M random finite `f64` bit patterns gives `d22129714af37cf9` on each —
two *different* inline expansions agreeing is evidence neither has a
platform-specific bug, though neither is a real CI leg's libcall path. Glibc's
and `msvcrt`'s actual `fmod` were not directly cross-executed from this
environment; the exactness argument above is what stands in for that, and it
is a stronger claim than a hash match would be, not a weaker one — a hash
match over any finite sample cannot rule out a divergence outside the sample,
where "any two correct implementations must agree" rules it out entirely,
provided both are in fact correct implementations.

The reason the fused expansion is a correct implementation is visible in the
assembly: it ends `movi.2d`/`fneg.2d`/`bit.16b`, which takes the magnitude from
the computed remainder and the **sign bit from the dividend** — reproducing
`fmod`'s sign-of-zero rule exactly, the one place a naive reimplementation
would most likely diverge from a real libm.

**Three more claims from earlier versions of this paragraph were wrong, and
each was corrected by measuring rather than by re-reasoning:**

- "The legs differ in 23.7% of cases" measured the *naive* model
  `x - trunc(x*0.5)*2`, without the sign fixup, against the real lowering:
  4,746,398 of 20,000,000, every one a sign of zero and none a value — exactly
  the count of negative zeros in the output. Naive-model versus reality, not
  leg versus leg.
- "The branch is safe because LLVM proved a range, and widening the argument
  type would make `-0.0` reachable" — `-0.0` is reachable **today** (`hue =
  -0.0` and `hue = -360.0` both give `sector` bits `8000000000000000`), and it
  does not matter: the next operation is `(t - 1.0)`, and `-0.0 - 1.0` and
  `+0.0 - 1.0` are both `bff0000000000000`. The sign is erased one operation
  later, before `.abs()` and before any `to_bits()` — but **not** before the
  six `sector < n` comparisons, which see `-0.0` and agree with `+0.0` on every
  one.
- "x86-64 inlines with no libcall at all, matching every CI leg" — false,
  corrected above; it was measured against a target that is not a CI leg.

Getting any of this wrong would invite a future reviewer either to reject a
legitimate `rem_euclid` in result-affecting code, or to wave through a
genuinely inexact operation on the grounds that we already ship one that
"differs 23.7% of the time" — neither of which is true.

That is safe, and the correct statement of why: **libm dispatches only on
operations IEEE-754 specifies exactly**, which is the same criterion
`clippy.toml` already applies. Verified rather than assumed — `cbrt` hashed over
2.4M inputs with and without the `arch` feature is bit-identical. The wrong
reason matters because it would wave through a future libm that added an `exp`
arch path; `exp` already has one for x86 without SSE2, invisible only because
the §13.4 matrix has no 32-bit leg.

Every `exp`, `ln`, `sin`, `cos` and `powf` in the workspace goes through it.
`sqrt` and the four arithmetic operations stay native — IEEE-754 specifies those
exactly. A direct `.exp()` or `.cos()` on an `f64` in library code is a
determinism bug, not a style choice.

The dev profile is deliberately `opt-level = 1` (deps at 2): the unit newtypes are
zero-cost only after inlining, and the beaker harness's value is that a
chemistry-tuning iteration takes seconds, not minutes.

## Architecture

Two seeds, and a strictly one-directional pipeline:

```
universe_seed → universe_gen  (elements, valences, bond energies, folding rules)
                    ↓          generated once, immutable, cached
world_seed  →  world_gen      (V1: terrain, vents, climate, abundances)
                    ↓
              simulation core  bulk 10⁶ patches ⇄ stochastic 10³ ⇄ molecular ~10  (V1)
                    ↓          V0 has none of this — one well-mixed beaker
              borbax-geometry → SVG (headless, golden-tested) | egui+wgpu (V1)
```

Crate dependency order — later depends on earlier, never the reverse:

`borbax-units` → `borbax-rng` → `borbax-universe` → `borbax-molecule` →
`borbax-reaction` → `borbax-beaker`, with `borbax-geometry` → `borbax-render`
and `borbax-cli` on top.

`borbax-molecule` carries the conceptual weight (graph → canonical form → 3D
embedding → signature → binding → fold → cavity) and is split into small files
because that is where nearly all iteration happens.

### The three ideas that everything else is a consequence of

**1. One mechanism, reused at every scale.** `affinity(A, B)` scores two shape
signatures over the 60-element icosahedral rotation group: bumps must meet
hollows, positive faces must meet negative ones. Molecules sticking together,
an enzyme recognising a substrate, a membrane self-assembling, and a pore gating
passage are all *that same function* at a different scale. If a feature seems to
need a second mechanism, question the feature (spec §8.3).

**2. Nothing about life is hardcoded.** No `struct Cell`, no replication rule, no
"if conditions are right, spawn a protocell". Catalysis is not implemented — a
folded chain that happens to have two well-placed cavities holds two reactants
adjacent, and the rate enhancement falls out of the geometry (§8.5). Protocells
and death are *detected* (via RAF closure), never declared. A special case for
life means the physics is wrong; fix the physics.

**3. Everything expensive is per-species, not per-molecule** (§8.6). Canonicalisation,
3D embedding, signature construction, folding, cavity extraction and per-bond decay
rates are pure functions of the species, computed once at intern time.
**No code reachable from a simulation step may call any of them.** Violating this
makes the performance budget unreachable by ~2 orders of magnitude, and the
regression usually arrives disguised as a small feature — most often a
per-molecule damage field, which would turn decay from one propensity channel
back into a sweep over every molecule.

### Design choices that are load-bearing, not incidental

- **Icosahedral sampling** (D ∈ {12, 42, 162} only) makes rotation an exact
  *permutation* of sample directions — precomputed, no interpolation, no
  resampling error. Intermediate resolutions destroy that and are not options.
  42 is the prior; V0 sweeps and **locks** the value (§22.2).
- **Stress majorization, not Laplacian eigenvectors**, for the 3D embedding.
  Eigenvector sign is arbitrary and degenerate eigenvalues are common in the
  symmetric graphs small molecules actually are — that arbitrariness would
  propagate straight into species identity (§8.2).
- **FCC lattice for folding** (12 neighbours, integer coords with even parity)
  and **3D throughout**: in 2D a chain has ~√n interior against √n boundary, so
  burial protects nothing and a cavity is an open notch. The catalysis and
  solvent-attack mechanisms both depend on real burial (§8.4).
- **Binding pairs each direction with its partner's *antipode*.** Two bodies in
  contact touch along opposite directions, so the kernel indexes `b` through
  the antipode of the rotated direction, never the rotated direction itself.
  Since Task 7 that composition is precomputed as `Geodesic::contact_perms(r)`,
  and the raw `rotation_perms(r)` is for Task 9's single-signature
  canonicalisation.

  **Two spellings break it and both are silent.** Reaching for
  `rotation_perms(r)` never applies the antipode. Applying `anti` a *second*
  time on top of `contact_perms` applies it twice, which cancels — and that is
  the one an earlier version of this bullet actively instructed, because it said
  the kernel must index through `ANTI[PERM[R][i]]` and `PERM` had been renamed
  under it. Measured at all three resolutions:
  `anti[contact_perms[r][i]] == rotation_perms[r][i]` in **100% of entries**
  (720/720, 2520/2520, 9720/9720). The apparently-defensive second `anti`
  produces exactly the defect.

  **Say which object is improper, because the two readings invert.** As
  *index permutations*, `anti` is `−I` and the correct composition
  `anti ∘ rotation` has determinant −1, so the correct table is the improper
  coset and both mistakes land on the 60 *proper* rotations. As *transformations
  of the partner body* it is the other way round: the correct kernel searches 60
  proper poses of B, and dropping the antipode searches B point-inverted, which
  is improper. An earlier version of this bullet asserted the second reading
  while the surrounding sentences were about the first, which is how it managed
  to be simultaneously right and backwards.

  **The consequence needs neither framing and is provable, not measured.** Since
  `mirror(B)[i] = B[anti[i]]`, pairing `a[i]` with `b[rotation_perms(r)[i]]` is
  *exactly* the correct kernel evaluated against `mirror(B)`. So the defect does
  not merely score wrongly — it computes A's affinity with B's enantiomer, which
  makes enantiomers interchangeable and homochirality (§22.8) impossible by
  construction rather than emergent. Nothing fails.

  **A complement built at the same index does not merely fail to test this — it
  inverts.** Measured at D=42 during Task 7's review, §8.3's kernel scored four
  ways: with `ANTI`, the same-index fixture gives **−4.602147** and the
  `ANTI`-built fixture gives **0.000000**; without `ANTI`, exactly the reverse.
  So an implementer who writes the obvious fixture sees the *correct* kernel
  score −4.6, "fixes" the kernel until it scores 0, and arrives at precisely
  this defect with a green test.

  **The zero is exact and structural, not a measurement.** With a same-index
  complement (`b.r[i] = gap − a.r[i]`, `b.a[i] = −a.a[i]`) the anti-less kernel
  gives `ds = a.r[i] + (gap − a.r[i]) − gap = 0` and `dc = a.a[i] − a.a[i] = 0`
  in every direction, so it scores exactly 0 at the identity — and 0 is the
  global maximum of a negated sum of squares, so no rotation can beat it. The
  correct kernel on that fixture gives `ds = a.r[i] − a.r[anti[i]] ≠ 0`. This
  holds for every universe seed, not just the one measured. Task 10 therefore carries a required probe:
  **substitute `g.rotation_perms(r)` for `g.contact_perms(r)`** in `affinity`
  and **both** `a_physical_complement_scores_zero` and
  `the_mirror_complement_does_not_fit` must fail. If only one does, the pair is
  not discriminating.

  The probe is stated as a *substitution* because the rename removed the thing
  its first wording told you to delete — there is no `anti[...]` in the index
  expression any more. A probe instruction naming an edit that does not exist
  is a probe nobody runs.
- **Binding searches rotations but not reflections** — 60 elements, not 120. This
  makes Borbax chemistry *handed*, so homochirality is an emergent result the
  simulation could produce rather than something excluded by construction (§22.8).
- **Decay is foundational, not a balancing mechanic.** Selection *is* differential
  persistence; without decay, fitness differences can never express themselves
  (§9.4). Every decay source is a physical process already being run for other
  reasons. The decay band is the most sensitive parameter in the system — suspect
  it first whenever a run produces nothing.

## Hard invariants

Violations of these are correctness bugs, not style preferences.

**Determinism** — same seeds → bit-identical output on macOS/aarch64,
Windows/x86-64 and Linux/x86-64, enforced by a CI golden matrix (§13.4, §13.6):

- All randomness comes from `borbax_rng::Stream` (counter-based, one stream per
  `Domain` × index). Never `rand::random`, never a thread-local RNG, never
  `RandomState`. `Stream` is deliberately not `Copy`.
- No `HashMap`/`HashSet` iteration in any result-affecting path. Use `BTreeMap`,
  or sort into a `Vec` first. Every float-key sort carries an ID tie-break.
- No wall-clock. No scheduler-ordered `.sum()`/`.reduce()` — fold into per-chunk
  partials indexed by chunk ID, combine in index order.
- No platform transcendentals. `exp`, `ln`, `sin`, `cos`, `powf` are not
  correctly-rounded by IEEE-754 and genuinely differ between libms; route them
  through a vendored pure-Rust implementation or a table. `+ - * /` and `sqrt`
  are exactly specified and may stay native.
- Float accumulation order is pinned and commented as load-bearing. Nobody
  "tidies" it later.
- `Mass` is fixed-point `i64`, so conservation is exact by construction. Every
  reaction conserves mass exactly, checked in tests.
- `f64` everywhere in simulation arithmetic; no `f32`.
- The `.prototools` Rust pin is part of the determinism contract. Bumping it
  invalidates golden hashes and requires deliberate regeneration, reviewed like
  a change to physics.

**Fiction guarantees** (§5). **G5 was withdrawn on 2026-08-06** — read §5 before
adding any enforcement here, because the burden of proof now sits with keeping a
constraint rather than with relaxing one. What remains, and what `cargo xtask`
still enforces:

- No real chemistry data anywhere — no element tables, reaction databases,
  molecular structures or sequence data. No data files in `borbax-universe`,
  `borbax-molecule`, `borbax-reaction`. (G1, `check_no_data_files`.)
- Generated element names are blocklist-checked against real symbols and names
  (`REAL_ELEMENT_SYMBOLS` in `borbax-universe/src/naming.rs`). (G2.)
- No real-world units. Temperature is `Thermal`, energy `Quanta`, distance
  `Span`, time `WorldYear`, mass `Mass` — distinct newtypes; cross-unit
  arithmetic must not compile. (G4 — and it earns its place on type-safety
  grounds regardless of the fiction argument.)

**Withdrawn, and do not reinstate without reopening §5:** the ban on real
chemical formats, on a Borbax↔real mapping table, and on real-world calibration
targets. Its two checks — `check_no_real_chemical_formats` and
`check_the_palette_reads_only_generated_properties` — were deleted with it.

**G3 and G6 were revised on 2026-08-06, same reasoning as G5.** They no longer
claim the physics is permanently non-isomorphic or permanently unpredictive —
both are now "not necessarily real", a property of the seed rather than a rule
about every seed. They no longer block #26; see the spec §5 for the exact
wording.

**Other:** no `unwrap()`, `expect()`, `panic!`, `todo!` or `unimplemented!` in
library code — all five are `deny` at workspace level.

**"Tests may unwrap freely" is true only with an attribute, and the attribute is
the point.** `clippy.toml` sets neither `allow-unwrap-in-tests` nor
`allow-expect-in-tests`, so the deny **reaches inside `#[cfg(test)]`**. Two
patterns are established and both are fine:

- a module-level opt-in — `borbax-units/src/lib.rs:751` carries
  `#[allow(clippy::unwrap_used, reason = "...")]` on its `mod tests`, and a
  bare `#[allow]` with no reason is a review finding;
- or no fallible spelling at all — `borbax-universe`, `borbax-rng` and
  `borbax-molecule` contain **zero** `.unwrap()` calls between them, using
  `assert!`/`assert_eq!` for conditions and
  `unwrap_or_else(|| unreachable!("why this cannot happen"))` where a fixture
  guarantees an `Option`. `clippy::unreachable` is deliberately not enabled: it
  is the spelling that names its own precondition.

Prefer the second in new code. `bonds.rs:1341` records why `unwrap` is rarely
the fix even when it is permitted — it turns "this constant stopped being valid"
into a different assertion firing with a message naming the wrong cause.

## Scope discipline

V0 is chemistry in a test tube: no world, no spatial simulation, no live viewer —
but SVG rendering from the start, because a shape-based chemistry cannot be
developed blind. V1 adds the region, multi-scale engine, keyframes, neutral
shadow and the Chronicle, targeting protocells. Genomes are **V2** and are
explicitly deferred however tempting they become; V1 succeeds or fails on §15.1
alone.

Two claims the design deliberately refuses to make, both from the prior-art
survey (§2.7): never assume self-maintaining units compose into higher-order
ones (the field's weakest result), and never declare a plateau by eyeballing a
flattening curve — fit competing models and compare, and read every novelty
metric against a neutral shadow run.

## Agents

Six project agents in `.claude/agents/`, most mapped to a named risk in §20
rather than to a generic role. They review; none of them writes code.

| Agent | Invoke when |
|---|---|
| `rust-developer-expert` | On any new module before it settles; when a type or trait is introduced; when code works but reads badly; as the last pass before a task is committed. Also owns **dependency review** — whether a well-established crate should be doing this instead of us. |
| `rust-performance-expert` | Before committing to a data structure in a hot path; when a task specifies performance targets; to audit that an optimisation preserved determinism. |
| `determinism-auditor` | Before adding parallelism or a dependency that touches results; when a golden hash moves; on any diff touching float arithmetic, collection iteration, sorting or RNG. Every hazard here is invisible on one machine. |
| `emergence-auditor` | Before adding a threshold, special case or per-molecule field; when a behaviour is not emerging and there is a temptation to help it along; on any chemistry or simulation diff. |
| `geometry-numerics-reviewer` | On canonicalisation, geodesic construction, the rotation table, stress majorization, signatures, binding, FCC folding, cavity detection — code that is either correct or silently wrong forever. |
| `alife-researcher` | When implementing an algorithm or metric from the literature (RAF, Gillespie, Bedau activity statistics, neutral networks, plateau fitting), or when a design decision rests on a cited result. Constrained never to return real chemistry (§5). |

The auditors are independent lenses on the same diff and can run in parallel.
Run the emergence and determinism auditors before any commit that touches
physics — those are the two invariants no test suite fully protects.

### Review precedence

The agents are told to argue with each other rather than defer, so genuine
conflicts will surface. This order breaks ties. **It ranks how expensive a
mistake is to discover late, not how important each concern is** — every one
of them matters.

1. **Fiction guarantees (§5).** G1, G2 and G4 stand unrevised and are treated
   as absolute in review — no result is worth breaching them, and there is
   nothing to trade against. Revising one of them is a §5 change in its own
   right, not a tiebreak call. G3 and G6 were revised 2026-08-06 to "not
   necessarily real", a property of the seed rather than a rule about every
   seed, and G5 was withdrawn outright the same day — a conflict against G3 or
   G6 is a design question at precedence 2, not an automatic veto. See §5 for
   the exact wording; this line is a pointer to it, not a restatement that can
   drift from it.
2. **Correctness of the physics.** Geometry and numerics that are silently
   wrong poison everything downstream and are the hardest thing here to
   detect.
3. **Determinism (§13.1, §13.4).** A result nobody can reproduce cannot be
   verified, shared, or debugged.
4. **Emergence invariants (§3, §8.6).** A special case for life means the
   physics is wrong.
5. **Performance (§17).** The budget is real, but a fast wrong answer is
   worthless.
6. **Idiom and readability.** Genuinely matters, and yields to the above.

Two notes on applying it. **2 above 3 is deliberate**: when a numerically
better formulation changes bit patterns, adopt it and regenerate goldens as a
deliberate physics change — do not preserve a worse formulation to protect a
hash. And **most conflicts are not about this order at all** — they are two
correct observations pulling opposite ways, where the resolution is a
measurement or a line of the spec, not a ruling.

### Probe the guard, and then probe the guard's bookkeeping

Every guard gets the mutation it exists to stop, run before the commit — the
question is never "does the suite pass" but "what would this catch that nothing
else does, and have I watched it catch that". That rule is well established
here and has repeatedly paid.

**It stops one level too early.** A guard has assertions *and* it has
instrumentation — counters, corpus filters, `examined > N` self-checks — and
the instrumentation is not covered by probing the assertion. Ask separately:
**does this counter count what its message says it counts?**

Measured on Task 8. `no_two_atoms_share_a_point` asserted `pairs > 2000` with
the message "pairs examined", and incremented `pairs` in the *outer* loop, so
it counted atoms: 4,818 against a true 19,258. The bar was four times weaker
than intended and reported a quantity it did not measure. Every assertion in
that test had been probed by mutation; its bookkeeping had not been read.

This is not testable, and that is the point — a test verifying its own
self-check is infinite regress. It is a reading job, and it belongs in the same
pass as the mutation probe. The same applies to a corpus filter that silently
selects nothing (`assert!(disconnected_seen > 0, ...)` is the shape that
catches it) and to any `assert!(x > N)` whose `N` was chosen before `x` was
measured.

## Dependencies

The best code is code someone else maintains. Reach for the settled crate:
`thiserror` for library errors, `anyhow` for binaries only, `clap` for CLI,
`proptest` for property tests, `insta` for snapshot and golden tests,
`criterion` for benchmarks, `libm` for portable transcendentals, `serde` for
serialisation, `rustc-hash` for deterministic hashing.

The asymmetry that governs the decision:

- **Dev-dependencies are cheap.** `proptest`, `insta`, `criterion`, `rstest`
  cannot affect simulation output. Low bar.
- **A runtime dependency in a result-affecting path is expensive.** It must be
  version-pinned, and upgrading it is a physics change requiring deliberate
  golden regeneration — exactly like a compiler bump (§18.1). Run the
  determinism auditor on any such addition.

Hand-rolled, and the reason matters: three-component vector maths is ~20 lines
using none of what `glam` or `nalgebra` sell. **Not** for determinism — both
would be bit-identical (`glam`'s f64 module has no arch dispatch and documents
cross-platform bit-identity; `nalgebra` special-cases dim 3 to `a + b + c`).
Claiming otherwise was a doctrine error worth naming, because a false
determinism argument gets used later to reject a dependency that is correct.

SMACOF stays hand-rolled on firmer ground: its fixed 240-iteration count is a
determinism requirement, and no solver crate offers "stop after exactly n
iterations regardless of convergence".

*Not* hand-rolled: transcendentals (`libm`), plateau model fitting
(`levenberg-marquardt`), snapshot tests (`insta` — do not hand-roll
`UPDATE_GOLDENS`), benchmarks (`criterion` — a timing assertion inside
`#[test]` is a flaky test that will get deleted, taking the guarantee with it).

Held in reserve: `nauty-Traces-sys` for canonicalisation if the hand-rolled
search ever becomes a bottleneck. It would be a result-affecting runtime
dependency *and* C FFI against a `forbid(unsafe_code)` workspace, so it stays
an escape hatch rather than a plan.
