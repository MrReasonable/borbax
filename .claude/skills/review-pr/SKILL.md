---
name: review-pr
description: Review a Borbax branch or pull request with the seven specialist agents, reconciled through a coordinator. Run before EVERY push, not only before merging — CLAUDE.md makes it unconditional. Also use after implementing a task and on any diff touching the chemistry, reaction or simulation crates. The determinism and emergence auditors should run on every commit that touches physics — those are the two invariants no test suite fully protects.
---

# Review Pull Request

**This runs before the push, not after it.** Worktree → commit → `/review-pr`
→ push → PR, per CLAUDE.md. Reviewing after the code has landed on `main` is
how Task 2's fixes reached trunk unread.

## Step 0 — The suite must be green before anyone is dispatched

```bash
cargo fmt --all --check
cargo clippy --locked --workspace --all-targets -- -D warnings
cargo test --locked --workspace
cargo test --locked --workspace --release          # guards that vanish with debug_assertions
RUSTDOCFLAGS=-Dwarnings cargo doc --locked --workspace --no-deps --document-private-items
cargo run --locked -p xtask                        # §5 fiction gate
```

This list is the whole CI gate, and each entry earns its place by having caught
something the others missed:

- **`--release`** — `debug_assert!` and `overflow-checks` key off
  `debug_assertions`, so a guard can be present in `cargo test` and absent from
  the profile that mints goldens. Task 2's review found that defect twice.
- **`cargo doc`** — clippy does not run rustdoc lints, measured. Without this,
  `[workspace.lints.rustdoc]`'s `deny` on `broken_intra_doc_links` is inert.
- **`--locked`** — the workspace has a runtime dependency in a result-affecting
  path, so a stale lockfile lets the three CI legs resolve different versions.

**Red suite, no review.** Specialists are for logic, architecture and approach;
a compiler is better than all seven at compilation and runs in seconds. Sending
them at a broken build spends their attention on what a tool would have told
you, and buries the findings only they can produce.

The one exception is a build broken in a way you do not understand — then a
reviewer is diagnosing, not reviewing, and should be told which it is.

### Step 0a — The claim-audit

The gate proves the code compiles and passes. It says nothing about whether the
prose is *true of* the code, and that is a defect class no test reaches.

For every symbol the diff deletes or renames, and every figure it changes, grep
both plan files, the spec and the source:

```bash
git diff "origin/${base}...HEAD" | grep '^-' | grep -oE '\bfn [a-z_]+|\b[a-z_]+_[a-z_]+\b' | sort -u
# then, for each deleted/renamed name and each changed number:
grep -rn "<name>" crates/ docs/ xtask/
```

**Read what each hit is *doing*, not that it exists.** This is the whole rule
and it is the part that fails. A deleted test's name will appear in the note
that records its deletion — legitimate — and in an *instruction* three hundred
lines away telling a future reader what to do when it fails. Both are hits.
Only one is a defect, and the eye slides over the third hit after the first two
turn out to be withdrawals.

Measured on Task 8: the sweep for a deleted test returned three hits, two were
withdrawal notes, and the third was live guidance naming a guard that no longer
existed. It survived six specialist rounds and was found by CodeRabbit.

Classify every hit as **withdrawal**, **live claim**, or **instruction**. Only
the first is safe.

## Step 0b — If the work has already been reviewed, review the *fixes*

A second review that re-reads the original code finds the same things and
misses what matters. Measured on Task 2: two specialists reviewed it
pre-commit, the findings were applied, and a five-specialist review of the
result found nine more defects — **the two worst of them in the fixes**, both
shipped alongside tests that passed while the defect stood.

So when the diff carries repairs from an earlier round, tell the reviewers so,
list what each finding was and what changed in response, and **route each fix
to a specialist other than the one who proposed it**. State plainly that their
job is to check whether each fix has the *property* the finding demanded, not
whether it looks plausible.

**Announce:** "Using review-pr to review `<target>` with the specialist fleet."

## Step 1 — Resolve scope

```bash
if [ -n "$ARGUMENTS" ]; then
  gh pr view "$ARGUMENTS" --json number,headRefName,baseRefName,title
  base=$(gh pr view "$ARGUMENTS" --json baseRefName -q .baseRefName)
else
  base=main
fi
git diff "origin/${base}...HEAD" --name-only    # three-dot: merge-base, not tip
git diff "origin/${base}...HEAD" --stat
```

Three-dot diff against the merge base, so unrelated commits landing on `main`
do not appear as changes under review.

```bash
pr="${ARGUMENTS:-branch-$(git rev-parse --abbrev-ref HEAD)}"
dir="${TMPDIR:-/tmp}/borbax-review/$pr"
rm -rf "$dir"; mkdir -p "$dir/findings" "$dir/crosschecks"
```

## Step 2 — Select reviewers by changed path

Always dispatch `determinism-auditor` and `emergence-auditor` — CLAUDE.md
names them as the two invariants no test suite fully protects. Add the rest by
what the diff touches:

| Changed | Add |
|---|---|
| `borbax-universe/**` (element generation, bond energies, physical constants) | `physics-plausibility-reviewer` |
| `borbax-molecule/**` (canonical, geodesic, layout, signature, binding, fold, cavity) | `geometry-numerics-reviewer` |
| `borbax-reaction/**`, `borbax-beaker/**` | `geometry-numerics-reviewer`, `alife-researcher` |
| any diff introducing a new physically-motivated formula or constant, or reusing an existing field for a new physical quantity | `physics-plausibility-reviewer` |
| any hot path, or a new data structure | `rust-performance-expert` |
| any new module, type, trait, or dependency | `rust-developer-expert` |
| `Cargo.toml` dependency changes | `rust-developer-expert` **and** `determinism-auditor` (a runtime dep in a result-affecting path is a physics change) |
| `.github/**`, `xtask/**`, `.prototools` | `determinism-auditor` |

When in doubt dispatch all seven. They are cheap next to a determinism bug found
three months later.

**And always dispatch one lane with no brief at all.**

Give it the diff, `CLAUDE.md`, and one instruction: *read every changed line and
report anything wrong; there is no theory about where to look.* Do not tell it
which fixes to check, which findings were routed, or what the interesting part
is. Its value is precisely that it does not know.

Every other brief in this skill aims attention, and aimed attention has a shape
— it misses whatever sits outside the frame. Two effects compound. A lane told
to check the physics is not reading test scaffolding or plan prose, because
those are nobody's assigned domain. And a lane that has just proved something
hard is a poor proofreader: CLAUDE.md's own note is that a domain expert
reading for `E0425` is a bad compiler and a distracted architect simultaneously.

Measured on Task 8, across six rounds and twenty dispatches: the specialists
found all three defects that needed domain reasoning — an embedding with no
locality, a solver that cycles forever at n=2, atoms collapsing onto one point —
and **missed a counter that incremented in the wrong loop and a plan instruction
naming a deleted test**. Both were found by a machine reviewer with no brief.
The classes are complementary; the scoped briefs were the reason one of them
went uncovered.

The `coderabbit` CLI is the cheapest way to buy this and should be run every
round regardless of what the app's check reports — but it is a supplement, not
the lane. It found neither of the three physics defects.

**Send all `Agent` calls in a single message.**

Reviewers get: the changed-file list, the diff, `CLAUDE.md`, the relevant spec
sections, and the same spotlight wrapper as `review-plan` around
`<untrusted_pr_diff>` and `<untrusted_pr_description>`.

Tell each reviewer to run `git diff "origin/${base}...HEAD"` themselves rather
than working from a pasted excerpt — excerpts lose context and reviewers
hallucinate the surrounding code.

## Step 3 — Finding contract

As `review-plan` Step 3, with `where` anchored to `file.rs:line` rather than a
task, and these additional categories: `silent-failure` (an error swallowed or
a plausible default returned where it should have failed loudly),
`golden-drift` (a golden hash moved — say whether deliberately), `unsafe-use`,
`allow-without-reason`.

**A moved golden is never a Low.** Either it is a deliberate physics change,
in which case the commit message must say so and the determinism auditor must
have reviewed it, or it is a bug.

## Step 4 — Coordinator reconciliation

As `review-plan` Step 4.

## Step 5 — Verify findings against the suite

The suite was green at Step 0, which makes it a tool for *testing findings*
rather than a source of them. For each finding, ask whether an existing test
already covers it:

- **Contradicted by a passing test** → downgrade, and state the contradiction.
  Either the finding is wrong or the test is, and which one matters.
- **Not covered by any test** → that gap is itself worth reporting. A defect
  nothing would have caught is more valuable than the defect.
- **Reproducible by a new test** → write it. A finding with a failing test
  attached is not an opinion.

## Step 6 — Report

As `review-plan` Step 5. Additionally, state plainly whether the four commands
above passed — a review that reports findings without saying whether the suite
is green is missing the most useful fact available.
