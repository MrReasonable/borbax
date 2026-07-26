---
name: review-pr
description: Review a Borbax branch or pull request with the six specialist agents, reconciled through a coordinator. Use before merging, after implementing a task, and on any diff touching the chemistry, reaction or simulation crates. The determinism and emergence auditors should run on every commit that touches physics — those are the two invariants no test suite fully protects.
---

# Review Pull Request

## Step 0 — The suite must be green before anyone is dispatched

```bash
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo run -p xtask                       # §5 fiction gate
```

**Red suite, no review.** Specialists are for logic, architecture and approach;
a compiler is better than all six at compilation and runs in seconds. Sending
them at a broken build spends their attention on what a tool would have told
you, and buries the findings only they can produce.

The one exception is a build broken in a way you do not understand — then a
reviewer is diagnosing, not reviewing, and should be told which it is.

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
| `borbax-molecule/**` (canonical, geodesic, layout, signature, binding, fold, cavity) | `geometry-numerics-reviewer` |
| `borbax-reaction/**`, `borbax-beaker/**` | `geometry-numerics-reviewer`, `alife-researcher` |
| any hot path, or a new data structure | `rust-performance-expert` |
| any new module, type, trait, or dependency | `rust-developer-expert` |
| `Cargo.toml` dependency changes | `rust-developer-expert` **and** `determinism-auditor` (a runtime dep in a result-affecting path is a physics change) |
| `.github/**`, `xtask/**`, `.prototools` | `determinism-auditor` |

When in doubt dispatch all six. They are cheap next to a determinism bug found
three months later.

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
