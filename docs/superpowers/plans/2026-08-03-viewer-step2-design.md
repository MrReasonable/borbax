# Viewer Step 2 — the periodic table panel: agreed design

Three-amigos session, 2026-08-03, two rounds. Product, developer and QA each
argued the step independently, then cross-examined. **Round 2 overturned four
round-1 positions, three of them the developer's own, and one of those was a
measurement he had reported wrongly.** That is recorded here because the
measurement is load-bearing for a test.

Supersedes the plan's Step 2 line where they disagree
(`2026-08-01-borbax-viewer.md:385-389`).

---

## The measurements everything rests on

Run against `Universe::generate` in this worktree, not recalled.

| Fact | Value | Consequence |
|---|---|---|
| Element count | **60..=120**, always | An empty or 1-element table is unreachable. Do not write that branch. |
| Periods (rows) | **2, 3 or 4** | Never assume 3. Never assume 7. |
| Groups (columns) | **0..=73** | Widest row measured **74 cells**. |
| Row lengths | seed 1: 16, 58, 6 · `emily`: 10, 34, 73 | Worst intra-universe ratio **74:1**. |
| Interior gaps | **0** over 1000 universes | No gap logic. No padding logic. |
| Duplicate `(period, group)` | **0** over 1000 universes | Placement is a bijection and is assertable. |
| `table.iter()` order | **exactly row-major**, 0 non-contiguous steps over 5000 universes | No 2D array, no sparse map. |
| Element **symbol** collisions | **0** in 2000 universes | Enforced by `taken` in `naming.rs:226`. |
| Element **name** collisions | **357 across 321 of 2000 universes (16.05%)** | First at seed 16: `"meax"` names both `Mx` and `Me`. |
| `Universe::generate` | **127 µs release, 836 µs debug** | 5.0% of a 16.7 ms frame in the profile `cargo run` gives you. |
| Non-finite property values | **0** across 180 864 elements × 6 fields | An `is_finite()` guard cannot fire and cannot be tested. |
| Worst-case grid extent | **3531 × 89 px**, independent of harness size | `Grid` columns are content-sized; there is no reflow to rescue a test. |

**The table is a left-aligned staircase with exactly one ragged edge.** Row
lengths are the shell capacities `k·n²+2` (`packing.rs:53-55`), so each row is
strictly longer than the one above until the table truncates mid-shell. Group
`g` is column `g` in every row.

**Aspect ratio is 3 × 74 where a real table is 7 × 18** — inverted, and by 4× on
the long axis. This is the largest gap between what "periodic table panel"
sounds like and what the data is. It is not a defect to be designed around; it
is what the generated physics produces, and folding it into a familiar block
would be a G6 breach (importing a real-world visual convention) *and* an
invention of structure Borbax does not have — there are no sub-shells.

---

## Scope

**In:** the grid (one cell per element at `(period, group)`, **symbol only**);
selection by click and by keyboard; a detail block below the grid; the two
`xtask` guards pulled forward from Step 7; the seam-check rework that must
precede them.

**Out, each with a reason:**

- **Colour by property** — Step 4 by name, and Step 7 assigns the
  real-chemistry-palette check to Step 4. Shipping colour here lands a palette
  ahead of its guard, which is this repository's recorded failure mode. A
  selection highlight is not a palette; it encodes no element property.
- **`abundance`** — it has no display scale, and any scale the viewer invents is
  arithmetic on a value that came out of a `borbax-*` call. (Measured: `{:.3}`
  renders `0.000` for **20.2%** of all elements, so the naive format is also
  wrong.)
- **`outer_fill_band`** — §7.1 says it is reserved and carries no chemistry.
- **Shell fill ("12 of 58")** — `packing::shell_size` is `pub(crate)`, so the
  only routes are reimplementing `k·n²+2` in the display layer or a cross-crate
  accessor. **`PeriodicTable::pattern()` must have zero call sites in
  `borbax-ui` at Step 2**, and that is checked.
- The solvent element, the depth dial, search, sort, filter, tooltips, 3D.

---

## Rulings

### The grid

Real 2D period × group grid; column index **is** `group`, row index **is**
`period`. Walk `table.iter()` and start a new row when `period` increments — no
`max`, no arithmetic, so the "is computing grid extents a physics violation?"
question does not arise.

Contiguity is a property of `generate_elements`, not a documented contract of
`PeriodicTable`. Since the panel's layout depends on it, one `borbax-ui` test
over ~200 seeds asserts `group == position_in_run`.

**Horizontal scroll. No reflow, no shrink-to-fit.** Reflow destroys group
alignment, and group alignment is where the families are — measured, group 0 →
valence 0 in every period, group 1 → valence 1, group 2 → valence 2. That
alignment is the educational payload. Shrink-to-fit puts 74 cells in ~1200 px,
i.e. 16 px each, which cannot hold a two-character symbol — and it would hide
the one thing Step 2 exists to show, that the table's *shape* differs per
universe.

**Cells show the symbol only**, against the plan's "symbol and name". Two
independent reasons: names are 5–9 characters and 74 of them per row is not a
grid; and **names collide in 16% of universes**, which would make
`egui_kittest`'s query API panic outright.

Row labels read `shell 1 / shell 2 / shell 3`, 1-based. `group` never reaches
the screen as a bare index — the word "group" imports the real table's family
meaning, which is measurably true at the left of the table and measurably false
at the right.

### Selection

**Lives inside `Outcome::Loaded { universe, selected }`, not as a
`ViewerState` field.** Product and developer reached this independently. It is
the second application of the argument in `Outcome`'s own doc: a
`selected` field beside `outcome` has four combinations and two are meaningless
— a selection in a universe that does not exist, and a selection **from a
previous universe**, which is the misattribution class `commit_typed_seed`
clears the phrase for one field along. Inside the variant, `reload`'s wholesale
`self.outcome = ...` destroys it by construction: there is no
`self.selected = None` line to forget.

`ElementId`, not an index. `select()` refuses an id this table has no slot for
rather than clamping — clamping "silently answers with a *different element's*
properties, which nothing downstream can detect" (`PeriodicTable::get`'s own
doc).

Detail block sits **below the grid, always present, fixed footprint**. Not a
tooltip (unreachable by keyboard, invisible in a screenshot, and screenshots are
how this gets shown to anyone). Not a side panel (horizontal space is the
scarcest resource on a 74-column screen). With nothing selected it shows the
same labels with `—` values, so the grid never reflows on first click.

### The nine property rows

Heading: `Mx · meax` — **symbol first**, because the symbol is what she just
clicked and it is the only identifier unique by construction.

| Label | Value | Gloss |
|---|---|---|
| `made of` | `41 base units` | |
| `shell` | `2 of 3` | |
| `outer shell` | `14 units` | |
| `bonding slots` | `4`, or `none — a closed shell` | |
| `mass` | `70.70` `mass units` | |
| `size` | `0.92` `spans` | |
| `surface` | `0.63` | how much of this element sits on the outside |
| `binding energy` | `2.60` `quanta per unit` | |
| `instability` | `0.075` | how far this sits from the most stable element in this universe |

Row shape is **four fields**, which is the reconciliation of a real conflict:

```rust
struct PropertyRow {
    label: &'static str,
    value: String,
    unit: &'static str,          // "" when dimensionless
    gloss: Option<&'static str>,
}
```

Product wanted the unit inside the value (`"0.92 spans"`) because that reads
better to a child. QA showed that destroys the guard: an allow-list would have
to parse sentences, and a unit assembled by `format!` would be indistinguishable
from a declared one. Four fields give product the reading (`panel.rs` emits
`value` and `unit` adjacently) **and** QA a set assertion over the attained
`unit` values. Glosses sit directly under the row they explain — a footnote
adjacent to nothing is a footnote nobody reads, and a fixed footnote position
reintroduces the index-into-order coupling `labelled_by` removed in Step 1b.

**Forbidden wording, and neither `xtask` guard can see it** (these are not
element names and not units): `electronegativity`, `half-life`, `radioactive`,
`decays at`, and anything with `per world-year`. §7.1 calls `instability`
"probability per world-year" while Task 15 consumes it as an unbounded Gillespie
propensity — a recorded, unresolved contradiction, and a UI that writes a rate
invents the resolution of an open physics question in the crate least entitled
to settle one. Held by `no_explanatory_line_borrows_a_word_from_real_chemistry`,
a unit test over the closed set of `&'static str` constants.

**`affinity` is `surface`, and there is no −1…+1 bar.** Measured range is
[0.140476, 1.0]; a bar on the documented scale puts every element in the top
half and teaches that this universe has no "negative" elements when the scale is
simply not reached. **`radius` is `size`, and no size trend is drawn** — it
rises monotonically where a real atomic radius falls, so a "radius across the
period" strip would manufacture exactly the comparison G3 exists to prevent.

**Formatting.** `{:.3}` fixed, never bare `{}` (shortest-round-trip `Display`
gives ragged widths). `Mass::to_f64()`, never `{:?}` — measured, `Debug` on
`Mass` prints raw 1/1024 sub-units, **off by 1024×**. No `is_finite()` guard:
zero non-finite values across 180 864 elements × 6 fields, so it is a branch
that cannot fire and therefore cannot be tested.

---

## The guards

### Seam check rework — lands FIRST, before any new file

Three holes, all measured on the shipped check
(`xtask/src/main.rs:2729-2816`):

1. It hard-codes exactly two filenames, so **any new `.rs` file is checked for
   nothing**.
2. **"Only `main.rs` names `eframe`" is enforced by nobody.** It holds today
   only because there are three files and two of them ban the string.
3. `if !viewer_src.exists() { return Ok(()) }` at line 2732 — **rename the crate
   and the whole seam check passes silently.** Found independently of the
   amigos; it is the same fail-open shape as the `crates/` rename that blinded
   `check_every_member_inherits_the_lints`.

The fix is a composition, and the developer conceded the important half: an
enumerated table alone is the hand-kept-list failure mode moved *into* the entry
(someone writes `("periodic.rs", &[])` and the hole is one word wide with a
green build), while a directory rule alone misses `state.rs` acquiring `egui`,
because it just changes category.

- **Tier 0** (default, any unnamed `.rs`): banned `["egui", "eframe"]`.
- **Tier 1** (named drawing files): banned `["eframe", "format!"]`.
- **Tier 2** (exactly `main.rs`): banned `[]`.
- **The count:** exactly one file names `eframe` in code, and it is `main.rs`.
- Missing directory ⇒ **failure**, not skip.

Plus `PeriodicTable::pattern()` at zero call sites. Its message must say "if you
are adding one, delete this check and say what it shows" — `peak` and `closures`
are legitimate things a later step will want, and a check whose only correct
response is deletion must say so or it gets worked around.

### Guard A — real element names and symbols

Lexed with `proc-macro2` like G5, **string literals only**, `#[doc]` bracket
groups skipped, pre-`#[cfg(test)]` cut. Vocabulary **duplicated in `xtask`**,
never imported: `xtask` depending on a chemistry crate to enforce §5 inverts the
gate, and a gate that goes down with the thing it guards is worse than one that
reports red. Held honest by a cross-check test in both directions with a floor
on each count — the `clippy.toml` ↔ `BANNED_CALLS` pattern.

**The rule is a split, and the split is measured**, not chosen:

| Tier | Match | Scope | False positives on this tree |
|---|---|---|---|
| Real **names** | whole segment, case-**in**sensitive | `crates/` | **0** |
| **Two-letter** symbols | whole segment, case-**sensitive** | `crates/` | **0** |
| **One-letter** symbols | whole **literal** | `crates/` | **0** |
| *(rejected)* any symbol, case-insensitive | | | **442** |

The case asymmetry contradicts the `elements.JSON` lesson at first glance and
**the numbers must be in the doc comment or a reviewer will "fix" it into 442
false positives**. It does not transfer: `elements.JSON` is a *filesystem*
lesson — on macOS and Windows those are the same file, so the rename is
invisible to whoever makes it. Two string literals differing in case are two
different tokens. And case is semantic in the symbol vocabulary (`Fe` is iron;
`fe`, `nO`, `aT`, `iN` are ordinary words — which is what the 442 counts) and
meaningless in the name vocabulary.

`experiments/` is a **declared exclusion** with a written reason, spelled the
way `UNSCANNED_CRATE_DIRS` declares `xtask` — it is measurement probes, ships to
nobody, prints to no child. It holds the only two hits (`"N"` in
`bin/{fusion,ptable}.rs`). Scope stays `crates/` rather than `borbax-ui` alone,
because a hard-coded `"C"` in `borbax-universe` is the G2 breach *at its
source*; the viewer is merely where it becomes visible.

Residual, stated so nobody thinks it is total: a hard-coded `"FE"` or `"fe"`
escapes. That is the spelling nobody writes.

### Guard B — real-world units

Same lexing and scope. Vocabulary: `Å`, `angstrom`, `kJ/mol`, `kcal`, `g/mol`,
`eV`, `nm`, `amu`, `dalton`, `kelvin`, `celsius`, `fahrenheit`, `joule`,
`newton`, `pascal`, `mole`, `gram`, `kilogram`, `metre`/`meter`, `litre`,
`hertz`, `°C`, `°F`, `µm`. Case-insensitive — no unit token collides with an
identifier in the tree.

**Not on the list, measured:** `second` (3 hits — `"{first} {second}"`, `"a
second layout pass"`), `pm` and `fm` (2 and 1 — and they are also *element
symbols*, so they collide with Guard A's fixtures), `bar`, `da`, bare `mol`.
And **`world-years` contains `year`**, so `year` must never be a bare entry.

**The `borbax-units` false positive, and why neither offered fix was taken.**
`unit!(Thermal, "Temperature, in thermals. Not Kelvin, not Celsius.")` puts four
G4 words in **bare string literals inside a parenthesised macro group** — the
`#[doc]` skip is a *Bracket* test (verified at `xtask:3090-3096`), so it does
not reach them, and neither does comment stripping. The crate that exists to
disclaim those units would fail the guard that enforces them.

- A **path exemption** is a hand-kept list of one, against the recorded *unknown
  ⇒ fail* rule.
- A **`Not `-prefix disclaimer rule** is worse than it looks: a guard satisfiable
  by editing the string it fired on trains people to edit strings, and the first
  thing anyone does when `"Å"` fires is try `"Not Å"`.

**Taken instead: change the `unit!` macro arm to `$(#[$doc:meta])*` and spell
the call sites as `///` comments.** The words then arrive as `#[doc = "…"]`
bracket groups and the **existing** skip covers them, with zero new exemption
surface. It *removes* an exemption rather than adding one. Cost is a three-line
macro-arm change plus four call sites in a physics crate; no formula moves, no
golden moves, and rendered rustdoc is unchanged. (The expansion is not quite
byte-identical — `///x` desugars to `#[doc = " x"]`, gaining a leading space.
Cosmetic, and stated rather than glossed.)

**Second leg, because a deny-list over an open vocabulary cannot invert:** no
unit newtype may gain a `Display` impl. Verified: none exists today. With no
`Display`, every unit word on screen is a literal in a crate the scan covers.

**Third leg, and it is the only fail-on-unknown half there can be:** the
acceptance test `every_unit_word_on_screen_is_one_this_universe_invented`, an
allow-list over the **attained set** of the `unit` field across every element of
many universes. It catches a real unit *assembled at runtime*, which no source
scan can see, and a *missing* unit word — a bare number attributed to nothing.
It must assert on the attained set, never against the declared constant: that is
a guard whose two sides are the same expression.

Allow-list, exactly: `""`, `base units`, `units`, `mass units`, `spans`,
`quanta per unit`, `bonding slots`. (**Seven, not the six this paragraph first
said** — the `made of` row carries `base units`, and the drift was caught by the
test failing on the first run rather than by the doc being re-read.) Plus two negative assertions — `thermals` and `world-years`
must **not** appear, since no Step 2 property carries a temperature or a time.

---

## The test list

Every test names the wrong implementation it kills that nothing else kills. A
test with no unique discriminator is cut, and the cuts are recorded.

### `state.rs` unit tests

- `every_element_sits_at_its_own_period_and_group_and_no_cell_holds_two` — kills
  `row = i / 18, col = i % 18` and every wrap-based layout that *looks* like a
  table. Both directions: each element reaches its cell, and the grid holds no
  id the table does not (the reverse arm kills padding cells carrying
  `ElementId::ZERO`).
- `the_grid_is_as_wide_and_as_tall_as_this_universe_needs` — kills an
  **over**-sized constant container. Placement alone cannot see this: a 7×18
  grid with correct placement fails only by being too narrow for group 73.
- `two_universes_do_not_have_to_agree_on_the_shape_of_the_table` — kills a shape
  computed once (`OnceLock`, a `static`, built in `new()`). **Assert on the
  `(rows, cols)` pair, never rows alone** — row counts collapse to `{2,3,4}` and
  `P(64 seeds all give 3 rows) = 5.5e-5`, a 1-in-18 000 flake.
- `a_selection_cannot_outlive_the_universe_it_was_made_in` — kills re-lifting
  `selected` to a `ViewerState` field (what someone does the first hour they
  want "keep the same element when I change the seed") **and**
  `Rejected { last_selected }`. Assert both paths: reload-to-valid and
  reload-to-rejected.
- `an_element_id_this_universe_has_no_slot_for_selects_nothing` — kills
  `table.get(id).unwrap()` and its silent cousin, clamping to the last slot.
- `nothing_is_selected_until_something_is_selected` — kills
  `selected = Some(ElementId::ZERO)` as a default, which opens the window
  painting element 0's properties under a heading nobody clicked.
- `the_cell_shows_the_symbol_and_the_heading_shows_the_name` — kills
  transposition at the display layer. Both are `String`, so
  `format!("{} {}", e.name, e.symbol)` compiles.
- `no_explanatory_line_borrows_a_word_from_real_chemistry` — the only thing
  standing between "surface, unlike electronegativity…" and a child's screen.

**Cut:** `every_element_in_the_universe_reaches_the_grid` (folded in as the
reverse arm above); `an_empty_table_paints_no_rows` (unreachable input —
`generate_elements` draws 60..=120 and `PeriodicTable::new` is `pub(crate)`).

### `tests/acceptance.rs`

- `the_properties_shown_are_the_selected_elements_and_not_a_neighbours` — kills
  an off-by-one at the read. Derive over **every** element of several universes.
- `every_unit_word_on_screen_is_one_this_universe_invented` — see Guard B's
  third leg.
- `no_element_reaches_the_screen_under_a_real_name_or_symbol` — kills
  display-side *derivation*: `&name[..2]`, `name.to_uppercase()`,
  `symbol.chars().next()` for a narrow cell. "helion" → `"He"` is one edit away,
  and **neither `xtask` guard can reach it** — there is no literal to scan.
- `selecting_an_element_does_not_regenerate_the_universe` — kills `select()`
  calling `reload()`, the plausible "refresh the view" implementation.
- `the_same_universe_gives_the_same_table_however_the_viewer_got_there` — kills
  a grid `push`ed to rather than rebuilt.

**Cut:** `an_element_selected_by_name_is_the_same_element_selected_by_seed` —
`the_boxes_agree` already pins name/seed identity and selection is downstream.

### `tests/panel.rs` — frame tests

Harness at **3600 × 400** for selection tests; real pointer `click()`
throughout. The far cell is picked **from `rows()`** (last cell of the widest
row), never a hard-coded label, and the test **asserts its bounds lie inside the
screen rect before clicking** — otherwise the day a symbol renders wider, the
failure reads "selection did not change" and sends someone hunting in
`select()`.

- `the_selected_elements_properties_appear_in_the_frame_the_cell_is_clicked` —
  this file's one defect, in Step 2's widget. Mutation: move the properties
  block above the grid.
- `clicking_a_cell_marks_that_cell_and_only_that_cell_as_chosen` — kills a cell
  drawn as `ui.label` inside a clickable region: works when a human clicks,
  carries `Role::Label`, no `Action::Click`, no `Toggled` — invisible to a
  screen reader and to every test.
- `the_grid_paints_as_many_rows_and_columns_as_this_universe_has` — kills
  `for p in 0..7 { for g in 0..18 }` in the paint body while the state is
  correct.
- `an_idle_table_does_not_rebuild_the_universe` — the twin of
  `an_unchanged_name_box_does_not_rebuild_the_universe`, whose own doc records
  that the state-level version failed **zero** tests because it never went
  through `draw`.
- `changing_the_seed_repaints_the_table_in_the_same_frame` — kills the grid
  built from a snapshot taken at the top of `draw`.
- `selecting_an_element_does_not_cost_a_second_layout_pass` — kills an
  **ungated** `request_discard` anywhere in the selection path. This is Step
  1b's recorded incident pre-empted: an ungated discard re-paints from updated
  state and *self-repairs* the frame-order defects the three tests above exist
  to catch — measured then, 4 failures collapsed to 1. Assert
  `num_completed_passes == 1`.
- `two_cells_never_answer_to_the_same_label` — kills labelling cells by group or
  period, which makes `get_by_label` ambiguous and every other frame test
  silently address the wrong cell.
- `every_cell_is_reachable_at_the_size_the_window_actually_opens` — harness at
  the **real** window size. This is the only test that means anything about the
  shipped program, and it forces the panel to ship a scroll rather than assume
  width. **The size is a `pub const` in `lib.rs` read by both `main.rs` and the
  test** — a retyped `800.0` is a number someone guessed.
- `the_grid_is_operable_from_the_keyboard` — a *separate concern* from
  reachability, and it is the one test that legitimately uses `focus()`.

**Forbidden by review, not by test** (and each must be named in a doc comment,
because nothing fails if someone writes them): `harness.run()` in place of
`step()`, which hides frame lag by definition; `click_accesskit()`, which is
documented to click invisible widgets and so makes a reachability test vacuous
by construction; `Harness::fit_contents()`, which calls `run_ok()` internally.

### `xtask` self-tests

Guard A: `a_hardcoded_real_element_name_in_a_display_string_is_caught`,
`a_hardcoded_real_symbol_is_caught_in_the_spellings_a_person_writes`,
`prose_may_state_the_prohibition`,
**`the_blocklist_the_guard_uses_is_the_blocklist_the_generator_uses`** (the
highest-value one — kills the duplicate silently narrowing while looking
maintained),
`the_generators_own_blocklist_file_is_the_one_exemption_and_it_must_exist`,
`a_case_flipped_real_name_is_still_caught_and_a_case_flipped_identifier_is_not`,
`code_behind_an_inactive_cfg_is_still_scanned`,
`a_file_that_cannot_be_lexed_is_reported_rather_than_skipped`,
`a_missing_scan_root_is_reported_rather_than_skipped`.

Guard B: `a_real_world_unit_in_a_label_is_caught_in_every_spelling_that_matters`,
**`the_disclaimers_are_doc_comments_and_a_bare_literal_is_not_exempt`** (fixture
carries both spellings side by side; under a path exemption or a disclaimer rule
the bare-literal line is silent, under the taken fix both fire),
`a_format_placeholder_is_not_a_unit`,
`a_new_invented_unit_type_fails_this_check_until_someone_names_its_word`.

Seam: `a_new_viewer_file_is_covered_without_being_named`,
`a_second_file_naming_eframe_is_reported`,
`a_missing_viewer_src_is_reported_rather_than_skipped`,
`a_pattern_call_in_the_viewer_is_reported`.

---

## Mutation probes

`--no-fail-fast`, count **named** tests, never the process exit code. Delete
**one line** of a guard, never the whole thing.

| One-line edit | Must fail |
|---|---|
| `table.get(id)` → `.unwrap_or(last)` | `an_element_id_this_universe_has_no_slot_for_selects_nothing` |
| `selected: None` → `Some(ElementId::ZERO)` | `nothing_is_selected_until_something_is_selected` |
| lift `selected` to a `ViewerState` field | `a_selection_cannot_outlive_the_universe_it_was_made_in` |
| `grid[(period, group)]` → `grid[i/18][i%18]` | `every_element_sits_at_its_own_period_and_group…` |
| `for p in 0..rows()` → `for p in 0..7` | `the_grid_paints_as_many_rows_and_columns…` |
| `selectable_value` → `ui.label` in a clickable frame | `clicking_a_cell_marks_that_cell_and_only_that_cell…` |
| label cells by `group` instead of `symbol` | `two_cells_never_answer_to_the_same_label` |
| move the properties block above the grid | `the_selected_elements_properties_appear_in_the_frame…` |
| add an unconditional `request_discard` in the selection path | `selecting_an_element_does_not_cost_a_second_layout_pass` |
| append `" Å"` to the size row | acceptance allow-list **and** Guard B |
| build the size unit as `format!("{}{}", "Å", "")` | the acceptance test **only** — Guard B is blind, which is why both exist |
| show `&e.name[..2]` in the cell | `no_element_reaches_the_screen_under_a_real_name_or_symbol` |
| delete one entry from `xtask`'s symbol list | `the_blocklist_the_guard_uses_is_the_blocklist…` |
| add `.to_ascii_lowercase()` to the symbol half | its negative arm — **and 442 hits on the real tree** |
| delete the `naming.rs` exemption line | 118+95 hits on the real tree |
| delete the `#[cfg(test)]` cut | fires on `"U-7F3A21C9"` at `state.rs:438` |
| add `second` or `pm` to Guard B's vocabulary | `a_format_placeholder_is_not_a_unit` |
| revert one `unit!` call site to a bare literal | `the_disclaimers_are_doc_comments_and_a_bare_literal_is_not_exempt` |
| add `("periodic.rs", &[])` to the seam table | the `eframe` count |
| `state.rs` acquires `egui` | Tier 0 |

Also, **after Step 2 lands**: move `ui.label(status_line)` to the top of `draw`
and record the failure count. If it *drops* from Step 1b's four, a
`request_discard` has been ungated somewhere.

---

## Shipped untested, stated rather than hidden

1. **Geometry.** Nothing asserts cells are *arranged* as a grid on screen. A
   table rendered as one long column with correct `(period, group)` metadata
   passes every test here. Pixel-position assertions were rejected: they pin a
   layout someone guessed.
2. **Per-frame layout cost** of emitting ~120 widgets. No benchmark; a timing
   assertion would be flaky.
3. **Readability.** A symbol visually truncated by cell width is invisible to
   the accessibility tree, which carries the full label.
4. **Colour** — Step 2 ships none deliberately; the CPK-palette breach is Step
   4's guard.
5. The `#[cfg(test)]` positional convention stays fail-**open** for a literal
   written *below* a test module.
6. **Guard A cannot see a real name assembled at runtime**; **Guard B cannot see
   a unit outside its vocabulary.** The acceptance allow-list closes only the
   second.
7. The alias hole is unchanged: `use borbax_universe::Universe as U;
   U::generate(..)` still passes the seam check.
8. Window creation, the event loop and GPU init remain untested, as at Step 1.

---

## Carried out of the session

- **The `Universe::generate` figure is wrong in six places.** They quote 146 µs
  (release) while making a **frame-budget** argument, and `cargo run` without
  `--release` is what Ian uses. Corrected to "127 µs release, 836 µs debug" at
  `state.rs:318`, `panel.rs:136`, `tests/panel.rs:84` and `:229`,
  `xtask:2770` and `:2809`. **`element.rs:321` is deliberately left** — it is
  about the cost of an `assert!` inside generation and has nothing to do with
  frames; harmonising it would read a frame argument into a crate that has none.
- **Element names collide in 16.05% of universes** and G2 does not require
  uniqueness. Harmless at Step 2 because only one name is on screen at a time,
  but a later step with a search box, or any formula display, will assume names
  are unique. Recorded as known-open.
- **`ScrollArea` does not cull**, so all cells stay in the accessibility tree.
  `show_rows`/`show_viewport` would remove off-screen cells from the tree, at
  which point no frame test and no screen reader can address them. Not taken at
  117 widgets; the tie-breaker if anyone wants it is a measurement of egui's
  per-frame layout cost, which nobody has run.

---

## Corrections found during implementation

Recorded here rather than silently applied, because in each case the design was
wrong and the tree said so.

- **`Na` inside `NaN` is not sodium.** The symbol tier was specified as
  whole-*segment* matching over `identifier_segments`, which also splits
  camelCase — so `NaN` becomes `["Na", "N"]`. Six real false positives, all §13.4
  float comments. Symbols now split on whole words only, and the residual
  (`"NaCl"` scores nothing) is written down.
- **`cArBoN` escaped the name tier** for the mirror-image reason: camelCase
  splitting shatters it into `c|Ar|Bo|N`. Names now match the union of both
  splittings.
- **The unit list missed `"\u{c5}"` written as an escape**, which renders as an
  angstrom sign but contains no such character for a substring match. Both
  codepoints are now listed in both spellings.
- **A `#[cfg(test)]` skip was set but never consulted**, so the module body it
  was meant to skip was recursed into anyway. Nothing on the tree fired — which
  is exactly the shape of a guard that is wrong and looks right.
- **`min_col_width(0.0)` would have shipped a real bug.** Tightening the grid let
  cells overflow their columns and overlap, so a click landed on the *neighbouring*
  element — silent, plausible, and caught only because three frame tests failed
  at once. Dropped; the spacing change alone gives the same benefit.
- **The reachability test asserts geometry, not a click, and that is forced.** A
  pointer `click()` does not reach widgets inside a *scrollable* `ScrollArea`
  under `egui_kittest`: measured, it works at 3600×900 (content fits) and fails
  at 900×600, 900×900 and 1400×600 with an identical widget rect. Not input
  timing (a hover frame, and one to four warm-up frames, change nothing) and not
  drag-to-scroll (excluding `ScrollSource::DRAG` changes nothing). The widget is
  correctly wired — `click_accesskit()` selects it, and a pointer click on
  "surprise me" outside the scroll area works at the same size. So the app is
  fine and the harness cannot express that one interaction.
- **The grid was unreadable and every test was green over it.** `Grid`'s default
  column spacing is sized for prose, so each two-character symbol got a ~55 px
  column: at 900 px only **eight** cells of a 58-cell row were on screen, and the
  one thing this step exists to show — that the shape changes with the seed — was
  invisible behind a scroll nobody would think to drag. Caught by *looking at the
  window*, which is the same way Step 1's opening-state defect was caught. With
  3 px spacing and a 1200 px window: seed 1 shows 16/16, 27/58, 6/6 and seed 7
  shows 8/8, 26/26, 27/51 — two of three rows complete in both.
- **`WINDOW_SIZE` is now a `pub const` in `lib.rs`**, read by both `main.rs` and
  the reachability test, for the reason the design gave: a retyped literal is a
  number somebody guessed.

**Step 1b's guard count held.** Moving `ui.label(state.status_line())` to the top
of `draw` still fails exactly **four** frame tests, as it did before this step.
That was checked because the design predicted the risk: a new `request_discard`
that is not gated self-repairs those defects and takes the count down.
