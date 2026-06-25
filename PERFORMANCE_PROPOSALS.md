# Performance & feature proposals — `particles` model

Scope: the live build is `MODEL_TYPE=lattice_particles`, which (per `src/models/CMakeLists.txt`) compiles `src/models/particles/` and
`Geometry`. So the hot path is `update_system` → `attempt_*` → `get_site_energy` / `measure_pair_energy` → `Geometry::get_neighbour`
/ `get_bond` / `get_interaction`. All proposals below target that path. The `fields`, `ising`, and `lattice_particles` directories are not
compiled and are ignored here.

---

## 1. Precompute a neighbour lookup table in `Geometry` (biggest win)

### Problem `Geometry::get_neighbour` (`src/geometry/geometry.cc:113`) recomputes lattice arithmetic on *every* call:
- `r_to_ijk` → 2 integer divisions,
- 3× `array_space::mod`, each `(a % b + b) % b` → 2 more `%` operations,
- `ijk_to_r`.

Integer division/modulo are among the slowest scalar integer ops. This function is called for every neighbour of every site touched, millions of times
per temperature step.

### Fix The neighbour relation is **static** for the whole run (the lattice never changes shape). Precompute it once at construction into a flat
array and turn `get_neighbour` into a single array read.

**`include/geometry/geometry.h`** — add private members/methods:

```cpp private: // ... existing members ... bond_struct bond_struct_m {bond_struct(chain)};

  // Precomputed neighbour table, flattened as // neighbour_table_m[site_ind * n_neighbours_m + bond_ind]. vec1i neighbour_table_m {};

  void set_lattice_properties(); // Lattice arithmetic, only used to populate neighbour_table_m. int compute_neighbour(const int site_ind, const int
  bond_ind) const; void build_neighbour_table(); ```

**`src/geometry/geometry.cc`** — rename the current body of `get_neighbour` to `compute_neighbour`, make `get_neighbour` a lookup, and build
the table at the end of `set_lattice_properties()` (both constructors already call it after `n_sites_m` is set):

```cpp int Geometry::compute_neighbour(const int site_ind, const int bond_ind) const { // ... exact current body of get_neighbour ... }

int Geometry::get_neighbour(const int site_ind, const int bond_ind) const { std::size_t idx {static_cast<std::size_t>(site_ind * n_neighbours_m
+ bond_ind)}; return neighbour_table_m[idx]; }

void Geometry::build_neighbour_table() { neighbour_table_m.resize( static_cast<std::size_t>(n_sites_m * n_neighbours_m)); for (int site {0}; site
< n_sites_m; ++site) { for (int bond {0}; bond < n_neighbours_m; ++bond) { std::size_t idx {static_cast<std::size_t>(site * n_neighbours_m + bond)};
neighbour_table_m[idx] = compute_neighbour(site, bond); } } } ```

Add `build_neighbour_table();` as the last statement of `set_lattice_properties()`.

### Cost / benefit
- Memory: `n_sites * n_neighbours` ints (e.g. 1M sites × 12 = 48 MB worst case; typically far less). Negligible for normal lattice sizes.
- Replaces ~4 integer divisions + 6 modulos per call with one load.
- Risk: very low. Pure refactor of an internal getter; behaviour identical.

---

## 2. Don't recompute `bond` that the caller already knows (algorithmic win)

### Problem `get_site_energy` (`src/models/particles/particles_interactions.cc:45`) does:

```cpp for (int bond {0}; bond < geometry.get_n_neighbours(); ++bond) { int neighbour_site {geometry.get_neighbour(site_index, bond)}; site_energy +=
get_contact_energy(state, site_index, neighbour_site, ...); } ```

`get_contact_energy` calls the 7-arg `geometry.get_interaction(..., site1, site2, ...)` overload (`geometry.h:166`), which calls
**`get_bond(site1, site2)`** — a loop over all `n_neighbours` that calls `get_neighbour` again to rediscover the `bond` we *already had* in
the loop above.

Net effect: `get_site_energy` is **O(n_neighbours²)** in `get_neighbour` calls instead of O(n_neighbours). For fcc (12 neighbours) that is ~144 vs
12 per site energy evaluation. (Proposal 1 makes each `get_neighbour` cheap, but this still multiplies the work; the two fixes compound.)

### Fix Thread the known `bond` through. Add a bond-aware overload of `get_contact_energy` and call it from `get_site_energy`.

**`include/models/particles/particles_interactions.h`** — add:

```cpp // Contact energy when the bond linking the two sites is already known. double get_contact_energy(state_struct& state, int site1, int site2,
int bond, interactions_struct& interactions, geometry_space::Geometry& geometry); ```

**`src/models/particles/particles_interactions.cc`** — implement it using the bond-taking `get_interaction` overload (`geometry.h:159`), and
have the existing `get_contact_energy` delegate via `get_bond` for callers that don't know the bond:

```cpp double get_contact_energy(state_struct& state, int site1, int site2, int bond, interactions_struct& interactions, geometry_space::Geometry&
geometry) { if (state.lattice_sites.is_empty(site1) || state.lattice_sites.is_empty(site2)) return 0.0; return geometry.get_interaction(
state.lattice_sites.get_orientation(site1), state.lattice_sites.get_type(site1), state.lattice_sites.get_orientation(site2),
state.lattice_sites.get_type(site2), bond, state.n_types, interactions.couplings); } ```

In `get_site_energy`, pass the loop `bond` directly:

```cpp site_energy += get_contact_energy( state, site_index, neighbour_site, bond, interactions, geometry); ```

Note: `measure_pair_energy` (`particles_update.cc:84`) already receives `bond` as a parameter and similarly recomputes it at line 110 — pass it
through there too.

### Cost / benefit
- Eliminates one full `get_bond` scan per neighbour in every energy evaluation (rotate, mutate, and both `get_site_energy` calls inside
  `measure_pair_energy`).
- Risk: low. Same numbers computed, just without rediscovering `bond`. Keep the no-bond `get_contact_energy` overload so other callers are
  unaffected.

---

## 3. Make `get_bond` O(1) (smaller win, optional)

### Problem After fixes 1–2, `get_bond` is only called once per swap attempt (in `measure_pair_energy` and
`attempt_rotate_and_swap_w_empty`). It still loops over `n_neighbours` calling `get_neighbour`. Minor now, but trivially removable.

### Fix `bond_struct_m.bond_index` already maps a displacement vector → bond index. Compute the displacement `(di, dj, dk)` (each wrapped to
the minimum-image range) and look it up in that map, returning `n_neighbours_m` on miss. This is already sketched in the commented-out code at
`geometry.cc:134-149`; it needs minimum-image wrapping of the displacement before the map lookup. Lower priority than 1 & 2.

---

## 4. Hoist per-call distribution objects (micro-optimization)

`pick_random_move` (`particles_update.cc:50`) and `is_move_accepted` (`particles_update.cc:277`) construct a `real_dist {0.0, 1.0}` on every
call. `std::uniform_real_distribution` construction is cheap but not free, and these are per-move-attempt hot calls. Make them `static
thread_local` (or a `constexpr` factory) so they are constructed once:

```cpp static thread_local real_dist proba_dist {0.0, 1.0}; ```

Risk: very low (the distribution is stateless for this engine). Marginal gain; do it only alongside 1 & 2.

---

## Recommended order
1. **Proposal 1** (neighbour table) — largest, lowest risk, self-contained.
2. **Proposal 2** (thread the bond) — algorithmic, compounds with 1.
3. Proposals 3–4 only if profiling still shows these as hot.

Validate after each: build with `make build`, run a short simulation and confirm the reported energies (`print_model_energy`) and saved averages
are **bit-for-bit unchanged** vs the current binary — all four changes are behaviour-preserving.

---

# Feature: write the MC acceptance (success) rate to the results files

### What "success rate" means here `update_system` (`particles_update.cc:10`) attempts `state.n_sites` moves per MC step. Each
`attempt_*` returns the energy change on acceptance and `0.0` on rejection. The acceptance rate is `accepted_attempts / total_attempts`,
optionally broken down per move type (very useful for tuning `move_probas` and diagnosing whether rotate/mutate/swap moves are effective at a given
T).

The current return-value convention is lossy: a genuinely accepted move with `delta_e == 0` is indistinguishable from a rejection. So acceptance
can't be inferred reliably from the returned energy — the `attempt_*`/`is_move_accepted` functions must report acceptance explicitly.

### Recommended implementation

**Step 1 — count acceptances at the single decision point.** Every move's accept/reject decision goes through `is_move_accepted`
(`particles_update.cc:277`). Add a lightweight counters struct and thread it (or a reference to it) into the update path. Minimal-churn option: add
the counters to `interactions_struct` or to a new `move_stats_struct` carried alongside `records`.

```cpp // new, e.g. in particles_records.h struct move_stats_struct { std::array<long, n_enum_moves> attempts {}; std::array<long, n_enum_moves>
accepts {}; void clear() { attempts.fill(0); accepts.fill(0); } }; ```

In `update_system`, record the chosen move and whether it was accepted. The cleanest signal is to have each `attempt_*` (or `is_move_accepted`)
bump the counters. Since all moves funnel through `is_move_accepted`, the lowest-churn approach is to increment `attempts[move]` in the
`update_system` switch and `accepts[move]` whenever the returned energy indicates acceptance — but because of the `delta_e == 0` ambiguity
above, prefer returning an explicit accepted flag.

Suggested signature change (small, local): have `update_system` track it directly, since it already knows `chosen_move`:

```cpp mc_moves chosen_move {pick_random_move(parameters)}; stats.attempts[chosen_move]++; double de {/* dispatch attempt_* as today */}; //
is_move_accepted is the authority; expose acceptance instead of inferring ```

To avoid the `delta_e==0` ambiguity, make `attempt_*` report acceptance. Two clean options:
- (a) return a small struct `{double de; bool accepted;}` from each `attempt_*`, or
- (b) keep returning `de` but increment `accepts[move]` inside the `else`/`if` branches of each `attempt_*` via a passed-in counter
  reference.

Option (b) is the least invasive: pass `move_stats_struct& stats` and the `chosen_move` into the `attempt_*` functions and do
`stats.accepts[move]++` in the acceptance branch (where they currently `return delta_e;`).

**Step 2 — accumulate during the averaging phase only.** In `mc::mc_simulate` (`mc_routines.cc:135`) the averaging loop is the meaningful window.
Call `stats.clear()` in `initialize_model_averages` and let `update_model_system` accumulate during the `mcs_av` loop. (Equilibration steps
should not count toward the reported rate.)

**Step 3 — write it out.** Mirror the existing `save_averages` pattern (`particles_averages.cc:31`). The averages file currently holds `{T, e_av,
e2_av}`. Append the rate(s):

```cpp // total attempts over the averaging phase = mcs_av * n_sites (per move type: // mcs_av * n_sites * move_probas[move] in expectation; use the
measured // attempts[move] as the denominator to be exact) vec1d output_vec = {T, averages.e_av, averages.e2_av}; for (int m {0}; m < n_enum_moves;
++m) { double rate {stats.attempts[m] > 0 ? static_cast<double>(stats.accepts[m]) / stats.attempts[m] : 0.0}; output_vec.push_back(rate);
} io_space::save_vector(output_vec, static_cast<int>(output_vec.size()), output_file); ```

Or, to keep the `esf_av` file format stable for existing Python analysis (`python/src/analysis/analyze_records.py` etc.), write a **separate
file** `acceptance_T_<T>.dat` with one row per move type (`move_name attempts accepts rate`). A separate file is the safer choice — it won't break
downstream parsers that assume 3 columns.

**Step 4 — gate behind an option** for consistency with the rest of the config: add `bool acceptance_option` + `std::string acceptance_output` to
`model_parameters_struct` (`particles_parameters.h:77`), read in its constructor, exactly like `e_av_option` / `e_av_output`.

### Summary of touch points
- `particles_parameters.h/.cc` — new option flag + output path.
- `particles_records.h` (or a new `move_stats.h`) — `move_stats_struct`.
- `particles_update.cc` — increment attempts/accepts (the `attempt_*` acceptance branches are the exact, unambiguous place).
- `particles_averages.cc` or `particles_records.cc` — write the rates, ideally to a dedicated file.
- `model.cc` / `mc_routines.cc` — clear stats at averaging start, save at end.

### Recommendation Use **option (b)** (pass a counter ref into `attempt_*`, increment in the acceptance branch) and write to a **dedicated
`acceptance_T_<T>.dat` file** gated behind a new option. This is unambiguous (no `delta_e==0` problem), keeps the existing `esf_av` format
intact, and follows the established option/output/save conventions.
