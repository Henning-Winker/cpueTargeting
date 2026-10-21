# Multispecies fleet system with fisher behaviour (dev-mmsy prototype)

Trial code for a multispecies effort-management simulator. It extends the shared-q
Schaefer model of the FAO RAP multispecies app (`mse.SEA`) with fleet segments,
métiers and a preference (expected-profit) function for effort allocation. The
system is conditioned on the landscape of `simulate_logbook()`. It sits outside the
package (`^dev$` in `.Rbuildignore`) and may move to its own repository.

Run the scripts from the repository root, with cpueTargeting installed.

| File | Content |
|---|---|
| `R/derive_landscape.R` | Rebuilds the `simulate_logbook()` landscape: carrying capacity by cell, depth, substrate and SST |
| `R/msfleet.R` | `ms_system()`, `allocate_effort()`, `ms_simulate()`, `ms_summary()` |
| `scripts/01_scan_global.R` | Total effort split equally across the segments, with fixed and behavioural allocation |
| `scripts/02_segment_grid.R`, `02b_report_grid.R` | 3D grid of segment caps (9,261 scenarios) and constrained optima |
| `scripts/03_checks.R` | Validation checks, sensitivity to γ, and a higher Dentex price |
| `scripts/04_hcr_demo.R` | Stochastic run with a segment-specific hockey-stick HCR |

## System

**Stocks.** One Schaefer pool per species, with r and K = B0 from `logbook_settings()`.

* Growth is annual and acts on the biomass at the start of the year.
* Catches are taken monthly as C = B(1 − e^−F).
* Seasonal availability of the pelagic species follows the logbook settings.
* Process error is optional: lognormal, with AR(1) autocorrelation.

**Métiers.** A métier is a tactic (Reef, Sediment, Pelagic) fished on a ground (shallow, or deep below 300 m, the depth of the logbook closure). That gives six métiers. The catchability of species i in métier m is

    q[m,i] = tac[tactic, i] · Σ_c w_mc K_ic / K_i

* `tac[tactic, i]` is the gear selectivity from the logbook settings.
* `K_ic / K_i` is species i's share of its carrying capacity in cell c, summed over the ground's cells.
* `w_mc` is how the métier's sets spread over those cells. It is uniform by default (`preference = 0`).

**Segments.** Demersal, pelagic and mixed, each a third of the fleet by default.

* Loyalty is 0.8 for the specialists and 0 for the mixed segment.
* The loyal share of a segment's effort stays on its home métiers.
* The rest chooses among all open métiers.

**Preference function.** Each month, segment g allocates its effort over métiers m as

    P(m) ∝ open_m · max(π_m, 0)^γ,   π_m = Σ_i p_i q_mi B_i avail_i − c_m

* π_m is the expected profit per unit of effort.
* Prices p are fixed, so the price ratios do not change.
* The cost per unit of effort c_m is `cost_share` × the mean unfished revenue rate, multiplied by 1.5 on the deep ground.
* γ = 0 means profit-blind allocation: equal shares among the allowed métiers.
* `behaviour = "fixed"` freezes the shares at the unfished allocation. This is the classical shared-q model.
* `info = "lagged"` uses last year's realised revenue per unit of effort instead of the current biomass.

**Management levers.**

* **Caps:** annual effort caps by segment (an n × G matrix, vectorised over scenarios and simulations).
* **Closures:** métier closures (`open`, a vector or a function of year).
* **Entry and exit:** optional, below the cap (`phi`).
* **HCR:** an optional `hcr(state)` returning multipliers on the caps.

**Reference points.** Single-species Bmsy = K/2, Bpa = 0.5 Bmsy and Blim = 0.25 Bmsy.

## Checks (seed 1 landscape)

* **Nesting.** At γ = 0, fixed and dynamic allocation give identical results (difference 0).
* **Analytic equilibrium.** The simulated B/K matches the analytic Schaefer B/K = 1 − F/r to within 0.005. The small gap comes from the monthly Baranov catch.
* **Speed:**
  * 9,261 equilibria (150 years × 12 months each) take about 50 s;
  * 200 stochastic simulations × 50 years take under 1 s.

## First findings

1. **Dentex is the choke species.**
   * Its habitat is narrow (the shallow reef only), which gives it the highest catchability.
   * It sets the effort limit everywhere.
   * Unconstrained MMSY sacrifices 3–4 species: yield is about 1,480 (fixed allocation) to 1,660 (dynamic), with several species driven to near zero.
   * With every species held at or above Bpa, yield falls to about 240–280.
2. **Segment caps release part of the choke.** Best yield with every species at or above Bpa (γ = 1.5, dynamic allocation):

   | Cap structure | Yield |
   |---|---|
   | Equal 1/3 split | 248 |
   | Each segment ≥ 20% of effort | 391 |
   | Unrestricted | 822 |

   The unrestricted optimum closes the demersal segment. Segment viability therefore has to be an explicit constraint (social or economic).
3. **Profit-seeking makes the choke worse.** Constrained yield at an equal split falls as γ rises (it is 248 at γ = 1.5):

   | γ | Constrained yield |
   |---|---|
   | 0 | 513 |
   | 0.5 | 324 |
   | 3 | 192 |

   Fishers concentrate on the valuable reef métier. Profit per unit of effort rises while yield falls. A higher Dentex price (12 → 30) lowers constrained yield further.
4. **Feedback control.** With caps at twice the static optimum and no HCR, the choke species collapses. A segment-specific hockey-stick prevents this:
   * each segment responds to the most depleted species that makes up ≥ 5% of its unfished catch;
   * Blim → Bpa;
   * result: P(B < Blim) = 0, mean yield 581.

## Design decisions still open

* **Biomass pools.** Use one pool per species, or sub-pools by ground with exchange? Sub-pools would add local depletion and hyperstability.
* **Within-ground effort concentration.** It is uniform by default. Freezing the unfished value^1.5 distribution makes the catchabilities span about 12 orders of magnitude.
* **Costs.** The level and structure of the cost per unit of effort, and whether costs differ by segment.
* **Expectations and inertia.** Current biomass or lagged CPUE, and whether to add an inertia term on last year's shares.
* **Constraints.** Every species at or above Bpa, or a rule such as 2/3 of species; and segment viability (minimum share or minimum profit).
* **Optimiser.** Replace the coarse log grid with `optim` or a finer grid near the optimum.
* **Fleet dynamics.** Effort creep and entry and exit below the cap (`phi > 0`).
