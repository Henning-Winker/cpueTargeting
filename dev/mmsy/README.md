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
| `scripts/06_optimise.R` | Constrained MMSY by lever set; demersal and pelagic caps with the mixed segment as a free agent; leakage |
| `scripts/07_viability_closure.R` | Constrained MMSY with a segment viability floor and a closure of the shallow reef métier |
| `scripts/08_fig_options.R` | Summary figure `figs/opt_options.png` |
| `scripts/11_feedback.R` | Feedback control from the tuned status quo: segment effort rule (trigger 0.8 Bmsy), choke-stock closure, change limits; `figs/fb_*.png` |
| `scripts/12_tune_trigger.R` | `ms_tune_hcr()`: tunes the effort and closure triggers to the risk limit; `figs/tune_trigger*.png` |
| `scripts/13_assessed_subset.R` | Management based on four assessed stocks |
| `scripts/14_yield_curves.R` | Equilibrium yield over effort for four conceptual treatments; `figs/yield_curves.png` |
| `scripts/09_tune.R` | `ms_tune()`: daily costs that make 35/35/30% of active days an open-access equilibrium, for two status targets |
| `R/plots.R`, `scripts/05_plots.R` | Figures of stock status, F, effort and yield, saved to `figs/`. Covers equilibrium against effort, Kobe plots of the status quo against the safeguarded optima, and projections from the status quo. Uses the balanced defaults |

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

## Equilibrium yield curves: conceptual treatments (`scripts/14_yield_curves.R`)

Yield over total days at the status-quo segment mix (tuned system, deterministic). Caps are scaled from the status-quo days, and the days fished are those of the vessels that break even.

| Treatment | MMSY: days, yield | Safeguarded (all ≥ Bpa): days, yield |
|---|---|---|
| Shared q (mse.SEA, analytic, calibrated to the status quo) | 1,157, 1,811 | the same (weakest 0.66) |
| Métiers, profit-blind (γ = 0) | 790, 1,628 (3 stocks < Bpa) | 212, 566 |
| Métiers, profit-seeking (γ = 1.5) | 1,993, 1,757 (2 stocks < Bpa) | 352, 945 |
| Profit-seeking + shallow reef closed | 1,952, 1,725 | the same |

For reference, the sum of single-species MSY is 1,884, and the status quo is 1,536 days with yield 1,745.

* **Shared q gives a dome with a false sense of safety.** All stocks are above Bpa at MMSY, and beyond it all of them collapse together.
* **Behaviour turns the dome into a plateau.** Fishers move away from depleted stocks, so yield is insensitive to effort over a wide range while the weak stocks erode underneath. Yield alone is a poor indicator.
* **The gap between MMSY and the safeguarded yield measures the technical interaction.** A métier measure (here the reef closure) closes it.

## Feedback control (`hcr_effort()`, `scripts/11_feedback.R`, `scripts/12_tune_trigger.R`)

**Uncertainty in the projections.**

* **Process error, JABBA form:** B_t+1 = (B_t + rB_t(1 − B_t/K) − C_t) · exp(η_t), with AR(1) 0.3 in time.
  * The SD of η for species i is 0.1 × r_i / median r (`proc_r_power = 1`), so more productive species are more variable. The range is 0.06 (blackspot seabream) to 0.25 (dolphinfish, chub mackerel).
  * An optional correlation across species (`proc_rho_sp`) is off by default.
  * A CV of 0.2 is far too large in this form. With CV 0.2 the stocks crash even under the rules: P(any stock < Blim) is 0.67 without control.
* **Assessment error:** lognormal CV 0.2 on the B/Bmsy the rule sees.
* **Implementation error:** lognormal CV 0.1 on the effort caps by segment (`impl_cv`, `impl_bias`).
* **Data lag:** `hcr_lag` = 2 by default. The rule sees the stocks as they were 2 years before the year its advice applies (data year y − 2, assessment y − 1, advice y). Before the projection, the data years are taken to be the status quo.

**The rule (`hcr_effort()`).**

* Each segment's days are a fraction of its status-quo days. The fraction ramps from 0 at Blim (0.25 Bmsy) to 1 at the trigger, applied to the most depleted stock the segment can catch substantially (catchability ≥ 20% of the stock's highest).
* **Optional closure.** A choke stock in `close` is protected by a closure instead. The open share (area or season) of the métiers that drive it (catchability ≥ 50% of its highest; for Dentex, the shallow reef) ramps with its own trigger (`Btrig_close`), and the stock is left out of the effort rule.
* `max_change` limits the annual change.
* A rule may return `list(mult, open)`. A partly closed métier fishes only its open share.

**Fixed trigger 0.8 Bmsy, 2-year lag.** 200 simulations; performance over years 31–50. The deterministic status quo has yield 1,745 and profit 1,489.

| Run | Yield | Profit | Days | P(any < Blim) | Worst stock P(B < Bpa) |
|---|---|---|---|---|---|
| No control | 1,667 | 1,500 | 1,529 | 0.066 | 0.92 (Dentex) |
| Effort rule | 1,463 | 2,621 | 695 | 0.004 | 0.57 (Dentex) |
| Effort rule + reef closure | 1,547 | 2,527 | 830 | 0.004 | 0.11 (Dentex) |
| Effort rule + closure, ±20% per year | 1,567 | 2,503 | 856 | 0.004 | 0.11 (Dentex) |
| Effort rule + closure, no assessment error | 1,610 | 2,390 | 937 | 0 | 0.15 (Dentex) |

**The trigger as a tuning parameter (`ms_tune_hcr()`).**

* A grid of triggers runs in one vectorised call with common random numbers (`crn`): every trigger sees the same process, assessment and implementation errors.
* The tuned rule has the highest yield among the triggers that meet the risk limit: every stock below Bpa at most 5% of the time (years 31–50).

| Rule | Tuned trigger(s), no lag | Tuned trigger(s), 2-year lag (default) | Yield (lag 2) | Profit (lag 2) | Days (lag 2) |
|---|---|---|---|---|---|
| Effort rule only | none up to 1.5 Bmsy | none (1.5 Bmsy gives 20%) | — | — | — |
| Effort rule + closure, one shared trigger | 0.85 (yield 1,538) | 1.00 | 1,467 | 2,673 | 723 |
| Separate triggers | effort 0.6, closure 3.0 (yield 1,606) | effort 0.6, closure 3.0 | 1,603 | 2,206 | 1,022 |

* **The lag hurts the shared-trigger rule but barely touches separate triggers.** The shared trigger has to rise from 0.85 to 1.00 Bmsy, costing about 5% of yield. With separate triggers the same values still pass (risk 3.0% → 4.4%), because the high closure trigger already gives Dentex a margin.

* **The closure decouples the choke species.** With separate triggers, the effort rule can be lenient (0.6), while the reef is partly closed until Dentex is well above Bmsy. A permanent full closure is close behind (1,572).
* **Yield is about 4% below no control, with profit up 47% and the risk limit met.** With one shared trigger, yield is lower and profit higher.

### Management based on four assessed stocks (`scripts/13_assessed_subset.R`)

Assessing every stock is often not feasible. `hcr_effort(assessed = ...)` lets the rule see only the assessed stocks; a closure needs its stock to be assessed.

* Tuning uses the risk to the assessed stocks only (`risk = c(worst_P_assessed = 0.05)`), which is all a manager could check.
* Performance reports the risk to all stocks.
* All set-ups use the 2-year lag and tuned triggers.

| Set-up | Tuned triggers | Yield | Profit | Days | Risk, assessed stocks | Risk, all stocks |
|---|---|---|---|---|---|---|
| All 10 assessed, effort + Dentex closure | effort 0.6, closure 3 | 1,603 | 2,206 | 1,022 | 0.044 | 0.044 (SBR) |
| 4 main commercial (HKE, SBR, DOL, VMA), effort only | 0.8 | 1,592 | 2,362 | 891 | 0.047 | **0.81 (DEC)** |
| 4 including the choke stock (DEC, SBR, HKE, VMA), effort + Dentex closure | effort 0.7, closure 3 | 1,615 | 2,195 | 1,016 | 0.028 | 0.028 (SBR) |
| 4 main + permanent shallow-reef closure (no Dentex data) | 0.6 | 1,609 | 2,117 | 1,086 | 0.047 | 0.072 (PAC) |

* **Assessing the main commercial stocks hides the choke.** The rule looks compliant but leaves Dentex below Bpa 81% of the time.
* **Four stocks are enough if the choke stock is among them.** With the closure, every stock is protected and yield is as high as with all ten assessed (1,615 against 1,603).
* **A permanent reef closure protects Dentex without data.** Common pandora (PAC), unassessed, then sits just above the limit (7.2%).

## Tuning the fleet to an open-access equilibrium (current defaults)

**Vessel heterogeneity.** Vessel efficiency within a segment is lognormal, with mean 1 and CV `vessel_cv` = 0.25, as in the random vessel effect of the logbook.

* Each year a segment's days move towards capacity × the share of its vessels that break even (`phi` is the adjustment rate).
* The active vessels' mean efficiency is the mean of that upper tail.

This replaces identical vessels, under which the composition of an open-access fleet is a knife-edge. A tiny change in cost hands the whole fishery to one segment, because every vessel in a segment breaks even at the same stock level. With heterogeneity, effort falls off smoothly as profit falls, and profit at equilibrium is positive: the more efficient vessels earn more than the marginal one.

**`ms_tune()`** solves the daily cost multipliers by segment (G unknowns) so that the open-access equilibrium hits G targets:

* the shares of active days (G − 1 values);
* one status target: median B/Bmsy, the weakest stock, or total days.

It uses Newton steps with a forward-difference Jacobian. The base point and the G perturbed points run in one vectorised call, warm-started and averaged over 30 years, and the best iterate is kept. Each tuning takes about 45 s. Results with capacity 1,000 days per segment and γ = 1.5:

| Status target | Daily costs (dem/pel/mix) | Active days (shares) | Open-access stocks |
|---|---|---|---|
| Median B/Bmsy = 0.9 | 1.00/1.69/1.31 | 370/364/311 (0.354/0.348/0.298) | median 0.90, Dentex lowest at 0.60 |
| Weakest stock = 0.4 (default in `05_plots.R`) | 0.63/1.08/0.84 | 549/528/459 (0.358/0.344/0.299) | yield 1,745, profit 1,489; Dentex 0.40, Seriola 0.47, blackspot seabream 0.50 |

**From the tuned status quo (weakest stock 0.4):**

| Option | Days (dem/pel/mix) | Yield | Profit |
|---|---|---|---|
| Safeguarded MMSY, one cap (same point as MMEY) | 125/126/125 | 991 | 2,386 |
| Safeguarded MMSY, segment caps each ≥ 20% | 132/380/136 | 1,521 | 2,735 |

* **Dentex is again the choke stock, and effort caps barely reach it.** When days are cut, fishers spend the remaining days on the most valuable métier, the shallow reef. Dentex keeps an F/Fmsy of about 1.65 at Bpa even with three-quarters of the days removed. Segment caps recover much of the yield, and a métier measure (reef closure, partial or seasonal) is the obvious next test.
* **The segment HCR must scale status-quo effort, not registered capacity.** Scaled to capacity, the rule never binds. Scaled to status-quo days, it cuts effort by about 40% with no loss of yield (about 1,720), and most stocks recover. Dentex stays at about 0.40.

### What makes Dentex a choke stock (`scripts/10_reef_scenarios.R`)

Each variant is retuned to the same targets (35/35/30% of days, weakest stock 0.4 at open access).

| Scenario | Safeguarded MMSY, one cap: days, yield, profit | Segment caps ≥ 20%: yield | Binding stock |
|---|---|---|---|
| Baseline | 125/126/125, 991, 2,386 | 1,521 | DEC |
| Dentex reef selectivity 1.0 → 0.6 | 364/424/397, 1,729, 2,089 | 1,793 | SBR / AMB |
| Reef fish 20% cheaper | 138/138/137, 1,056, 2,522 | 1,593 | DEC |
| Both | 370/454/418, 1,743, 1,804 | 1,797 | SBR / AMB |

* **Catchability decides the choke; value hardly matters.** With a lower Dentex catchability, one global cap keeps all stocks at or above Bpa with about 25% fewer days than the status quo and no loss of yield. Lower reef prices change little.
* **Dentex chokes with very little effort.** The shallow reef takes only 3–10% of all days. Its catchability relative to its growth rate is what makes it the choke stock.
* **Catchability (habitat breadth) and price belong on the uncertainty axes**, not among the settings to tune away.

## Specialist efficiency and the capacity-limited status quo

**Type efficiency.** Specialists fish more efficiently than generalists.

* `simulate_logbook()` has a new setting, `type_efficiency`, which multiplies the random vessel effect by vessel type. The default is 1 for all types, so existing results are unchanged.
* `ms_settings()` sets mixed vessels to 0.75.
* `ms_system()` uses the same values as the segment efficiency (`eff_home`, `eff_away`). Efficiency multiplies catchability, so F comes from effective effort, and each segment chooses métiers by its own profit, eff × revenue − cost.

**Daily cost follows efficiency** (`cost_seg`, default = efficiency): smaller, cheaper vessels. With equal daily costs, a mixed segment 25% less efficient is excluded:

* under open access it disappears;
* with capacity limits it is squeezed out once the fishery is overexploited, because it can't break even where the specialists do.

At a daily cost of 0.85 it is excluded from 1,000–1,500 days per segment upwards.

**Status quo** (`ms_status_quo()`): capacity-limited, with 1,000 days per segment (a third of the fleet each). Profit decides how many days are used.

* All segments use their full capacity.
* Yield 1,815, profit 535.
* Seriola 0.42 (below Bpa), albacore 0.67, median B/Bmsy about 0.9.

**Safeguarded optima** (γ = 1.5, caps limited to capacity):

| Option | Days (dem/pel/mix) | Yield | Profit |
|---|---|---|---|
| Safeguarded MMSY, one cap | 672 per segment | 1,846 | 920 |
| Safeguarded MMSY, caps each ≥ 20% | 558/469/857 | 1,847 | 1,084 |
| MMEY, one cap | 196 per segment | 1,352 | 2,322 |

**Segment HCR** (responsibility by catchability): all three segments answer for Seriola. Each cuts back to about 700–800 days, Seriola moves up to just below Bpa, and yield is unchanged.

## Balanced reference system

The first version was unbalanced in two ways:

* **Cheap fishing.** With a common cost level, the reef paid long after Dentex was gone.
* **One choke stock.** Dentex choked every option. Low-value métiers stopped paying long before high-value ones, so stocks were either collapsed or unfished.

The defaults now give a plausible unmanaged state: the **open-access equilibrium** (`ms_open_access()`), in which caps don't bind and effort enters and leaves with each segment's profit. Three settings achieve this:

* **Relative costs.** `cost_mode = "relative"`, `cost_share = 0.4`: each métier's cost per unit of effort is 40% of its own unfished revenue rate (with availability). Every métier then breaks even at the same relative revenue loss.
* **Seriola troll catchability 0.4 instead of 0.8** (`ms_settings()`). With 0.8, Seriola is a large by-catch in both reef and troll métiers but the main revenue in neither, and it collapses under open access.
* **Economic exit on by default** (`phi = 0.5`).

**Status quo (open access).**

* Effort 848/1,425/1,161 (demersal/pelagic/mixed), yield 1,803, profit about 0.
* Median B/Bmsy 0.90.
* Several stocks below Bmsy: Dentex 0.80, albacore 0.60.
* One choke stock: Seriola 0.35, between Blim and Bpa.

**Safeguarded optima (every stock ≥ Bpa, γ = 1.5).**

| Option | Effort | Yield | Profit |
|---|---|---|---|
| Status quo, open access | 848/1,425/1,161 | 1,803 | about 0 |
| Safeguarded MMSY, one cap | 565 per segment | 1,846 | 1,054 |
| Safeguarded MMSY, segment caps each ≥ 20% | 767/429/680 | 1,855 | 1,036 |
| MMEY, one cap | 183 per segment | 1,369 | 2,358 |
| Demersal + pelagic caps, mixed free at 30% of open-access effort | — | 1,874 | 897 |

* The one-cap safeguarded MMSY halves effort, restores Seriola to Bpa and turns profit positive.
* The optimum is flat, so a viability floor costs almost nothing.
* With caps unrestricted, the optimiser still drifts towards closing a segment, but for a gain of only about 3%.
* Behaviour now helps. Along an equal-split scan, the highest effort that keeps every stock ≥ Bpa is:
  * about 1,740 for a profit-seeking fleet;
  * about 190 when allocation is profit-blind;
  * about 80 when allocation is frozen.

  Fishers leave depleted métiers.

**A segment HCR depends on who answers for a stock.** If a segment answers only for stocks that are ≥ 5% of its own catch, nobody answers for Seriola (a small share of every segment's catch). The demersal segment is closed, the mixed segment moves in, and Seriola stays at 0.33. With responsibility defined by catchability (any usable métier has ≥ 20% of the stock's highest q), Seriola recovers to about 0.42. A ramp between Blim and Bpa still settles inside the ramp, so the trigger should sit above Bpa.

Results under "Safeguarded MMSY: optimisation results" below come from the first, unbalanced defaults (`cost_mode = "common"`, `cost_share = 0.3`, Seriola troll 0.8, `phi = 0` unless stated).

## Safeguarded MMSY: optimisation results (initial, unbalanced system)

`ms_optimise()` maximises equilibrium yield, revenue or profit over effort caps, subject to every stock staying at or above Bpa. Two more constraints are optional:

* a viability floor (`min_share`);
* métier closures (`open`).

It uses a vectorised zooming grid in log(1 + effort). Gradient methods fail on the kinks and small cycles of a behavioural fleet.

All runs below use a profit-seeking fleet (γ = 1.5) with entry and exit (`phi = 0.5`).

* **Economic exit is not enough.** With entry and exit, effort leaves unprofitable segments and profit settles at about zero. Dentex still collapses, because the other reef species keep the reef métier profitable.
* **Effort caps alone are choke-limited.** Safeguarded yield:
  * one global cap: 288;
  * segment caps with each segment ≥ 20%: 404;
  * otherwise the optimiser closes the demersal fleet (yield 841–863).

  Yield and revenue give the same optimum.
* **The mixed segment as a free agent.** Even uncontrolled, about 87% of its effort goes to demersal métiers. At 20 effort units it uses up Dentex's whole safe budget, so the demersal cap must be near zero. No caps on the specialists meet Bpa once mixed effort is somewhere between 40 and 60. Tightening the demersal cap moves only a little more mixed effort onto demersal métiers (about 15.5 → 17.6).
* **A métier measure decouples the choke species.** Closing the shallow reef métier raises safeguarded yield to about 1,270, 4.4 times the global-cap value. Greater amberjack (AMB) becomes the binding stock, and the split between segments then hardly matters. Dentex is unexploited (B/Bmsy 1.98), so a partial or seasonal closure should do better still.

## Design decisions still open

* **Biomass pools.** Use one pool per species, or sub-pools by ground with exchange? Sub-pools would add local depletion and hyperstability.
* **Within-ground effort concentration.** It is uniform by default. Freezing the unfished value^1.5 distribution makes the catchabilities span about 12 orders of magnitude.
* **Costs.** The level and structure of the cost per unit of effort, and whether costs differ by segment.
* **Expectations and inertia.** Current biomass or lagged CPUE, and whether to add an inertia term on last year's shares.
* **Constraints.** Every species at or above Bpa, or a rule such as 2/3 of species; and segment viability (minimum share or minimum profit).
* **Optimiser.** Replace the coarse log grid with `optim` or a finer grid near the optimum.
* **Fleet dynamics.** Effort creep and entry and exit below the cap (`phi > 0`).

## Towards a package and Shiny app

* **Engine (pure R, no Shiny):**
  * `ms_system()`, built from defaults (the logbook landscape) or from user tables: species (r, K, price, group), métier × species catchability, segments (share, home métiers, loyalty), costs;
  * `ms_simulate()`, `ms_optimise()`, HCR constructors, `ms_summary()`, plots.
* **Shiny:**
  1. **Set-up:** defaults, or uploaded or edited tables. The catchability matrix could come from catch composition by cluster (cpueTargeting step 6).
  2. **Equilibrium explorer:** effort sliders by segment, closures, γ and loyalty, with live status, F, yield and the métier allocation.
  3. **Optimiser:** objective, Bpa, viability floor, levers, and which segments are free agents.
  4. **Projections:** HCR choice, uncertainty, trade-off plots.
* **Speed:** an equilibrium takes milliseconds; an optimisation takes 7–50 s, depending on the number of levers.
