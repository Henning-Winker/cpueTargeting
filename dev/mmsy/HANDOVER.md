# From logbook simulator to multispecies harvest-strategy simulator: hand-over (21 October 2026)

Prototype: branch `dev-mmsy` of https://github.com/Henning-Winker/cpueTargeting, folder `dev/mmsy/` (excluded from the package build). The results log is `dev/mmsy/README.md` (also in the project as `claude/mmsy_system_design.md`). This note is the hand-over for setting up a new, stand-alone package with a Shiny app.

## 1. Purpose

The tool explores **effort-based harvest strategies in mixed fisheries** where fishers respond to profit. It answers three questions:

* What multispecies MSY (MMSY) can be achieved when every stock is kept at or above a safeguard (Bpa = 0.5 Bmsy, Blim = 0.25 Bmsy)?
* How much does it matter how effort is split across fleet segments (pelagic, mixed, demersal) and whether there are métier closures?
* Which feedback rules (effort, closures) meet a risk limit under realistic uncertainty, and at what cost?

It replaces the shared-q model of the FAO RAP app (`mse.SEA`), in which every effort rule only moves along one line (F_i / F_j = q_i / q_j is fixed), with a fleet that chooses where to fish.

## 2. How it evolved from the logbook simulator

### Kept from `simulate_logbook()`

* **Species.** The ten species, with group, r, B0 = K, price and seasonal availability (pelagic migrants).
* **Habitat landscape.** Depth, reef patches and SST, giving the species' spatial distributions. It is used only once, to derive the catchability matrix.
* **Tactics and selectivity.** Reef (handline), Sediment (longline) and Pelagic (troll), with the gear selectivity matrix.
* **Fleet types.** Demersal, pelagic and mixed vessels, with home tactics and loyalty (0.8 / 0.8 / 0).
* **Choice by expected value.** The idea that trips choose tactics by expected value (price × catchability × biomass).
* **Vessel heterogeneity.** The random vessel effect (lognormal, CV 0.25), plus a new type multiplier (`type_efficiency` in `logbook_settings()`, default 1). With mixed = 0.75, generalists are less efficient.

### Simplified

| Logbook simulator | Fleet system |
|---|---|
| 40 × 40 grid, cell-level Schaefer and movement | One Schaefer pool per stock (fast mixing) |
| Cells and fishing positions | Grounds: shallow/deep (split at the 300 m closure depth) × reef/sediment, giving six métiers (tactic × ground) |
| Catchability by cell | q[m, i] = selectivity × mean share of K per cell of the ground. The within-ground set distribution is uniform: freezing the unfished value^1.5 made the q's span 12 orders of magnitude |
| Trips, sets, Tweedie catch noise | Effort in days by segment × métier; Baranov-type monthly catch |
| Logbook reporting layer (positions, rounding, under-reporting, synonyms) | Dropped (it remains in the logbook simulator as an optional observation model) |
| Individual vessels, entry/exit years, creep | A lognormal efficiency distribution per segment; days follow the share of vessels that break even |
| Price periods | Fixed price ratios |

### Added

* **Economics.** A daily cost per métier, which by default is 40% of the métier's own unfished revenue rate, so every métier breaks even at the same relative revenue loss. On top of that come a daily cost multiplier and an efficiency per segment.
* **Preference function.** Each segment allocates its days over métiers with P(m) ∝ max(eff_gm · value_m − cost_gm, 0)^γ, split into a loyal part (home métiers) and a free part (all métiers). γ = 0 is profit-blind; γ = 1.5 is the default.
* **Effort supply with vessel heterogeneity.** Each segment's days move towards capacity × the share of vessels that break even. The active vessels' mean efficiency is the mean of that upper tail. Without heterogeneity the open-access composition is a knife-edge: one segment takes everything.
* **Tuning the fleet (`ms_tune()`).** Solves the daily costs by segment so that a target fleet composition (35/35/30% of days) and a stock-status target (e.g. weakest stock 0.4 Bmsy) form the open-access equilibrium, the status quo. It uses Newton steps with a vectorised Jacobian.
* **Constrained optimisation (`ms_optimise()`).** Maximises equilibrium yield, revenue or profit over effort caps subject to every stock ≥ Bpa. Options: a viability floor per segment, métier closures, and free-agent segments.
* **Feedback control (`hcr_effort()`).** Segment days ramp from 0 at Blim to full at a trigger, applied to the most depleted stock each segment can catch substantially. Options:
  * a choke-stock closure (the open share of the métiers that drive it), with its own trigger;
  * change limits;
  * a subset of assessed stocks.
* **Uncertainty.**
  * Process error in JABBA form on the whole process equation, with SD 0.1 × r_i / median r and AR(1) 0.3 (an optional correlation across species).
  * Assessment error, CV 0.2.
  * Implementation error on the caps, CV 0.1.
  * A data lag of 2 years.
* **Tuning the trigger (`ms_tune_hcr()`).** A grid of triggers runs in one vectorised call with common random numbers. The tuned rule maximises yield subject to a risk limit (every stock below Bpa at most 5% of the time).
* **Diagnostics.** `ms_performance()`, equilibrium yield curves for different conceptual treatments, Kobe plots and projection plots.

## 3. Structure (equations)

```
Stocks (S):   B_i,t+1 = (B_i,t + r_i B_i,t (1 − B_i,t/K_i) − C_i,t) · exp(η_i,t)
              η_i,t ~ AR(1), SD_i = σ · r_i / median(r)        (σ = 0.1)
Catch (monthly):  C_i = B_i (1 − exp(−F_i)),  F_i = Σ_m E^eff_m · open_m · q_m,i · avail_i,month
Effective effort: E^eff_m = Σ_g E_g,m · eff_g,m · ē_g    (ē_g = mean efficiency of active vessels)
Allocation:   E_g,m = E_g [λ_g P(m | home_g) + (1 − λ_g) P(m | all)]
              P(m | set) ∝ open_m · max(eff_g,m ē_g · v_m − c_m · c_g, 0)^γ,  v_m = Σ_i p_i q_m,i B_i avail_i
Effort supply (yearly):  E_g ← E_g + φ (min(cap_g, capacity_g · P(e > e*_g)) − E_g)
              e*_g = cost per day / revenue per day of a mean vessel; e ~ lognormal(mean 1, CV 0.25)
Rule (year y, data year y − 2):  cap_g = E_ref,g · ramp(min_{i ∈ own_g} B̂_i/Bmsy_i; Blim, Btrig) · exp(ε_impl)
              open_m = ramp(B̂_choke/Bmsy; Blim, Btrig_close) for the métiers that drive the choke stock
```

**Defaults (balanced reference system).**

* Relative costs, `cost_share` 0.4.
* Seriola troll selectivity 0.4 instead of 0.8. With 0.8, the reef + troll by-catch collapses it under open access.
* Mixed segment efficiency 0.75; daily cost multiplier tuned.
* Capacity 1,000 days per segment.
* Status quo tuned to 35/35/30% of days, with the weakest stock at 0.4 Bmsy:
  * yield 1,745, profit 1,489;
  * Dentex 0.40, Seriola 0.47, blackspot seabream 0.50.

## 4. Key findings (Atlantic example, tuned system)

1. **Effort caps alone are choke-limited.**
   * Dentex (narrow reef habitat, the highest catchability relative to its growth rate) chokes the fishery even though the shallow reef takes only 3–10% of the days.
   * When days are cut, the remaining days go to the most valuable métier, so Dentex keeps a high F.
   * Lower Dentex catchability removes the choke; cheaper reef fish barely help.
2. **Economic exit protects no stock that is a by-catch of a profitable métier.** Profit settles near zero while Dentex collapses.
3. **Equilibrium yield curves** (`figs/yield_curves.png`):
   * Shared q gives a dome with a false sense of safety: all stocks are fine at MMSY and collapse together beyond it.
   * Behaviour turns the dome into a plateau: yield is insensitive to effort while weak stocks erode.
   * The gap between MMSY and the safeguarded yield measures the technical interaction. A métier measure closes it.
4. **Feedback control** (2-year lag; risk limit P(B < Bpa) ≤ 5% per stock):
   * The effort rule alone never meets the limit.
   * Effort rule + reef closure meets it. With separate triggers (effort 0.6, closure 3.0 Bmsy), yield is 1,603 against 1,667 without control, and profit is +47%.
   * The lag forces a shared trigger up from 0.85 to 1.0 Bmsy (−5% yield) but barely affects separate triggers.
5. **Assessing four stocks is enough if the choke stock is one of them.**
   * Assessing the four main commercial stocks hides Dentex (below Bpa 81% of the time).
   * A permanent reef closure protects Dentex without data.

## 5. Proposed new package

Working name: `mmsyfleet` (alternatives: `mixfishHS`, `fleetMMSY`). It should be stand-alone and **not depend on cpueTargeting**. The Atlantic example ships as tables (`dev/mmsy/extdata/atlantic_*.csv`, exported by `scripts/15_export_tables.R`). `derive_landscape()` stays in cpueTargeting as the way the tables were made.

### Layout

```
mmsyfleet/
├── DESCRIPTION            Imports: stats, ggplot2; Suggests: shiny, bslib, DT, patchwork, ggrepel, testthat, knitr, rmarkdown
├── NAMESPACE              (roxygen2)
├── R/
│   ├── system.R           ms_system() from tables: species, métiers, catchability, segments; validation; print method
│   ├── presets.R          atlantic_system(), gulf_of_thailand_system() (read inst/extdata)
│   ├── allocation.R       allocate_effort()
│   ├── simulate.R         ms_simulate() (dynamics, effort supply, rules, errors, lag, crn)
│   ├── equilibrium.R      ms_open_access(), ms_status_quo(), ms_tune(), ms_yield_curve()
│   ├── optimise.R         ms_optimise()
│   ├── rules.R            hcr_effort() (+ future: CPUE-based rules, group-aggregated assessments)
│   ├── tune_rules.R       ms_tune_hcr()
│   ├── performance.R      ms_summary(), perf_rows(), ms_performance()
│   ├── plots.R            theme, palette, equilibrium, Kobe, projection, yield-curve and tuning plots
│   └── utils.R
├── inst/
│   ├── extdata/           atlantic_*.csv, gulf_of_thailand_*.csv
│   └── shiny/mmsyfleet/   app.R + modules (set-up, equilibrium explorer, optimiser, feedback/MSE, tuning)
├── tests/testthat/        see the checks below
├── vignettes/             intro, conceptual yield curves, harvest strategies, Gulf of Thailand case study
├── _pkgdown.yml           site in docs/
└── .github/workflows/     R-CMD-check.yaml, pkgdown.yaml
```

`run_app()` exported: `shiny::runApp(system.file("shiny/mmsyfleet", package = "mmsyfleet"))`.

### Shiny tabs

1. **Set-up.** Choose a preset or upload or edit the tables. Tune the status quo (target shares and status). Show the catchability heatmap.
2. **Equilibrium explorer.** Sliders for days by segment, closures, γ and loyalty. Live status, F, yield, profit and the métier allocation.
3. **Yield curves.** The conceptual treatments compared (shared q against behaviour against closures).
4. **Optimiser.** Objective, Bpa, viability floor, free-agent segments, closures.
5. **Harvest strategies.** Rule type, triggers, assessed stocks, lag, errors. Projections and the performance table. Trade-off plot across runs.
6. **Tuning.** A risk limit gives the tuned triggers (1D curve or 2D grid).

The grid computations run in about 10–60 s, so use `shiny::ExtendedTask` or a progress bar.

### Tests (port from the prototype checks)

* **Nesting:** at γ = 0, fixed and dynamic allocation give identical results.
* **Analytic equilibrium:** with fixed allocation, B/K matches 1 − F/r to within 0.01.
* **Efficiency:** with all efficiencies 1, the earlier status quo is reproduced.
* **Common random numbers:** identical rows give identical trajectories, and a grid with one trigger equals a single run.
* **Lag:** `hcr_lag = 0` reproduces the no-lag tuning.
* **`ms_tune()`:** reaches the target shares to within 0.01 on the preset.
* **Input validation:** `ms_system()` rejects tables that don't match.

### Steps for tomorrow

1. Create the repository (e.g. `Henning-Winker/mmsyfleet`) and run `usethis::create_package()`, `use_mit_license()`, `use_testthat()`, `use_github_action("check-standard")`, `use_pkgdown_github_pages()`.
2. Copy `dev/mmsy/R/msfleet.R` and `plots.R`, split into the files above, and add roxygen headers. The functions already have roxygen-style comments.
3. Write `ms_system()` from tables, the main refactor: it currently builds from the logbook landscape. Add `atlantic_system()` reading `inst/extdata/atlantic_*.csv`.
4. Port the checks to `tests/testthat`.
5. Port scripts 05, 11, 12 and 14 into the first two vignettes.
6. Write the Shiny skeleton (tabs 1, 2, 3 first).
7. Add the Gulf of Thailand preset (section 6).

## 6. Gulf of Thailand example: design

**Why.** This is the FAO RAP context of `mse.SEA` and the GEF/FAO GoTFish project (Cambodia, Malaysia, Thailand, Viet Nam; 2026–2030). GoTFish targets the demersal trawl fishery, the pelagic purse-seine fishery and the small-scale coastal fishery, and states that more than 40% of stocks are overexploited ([SFP](https://sustainablefish.org/strengthening-fisheries-governance-in-the-gulf-of-thailand/)).

### What the literature gives

* **Depletion history.**
  * Trawl survey biomass in 1995 was about 8% of 1961 ([Kongprom et al. 2003](https://digitalarchive.worldfishcenter.org/handle/20.500.12348/2126)).
  * In 1995, demersal effort was about twice the level for MSY (Fox model), and most of 23 species had E ≥ 0.79.
  * Survey CPUE fell from about 300 kg/h (1961) to about 20–24 kg/h ([FAO country profile](https://www.fao.org/fishery/docs/DOCUMENT/fcp/en/FI_CP_TH.pdf); [SEAFDEC refugia](https://seafdec.or.th/fisheries-refugia/home/country-activities/thailand/background/status-and-trends.html)).
* **Catch composition.**
  * Trash fish are about 60% of the trawl catch, and 18–32% of it is juveniles of commercial species.
  * Indo-Pacific mackerel is fully exploited; sardines are overexploited.
  * Squid (*Photololigo duvaucelii*) has been fully exploited since the mid-1980s; penaeid shrimp is overexploited ([FAO](https://www.fao.org/fishery/docs/DOCUMENT/fcp/en/FI_CP_TH.pdf); [SEAFDEC](https://seafdec.or.th/fisheries-refugia/home/country-activities/thailand/background/status-and-trends.html)).
* **Management (Thailand).** Under the Royal Ordinance on Fisheries 2015:
  * three assessed groups (demersal fish, pelagic fish, anchovy);
  * TAC at 95% of MSY;
  * licences limited to a total allowable effort, with fishing-day allocations for high-efficiency gear;
  * seasonal closures for Indo-Pacific mackerel spawning;
  * small-scale demarcation zones ([Thai DoF 2023](https://infofish.org/v4/media/attachments/2023/10/19/tha---fao-workshop-fisheries-support-measures-oct-2023---thai-dof-final--9-oct-66.pdf); [DoF](https://www4.fisheries.go.th/dof_en/view_news/426)).

  2016 limits for the Gulf, which give the K scale:

  | Group | Limit |
  |---|---|
  | Fish | 230,803 t/yr |
  | Marine benthos | 715,294 t |
  | Anchovy | 172,607 t |

  Source: [Bangkok Post](https://www.bangkokpost.com/thailand/general/859528/new-fishing-restrictions-to-take-effect).
* **Fleet.** In 2022, 9,608 commercial vessels (≥ 10 GT) took 21% of the catch and 50,639 small-scale vessels took 79% (national, [DoF 2023](https://infofish.org/v4/media/attachments/2023/10/19/tha---fao-workshop-fisheries-support-measures-oct-2023---thai-dof-final--9-oct-66.pdf)). Since 2015, about 80% of effort has been controlled through vessel numbers and monitoring ([ESL](https://www.esl.co.jp/products/menu/thailand_2ndworkshop/(0)%20Fisheries%20regime%20shift%20(q)%20(GOT).pdf)).
* **Mackerel closure.** In the upper southern Gulf, trawl and purse seine with mesh < 4.7 cm are banned from 1 Feb to 15 May (since 1984; the area was negotiated down to about 50,000 km²) ([FAO](https://www.fao.org/fishery/docs/DOCUMENT/fcp/en/FI_CP_TH.pdf); [Saikliang](https://www.ais.unwater.org/aisresources/projectdoc/rpt/885/IWLEARN/Saikliang-fisheries-refugias.html)).

### Species groups from the mse.SEA input (`extdata/got_species_mseSEA.csv`)

The uploaded `GOT_full.csv` has 18 species groups with K, r, 2024 depletion (B/K), q and a group label. These replace the indicative r and K of the table below as the starting point.

* **Status.** 13 of 18 groups are below Bpa (B/Bmsy < 0.5), and 4 are below Blim: bigeyes 0.06, croakers 0.12, threadfin bream 0.14, Indo-Pacific mackerel 0.16.
* **Slow stocks (r ≈ 0.05–0.06).** Snapper, catfish, croakers, rays and sharks are the choke candidates.
* **Yield.** The sum of MSY is 475 kt; equilibrium yield at the 2024 depletions is 184 kt.
* **Shrimp looks anomalous.** Depletion 0.97 (nearly unfished), r 1.44, and MSY 197 kt, which is 42% of the total. This conflicts with FAO's "overexploited". Check the assessment.
* **The group labels mix habitat and gear.** Rays, sharks, lizardfish, threadfin bream, bigeyes, hairtail and pomfret are labelled "Pelagic". For segments and métiers, regroup them by gear (trawl / purse seine / small-scale).
* **The q values are placeholders,** and the shared-q model can't reproduce the depletions. Fitting one shared effort to them (E ≈ 617) gives 0 for snapper, catfish and croakers, but 0.97–0.99 for bigeyes, Indo-Pacific mackerel and squid. The 2024 pattern needs group-specific fishing mortality, F_i = r_i (1 − d_i) at equilibrium.
* **Conditioning in the new system.** Tune catchability scalers by stock (or by métier × group) together with the segment costs, so that the status quo reproduces the 2024 depletions. This extends `ms_tune()` from one status target to a depletion vector. The bigeyes and croakers depletions (0.03–0.06) are probably not equilibrium states, so condition on the trend, or start the projection from the 2024 state.

### Proposed system (r, K split and prices are indicative, to be refined)

**Stocks**

| Stock | Group | Main gears | r (indicative) | Price (relative) | Role |
|---|---|---|---|---|---|
| Threadfin breams (*Nemipterus*) | Demersal | Trawl | 0.9 | 2 | Main demersal food fish |
| Lizardfishes (*Saurida*) | Demersal | Trawl | 0.8 | 1.5 | |
| Bigeyes (*Priacanthus*) | Demersal | Trawl | 0.8 | 1.5 | |
| Squids (*Photololigo*) | Demersal | Trawl, light lure | 1.2 | 6 | High value, fast |
| Penaeid shrimp | Demersal | Push net, trawl, small-scale | 1.5 | 10 | High value, coastal |
| Blue swimming crab | Demersal | Gillnet, trap (small-scale) | 1.0 | 8 | Small-scale target |
| Trash fish (low value, juveniles) | Demersal | Trawl, push net | 1.0 | 0.3 | Large biomass, fishmeal |
| Rays and sharks | Demersal by-catch | Trawl | 0.12 | 1 | **Low-productivity, low-value choke** (ETP angle of GoTFish) |
| Indo-Pacific mackerel | Pelagic | Purse seine, encircling gillnet | 1.0 | 2.5 | Seasonal spawning closure |
| Sardines and round scads | Pelagic | Purse seine | 1.0 | 1 | |
| Anchovy (*Stolephorus*) | Anchovy | Light-luring anchovy purse seine | 1.5 | 1 | Fish sauce chain |

* K by group from the 2016 limits read as MSY, with K = 4 MSY / r, split within groups by share.

**Segments**

| Segment | Fleet | Capacity, efficiency, cost | Home métiers |
|---|---|---|---|
| Trawl | Otter-board and pair trawl, commercial | Capacity from licences | Inshore and offshore trawl |
| Purse seine | Commercial | | Mackerel/scad purse seine (+ anchovy light purse seine, or a separate anchovy segment) |
| Small-scale | Gillnets, traps, push nets, multi-gear | Most vessels, efficiency about 0.3–0.5 per day, low cost | Coastal métiers; the free agent |

**Métiers**

* Trawl inshore and trawl offshore. The 3 km inshore trawl ban acts as a closure.
* Purse seine (mackerel/scads), anchovy light purse seine.
* Coastal gillnet/trap, push net.

**Contrast with the Atlantic example**

* The choke stock is a **low-value, low-productivity by-catch** (rays and sharks). Profit-seeking does not protect it, because its value does not attract effort, but nothing makes the trawl unprofitable either.
* **Assessments are by aggregated group** (demersal, pelagic, anchovy), so depletion within a group can hide.
* **The closure is seasonal** (the mackerel spawning closure), not spatial.

### Engine changes needed for the Gulf of Thailand

1. **`ms_system()` from tables** (also needed for the package).
2. **Seasonal closures:** `open` as a months × métiers matrix, or a function of (year, month). It is yearly only now.
3. **Group-aggregated assessments in `hcr_effort()`:** B̂ of a group = Σ B over its stocks against Σ Bmsy, with an option to tune on the group risk and report stock risk.
4. **Four segments:** the engine is generic in G and M, so this just needs tests.
5. Optional: a **juvenile/trash-fish link** (growth overfishing of commercial demersals through the trawl and push-net take of juveniles) as a simple loss term. This should be a later refinement.

## 7. Open points

* **CPUE indices as the rule's signal:** nominal CPUE carries targeting bias; standardised CPUE from the cpueTargeting toolbox could be used instead. Ideally the logbook simulator would serve as the observation model.
* **Robustness set:** Dentex catchability, γ, `vessel_cv`, effort creep (q rising over time), implementation bias, and price shocks.
* **Partial or seasonal reef closure** instead of a ramp on the open share.
* **Trigger above Bpa:** a ramp between Blim and Bpa settles inside the ramp. The tuned triggers (≥ 0.6–1.0 Bmsy) handle this.
