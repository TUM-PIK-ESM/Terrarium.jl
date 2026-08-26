# Investigation: insufficient high-latitude winter snow accumulation (SEB/snow-fraction coupling)

> Status: **planned**. Investigation + candidate fixes for the reported bug (SWE < 0.001 m at northern
> latitudes in winter, likely from excessive melt).

Date of initial draft: 2026-08-27

Base revision: 836bbb44b2e04c8e2ff0882be7d61ae1d2e6094f

## Originating prompt

> There currently appears to be a bug where not enough snow is accumulating in the northern latitudes
> during winter (<0.001 m SWE), likely due to it being melted too fast. Please review the current
> implementation as well as the text from the GitHub issue [pasted below] and prepare a plan for an
> investigation + potential fixes.

GitHub issue (verbatim, condensed): the surface energy balance (SEB) uses `ground_heat_flux` to force
the top of the snowpack. With fractional snow cover and no tiling, the SEB blends ground/snow inputs
(temperature, conductivity) into a single bulk `ground_heat_flux`, which is then applied unmediated as
the snowpack's top boundary flux. This likely produces unrealistic snow-top fluxes. The ideal fix is
tiling (solve the SEB separately for snow-covered and snow-free fractions); an interim alternative is to
force the *actual ground* with `ground_heat_flux` and partition it in the SEB into `G_g` (snow-free
conductive flux) and `G_s` (snow-covered conductive flux):
`R_net(Ts) + H(Ts) + LE(Ts) - G_g(Ts) - G_s(Ts) = 0`. A remaining concern even under partitioning: the
skin temperature is still solved with a single aggregated albedo in `R_net`, so the snow-covered fraction
still sees an unrealistic radiation budget.

## Revision log

- Rev 1 (2026-08-27): initial draft from code review of `skin_temperature.jl`, `surface_energy_balance.jl`,
  `snow_energy.jl`, `snow_interfaces.jl`, `albedo.jl`, and `land_model.jl`.

## Problem description

Northern-latitude winter SWE stays below 0.001 m — snow that does accumulate appears to melt (or
sublimate/never accumulate) almost immediately, rather than persisting through the cold season.

## Background: how the bulk-blend forcing actually reaches the snowpack today

Confirmed by reading the current implementation (not yet run/instrumented):

1. **Skin temperature solve blends the conduction target, not just diagnostics.**
   [`ground_thermal_interface`](../../src/processes/surface/skin_temperature.jl#L121-L132) (snow-aware
   overload) computes `Tg_eff = (1-f)·T_ground + f·T_snow`, `κ_eff = (1-f)·κₛ + f·κ_snow`,
   `Δz_eff = (1-f)·Δz_ground + f·d_snow`, and the implicit solve uses these to get
   `Ts = Tg_eff - G·Δz_eff/(2·κ_eff)` ([`skin_temperature.jl:94-97`](../../src/processes/surface/skin_temperature.jl#L94-L97)).
   `G` itself closes the SEB as `G = R_net(Ts) + H_s(Ts) + H_l(Ts)`
   ([`skin_temperature.jl:141-144`](../../src/processes/surface/skin_temperature.jl#L141-L144)). So there is
   a *single* skin temperature and a *single* bulk `G` for the whole grid cell, exactly as the issue
   describes — this is not tiled.
2. **That same bulk `G` is fed directly into the snowpack's top-of-pack energy tendency, unmediated by
   `f_snow`.** In `LandModel`, `snow_flux_aliases = (; surface_heat_flux = ground_heat_flux, ...)`
   ([`land_model.jl:81`](../../src/models/coupled/land_model.jl#L81)) — no `f_snow` weighting at all — and
   [`compute_snow_energy_tendency`](../../src/processes/snow/energy/snow_energy.jl#L79-L100) reads
   `Q_top = fields.surface_heat_flux[i,j]` directly as `dŪ_snow/dt = Q_base - Q_top + ...`. So a grid cell
   with, say, `f_snow = 0.1` still has the *entire* bulk `G` (computed from a 90%-bare-ground-weighted
   conduction target and a blended albedo) driving 100% of the snowpack's top-of-pack tendency. This is
   the mechanism most directly implicated by the issue and the most likely single cause of excess melt:
   thin/patchy snow gets hit with a flux sized for a much larger effective bare-ground fraction.
3. **The soil-top flux is already properly `f_snow`-blended** —
   [`compute_snow_soil_heat_flux`](../../src/processes/snow/snow_interfaces.jl#L37-L52) computes
   `f·Q_base + (1-f)·G` for `soil_heat_flux` — but this only fixes the *soil's* boundary condition, not the
   snowpack's own top-of-pack forcing from step 2.
4. **Albedo is also bulk-blended before the skin-temperature solve.**
   [`compute_albedo` (`DiagnosticAlbedo`)](../../src/processes/surface/albedo.jl#L94-L118) computes a single
   `α_eff = (1-f_snow)(1-f_veg)·α₀ + (1-f_snow)·f_veg·α_veg + f_snow·α_snow` used for the one skin-temperature
   solve. For small `f_snow` this washes out snow's high albedo, inflating absorbed shortwave (and hence
   `G`, `Ts`) relative to what the snow-covered patch would actually see — compounding point 2.
5. A prior design doc ([`2026-07-10_PLAN_01B_surface_ground_heat_flux_split.md`](../2026-07/2026-07-10_PLAN_01B_surface_ground_heat_flux_split.md))
   already flagged the semantic ambiguity between the SEB closure flux and the soil BC and split them by
   name (`ground_heat_flux` vs `soil_heat_flux`), fixing the soil side. It explicitly deferred the
   *snow-top* forcing question — this is that deferred problem.

This is consistent with, and sharpens, the issue text: the reporter suspected the bulk blend was the
culprit but had not yet traced it all the way to `surface_heat_flux = ground_heat_flux` with **zero**
`f_snow` weighting at the snow-tendency boundary. That zero-weighting (not just "some blending exists")
is likely the dominant term, especially for the thin/patchy snow (`f_snow` small) typical of snow onset
and ablation in shoulder seasons — exactly when spurious melt would show up as chronically-suppressed
peak SWE.

## Investigation plan (before committing to a fix)

1. **Instrument and confirm the mechanism.** Add a diagnostic run (single grid cell or small regional
   subset at a northern-latitude site with known snow cover, e.g. from ERA5/CRU forcing already used in
   existing examples/tests) that logs `f_snow`, `ground_heat_flux` (bulk `G`), `compute_snow_soil_heat_flux`
   (`Q_base` blend), `snow_temperature`, `skin_temperature`, `albedo`, and the snow energy tendency terms
   (`Q_top`, `Q_base`, `Q_prcp`, `Q_subl`) through a full autumn-onset → winter → melt cycle. Confirm that
   (a) `Q_top` tracks the *bulk* `G` rather than a snow-covered-fraction-appropriate flux, and (b) melt
   events correlate with periods of small `f_snow` (patchy cover) rather than genuine warm spells.
2. **Quantify the albedo effect in isolation.** With `f_snow` diagnostics from (1), estimate how much
   absorbed shortwave increases from blended vs. snow-only albedo at representative small-`f_snow` values,
   to judge whether fixing point 2 alone (without touching albedo) would be enough, or whether point 4
   also needs addressing per the issue's closing concern.
3. **Check `snow_cover_fraction`'s functional form** ([`snow_cover.jl`](../../src/processes/snow/mass/snow_cover.jl))
   — if it saturates slowly (small `f_snow` for a wide range of low SWE), the zero-weighting bug in point 2
   would bite hardest exactly during accumulation onset, which matches "not enough snow accumulating."
4. **Rule out independent causes** before attributing the whole deficit to SEB coupling: verify
   snowfall-phase partitioning (rain/snow split by `T_air`), sublimation magnitude
   ([`compute_snow_sublimation_flux`](../../src/processes/snow/mass/snow_mass.jl#L13-L29)), and the
   `min_conduction_thickness` floor behavior for very thin packs (already flagged as thermally transient by
   design — worth confirming it isn't itself injecting spurious heat during pack initiation).

## Candidate fixes (in increasing order of scope)

**Fix A — weight the snowpack's top-of-pack flux by `f_snow` (minimal, addresses point 2 directly).**
Instead of aliasing `surface_heat_flux = ground_heat_flux` unweighted, apply the same `f_snow` blend
already used for the soil side: derive a snow-only top flux and multiply, or (more physically) recognize
that with `f_snow` already baked into the *conduction target* (`Tg_eff`, `κ_eff`, `Δz_eff` blend) `G` is
not simply separable into "the snow's share" — this needs care, not just multiplying `G` by `f_snow`
(that alone would not correct the albedo/conduction-target distortion in point 4, only rescale the
already-distorted flux). This is the smallest change but only partially matches the issue's own proposed
interim fix.

**Fix B — issue's proposed interim fix: fix `ground_heat_flux` to the true ground interface and
explicitly partition `G_g`/`G_s` in the SEB residual.** Redefine the skin-temperature residual as
`R_net(Ts) + H(Ts) + LE(Ts) - G_g(Ts) - G_s(Ts) = 0`, where `G_g` uses the snow-free conduction target
(pure ground) and `G_s` uses the pure-snow conduction target, each weighted by `(1-f_snow)`/`f_snow`
respectively (no blending of `Tg`/`κ`/`Δz` themselves — only the two resulting fluxes are combined).
`G_s` (or `G_s/f_snow`, the per-area snow flux) then drives the snowpack tendency directly, replacing the
current unweighted alias. This directly fixes point 2 and is a bounded, well-specified change to
`skin_temperature.jl` (new `compute_ground_heat_flux` partition), `land_model.jl` (wire `surface_heat_flux`
to the new `G_s`-derived field instead of raw `ground_heat_flux`), and `snow_interfaces.jl` (`Q_base` blend
becomes consistent with `G_s`'s ground-facing counterpart, avoiding double-counting). Leaves the
single-albedo/single-`Ts` limitation (point 4) unresolved, as the issue itself notes.

**Fix C — full tiling.** Solve the SEB independently for the snow-covered and snow-free sub-fractions
(two skin temperatures, two albedos, two turbulent-flux solves), then area-weight only the *aggregate*
diagnostics (mean skin temperature for output, total sensible/latent flux for atmosphere coupling). This
is the "ideal" fix per the issue and resolves both point 2 and point 4, but is a substantially larger
change: doubles the per-cell nonlinear solve cost, touches `SurfaceEnergyBalance`, `skin_temperature.jl`,
`turbulent_fluxes.jl`, `radiative_fluxes.jl`, `albedo.jl`, and the snow/soil coupling in `land_model.jl`
and `snow_interfaces.jl`. Also raises open questions about how vegetation fraction composes with snow
tiling (currently vegetation and snow fractions are both blended into the same single albedo formula) and
about GPU/Reactant kernel cost of two solves per cell.

## Recommendation

Run the investigation steps first (mechanism confirmation is cheap and will validate whether Fix B alone
is sufficient before committing to the much larger Fix C). Given the issue author's own framing ("another
interim solution" for B, deferred judgement on tiling pending their own investigation), propose starting
with **Fix B** as a scoped, physically-motivated interim PR, explicitly flagging the residual
single-albedo limitation as known-limitation / future-work (candidate follow-up: apply the same `G_g`/`G_s`
partition idea to albedo, i.e. compute `R_net` twice with `α_ground`/`α_snow` and combine post-hoc, without
going all the way to full turbulent-flux tiling — worth scoping as **Fix B+** if the investigation shows
albedo blending is still a significant residual error after B).

## Summary of changes

(To be filled in once a fix is selected and approved — this document currently proposes B as the
scoped candidate. No code changes have been made yet.)

## Testing and verification

- Regression test: no-snow behavior must be bit-for-bit unchanged (same pattern as PR-B in
  `2026-07-10_PLAN_01B_surface_ground_heat_flux_split.md`).
- New unit test(s) in `test/snow/` (or `test/surface/`) exercising the SEB residual with partitioned
  `G_g`/`G_s` at a few `f_snow` values, checking energy conservation: `G_g + G_s` (area-weighted) should
  recover the same total ground-column energy input as today's bulk `G` in the `f_snow → 0` and
  `f_snow → 1` limits.
- End-to-end site-level regression against the northern-latitude case used in the investigation step,
  checking peak SWE and snow-season length against a reasonable reference (reanalysis snow cover or a
  previous known-good model version, if available).
- Differentiability: re-run `test/differentiability` and `test/reactant/autodiff.jl` relevant snow/SEB
  cases after the change, since `compute_ground_heat_flux` and `ground_thermal_interface` sit inside the
  implicit skin-temperature solve's differentiated call graph.

## Documentation changes

- Update the SEB residual equation and prose in the `ImplicitSkinTemperature` docstring
  ([`skin_temperature.jl:42-50`](../../src/processes/surface/skin_temperature.jl#L42-L50)) to reflect the
  `G_g`/`G_s` partition.
- Update `docs/src` process pages for surface energy balance / snow coupling (per `AGENTS.md` doc-page
  conventions: Overview/Implementations/Methods/Kernel-functions sections) once the concrete
  implementation is decided.

## Known limitations

- Fix B does not address the aggregated-albedo distortion of `R_net` (issue's closing concern); this
  remains until Fix B+ or full tiling (Fix C) is implemented.
- Fix B still solves a single skin temperature, so turbulent fluxes (`H_s`, `H_l`) are not tiled either —
  only the conductive term is partitioned.

## Future work

- Full tiling (Fix C), if the investigation or Fix B's post-hoc evaluation shows residual error from the
  shared skin temperature/albedo is still significant.
- Revisit `snow_cover_fraction`'s functional form if step 3 of the investigation shows pathological
  behavior at low SWE.
