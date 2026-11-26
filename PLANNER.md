# Spatial sampling calculator and protocol planner

Open `index.html` in a modern browser. Inputs and results stay on the device. Enter a measured single-request probability or multiply the placement probability, full-budget probability given placement, and dual recovery given full budget. The calculator gives the integer number of entire requests needed for 80% or 95% probability of at least one success. A request can comprise multiple sites. Positive intraclass correlation uses a shared beta propensity. Counts above one trillion are reported separately from a true zero-probability event. Set site diameter or square side to calculate the geometric footprint.

With Python and the dependencies of the v0.5.0 core installed, run these commands from the repository root.

```sh
python protocol_plan.py --input protocol_demo.csv --section demo --targets TargetA TargetB --budget 600 --sites 6 12 20 --area-mm2 0.04 --size-um 55 --um-per-unit 1 --draws 2000 --out demo_output
python target_pair_lookup.py --input protocol_demo.csv --targets TargetA TargetB --reference-K 50 --out pair_lookup.csv
python extraction_bias_sim.py --mode generative --out retention_model
```

For the layout planner, each unique site requires its section, coordinates, measured total and target counts, predefined dark/light stratum and a binary `geometry_ok` flag calculated for the requested footprint. The program checks pairwise overlap and area after generating nominal layouts. For another diameter or shape, reassess the clearance flags. Outputs include feasible example coordinates, geometrical placement, budget attainment, the two molecular endpoints and the three-objective Pareto set. RMSE is evaluated for full layouts containing both strata. The uniform target-pair lookup uses all twelve-unit requests as its denominator and tests the eligible integer depths until success or a finite-population bound. Use independent biological units for inference.

The numerical extraction model compares across-region and lateral placements. Empirical mode accepts a separate pre-registered table of legal rectangles as documented in the script. Native data need to come from the user's own pilot or an authorized source. All examples supplied here are synthetic. These extensions accompany the versioned core at https://doi.org/10.5281/zenodo.24527452 .
