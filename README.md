# Reproducible spatial sampling, version 0.5.0

## Human disc suspensions, from official files to both endpoints

Run from the extracted software directory using Python 3.12 and the pinned dependencies in requirements-lock.txt.

```text
python reproduce_omix.py --archives public_inputs --out public_results
```

This command downloads any missing archive from NGDC, verifies its frozen size and SHA-256, reads each designated matrix exactly once, and checks the complete result independently by Newton identities. To reuse a verified local mirror with networking disabled, add `--offline`. Rerun only the independent calculation with `python check_public.py public_results`.

The catalogue labels archive 02 as Patient 1. Inspection of the November 25, 2026 downloads found five matrix directories inside file 02 and one in each other file. The selected inputs follow the five official catalogue assignments below. Across the other four libraries, all twelve duplicated components (matrix, features, barcodes) match by their decompressed SHA-256. File 02 alone remains a supported legacy input when its entire digest matches. The release command uses all five official files. These SHA-256 values were calculated from downloaded bytes and are supplied for independent verification. [NGDC record](https://ngdc.cncb.ac.cn/omix/release/OMIX017533), accessed 2026-11-25.

| Official ZIP | NGDC catalogue label | Directory scored | Matrix members in ZIP | SHA-256 (downloaded ZIP) |
|---|---|---|---:|---|
| [OMIX017533-01.zip](https://download.cncb.ac.cn/OMIX/OMIX017533/OMIX017533-01.zip) | Patient 4 | IVD4/matrix | 1 | `fc3a6e4302556a7f7eb7d5fe4a79487ea651dbf58a7f0b9fcf0f0724d5425cd0` |
| [OMIX017533-02.zip](https://download.cncb.ac.cn/OMIX/OMIX017533/OMIX017533-02.zip) | Patient 1 | IVD1/matrix | 5 | `f0fd52d8911bda8311b9e448c22e6fb6f5303d4e53e908cb1e81e006be911fc9` |
| [OMIX017533-03.zip](https://download.cncb.ac.cn/OMIX/OMIX017533/OMIX017533-03.zip) | Patient 5 | IVD5/matrix | 1 | `5d91e5b593631edfc1c59a60dd153028df439ff1688f81866f5686cc8a6935a0` |
| [OMIX017533-04.zip](https://download.cncb.ac.cn/OMIX/OMIX017533/OMIX017533-04.zip) | Patient 2 | IVD2/matrix | 1 | `1605fdaa128fb6ed272746ceb30b3e6c99c520594ed3312a6bd59e3bc2a37c5c` |
| [OMIX017533-05.zip](https://download.cncb.ac.cn/OMIX/OMIX017533/OMIX017533-05.zip) | Patient 3 | IVD3/matrix | 1 | `a2d398c2dbf192255b6483655c8c79d5ec0aa2c48ab6eab85180bf2258b27717` |

`strict_dual` is the probability of a full twelve-barcode request containing both assay targets anywhere. `within_barcode_strict` additionally requires at least one sampled barcode to detect both targets. Both denominators include requests that fail the 521-nontarget-UMI eligibility gate. The latter readout measures barcode-level co-detection and is sensitive to capture and multiplets. The two columns are probabilities on [0,1]. Do not interpret them as a person-level phenotype.

## Compare layouts on a synthetic or investigator-supplied frame

```text
python -m pip install -r requirements-lock.txt
python -m unittest discover -v
python run_sections.py --input synthetic_pilot.csv --out results
python resource_analysis.py --input synthetic_pilot.csv --fixed600 results/fixed600.csv --out results
python plot_matched.py --input synthetic_pilot.csv --results results --out figures
```

The included grid is simulated. To use your own registered measurements, provide a compatible site table and a fresh output directory. The study rat site-level matrix is not distributed. The notebook runs on the synthetic grid. The OMIX commands above retrieve the public data directly from the originating repository. Set OMIX_ARCHIVES to a local mirror to enable the optional real-file test.

## Compare choices for a measured pilot

```text
python choose_design.py --input synthetic_pilot.csv --section Synthetic-A --targets Acan Col2a1 --budget 600 --sites 6 12 20 30 --max-sites 12 --spot-diameter-um 55 --out pilot_plans.csv
```

Supply an eligible, registered table with finite nonnegative integer target counts and coordinates, dark/light strata and site IDs. All candidate plans in this example use the same per-target eligibility frame. Both the error-minimising and detection-maximising designs are reported. A pilot frontier minimises contrast RMSE and summed spot footprint and maximises request-level dual recovery. Cost, physical extraction and actual library yield need separate measurements. Per-site footprint assumes non-overlapping circular sites; it excludes intervening tissue and is not the bounding envelope. Simulations condition on captured molecules.

The MIT licence applies to original software. Data and image source terms are listed separately. See DATA_SOURCES.md for provenance and Worked_example.ipynb for a short example.
