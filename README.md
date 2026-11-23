# Guarantee Bench

A synthetic benchmark for checking implementations of statistical guarantees on scored subclaims. Research and original reference implementation by Aaron Michael Weis, Oneiron Holdings LLC. [Author ORCID](https://orcid.org/0009-0002-1906-148X)

[Read the synthetic-edition paper](PAPER.md) · [Simulation results](results/RESULTS.md) · [Verification record](checks.json)

## Run the benchmark

The tools use the Python standard library. The full run uses seed 20261117 by default. Run these commands from the folder containing `guarantee_bench.py` with a new output directory and a new mutation report filename. The public edition's tests and measurements are documented in the paper.

```sh
python -B -m unittest -v test_guarantee_bench
python -B guarantee_bench.py bench --out repeat-results --seed 20261117
python -B mutate.py --out repeat-mutations.txt
```

For a shorter simulation use `python -B guarantee_bench.py bench --quick --out quick-results`. The quick option changes the number of draws. The numerical results in the paper use the full run. The only non-deterministic field in the full seeded simulation's JSON is elapsed runtime. Identical numeric output across Python versions is not guaranteed.

## Try the file filter

The example files contain invented claim IDs and scores. Each calibration row needs a list of claims with a finite numeric `score` and a boolean `correct` label. Each test claim needs a finite score. A higher score favors keeping the claim. An optional string `group` permits calibration per group with `--by-group`. The threshold for an unseen group keeps no claims. The command writes a new output file.

```sh
python -B guarantee_bench.py filter --cal example/calibration.jsonl --test example/test.jsonl --alpha 0.25 --out kept.jsonl
```

The filter retains the claims strictly above the threshold. In its JSON, a null global threshold indicates that no finite threshold can be obtained, and the tool keeps no claims. Under group calibration the summary’s global threshold is also null and each output row records its own threshold. A threshold of `"-inf"` indicates that all finite scores clear. An empty kept set counts as all kept claims correct, so inspect the retained count alongside the coverage number.

## Source and reuse

Aaron Michael Weis is the author of this research and its benchmark design. The manuscript distinguishes his implementation work from the authors of the statistical methods. The paper and the original documentation and result tables are licensed under [CC BY 4.0](LICENSE-PAPER.txt). The original Python code and synthetic example fixtures are licensed under the [MIT License](LICENSE-CODE.txt).
