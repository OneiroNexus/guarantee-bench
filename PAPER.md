# Guarantee Bench

## Checking statistical guarantees on scored language-model answers with synthetic controls

**Aaron Michael Weis · Oneiron Holdings LLC · Synthetic edition 0.1.0 · 22 November 2026**  
[ORCID 0009-0002-1906-148X](https://orcid.org/0009-0002-1906-148X)

### Abstract

A filter can make every statement it keeps correct by keeping nothing. That makes the probability of a fully correct retained answer difficult to interpret without a second number, the proportion retained. This paper presents a runnable reference benchmark for five published statistical procedures. It checks how an implementation handles finite-sample ranks, tied scores, subgroup calibration, abstention and a tolerance chosen after scoring. Every experiment uses generated claims or confidence values. In a mixture of two synthetic populations, pooled calibration produced mean whole-answer coverage of 0.8998 across both groups and 0.7478 in the harder group. Calibrating the groups separately raised harder-group coverage to 0.9049, with its mean retained share falling from 0.5668 to 0.4435. An independent arithmetic recount found 222 cases in a grid of 100,050 in which a floating-point ceiling gave a different rank from exact arithmetic. The study demonstrates how to test an implementation and read a guarantee at its stated scope.

### 1. The question

Suppose a system breaks one generated answer into several claims. Each gets a score and, for calibration, a correctness label. Which statements may the system keep, what probability does a theorem bound, and how much of the answer survives the filter? The answer depends on the unit being judged. An error per incoming question, an error among questions that received an answer, and a retained answer containing any error have different denominators.

The contribution here is the assembly of exact arithmetic checks, generated populations with known behavior, seeded simulations and deliberately altered program variants into one benchmark. The statistical procedures have their own authors. Their results provide targets against which an implementation can be checked. [Conformal factuality](https://arxiv.org/abs/2402.10978) · [Conformal prediction and group calibration](https://arxiv.org/abs/2107.07511) · [Conformal risk control](https://arxiv.org/abs/2208.02814) · [Learn then Test](https://arxiv.org/abs/2110.01052) · [E-scores](https://arxiv.org/abs/2510.25770)

### 2. How one filter works

For each calibration answer, the benchmark records the largest score among its incorrect claims. If all claims in that answer are correct, its score for calibration is negative infinity. With $n$ calibration answers and error tolerance $\alpha$, the implementation selects the $k$-th smallest calibration score, where $k=\lceil(n+1)(1-\alpha)\rceil$. If $k>n$, the threshold is positive infinity. A new answer retains the claims whose scores strictly exceed this threshold.

With exchangeable calibration and test answers, this construction targets marginal coverage of at least $1-\alpha$. Here, coverage means that *all retained claims in a new answer are correct*. Exchangeability means that the joint distribution of these answers is invariant to permuting their positions. The probability averages over the calibration sample and the new answer. It offers no separate guarantee for a particular calibration sample or an arbitrary subgroup. Under this definition an empty answer counts as covered. Retention measures what that convention can hide. [Mohri and Hashimoto](https://arxiv.org/abs/2402.10978) · [Angelopoulos and Bates](https://arxiv.org/abs/2107.07511)

The small accompanying example makes the tie rule visible. Four calibration answers have incorrect-claim scores of 0.1, 0.2, 0.3 and 0.4. At $\alpha=0.25$, the rank is four and the threshold is 0.4. One of the two test answers retains a single claim scored 0.6. The other retains none. Both test answers satisfy the all-retained-correct criterion, while only one of the four test claims survives. [Example inputs and reproduction instructions](README.md#try-the-file-filter)

### 3. Experiments and controls

The exact checks rotate each element of a fixed pool through the test position and calibrate on the remaining elements. Enumerating those positions makes finite-sample rank identities directly testable. Separate cases exercise ties, the empty-retention rule, group-specific thresholds, risk-control loss, fixed-order testing and score-based error bounds.

The seeded experiments use synthetic answer populations. In the claim populations, correctness labels and scores are generated from declared distributions. For each calibrated threshold, the program can calculate the implied coverage and retention in its own population. In the abstention experiment, the policy can decline an entire question. The suite records both the chance of answering incorrectly on an incoming question and the error rate conditional on answering. The e-score experiment varies how the error tolerance is chosen. All five procedures operate on scores and labels. No model inference was needed.

The full simulation uses seed 20261117. The source and tables record the draw counts and parameters for each experiment. The 29 unit tests include deterministic identities, fixed-seed simulation checks and a command-line fixture. Separately, the mutation runner applies 28 specified changes to scratch copies of the implementation and reruns the suite on each copy. A detected change is one for which the tests actually report a failure. [Source](guarantee_bench.py) · [Tests](test_guarantee_bench.py) · [Mutation definitions](mutate.py)

### 4. Results

| Check | Expected behavior or comparison | Observed in this run |
|---|---|---|
| Continuous scores, 20 calibration answers, tolerance 0.10 | The exact rotation target is 19/21, or 0.90476. | Across 4,000 draws, mean coverage was 0.90463. Monte Carlo standard error was 0.00097. |
| Mixed population, 300 calibration answers, tolerance 0.10 | Pooled calibration has a marginal target. | Across 2,000 draws, pooled coverage was 0.8998 overall, 0.9649 in the easier group and 0.7478 in the harder group. |
| The same population with group-specific calibration | Compare retention at both groups’ coverage levels. | Easier-group coverage and retention were 0.9020 and 0.9380. Harder-group coverage and retention were 0.9049 and 0.4435. |
| Distribution shift | Change the generator’s per-claim error probability from 0.10 to 0.25 after calibration. | Across 400 repetitions, mean coverage was 0.8991 under matched conditions and 0.7645 after the change. |
| Abstention in the population where correctness probability equals confidence | Joint and selective error use different denominators. | Conformal risk control yielded mean joint error of 0.0998, selective error of 0.2233 and an answer rate of 0.4466. |
| Numerical boundary | Compare a binary-floating ceiling with an integer ceiling on the same rational inputs. | The ranks differed on 222 of 100,050 parameter pairs. At $n=24$ and $\alpha=0.44$, the exact rank is 14 and the floating computation returned 15. |

The 0.8998 pooled estimate is the result of a finite simulation. A small excursion around the marginal target does not resolve whether a theorem holds. The harder-group result answers a different question that the pooled target leaves open. In the shift example, the calibration population itself has changed. Each observation therefore carries its sampling unit, generator and denominator.

All 29 unit tests passed. The suite detected all 28 prespecified program changes, with none skipped or timed out. An independent integer-arithmetic calculation reproduced the 222 rank discrepancies. The file filter also reproduced the two fictional examples in Section 2, refused an existing output path, and rejected a nonfinite test score before creating an output file.

The complete results include the Learn then Test decision grid and the e-score and p-score comparisons under the chosen tolerance rules. [Full tables](results/RESULTS.md) · [Machine-readable results](results/RESULTS.json) · [Checks from this run](checks.json)

### 5. What follows for readers

Before relying on a statistical guarantee, write down the unit counted as a success, the population over which the bound applies, and the share of material returned. A system can meet a whole-answer target after removing information its users value. Report retention and the frequency of empty answers alongside correctness. A subgroup with few calibration examples may require its own analysis. A change in prompts, scoring behavior or user population can change the distribution on which calibration depended.

These experiments use invented correctness labels and specified scoring distributions. They measure implementation behavior on those fixtures. They establish no error rate for a real language model, human annotator, deployment or commercial product. Mutation testing measures the response to the particular changes made. The paper does not estimate how often those defects occur elsewhere.

### Authorship and reuse

Aaron Michael Weis conducted and directed the research and designed the benchmark study. Claude assisted with drafting and coding the original reference implementation. Codex prepared this synthetic edition and ran the documented checks. Original source code and fictional example fixtures are MIT licensed. The paper and its original tables use CC BY 4.0. Published statistical methods and cited works retain their respective authorship and terms. [Code license](LICENSE-CODE.txt) · [Paper license](LICENSE-PAPER.txt)

### References

1. Christopher Mohri and Tatsunori Hashimoto. *Language Models with Conformal Factuality Guarantees* (2024). [Paper](https://arxiv.org/abs/2402.10978)
2. Anastasios N. Angelopoulos and Stephen Bates. *A Gentle Introduction to Conformal Prediction and Distribution-Free Uncertainty Quantification* (2021, revised 2022). [Paper](https://arxiv.org/abs/2107.07511)
3. Anastasios N. Angelopoulos, Stephen Bates, Adam Fisch, Lihua Lei and Tal Schuster. *Conformal Risk Control* (2022, revised 2025). [Paper](https://arxiv.org/abs/2208.02814)
4. Anastasios N. Angelopoulos, Stephen Bates, Emmanuel J. Candès, Michael I. Jordan and Lihua Lei. *Learn then Test: Calibrating Predictive Algorithms to Achieve Risk Control* (2021, revised 2022). [Paper](https://arxiv.org/abs/2110.01052)
5. Guneet S. Dhillon, Javier González, Teodora Pandeva and Alicia Curth. *E-Scores for (In)Correctness Assessment of Generative Model Outputs* (2025 preprint). [Paper](https://arxiv.org/abs/2510.25770)
