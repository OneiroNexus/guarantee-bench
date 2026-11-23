#!/usr/bin/env python3
"""Run the 28 prespecified program variants against the local synthetic tests.

Use a fresh output path. A timed-out or interrupted test run is reported as an
incomplete check. The original source files are read-only inputs.
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SRC, TESTS = "guarantee_bench.py", "test_guarantee_bench.py"

MUTANTS = [
    ("M01", "rank: drop the finite-sample correction, n(1-a) for (n+1)(1-a)",
     "return math.ceil((n + 1) * (1 - a))", "return math.ceil(n * (1 - a))"),
    ("M02", "rank: floor instead of ceiling",
     "return math.ceil((n + 1) * (1 - a))", "return math.floor((n + 1) * (1 - a))"),
    ("M03", "rank: one rank too conservative",
     "return math.ceil((n + 1) * (1 - a))", "return math.ceil((n + 1) * (1 - a)) + 1"),
    ("M04", "threshold: clamp k to n instead of keeping nothing when k > n",
     "    if k > n:\n        return INF\n    return sorted(scores)[k - 1]",
     "    if k > n:\n        k = n\n    return sorted(scores)[k - 1]"),
    ("M05", "threshold: index one rank too low",
     "return sorted(scores)[k - 1]", "return sorted(scores)[k - 2]"),
    ("M06", "levels in floating point instead of exact rationals",
     "return Fraction(str(x))", "return x"),
    ("M07", "nonconformity from correct claims instead of incorrect ones",
     "if not ok and score > r:", "if ok and score > r:"),
    ("M08", "nonconformity never updated (minimum instead of maximum)",
     "if not ok and score > r:", "if not ok and score < r:"),
    ("M09", "keep claims with score >= tau (tie convention)",
     "return [c for c in claims if c[0] > tau]", "return [c for c in claims if c[0] >= tau]"),
    ("M10", "output judged correct if any claim is correct",
     "return all(ok for _, ok in claims)", "return any(ok for _, ok in claims)"),
    ("M11", "group-balanced: pool all groups into one calibration set",
     "pools.setdefault(g, []).append(max_false_score(claims))",
     "pools.setdefault(\"pooled\", []).append(max_false_score(claims))"),
    ("M12", "risk control: drop the B/(n+1) term",
     "budget = as_fraction(alpha) * (n + 1) - 1", "budget = as_fraction(alpha) * (n + 1)"),
    ("M13", "risk control: ceiling instead of floor on the error budget",
     "m = math.floor(budget)", "m = math.ceil(budget)"),
    ("M14", "risk control: threshold one wrong item too lenient",
     "return wrong[m]", "return wrong[m - 1]"),
    ("M15", "answer when confidence >= lam (breaks right-continuity)",
     "    return conf > lam", "    return conf >= lam"),
    ("M16", "Learn then Test: keep testing after the first non-rejection",
     "        if p <= delta:\n            lam_hat = lam\n        else:\n            break",
     "        if p <= delta:\n            lam_hat = lam"),
    ("M17", "Learn then Test: p-value off by one (P(X < e) instead of P(X <= e))",
     "p = binom_cdf(e_j, n_j, level) if n_j else 1.0", "p = binom_cdf(e_j - 1, n_j, level) if n_j else 1.0"),
    ("M18", "binomial CDF off by one",
     "return _binom_range(0, k, n, p)", "return _binom_range(0, k - 1, n, p)"),
    ("M19", "e-value: n instead of n + 1",
     "out.append((n + 1) * f / den if den > 0 else 0.0)", "out.append(n * f / den if den > 0 else 0.0)"),
    ("M20", "e-value: test score left out of the denominator",
     "den = f + s", "den = s"),
    ("M21", "e-scores: f* taken over all responses, not only incorrect ones",
     "if not ok and f > m:", "if f > m:"),
    ("M22", "p-scores: drop the +1",
     "return [(1 + len(srt) - bisect_left(srt, f)) / (n + 1) for f in test_fs]",
     "return [(len(srt) - bisect_left(srt, f)) / (n + 1) for f in test_fs]"),
    ("M23", "size distortion: error counted only strictly below the level",
     "err = any(s <= level and not ok for s, ok in zip(scores, oks))",
     "err = any(s < level and not ok for s, ok in zip(scores, oks))"),
    ("M24", "rank histogram: ties not broken at random",
     "r = lo + rng.randint(0, hi - lo)", "r = lo"),
    ("M25", "permutation p-value without the +1",
     "\"p_value\": (1 + ge) / (1 + n_perm)}", "\"p_value\": ge / n_perm}"),
    ("M26", "ground truth: wrong coverage formula in the claim world",
     "c = (1.0 - g[\"p\"] * (1.0 - self._gf(tau, g[\"cF\"]))) ** g[\"K\"]",
     "c = (1.0 - g[\"p\"] * self._gf(tau, g[\"cF\"])) ** g[\"K\"]"),
    ("M27", "ground truth: Beta CDF parameters swapped",
     "return binom_sf(a, a + b - 1, x)", "return binom_sf(b, a + b - 1, x)"),
    ("M28", "ground truth: wrong joint-error formula in the answer world",
     "return (1.0 - x) - (1.0 - x ** (g + 1)) / (g + 1)", "return (1.0 - x) - (1.0 - x ** g) / (g + 1)"),
]

# Mutants whose survival is expected, with the reason (judged before the run; checked after it).
# First run (17 Nov 2026): 24 of 28 killed. M16 and M17 survived (gaps, closed by two known-answer
# tests); M24 and M25 survived (closed by a flat-picture test and a p-value floor test).
EXPECTED_SURVIVORS = {}


def run_suite(folder, timeout):
    t0 = time.time()
    try:
        p = subprocess.run([sys.executable, "-m", "unittest", "-q", "test_guarantee_bench"], cwd=folder,
                           capture_output=True, text=True, timeout=timeout)
        rc, tail = p.returncode, (p.stderr.strip().splitlines() or [""])[-1]
    except subprocess.TimeoutExpired:
        rc, tail = -9, "timeout"
    return rc, time.time() - t0, tail


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, help="new result file")
    options = parser.parse_args()
    if os.path.exists(options.out):
        raise SystemExit("Choose a new mutation report path.")
    timeout = 600
    src = open(os.path.join(HERE, SRC), encoding="utf-8").read()
    with tempfile.TemporaryDirectory() as d:
        shutil.copy(os.path.join(HERE, SRC), d)
        shutil.copy(os.path.join(HERE, TESTS), d)
        rc, secs, tail = run_suite(d, timeout)
    if rc != 0:
        sys.exit("the unmutated suite must pass first: " + tail)
    lines = ["Guarantee Bench mutation testing", "Unmutated suite: PASS in %.1f s" % secs, ""]
    killed = survived = skipped = 0
    for mid, desc, old, new in MUTANTS:
        if src.count(old) != 1:
            lines.append("%s NOT APPLIED (%d matches)  %s" % (mid, src.count(old), desc))
            skipped += 1
            continue
        with tempfile.TemporaryDirectory() as d:
            shutil.copy(os.path.join(HERE, TESTS), d)
            with open(os.path.join(d, SRC), "w", encoding="utf-8") as fh:
                fh.write(src.replace(old, new))
            rc, secs, tail = run_suite(d, timeout)
        if rc == -9:
            lines.append("%s INCOMPLETE: test process timed out" % mid)
            skipped += 1
        elif rc != 0 and tail.startswith("FAILED"):
            killed += 1
            lines.append("%s KILLED    %5.1f s  %s  [%s]" % (mid, secs, desc, tail))
        elif rc != 0:
            lines.append("%s INCOMPLETE: test process failed (%s)" % (mid, tail))
            skipped += 1
        else:
            survived += 1
            why = EXPECTED_SURVIVORS.get(mid, "NOT EXPECTED: the suite has a gap here")
            lines.append("%s SURVIVED  %5.1f s  %s  -- %s" % (mid, secs, desc, why))
        print(lines[-1], flush=True)
    total = killed + survived
    unexpected = [l for l in lines if "NOT EXPECTED" in l]
    lines += ["", "Killed %d of %d applied mutants (%.0f%%); %d not applied." % (killed, total, (100.0 * killed / total if total else 0.0), skipped),
              "Unexpected survivors: %d" % len(unexpected)]
    with open(options.out, "x", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines[-2:]))
    return int(bool(unexpected or skipped))


if __name__ == "__main__":
    sys.exit(main())
