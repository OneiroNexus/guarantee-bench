#!/usr/bin/env python3
"""Tests for Guarantee Bench.  Run:  python -m unittest -v test_guarantee_bench

Three kinds of check:
  * exact arithmetic and the rotation identities that the published proofs rest on
    (deterministic: no tolerance, no randomness in the verdict);
  * simulation against worlds whose coverage, risk and retention are known exactly,
    judged by exact binomial tests with a declared false-fail rate of 1e-6 per assertion;
  * the command-line filter on small files.
Seeds are fixed, so every run gives the same verdicts."""
import json
import math
import os
import random
import tempfile
import unittest
from fractions import Fraction

import guarantee_bench as gb

DELTA = 1e-6


def at_least(x, r, level):
    """Is x successes in r trials consistent with a success rate >= level?"""
    return gb.binom_cdf(x, r, level) >= DELTA


def at_most(x, r, level):
    return gb.binom_sf(x, r, level) >= DELTA


def equal(x, r, level):
    return at_least(x, r, level) and at_most(x, r, level)


def random_outputs(rng, count, tied=False, force_false=False):
    outs = []
    for _ in range(count):
        k = rng.randint(1, 6)
        if tied:
            claims = [(rng.randint(0, 5), rng.random() < 0.7) for _ in range(k)]
        else:
            claims = [(rng.random(), rng.random() < 0.7) for _ in range(k)]
        if force_false:
            claims.append((rng.random(), False))
        outs.append(claims)
    return outs


def random_prompts(rng, count):
    return [[(rng.expovariate(1.0), rng.random() < 0.5) for _ in range(rng.randint(1, 4))]
            for _ in range(count)]


class Arithmetic(unittest.TestCase):
    def test_rank_known_values(self):
        self.assertEqual(gb.conformal_rank(20, 0.1), 19)
        self.assertEqual(gb.conformal_rank(9, 0.1), 9)
        self.assertEqual(gb.conformal_rank(4, 0.1), 5)
        self.assertEqual(gb.conformal_rank(25, 0.2), 21)
        self.assertEqual(gb.conformal_rank(49, 0.1), 45)

    def test_rank_matches_exact_rationals(self):
        for n in range(0, 301):
            for a100 in range(1, 100):
                exact = math.ceil((n + 1) * (1 - Fraction(a100, 100)))
                self.assertEqual(gb.conformal_rank(n, a100 / 100), exact, (n, a100))

    def test_alpha_must_be_inside_zero_one(self):
        for bad in (0, 1, -0.1, 1.5):
            with self.assertRaises(ValueError):
                gb.conformal_rank(10, bad)

    def test_threshold(self):
        self.assertEqual(gb.conformal_threshold([5, 1, 4, 2, 3, 9, 8, 7, 6, 0], 0.2), 8)
        self.assertEqual(gb.conformal_threshold([0.3, 0.7, 0.1], 0.1), math.inf)
        self.assertEqual(gb.filter_claims([(0.9, False), (0.8, True)], math.inf), [])
        self.assertEqual(gb.max_false_score([(0.2, True), (0.6, False), (0.4, False)]), 0.6)
        self.assertEqual(gb.max_false_score([(0.2, True)]), -math.inf)

    def test_binomial_against_exact(self):
        for n in range(0, 13):
            for p in (Fraction(1, 10), Fraction(1, 3), Fraction(7, 10)):
                for k in range(-1, n + 2):
                    exact = sum(Fraction(math.comb(n, i)) * p ** i * (1 - p) ** (n - i)
                                for i in range(0, min(k, n) + 1)) if k >= 0 else Fraction(0)
                    self.assertAlmostEqual(gb.binom_cdf(k, n, float(p)), float(exact), places=12)
                    exact_sf = sum(Fraction(math.comb(n, i)) * p ** i * (1 - p) ** (n - i)
                                   for i in range(max(k, 0), n + 1))
                    self.assertAlmostEqual(gb.binom_sf(k, n, float(p)), float(exact_sf), places=12)

    def test_beta_cdf_integer_parameters(self):
        for x in (0.1, 0.35, 0.8):
            self.assertAlmostEqual(gb.beta_cdf_int(x, 1, 1), x, places=12)
            self.assertAlmostEqual(gb.beta_cdf_int(x, 2, 1), x * x, places=12)
            self.assertAlmostEqual(gb.beta_cdf_int(x, 1, 2), 1 - (1 - x) ** 2, places=12)
            self.assertAlmostEqual(gb.beta_cdf_int(x, 3, 2), 4 * x ** 3 - 3 * x ** 4, places=12)


class RotationIdentities(unittest.TestCase):
    """Each published guarantee rests on an exchangeability argument that can be run as an exact
    identity on any multiset of n + 1 items. A correct implementation satisfies it every time."""

    def test_factuality_distinct_scores_gives_exactly_k(self):
        rng = random.Random(1)
        for _ in range(300):
            n = rng.randint(1, 40)
            alpha = rng.choice([0.05, 0.1, 0.2, 0.3, 0.5])
            outs = random_outputs(rng, n + 1, force_false=True)
            self.assertEqual(gb.rotation_coverage_count(outs, alpha), gb.conformal_rank(n, alpha))

    def test_factuality_tied_scores_gives_at_least_k(self):
        rng = random.Random(2)
        for _ in range(400):
            n = rng.randint(1, 40)
            alpha = rng.choice([0.05, 0.1, 0.2, 0.3, 0.5])
            outs = random_outputs(rng, n + 1, tied=True)
            self.assertGreaterEqual(gb.rotation_coverage_count(outs, alpha), gb.conformal_rank(n, alpha))

    def test_at_least_rule_breaks_with_ties(self):
        """Finding, not a requirement of the code: with tied scores the 'score >= t' rule can fall below k."""
        rng = random.Random(3)
        short = 0
        for _ in range(300):
            n = rng.randint(5, 40)
            outs = random_outputs(rng, n + 1, tied=True)
            if gb.rotation_coverage_count(outs, 0.1, keep=gb.filter_claims_at_least) < gb.conformal_rank(n, 0.1):
                short += 1
        self.assertGreater(short, 0)

    def test_crc_rotation_loss_at_most_alpha(self):
        rng = random.Random(4)
        for _ in range(400):
            n = rng.randint(1, 60)
            alpha = rng.choice([0.05, 0.1, 0.2, 0.3])
            items = [(rng.randint(0, 9) / 10, rng.random() < 0.6) for _ in range(n + 1)]
            self.assertLessEqual(gb.rotation_crc_loss(items, alpha), gb.as_fraction(alpha))

    def test_e_value_mean_is_exactly_one(self):
        rng = random.Random(5)
        checked = 0
        for _ in range(200):
            prompts = random_prompts(rng, rng.randint(2, 30))
            if sum(gb.worst_incorrect(p) for p in prompts) == 0:
                continue
            self.assertAlmostEqual(gb.rotation_e_mean(prompts), 1.0, places=9)
            checked += 1
        self.assertGreater(checked, 150)

    def test_e_scores_size_distortion_at_most_one(self):
        rng = random.Random(6)
        for _ in range(200):
            prompts = random_prompts(rng, rng.randint(2, 30))
            for rule, kw in (("fixed", {}), ("best", {}), ("fraction", {"fraction": 0.5})):
                self.assertLessEqual(gb.rotation_size_distortion(prompts, rule, "e", **kw), 1 + 1e-12)

    def test_p_scores_valid_at_a_fixed_level(self):
        rng = random.Random(7)
        for _ in range(300):
            n = rng.randint(1, 40)
            prompts = random_prompts(rng, n + 1)
            fst = [gb.worst_incorrect(p) for p in prompts]
            errs = 0
            for i, p in enumerate(prompts):
                sc = gb.p_scores([f for f, _ in p], fst[:i] + fst[i + 1:])
                errs += any(s <= 0.1 and not ok for s, (_, ok) in zip(sc, p))
            self.assertLessEqual(errs, math.floor(0.1 * (n + 1)))

    def test_p_scores_can_exceed_one_after_the_fact(self):
        """Finding: choosing the level after seeing p-scores can push size distortion above 1."""
        rng = random.Random(8)
        worst = max(gb.rotation_size_distortion(random_prompts(rng, 30), "best", "p") for _ in range(60))
        self.assertGreater(worst, 1.0)


class Worlds(unittest.TestCase):
    """The ground-truth formulas themselves are checked against sampling."""

    def test_claim_world_formulas(self):
        rng = random.Random(9)
        w = gb.ClaimWorld([dict(name="a", weight=1, K=5, p=0.2, cF=2.0, cT=2.0)])
        for tau in (0.2, 0.45):
            r, hits, kept, claims = 4000, 0, 0, 0
            for _ in range(r):
                _, c = w.sample(rng)
                hits += gb.max_false_score(c) <= tau
                kept += len(gb.filter_claims(c, tau))
                claims += len(c)
            self.assertTrue(equal(hits, r, w.coverage(tau)))
            self.assertTrue(equal(kept, claims, w.retention(tau)))

    def test_answer_world_formulas(self):
        rng = random.Random(10)
        for gamma in (1.0, 3.0):
            w = gb.AnswerWorld(gamma)
            r, lam = 6000, 0.4
            joint = sum(1 for _ in range(r) for c, ok in [w.sample(rng)] if c > lam and not ok)
            self.assertTrue(equal(joint, r, w.joint_error(lam)))


class Simulation(unittest.TestCase):
    def test_coverage_lower_bound(self):
        rng = random.Random(11)
        w = gb.ClaimWorld([dict(name="a", weight=1, K=5, p=0.2, cF=2.0, cT=2.0)])
        n, alpha, r, hits = 20, 0.1, 3000, 0
        for _ in range(r):
            tau = gb.calibrate_factuality([w.sample(rng)[1] for _ in range(n)], alpha)
            hits += gb.all_correct(gb.filter_claims(w.sample(rng)[1], tau))
        self.assertTrue(at_least(hits, r, 1 - alpha))

    def test_exact_coverage_with_continuous_scores(self):
        rng = random.Random(12)
        w = gb.ClaimWorld([dict(name="c", weight=1, K=4, p=0.0, cF=2.0, cT=2.0)], one_false=True)
        n, alpha, r, hits = 20, 0.1, 3000, 0
        for _ in range(r):
            tau = gb.calibrate_factuality([w.sample(rng)[1] for _ in range(n)], alpha)
            hits += gb.all_correct(gb.filter_claims(w.sample(rng)[1], tau))
        self.assertTrue(equal(hits, r, gb.conformal_rank(n, alpha) / (n + 1)))

    def test_too_few_calibration_outputs_keeps_nothing(self):
        rng = random.Random(13)
        w = gb.ClaimWorld([dict(name="a", weight=1, K=5, p=0.5, cF=2.0, cT=2.0)])
        tau = gb.calibrate_factuality([w.sample(rng)[1] for _ in range(5)], 0.1)
        self.assertEqual(tau, math.inf)

    def test_group_balanced_exact_per_group(self):
        rng = random.Random(14)
        w = gb.ClaimWorld([dict(name="easy", weight=1, K=4, p=0.0, cF=4.0, cT=2.0),
                           dict(name="hard", weight=1, K=4, p=0.0, cF=1.0, cT=2.0)], one_false=True)
        ng, alpha, r = 30, 0.1, 2000
        target = gb.conformal_rank(ng, alpha) / (ng + 1)
        hits = {"easy": 0, "hard": 0}
        for _ in range(r):
            cal, groups = [], []
            for g in ("easy", "hard"):
                for _ in range(ng):
                    cal.append(w.sample(rng, g)[1])
                    groups.append(g)
            taus = gb.calibrate_by_group(cal, groups, alpha)
            for g in ("easy", "hard"):
                hits[g] += gb.all_correct(gb.filter_claims(w.sample(rng, g)[1], gb.group_threshold(taus, g)))
        for g in hits:
            self.assertTrue(equal(hits[g], r, target), (g, hits[g]))

    def test_crc_upper_and_lower_bounds(self):
        rng = random.Random(15)
        w = gb.AnswerWorld(1.0)
        n, alpha, r, loss = 100, 0.1, 3000, 0
        for _ in range(r):
            lam = gb.crc_threshold([w.sample(rng) for _ in range(n)], alpha)
            c, ok = w.sample(rng)
            loss += gb.answered(c, lam) and not ok
        self.assertTrue(at_most(loss, r, alpha))
        self.assertTrue(at_least(loss, r, alpha - 2 / (n + 1)))

    def test_ltt_validity_and_power(self):
        rng = random.Random(16)
        w = gb.AnswerWorld(1.0)
        grid = [i / 100 for i in range(94, -1, -1)]   # fixed in advance; starts where answers exist
        n, alpha, delta, r = 2000, 0.1, 0.1, 400
        bad = answering = 0
        for _ in range(r):
            lam = gb.ltt_selective_threshold([w.sample(rng) for _ in range(n)], alpha, delta, grid)
            bad += w.selective_error(lam) > alpha
            answering += lam < 1
        self.assertTrue(at_most(bad, r, delta))
        self.assertGreater(answering, r // 2)

    def test_ltt_stops_at_the_first_acceptance(self):
        """Known answer: the published fixed-sequence procedure stops at the first non-rejection."""
        items = [(0.95, True), (0.95, False), (0.95, True)] + [(0.5, True)] * 100
        self.assertEqual(gb.ltt_selective_threshold(items, 0.1, 0.1, [0.9, 0.4]), math.inf)

    def test_ltt_exact_p_value_boundary(self):
        """Known answer: P(Bin(22, 0.1) <= 0) = 0.0985 rejects at delta = 0.1; P(Bin(22, 0.1) <= 1) = 0.339 does not."""
        clean = [(0.9, True)] * 22 + [(0.1, True)] * 10
        one_wrong = [(0.9, True)] * 21 + [(0.9, False)] + [(0.1, True)] * 10
        self.assertEqual(gb.ltt_selective_threshold(clean, 0.1, 0.1, [0.5]), 0.5)
        self.assertEqual(gb.ltt_selective_threshold(one_wrong, 0.1, 0.1, [0.5]), math.inf)

    def test_ltt_grid_must_run_strict_to_lenient(self):
        with self.assertRaises(ValueError):
            gb.ltt_selective_threshold([(0.5, True)], 0.1, 0.1, [0.1, 0.9])

    def test_e_scores_size_distortion_and_power(self):
        rng = random.Random(17)
        w = gb.ResponseWorld(m=4, beta=2.0, sigma=0.5)
        n, r = 100, 2000
        terms, kept, correct = [], 0, 0
        for _ in range(r):
            cal = [gb.worst_incorrect(w.sample(rng)) for _ in range(n)]
            test = w.sample(rng)
            sc = gb.e_scores([f for f, _ in test], cal)
            oks = [ok for _, ok in test]
            terms.append(gb.distortion_term(sc, oks, gb.posthoc_level(sc, "best")))
            kept += sum(1 for s, ok in zip(sc, oks) if ok and s <= 0.2)
            correct += sum(oks)
        m = math.fsum(terms) / r
        se = math.sqrt(math.fsum((t - m) ** 2 for t in terms) / (r - 1) / r)
        self.assertLessEqual(m - 5 * se, 1.0)
        self.assertGreater(kept / correct, POWER_FLOOR)


POWER_FLOOR = 0.9   # the unmutated release keeps 100% here (seed 17); a filter that keeps
                    # less than 90% of correct responses in this easy world has lost its power


class Exchangeability(unittest.TestCase):
    def test_null_and_shift(self):
        rng = random.Random(18)
        a = gb.ClaimWorld([dict(name="a", weight=1, K=5, p=0.15, cF=2.0, cT=2.0)])
        b = gb.ClaimWorld([dict(name="b", weight=1, K=5, p=0.45, cF=2.0, cT=2.0)])
        cal = [gb.max_false_score(a.sample(rng)[1]) for _ in range(150)]
        same = [gb.max_false_score(a.sample(rng)[1]) for _ in range(150)]
        moved = [gb.max_false_score(b.sample(rng)[1]) for _ in range(150)]
        self.assertGreater(gb.exchangeability_check(cal, same, n_perm=499, seed=1)["p_value"], 0.001)
        p_moved = gb.exchangeability_check(cal, moved, n_perm=499, seed=1)["p_value"]
        self.assertLess(p_moved, 0.01)
        self.assertGreaterEqual(p_moved, 1 / 500)   # a permutation p-value is never below 1/(n_perm + 1)

    def test_flat_picture_with_ties(self):
        """With many tied scores (outputs with no incorrect claim), random tie-breaking keeps the
        rank histogram flat under exchangeability; without it, the lowest bin overflows."""
        rng = random.Random(20)
        a = gb.ClaimWorld([dict(name="a", weight=1, K=5, p=0.15, cF=2.0, cT=2.0)])
        cal = [gb.max_false_score(a.sample(rng)[1]) for _ in range(150)]
        same = [gb.max_false_score(a.sample(rng)[1]) for _ in range(300)]
        res = gb.exchangeability_check(cal, same, n_perm=1, seed=2)
        self.assertLess(res["counts"][0] / res["expected"][0], 2.0)


class CommandLine(unittest.TestCase):
    def test_filter_files(self):
        rng = random.Random(19)
        w = gb.ClaimWorld([dict(name="a", weight=1, K=4, p=0.2, cF=2.0, cT=2.0)])
        with tempfile.TemporaryDirectory() as d:
            cal_p, test_p, out_p = (os.path.join(d, x) for x in ("cal.jsonl", "test.jsonl", "kept.jsonl"))
            with open(cal_p, "w", encoding="utf-8") as fh:
                for i in range(40):
                    fh.write(json.dumps({"id": i, "claims": [{"text": "c", "score": s, "correct": ok}
                                                             for s, ok in w.sample(rng)[1]]}) + "\n")
            tests = [w.sample(rng)[1] for _ in range(10)]
            with open(test_p, "w", encoding="utf-8") as fh:
                for i, c in enumerate(tests):
                    fh.write(json.dumps({"id": i, "claims": [{"text": "c", "score": s} for s, _ in c]}) + "\n")
            s = gb.run_filter(cal_p, test_p, 0.1, out_p)
            self.assertEqual(s["n_calibration"], 40)
            self.assertEqual(s["k"], 37)
            with open(out_p, encoding="utf-8") as fh:
                rows = [json.loads(x) for x in fh]
            self.assertEqual(len(rows), 10)
            for row, c in zip(rows, tests):
                self.assertTrue(all(k["score"] > s["threshold"] for k in row["kept"]))
                self.assertEqual(row["removed"] + len(row["kept"]), len(c))
            with self.assertRaises(SystemExit):
                gb.run_filter(cal_p, test_p, 0.1, out_p)


if __name__ == "__main__":
    unittest.main()
