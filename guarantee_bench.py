#!/usr/bin/env python3
# Guarantee Bench 0.1.0
# Aaron Michael Weis (Oneironexus · Oneiron Holdings LLC · ORCID 0009-0002-1906-148X), 2026.
# Standard library only. No network access, no model calls, no telemetry.
"""Guarantee Bench. Check published statistical procedures on fictional, scored answers.

Original reference implementation by Aaron Michael Weis, Oneiron Holdings LLC.
The synthetic edition contains exact finite-sample identities, generated control
worlds and a file filter. Method references and assumptions appear in PAPER.md.
"""

import argparse
import json
import math
import os
import random
import sys
import time
from bisect import bisect_left, bisect_right
from fractions import Fraction

__version__ = "0.1.0"
INF = math.inf

# ------------------------------------------------------------------ exact arithmetic


def as_fraction(x):
    """Exact rational value of a level such as 0.1, taken from its shortest decimal form."""
    if isinstance(x, Fraction):
        return x
    if isinstance(x, int):
        return Fraction(x)
    return Fraction(str(x))


def conformal_rank(n, alpha):
    """k = ceil((n + 1)(1 - alpha)), computed with exact rationals."""
    a = as_fraction(alpha)
    if not 0 < a < 1:
        raise ValueError("alpha must lie strictly between 0 and 1")
    if n < 0:
        raise ValueError("n must be non-negative")
    return math.ceil((n + 1) * (1 - a))


def conformal_threshold(scores, alpha):
    """The k-th smallest calibration score; +inf when k > n (no finite threshold exists)."""
    n = len(scores)
    k = conformal_rank(n, alpha)
    if k > n:
        return INF
    return sorted(scores)[k - 1]


# ------------------------------------------------------------------ conformal factuality
# A claim is a pair (score, ok): ok is True for a correct claim, False for an incorrect one.


def max_false_score(claims):
    """Nonconformity of one output: the highest score held by an incorrect claim (-inf if none)."""
    r = -INF
    for score, ok in claims:
        if not ok and score > r:
            r = score
    return r


def filter_claims(claims, tau):
    """Keep the claims whose score is strictly above tau."""
    return [c for c in claims if c[0] > tau]


def filter_claims_at_least(claims, tau):
    """The 'score >= tau' convention, kept for comparison only (it needs distinct scores)."""
    return [c for c in claims if c[0] >= tau]


def all_correct(claims):
    return all(ok for _, ok in claims)


def calibrate_factuality(cal_outputs, alpha):
    return conformal_threshold([max_false_score(o) for o in cal_outputs], alpha)


def calibrate_by_group(cal_outputs, cal_groups, alpha):
    pools = {}
    for claims, g in zip(cal_outputs, cal_groups):
        pools.setdefault(g, []).append(max_false_score(claims))
    return {g: conformal_threshold(v, alpha) for g, v in pools.items()}


def group_threshold(taus, g):
    """A group never seen in calibration keeps nothing."""
    return taus.get(g, INF)


# ------------------------------------------------------------------ conformal risk control
# An item is (confidence, ok). The policy answers when confidence > lam and abstains otherwise.
# Loss L(lam) = 1{answered and wrong}: non-increasing and right-continuous in lam, bounded by B = 1.


def answered(conf, lam):
    return conf > lam


def crc_threshold(cal_items, alpha):
    """lam_hat = inf{lam : (sum_i L_i(lam) + B) / (n + 1) <= alpha}, B = 1.
    Returns -inf (answer everything) or +inf (answer nothing) at the ends."""
    n = len(cal_items)
    budget = as_fraction(alpha) * (n + 1) - 1
    if budget < 0:
        return INF
    m = math.floor(budget)
    wrong = sorted((c for c, ok in cal_items if not ok), reverse=True)
    if len(wrong) <= m:
        return -INF
    return wrong[m]


# ------------------------------------------------------------------ binomial and beta


def _binom_range(lo, hi, n, p):
    """P(lo <= Bin(n, p) <= hi), summed in log space."""
    lo, hi = max(lo, 0), min(hi, n)
    if lo > hi:
        return 0.0
    if p <= 0.0:
        return 1.0 if lo == 0 else 0.0
    if p >= 1.0:
        return 1.0 if hi == n else 0.0
    lp, lq, lc = math.log(p), math.log1p(-p), math.lgamma(n + 1)
    ts = [lc - math.lgamma(i + 1) - math.lgamma(n - i + 1) + i * lp + (n - i) * lq
          for i in range(lo, hi + 1)]
    mx = max(ts)
    return min(1.0, math.exp(mx) * math.fsum(math.exp(t - mx) for t in ts))


def binom_cdf(k, n, p):
    """P(Bin(n, p) <= k)."""
    return _binom_range(0, k, n, p)


def binom_sf(k, n, p):
    """P(Bin(n, p) >= k)."""
    return _binom_range(k, n, n, p)


def beta_cdf_int(x, a, b):
    """CDF of Beta(a, b) at x for positive integers a, b: P(Bin(a + b - 1, x) >= a)."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    return binom_sf(a, a + b - 1, x)


# ------------------------------------------------------------------ Learn then Test


def ltt_selective_threshold(cal_items, alpha, delta, grid):
    """Certify P(wrong | answered) <= alpha with probability >= 1 - delta.
    grid: thresholds fixed before the data are seen, from strict (few answers) to lenient.
    H_j: selective risk at grid[j] > alpha, tested in that fixed order at level delta with
    p_j = P(Bin(n_j, alpha) <= e_j) on the n_j answered calibration items (e_j of them wrong).
    Returns the last threshold rejected, or +inf (answer nothing)."""
    if any(b > a for a, b in zip(grid, grid[1:])):
        raise ValueError("grid must run from strict to lenient (descending)")
    items = sorted(cal_items, key=lambda t: -t[0])
    level = float(as_fraction(alpha))
    lam_hat = INF
    i = n_j = e_j = 0
    for lam in grid:
        while i < len(items) and answered(items[i][0], lam):
            n_j += 1
            if not items[i][1]:
                e_j += 1
            i += 1
        p = binom_cdf(e_j, n_j, level) if n_j else 1.0
        if p <= delta:
            lam_hat = lam
        else:
            break
    return lam_hat


# ------------------------------------------------------------------ e-scores
# A prompt carries responses [(f, ok)], f >= 0 from a fixed scorer (larger means more likely correct).


def worst_incorrect(responses):
    """f* for a calibration prompt: the largest f among its incorrect responses (0 if none)."""
    m = 0.0
    for f, ok in responses:
        if not ok and f > m:
            m = f
    return m


def e_values(test_fs, cal_fstars):
    """E(y) = (n + 1) f(y) / (f(y) + sum_i f*_i). The e-score is 1 / E(y)."""
    n = len(cal_fstars)
    s = math.fsum(cal_fstars)
    out = []
    for f in test_fs:
        den = f + s
        out.append((n + 1) * f / den if den > 0 else 0.0)
    return out


def e_scores(test_fs, cal_fstars):
    return [1.0 / e if e > 0 else INF for e in e_values(test_fs, cal_fstars)]


def p_scores(test_fs, cal_fstars):
    """Conformal p-scores, for comparison: (1 + #{i : f*_i >= f}) / (n + 1)."""
    n = len(cal_fstars)
    srt = sorted(cal_fstars)
    return [(1 + len(srt) - bisect_left(srt, f)) / (n + 1) for f in test_fs]


def posthoc_level(scores, rule, fraction=0.5, fixed=0.1):
    """Tolerance chosen after seeing the scores (rules of this release, for testing)."""
    if rule == "fixed":
        return fixed
    srt = sorted(scores)
    if rule == "best":          # keep the single best response
        return srt[0]
    if rule == "fraction":      # keep the best fraction of the responses
        return srt[max(0, math.ceil(fraction * len(srt)) - 1)]
    raise ValueError("unknown rule: %s" % rule)


def distortion_term(scores, oks, level):
    """1{an incorrect response is kept at this level} / level."""
    err = any(s <= level and not ok for s, ok in zip(scores, oks))
    if not err or level == INF:
        return 0.0
    return 1.0 / level


# ------------------------------------------------------------------ rotation identities


def rotation_coverage_count(outputs, alpha, keep=filter_claims):
    """Over n + 1 outputs, put each one in the test seat in turn and calibrate on the rest.
    Returns how many come out fully correct. The exchangeability argument says the count is
    at least k = ceil((n + 1)(1 - alpha)), and exactly k when the nonconformity scores are distinct."""
    rs = [max_false_score(o) for o in outputs]
    count = 0
    for i, o in enumerate(outputs):
        tau = conformal_threshold(rs[:i] + rs[i + 1:], alpha)
        if all_correct(keep(o, tau)):
            count += 1
    return count


def rotation_crc_loss(items, alpha):
    """Mean test loss over the n + 1 rotations; the conformal risk control argument says <= alpha."""
    total = 0
    for i, (c, ok) in enumerate(items):
        lam = crc_threshold(items[:i] + items[i + 1:], alpha)
        if answered(c, lam) and not ok:
            total += 1
    return Fraction(total, len(items))


def rotation_e_mean(prompts):
    """Mean over rotations of the e-value of each prompt's worst incorrect response (1 when sum f* > 0)."""
    fst = [worst_incorrect(p) for p in prompts]
    es = [e_values([fst[i]], fst[:i] + fst[i + 1:])[0] for i in range(len(prompts))]
    return math.fsum(es) / len(es)


def rotation_size_distortion(prompts, rule, kind="e", **kw):
    """Mean over rotations of 1{error} / level for a post-hoc rule (e-scores: <= 1)."""
    fst = [worst_incorrect(p) for p in prompts]
    vals = []
    for i, p in enumerate(prompts):
        cal = fst[:i] + fst[i + 1:]
        fs = [f for f, _ in p]
        oks = [ok for _, ok in p]
        sc = e_scores(fs, cal) if kind == "e" else p_scores(fs, cal)
        vals.append(distortion_term(sc, oks, posthoc_level(sc, rule, **kw)))
    return math.fsum(vals) / len(vals)


# ------------------------------------------------------------------ exchangeability check


def _rank_counts(test, ref_sorted, bins, rng):
    n = len(ref_sorted)
    counts = [0] * bins
    for x in test:
        lo, hi = bisect_left(ref_sorted, x), bisect_right(ref_sorted, x)
        r = lo + rng.randint(0, hi - lo)          # ties broken at random: r in 0..n
        counts[r * bins // (n + 1)] += 1
    return counts


def exchangeability_check(cal, test, bins=10, n_perm=999, seed=0):
    """One batch: rank histogram of the test scores against the calibration scores, and a
    permutation p-value for its chi-square distance from the shape expected under exchangeability.
    Valid when the n + T scores are exchangeable. It says nothing about batches it did not see,
    and a large p-value does not establish exchangeability."""
    rng = random.Random(seed)
    n, t = len(cal), len(test)
    width = [0] * bins
    for r in range(n + 1):
        width[r * bins // (n + 1)] += 1
    expected = [t * w / (n + 1) for w in width]

    def stat(c, s):
        cnt = _rank_counts(s, sorted(c), bins, rng)
        return math.fsum((o - e) ** 2 / e for o, e in zip(cnt, expected) if e > 0), cnt

    observed, counts = stat(cal, test)
    pooled = list(cal) + list(test)
    ge = 0
    for _ in range(n_perm):
        rng.shuffle(pooled)
        s, _ = stat(pooled[:n], pooled[n:])
        if s >= observed:
            ge += 1
    return {"counts": counts, "expected": [round(e, 3) for e in expected],
            "statistic": observed, "p_value": (1 + ge) / (1 + n_perm)}


def ks_uniform(us):
    """Kolmogorov-Smirnov distance from Uniform(0, 1) and its asymptotic p-value."""
    xs = sorted(us)
    n = len(xs)
    d = max(max((i + 1) / n - x, x - i / n) for i, x in enumerate(xs))
    lam = (math.sqrt(n) + 0.12 + 0.11 / math.sqrt(n)) * d
    p = 2 * sum((-1) ** (j - 1) * math.exp(-2 * j * j * lam * lam) for j in range(1, 101))
    return d, max(0.0, min(1.0, p))


# ------------------------------------------------------------------ worlds with exact ground truth


class ClaimWorld:
    """Outputs of K claims. Each claim is incorrect with probability p, independently, or
    exactly one claim per output is incorrect (one_false=True, which makes the nonconformity
    score continuous). Incorrect-claim scores ~ Beta(1, cF); correct-claim scores ~ Beta(cT, 1).
    Groups mix with weights. Coverage and retention are exact for any threshold."""

    def __init__(self, groups, one_false=False):
        self.groups = [dict(g) for g in groups]
        self.one_false = one_false
        total = float(sum(g["weight"] for g in groups))
        acc, self.cum = 0.0, []
        for g in self.groups:
            acc += g["weight"] / total
            self.cum.append(acc)

    def _pick(self, rng, name=None):
        if name is not None:
            return next(g for g in self.groups if g["name"] == name)
        u = rng.random()
        for g, c in zip(self.groups, self.cum):
            if u <= c:
                return g
        return self.groups[-1]

    def sample(self, rng, group=None):
        g = self._pick(rng, group)
        k = g["K"]
        if self.one_false:
            j = rng.randrange(k)
            oks = [i != j for i in range(k)]
        else:
            oks = [rng.random() >= g["p"] for _ in range(k)]
        claims = []
        for ok in oks:
            u = rng.random()
            s = u ** (1.0 / g["cT"]) if ok else 1.0 - (1.0 - u) ** (1.0 / g["cF"])
            claims.append((s, ok))
        return g["name"], claims

    @staticmethod
    def _gf(t, c):
        if t <= 0:
            return 0.0
        if t >= 1:
            return 1.0
        return 1.0 - (1.0 - t) ** c

    @staticmethod
    def _gt(t, c):
        if t <= 0:
            return 0.0
        if t >= 1:
            return 1.0
        return t ** c

    def coverage(self, tau, group=None):
        """Exact P(all kept claims correct) when claims scoring above tau are kept."""
        num = den = 0.0
        for g in self.groups:
            if group is not None and g["name"] != group:
                continue
            if tau == INF:
                c = 1.0
            elif self.one_false:
                c = self._gf(tau, g["cF"])
            else:
                c = (1.0 - g["p"] * (1.0 - self._gf(tau, g["cF"]))) ** g["K"]
            num += g["weight"] * c
            den += g["weight"]
        return num / den

    def retention(self, tau, group=None):
        """Exact expected share of claims kept (claim-weighted)."""
        num = den = 0.0
        for g in self.groups:
            if group is not None and g["name"] != group:
                continue
            if tau == INF:
                kt = kf = 0.0
            else:
                kt, kf = 1.0 - self._gt(tau, g["cT"]), 1.0 - self._gf(tau, g["cF"])
            q = 1.0 / g["K"] if self.one_false else g["p"]
            num += g["weight"] * g["K"] * ((1 - q) * kt + q * kf)
            den += g["weight"] * g["K"]
        return num / den


class AnswerWorld:
    """Confidence c ~ Uniform(0, 1); P(correct | c) = c ** gamma. Exact answer rate, joint error
    P(answered and wrong) and selective error P(wrong | answered) for 'answer when c > lam'."""

    def __init__(self, gamma=1.0):
        self.gamma = gamma

    def sample(self, rng):
        c = rng.random()
        return c, rng.random() < c ** self.gamma

    def answer_rate(self, lam):
        if lam >= 1:
            return 0.0
        return 1.0 - max(lam, 0.0)

    def joint_error(self, lam):
        if lam >= 1:
            return 0.0
        x, g = max(lam, 0.0), self.gamma
        return (1.0 - x) - (1.0 - x ** (g + 1)) / (g + 1)

    def selective_error(self, lam):
        a = self.answer_rate(lam)
        return self.joint_error(lam) / a if a > 0 else 0.0


class ResponseWorld:
    """m responses per prompt. Prompt difficulty d ~ Uniform(0, 1); each response is correct with
    probability 1 - d. A fixed, imperfect scorer gives f = exp(beta * (+1 or -1) + sigma * Z)."""

    def __init__(self, m=4, beta=1.0, sigma=1.0):
        self.m, self.beta, self.sigma = m, beta, sigma

    def sample(self, rng):
        d = rng.random()
        out = []
        for _ in range(self.m):
            ok = rng.random() >= d
            out.append((math.exp(self.beta * (1 if ok else -1) + self.sigma * rng.gauss(0.0, 1.0)), ok))
        return out


# ------------------------------------------------------------------ reporting helpers


def _mean(xs):
    return math.fsum(xs) / len(xs) if xs else float("nan")


def _sd(xs):
    if len(xs) < 2:
        return 0.0
    m = _mean(xs)
    return math.sqrt(math.fsum((x - m) ** 2 for x in xs) / (len(xs) - 1))


def _se(xs):
    return _sd(xs) / math.sqrt(len(xs)) if xs else float("nan")


def _md_table(headers, rows):
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(x) for x in r) + " |")
    return "\n".join(out)


def _f(x, d=4):
    if isinstance(x, float):
        if x == INF:
            return "inf"
        if x == -INF:
            return "-inf"
        return ("%." + str(d) + "f") % x
    return str(x)


def _new_dir(path):
    try:
        os.makedirs(path, exist_ok=False)
    except FileExistsError:
        raise SystemExit("use a new output directory: %s" % path) from None


# ------------------------------------------------------------------ bench


def run_bench(out_dir, quick=False, seed=20261117):
    """Run the synthetic experiments for this edition with a fixed seed."""
    t0 = time.time()
    rng = random.Random(seed)
    big = 600 if quick else 4000
    res = {"version": __version__, "seed": seed, "quick": quick}
    md = ["# Guarantee Bench %s results (seed %d%s)" % (__version__, seed, ", quick" if quick else "")]

    # B1 conformal factuality in a generated claim world
    world = ClaimWorld([dict(name="all", weight=1, K=6, p=0.15, cF=3.0, cT=3.0)])
    alphas = (0.05, 0.1, 0.2)
    rows, out = [], []
    for n in (20, 50, 200):
        cov = {a: [] for a in alphas}
        ret = {a: [] for a in alphas}
        none = {a: 0 for a in alphas}
        for _ in range(big):
            rs = [max_false_score(world.sample(rng)[1]) for _ in range(n)]
            for a in alphas:
                tau = conformal_threshold(rs, a)
                cov[a].append(world.coverage(tau))
                ret[a].append(world.retention(tau))
                none[a] += tau == INF
        for a in alphas:
            k = conformal_rank(n, a)
            d = dict(n=n, alpha=a, k=k, k_over_n1=k / (n + 1), mean_coverage=_mean(cov[a]),
                     se=_se(cov[a]), mean_retention=_mean(ret[a]), keep_nothing_share=none[a] / big,
                     base_factuality=world.coverage(-INF))
            out.append(d)
            rows.append([n, a, k, _f(d["k_over_n1"]), _f(d["mean_coverage"]), _f(d["se"], 5),
                         _f(d["mean_retention"]), _f(d["keep_nothing_share"], 3)])
    res["B1_factuality"] = out
    md += ["", "## B1 Conformal factuality, K=6 claims, p=0.15, cF=cT=3 (base factuality %.4f), %d draws"
           % (world.coverage(-INF), big),
           _md_table(["n", "alpha", "k", "k/(n+1)", "mean coverage", "s.e.", "mean retention",
                      "keep-nothing share"], rows)]

    # B2 continuous world: exact value k/(n+1) and the Beta law of conditional coverage
    cw = ClaimWorld([dict(name="c", weight=1, K=6, p=0.0, cF=3.0, cT=3.0)], one_false=True)
    rows, out = [], []
    for n in (20, 100):
        a = 0.1
        k = conformal_rank(n, a)
        covs, pits = [], []
        for _ in range(big):
            tau = conformal_threshold([max_false_score(cw.sample(rng)[1]) for _ in range(n)], a)
            c = cw.coverage(tau)
            covs.append(c)
            pits.append(beta_cdf_int(c, k, n + 1 - k))
        dks, pks = ks_uniform(pits)
        d = dict(n=n, alpha=a, k=k, exact=k / (n + 1), upper=1 - a + 1 / (n + 1), mean=_mean(covs),
                 se=_se(covs), sd=_sd(covs), beta_sd=math.sqrt(k * (n + 1 - k) / ((n + 1) ** 2 * (n + 2))),
                 ks_d=dks, ks_p=pks)
        out.append(d)
        rows.append([n, k, _f(d["exact"]), _f(d["mean"]), _f(d["se"], 5), _f(d["upper"]), _f(d["sd"]),
                     _f(d["beta_sd"]), _f(dks), _f(pks, 3)])
    res["B2_beta_law"] = out
    md += ["", "## B2 Continuous world (one incorrect claim per output), alpha = 0.1, %d draws" % big,
           _md_table(["n", "k", "k/(n+1)", "mean coverage", "s.e.", "upper bound", "sd", "Beta sd",
                      "KS D", "KS p"], rows)]

    # B3 subgroup gap: marginal versus group-balanced calibration
    gw = ClaimWorld([dict(name="easy", weight=0.7, K=6, p=0.05, cF=3.0, cT=3.0),
                     dict(name="hard", weight=0.3, K=6, p=0.40, cF=3.0, cT=3.0)])
    n, a = 300, 0.1
    acc = {k_: [] for k_ in ("m_easy_cov", "m_hard_cov", "m_easy_ret", "m_hard_ret",
                             "g_easy_cov", "g_hard_cov", "g_easy_ret", "g_hard_ret", "m_all_cov")}
    for _ in range(big // 2):
        samp = [gw.sample(rng) for _ in range(n)]
        tau = calibrate_factuality([c for _, c in samp], a)
        taus = calibrate_by_group([c for _, c in samp], [g for g, _ in samp], a)
        acc["m_all_cov"].append(gw.coverage(tau))
        for g in ("easy", "hard"):
            acc["m_%s_cov" % g].append(gw.coverage(tau, g))
            acc["m_%s_ret" % g].append(gw.retention(tau, g))
            tg = group_threshold(taus, g)
            acc["g_%s_cov" % g].append(gw.coverage(tg, g))
            acc["g_%s_ret" % g].append(gw.retention(tg, g))
    res["B3_subgroups"] = {k_: _mean(v) for k_, v in acc.items()}
    res["B3_subgroups"]["base_easy"] = gw.coverage(-INF, "easy")
    res["B3_subgroups"]["base_hard"] = gw.coverage(-INF, "hard")
    s = res["B3_subgroups"]
    md += ["", "## B3 Subgroups (easy 70%%: p=0.05; hard 30%%: p=0.40; K=6), n=300, alpha=0.1, %d draws"
           % (big // 2),
           _md_table(["calibration", "overall coverage", "easy coverage", "hard coverage",
                      "easy retention", "hard retention"],
                     [["marginal", _f(s["m_all_cov"]), _f(s["m_easy_cov"]), _f(s["m_hard_cov"]),
                       _f(s["m_easy_ret"]), _f(s["m_hard_ret"])],
                      ["group-balanced", "-", _f(s["g_easy_cov"]), _f(s["g_hard_cov"]),
                       _f(s["g_easy_ret"]), _f(s["g_hard_ret"])],
                      ["no filter (base)", "-", _f(s["base_easy"]), _f(s["base_hard"]), "1", "1"]])]

    # B4 shift, and the one-batch check
    wa = ClaimWorld([dict(name="a", weight=1, K=6, p=0.10, cF=3.0, cT=3.0)])
    wb = ClaimWorld([dict(name="b", weight=1, K=6, p=0.25, cF=3.0, cT=3.0)])
    n, a, t_batch = 200, 0.1, 200
    reps = 100 if quick else 400
    cov_a, cov_b, rej_null, rej_shift = [], [], 0, 0
    for r_ in range(reps):
        cal = [max_false_score(wa.sample(rng)[1]) for _ in range(n)]
        tau = conformal_threshold(cal, a)
        cov_a.append(wa.coverage(tau))
        cov_b.append(wb.coverage(tau))
        same = [max_false_score(wa.sample(rng)[1]) for _ in range(t_batch)]
        moved = [max_false_score(wb.sample(rng)[1]) for _ in range(t_batch)]
        rej_null += exchangeability_check(cal, same, n_perm=199, seed=r_)["p_value"] <= 0.05
        rej_shift += exchangeability_check(cal, moved, n_perm=199, seed=r_)["p_value"] <= 0.05
    res["B4_shift"] = dict(n=n, alpha=a, reps=reps, coverage_same=_mean(cov_a), coverage_shifted=_mean(cov_b),
                           base_a=wa.coverage(-INF), base_b=wb.coverage(-INF),
                           check_rejects_null=rej_null / reps, check_rejects_shift=rej_shift / reps)
    s = res["B4_shift"]
    md += ["", "## B4 Shift: calibrate on p=0.10, test on p=0.25 (K=6), n=200, alpha=0.1, %d repetitions" % reps,
           _md_table(["quantity", "value"],
                     [["coverage, same world", _f(s["coverage_same"])],
                      ["coverage, shifted world", _f(s["coverage_shifted"])],
                      ["base factuality, same / shifted", "%s / %s" % (_f(s["base_a"]), _f(s["base_b"]))],
                      ["one-batch check rejects at 0.05, same world (batch of 200)", _f(s["check_rejects_null"], 3)],
                      ["one-batch check rejects at 0.05, shifted world (batch of 200)", _f(s["check_rejects_shift"], 3)]])]

    # B5 abstention: conformal risk control versus the selective-risk certificate
    rows, out = [], []
    grid = [i / 200 for i in range(188, -1, -1)]   # 0.94 down to 0, fixed before any data
    n, a, delta = 2000, 0.1, 0.1
    for gamma in (1.0, 3.0):
        aw = AnswerWorld(gamma)
        jr, sr, ar, viol, lj, ls_ = [], [], [], 0, [], []
        for _ in range(big // 2):
            cal = [aw.sample(rng) for _ in range(n)]
            lam = crc_threshold(cal, a)
            jr.append(aw.joint_error(lam))
            sr.append(aw.selective_error(lam))
            ar.append(aw.answer_rate(lam))
            lam2 = ltt_selective_threshold(cal, a, delta, grid)
            lj.append(aw.answer_rate(lam2))
            ls_.append(aw.selective_error(lam2))
            viol += aw.selective_error(lam2) > a
        d = dict(gamma=gamma, crc_joint=_mean(jr), crc_selective=_mean(sr), crc_answer=_mean(ar),
                 crc_lower=a - 2 / (n + 1), ltt_answer=_mean(lj), ltt_selective=_mean(ls_),
                 ltt_violation=viol / (big // 2), base_error=aw.joint_error(-INF))
        out.append(d)
        rows.append([gamma, _f(d["base_error"]), _f(d["crc_joint"]), _f(d["crc_selective"]), _f(d["crc_answer"]),
                     _f(d["ltt_selective"]), _f(d["ltt_answer"]), _f(d["ltt_violation"], 3)])
    res["B5_abstention"] = out
    md += ["", "## B5 Abstention, n=2000, alpha=0.1, delta=0.1, LTT grid 0.94 to 0 by 0.005, %d draws "
               "(P(correct | c) = c^gamma)" % (big // 2),
           _md_table(["gamma", "error, no abstention", "CRC joint error", "CRC selective error", "CRC answer rate",
                      "LTT selective error", "LTT answer rate", "LTT share above alpha"], rows)]

    # B6 e-scores versus p-scores under post-hoc tolerance
    rw = ResponseWorld(m=4, beta=1.0, sigma=1.0)
    n = 200
    rules = [("fixed", {}), ("best", {}), ("fraction", {"fraction": 0.5})]
    terms = {(kind, r): [] for kind in ("e", "p") for r, _ in rules}
    kept_ok = {kind: [] for kind in ("e", "p")}
    for _ in range(big):
        cal = [worst_incorrect(rw.sample(rng)) for _ in range(n)]
        test = rw.sample(rng)
        fs = [f for f, _ in test]
        oks = [ok for _, ok in test]
        for kind in ("e", "p"):
            sc = e_scores(fs, cal) if kind == "e" else p_scores(fs, cal)
            for r, kw in rules:
                terms[(kind, r)].append(distortion_term(sc, oks, posthoc_level(sc, r, **kw)))
            n_ok = sum(oks)
            if n_ok:
                kept_ok[kind].append(sum(1 for s_, ok in zip(sc, oks) if ok and s_ <= 0.1) / n_ok)
    out, rows = [], []
    for r, _ in rules:
        d = dict(rule=r, e_sd=_mean(terms[("e", r)]), e_se=_se(terms[("e", r)]),
                 p_sd=_mean(terms[("p", r)]), p_se=_se(terms[("p", r)]))
        out.append(d)
        rows.append([r, _f(d["e_sd"]), _f(d["e_se"]), _f(d["p_sd"]), _f(d["p_se"])])
    res["B6_escores"] = dict(rows=out, correct_kept_at_0_1_e=_mean(kept_ok["e"]),
                             correct_kept_at_0_1_p=_mean(kept_ok["p"]), n=n)
    hn = math.fsum(1.0 / j for j in range(1, n + 2))
    res["B6_escores"]["harmonic_H_n1"] = hn
    md += ["", "## B6 Size distortion E[1{error} / level], n=200 calibration prompts, 4 responses each, %d draws" % big,
           _md_table(["tolerance rule", "e-scores", "s.e.", "p-scores", "s.e."], rows),
           "", "Share of correct responses kept at a fixed level 0.1: e-scores %s, p-scores %s. H_(n+1) = %s."
           % (_f(res["B6_escores"]["correct_kept_at_0_1_e"]), _f(res["B6_escores"]["correct_kept_at_0_1_p"]), _f(hn))]

    # B7 the floating-point ceiling
    mism, total, examples = 0, 0, []
    for n_ in range(0, 2001):
        for a100 in range(1, 51):
            total += 1
            flt = math.ceil((n_ + 1) * (1 - a100 / 100))
            if flt != conformal_rank(n_, Fraction(a100, 100)):
                mism += 1
                if len(examples) < 8:
                    examples.append((n_, a100 / 100, flt, conformal_rank(n_, Fraction(a100, 100))))
    res["B7_float"] = dict(pairs=total, mismatches=mism, examples=examples)
    md += ["", "## B7 Floating-point ceiling: (n, alpha) pairs with n = 0..2000, alpha = 0.01..0.50",
           "%d of %d pairs give a different k in float arithmetic. Examples (n, alpha, float k, exact k): %s"
           % (mism, total, examples)]

    res["seconds"] = round(time.time() - t0, 1)
    md += ["", "Run time: %.1f s." % res["seconds"]]
    _new_dir(out_dir)
    with open(os.path.join(out_dir, "RESULTS.json"), "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1)
    with open(os.path.join(out_dir, "RESULTS.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(md) + "\n")
    return res


# ------------------------------------------------------------------ filter your own outputs


def _load_jsonl(path):
    rows = []
    with open(path, encoding="utf-8") as fh:
        for ln, line in enumerate(fh, 1):
            line = line.strip()
            if not line:
                continue
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError as e:
                raise SystemExit("%s:%d is not valid JSON (%s)" % (path, ln, e.msg))
    return rows


def _claims_of(row, need_labels):
    out = []
    for c in row.get("claims", []):
        s = float(c["score"])
        if not math.isfinite(s):
            raise SystemExit("a claim score must be finite in row %r" % row.get("id"))
        lab = c.get("correct")
        if need_labels and not isinstance(lab, bool):
            raise SystemExit("calibration row %r has a claim without a true/false 'correct' field" % row.get("id"))
        out.append((s, lab))
    return out


def run_filter(cal_path, test_path, alpha, out_path, by_group=False):
    """Calibrate on labeled outputs and filter new ones. Writes one JSON line per test output."""
    if os.path.exists(out_path):
        raise SystemExit("%s exists; choose a new file name" % out_path)
    cal_rows, test_rows = _load_jsonl(cal_path), _load_jsonl(test_path)
    cal = [_claims_of(r, True) for r in cal_rows]
    for row in test_rows:
        _claims_of(row, False)
    if by_group:
        taus = calibrate_by_group(cal, [r.get("group") for r in cal_rows], alpha)
    else:
        tau0 = calibrate_factuality(cal, alpha)
    n = len(cal)
    k = conformal_rank(n, alpha)
    realized = []
    with open(out_path, "x", encoding="utf-8") as fh:
        for row in test_rows:
            tau = group_threshold(taus, row.get("group")) if by_group else tau0
            kept_idx = [i for i, c in enumerate(row.get("claims", [])) if float(c["score"]) > tau]
            kept = [row["claims"][i] for i in kept_idx]
            labs = [c.get("correct") for c in kept]
            ok = all(labs) if all(isinstance(x, bool) for x in labs) else None
            if ok is not None:
                realized.append(ok)
            fh.write(json.dumps({"id": row.get("id"), "threshold": tau if math.isfinite(tau) else ("-inf" if tau < 0 else None),
                                 "kept": kept, "removed": len(row.get("claims", [])) - len(kept),
                                 "all_kept_correct": ok}) + "\n")
    summary = dict(n_calibration=n, alpha=alpha, k=k, by_group=by_group,
                   threshold=None if by_group else (tau0 if math.isfinite(tau0) else ("-inf" if tau0 < 0 else None)),
                   test_outputs=len(test_rows), realized_coverage=_mean(realized) if realized else None)
    return summary


GUARANTEE_TEXT = (
    "Guarantee (Mohri & Hashimoto 2024, Thm 4.1, in this release's convention): if the calibration "
    "outputs and a new output are exchangeable, P(every kept claim of the new output is correct) >= 1 - alpha, "
    "averaged over calibration sets and new outputs. It is not a statement about any one output, any one "
    "subgroup or any one calibration set, and it does not survive a change in the kind of prompts.")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="guarantee_bench", description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd")
    b = sub.add_parser("bench", help="run the simulation results")
    b.add_argument("--out", required=True)
    b.add_argument("--quick", action="store_true")
    b.add_argument("--seed", type=int, default=20261117)
    f = sub.add_parser("filter", help="calibrate on labeled outputs and filter new ones")
    f.add_argument("--cal", required=True)
    f.add_argument("--test", required=True)
    f.add_argument("--alpha", type=float, required=True)
    f.add_argument("--out", required=True)
    f.add_argument("--by-group", action="store_true")
    a = ap.parse_args(argv)
    if a.cmd == "bench":
        res = run_bench(a.out, quick=a.quick, seed=a.seed)
        print("bench done in %.1f s; see %s" % (res["seconds"], os.path.join(a.out, "RESULTS.md")))
    elif a.cmd == "filter":
        s = run_filter(a.cal, a.test, a.alpha, a.out, by_group=a.by_group)
        print(json.dumps(s, indent=1))
        print(GUARANTEE_TEXT)
    else:
        ap.print_help()
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
