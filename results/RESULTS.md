# Guarantee Bench 0.1.0 results (seed 20261117)

## B1 Conformal factuality, K=6 claims, p=0.15, cF=cT=3 (base factuality 0.3771), 4000 draws
| n | alpha | k | k/(n+1) | mean coverage | s.e. | mean retention | keep-nothing share |
|---|---|---|---|---|---|---|---|
| 20 | 0.05 | 20 | 0.9524 | 0.9531 | 0.00072 | 0.5846 | 0.000 |
| 20 | 0.1 | 19 | 0.9048 | 0.9052 | 0.00099 | 0.7128 | 0.000 |
| 20 | 0.2 | 17 | 0.8095 | 0.8104 | 0.00132 | 0.8204 | 0.000 |
| 50 | 0.05 | 49 | 0.9608 | 0.9607 | 0.00043 | 0.5952 | 0.000 |
| 50 | 0.1 | 46 | 0.9020 | 0.9009 | 0.00067 | 0.7386 | 0.000 |
| 50 | 0.2 | 41 | 0.8039 | 0.8010 | 0.00089 | 0.8352 | 0.000 |
| 200 | 0.05 | 191 | 0.9502 | 0.9503 | 0.00024 | 0.6528 | 0.000 |
| 200 | 0.1 | 181 | 0.9005 | 0.9003 | 0.00033 | 0.7496 | 0.000 |
| 200 | 0.2 | 161 | 0.8010 | 0.8011 | 0.00044 | 0.8394 | 0.000 |

## B2 Continuous world (one incorrect claim per output), alpha = 0.1, 4000 draws
| n | k | k/(n+1) | mean coverage | s.e. | upper bound | sd | Beta sd | KS D | KS p |
|---|---|---|---|---|---|---|---|---|---|
| 20 | 19 | 0.9048 | 0.9046 | 0.00097 | 0.9476 | 0.0614 | 0.0626 | 0.0097 | 0.843 |
| 100 | 91 | 0.9010 | 0.9010 | 0.00047 | 0.9099 | 0.0298 | 0.0296 | 0.0126 | 0.551 |

## B3 Subgroups (easy 70%: p=0.05; hard 30%: p=0.40; K=6), n=300, alpha=0.1, 2000 draws
| calibration | overall coverage | easy coverage | hard coverage | easy retention | hard retention |
|---|---|---|---|---|---|
| marginal | 0.8998 | 0.9649 | 0.7478 | 0.8282 | 0.5668 |
| group-balanced | - | 0.9020 | 0.9049 | 0.9380 | 0.4435 |
| no filter (base) | - | 0.7351 | 0.0467 | 1 | 1 |

## B4 Shift: calibrate on p=0.10, test on p=0.25 (K=6), n=200, alpha=0.1, 400 repetitions
| quantity | value |
|---|---|
| coverage, same world | 0.8991 |
| coverage, shifted world | 0.7645 |
| base factuality, same / shifted | 0.5314 / 0.1780 |
| one-batch check rejects at 0.05, same world (batch of 200) | 0.048 |
| one-batch check rejects at 0.05, shifted world (batch of 200) | 1.000 |

## B5 Abstention, n=2000, alpha=0.1, delta=0.1, LTT grid 0.94 to 0 by 0.005, 2000 draws (P(correct | c) = c^gamma)
| gamma | error, no abstention | CRC joint error | CRC selective error | CRC answer rate | LTT selective error | LTT answer rate | LTT share above alpha |
|---|---|---|---|---|---|---|---|
| 1.0 | 0.5000 | 0.0998 | 0.2233 | 0.4466 | 0.0708 | 0.1415 | 0.077 |
| 3.0 | 0.7500 | 0.1000 | 0.3513 | 0.2843 | 0.0161 | 0.0112 | 0.065 |

## B6 Size distortion E[1{error} / level], n=200 calibration prompts, 4 responses each, 4000 draws
| tolerance rule | e-scores | s.e. | p-scores | s.e. |
|---|---|---|---|---|
| fixed | 0.0275 | 0.0083 | 1.0225 | 0.0479 |
| best | 0.4810 | 0.0201 | 4.0134 | 0.3039 |
| fraction | 0.5260 | 0.0145 | 2.7460 | 0.1863 |

Share of correct responses kept at a fixed level 0.1: e-scores 0.1376, p-scores 0.6497. H_(n+1) = 5.8830.

## B7 Floating-point ceiling: (n, alpha) pairs with n = 0..2000, alpha = 0.01..0.50
222 of 100050 pairs give a different k in float arithmetic. Examples (n, alpha, float k, exact k): [(24, 0.44, 15, 14), (49, 0.42, 30, 29), (49, 0.44, 29, 28), (74, 0.44, 43, 42), (99, 0.41, 60, 59), (99, 0.42, 59, 58), (99, 0.43, 58, 57), (99, 0.44, 57, 56)]

Run time: 121.8 s.
