"""
This is not an implementation of SQISign

The purpose of this file is simply a timing of the
costs of computing and evaluating isogenies using the
available torsion from the SQISign parameters
"""

# Local imports
from kummer_line import KummerLine
from utilities import compute_quadratic_twist, fix_even_torsion
from benchmark_utils import compare_isogeny

proof.all(False)

PARAM = "p3923"
for arg in sys.argv[1:]:
    if arg.lower() in ["--p6983"]:
        PARAM = "p6983"
    if arg.lower() in ["--p3923"]:
        PARAM = "p3923"

if PARAM == "p3923":
    # p6983
    # p - 1 = 2 * 3^53 * 43 * 103^2 * 109 * 199 * 227 * 419 * 491 * 569 * 631 * 677 * 857 * 859 * 883 * 1019 * 1171 * 1879 * 2713 * 4283
    # p + 1 = 2^33 * 5^21 * 7^2 * 11 * 31 * 83 * 107 * 137 * 751 * 827 * 3691 * 4019 * 6983 * 517434778561 * 26602537156291
    p = 73743043621499797449074820543863456997944695372324032511999999999999999999999

    A_torsion = 2^33 * 5^21 * 7^2 * 11 * 31 * 83 * 107 * 137 * 751 * 827 * 3691 * 4019 * 6983
    A_cofactor = (p+1) // A_torsion
    B_torsion = 2 * 3^53 * 43 * 103^2 * 109 * 199 * 227 * 419 * 491 * 569 * 631 * 677 * 857 * 859 * 883 * 1019 * 1171 * 1879 * 2713 * 4283
    B_cofactor = (p-1) // B_torsion

else:
    # p - 1 = 2 * 3^65 * 13 * 17 * 43 * 79 * 157 * 239 * 271 * 283 * 307 * 563 * 599 * 607 * 619 * 743 * 827 * 941 * 2357 * 10069
    # p + 1 = 2^65 * 5^2 * 7 * 11 * 19 * 29^2 * 37^2 * 47 * 197 * 263 * 281 * 461 * 521 * 3923 * 62731 * 96362257 * 3924006112952623
    p = 23759399264157352358673788613307970528646815114090876784643387662192449945599

    A_torsion = 2^65 * 5^2 * 7 * 11 * 19 * 29^2 * 37^2 * 47 * 197 * 263 * 281 * 461 * 521 * 3923
    A_cofactor = (p+1) // A_torsion
    B_torsion = 2 * 3^65 * 13 * 17 * 43 * 79 * 157 * 239 * 271 * 283 * 307 * 563 * 599 * 607 * 619 * 743 * 827 * 941 * 2357
    B_cofactor = (p-1) // B_torsion

# Compute base field and elliptic curve
Fp2.<i> = GF(p^2)
E0 = EllipticCurve(Fp2, [1, 0])

# Compute the quadratic twist of E0 for torsion
# generation
D_twist, E0t = compute_quadratic_twist(E0)

# Compute points on E0 and E0t with order 
# A_torsion and B_torsion
PA, QA = [A_cofactor*X for X in E0.gens()]
PB, QB = [B_cofactor*X for X in E0t.gens()]

# Ensure that 2^(b-1) * Q = (0,0) to 
# stop the kernel ever being the point (0,0)
PA, QA = fix_even_torsion(PA, QA)
PB, QB = fix_even_torsion(PB, QB, twist=True)

# Represent x-coordinate on the Kummer line of
# E0
L = KummerLine(E0)
xPA, xQA, xPQA = [L(X) for X in (PA, QA, PA-QA)]

# Points in the B torsion are on E0_twist, so we map
# the x-coords back to E0 
xPB, xQB, xPQB = [L(D_twist * X[0]) for X in (PB, QB, PB-QB)]

# Compare isogeny timings for p+1 torsion and p-1 torsion
compare_isogeny(PA, QA, xPA, xQA, A_torsion)
compare_isogeny(PB, QB, xPB, xQB, B_torsion)
