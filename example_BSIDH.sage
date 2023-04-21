"""
Implementation of BSIDH https://ia.cr/2019/1145
to test the Kummer Line isogenies and compare timings
"""

# Python imports
import time 
import sys

# Local imports

from kummer_line import KummerLine
from kummer_isogeny import KummerLineIsogeny
from utilities import compute_quadratic_twist, fix_even_torsion, compare_isogeny

# ============================== #
#     B-SIDH parameter setup     #
# ============================== #

# Set which parameter set to use
# Default is example three
PARAM = "THREE"
for arg in sys.argv[1:]:
    if arg.lower() in ["--two"]:
        PARAM = "TWO"

if PARAM == "TWO":
    # p + 1 = 2^4 * 3 * 7^16 * 17^9 * 31^8 * 311 * 571 * 1321 * 5119 * 6011 * 14207 * 28477 * 76667 * 315668179
    # p - 1 = 2 * 11^18 * 19 * 23^13 * 47 * 79 * 83 * 89 * 151 * 3347 * 17449 * 33461 * 51193 * 258434945441
    p = 0x1935BECE108DC6C0AAD0712181BB1A414E6A8AAA6B510FC29826190FE7EDA80F
    A_torsion = 2^4 * 3 * 7^16 * 17^9 * 31^8 * 311 * 571 * 1321 * 5119 * 6011 * 14207 * 28477 * 76667
    A_cofactor = (p+1) // A_torsion
    B_torsion = 11^18 * 19 * 23^13 * 47 * 79 * 83 * 89 * 151 * 3347 * 17449 * 33461 * 51193
    B_cofactor = (p-1) // B_torsion
else:
    # Example 3 from BSIDH
    # p + 1 = 2^110 * 5 * 7^2 * 67 * 223 * 4229 * 9787 * 13399 * 21521 * 32257 * 47353 * 616228535059
    # p - 1 = 2 * 3^34 * 11 * 17 * 19^2 * 29 * 37 * 53^2 * 97 * 107 * 109 * 131 * 137 * 197 * 199 * 227 * 251 * 5519 * 9091 * 33997 * 38201 * 460409 * 5781011
    p = 0x76042798BBFB78AEBD02490BD2635DEC131ABFFFFFFFFFFFFFFFFFFFFFFFFFFF
    A_torsion = 2^110 * 5 * 7^2 * 67 * 223 * 4229 * 9787 * 13399 * 21521 * 32257 * 47353
    A_cofactor = (p+1) // A_torsion
    B_torsion = 3^34 * 11 * 17 * 19^2 * 29 * 37 * 53^2 * 97 * 107 * 109 * 131 * 137 * 197 * 199 * 227 * 251 * 5519 * 9091 * 33997 * 38201
    B_cofactor = (p-1) // B_torsion

# We need our isogenies to be coprime
# for the protocol to work
assert gcd(A_torsion, B_torsion) == 1

# Construct the base-field and elliptic curve
Fp2.<z2> = GF(p^2)
E0 = EllipticCurve(Fp2, [1,0])
E0.set_order((p + 1)^2)
assert E0.is_supersingular()

# Compute the quadratic twist of E0 for torsion
# generation
D_twist, E0t = compute_quadratic_twist(E0)

# ============================== #
#      B-SIDH parameter gen.     #
# ============================== #

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

# Points in the A torsion are on E0
A_torsion_data = [L(X) for X in (PA, QA, PA-QA)]

# Points in the B torsion are on E0_twist, so we map
# the x-coords back to E0 (y will be irrational but
# Kummer Line doesn't care).
B_torsion_data = [L(D_twist * X[0]) for X in (PB, QB, PB-QB)]

# Check orders of x-only points
assert (A_torsion * A_torsion_data[0]).is_zero()
assert (A_torsion * A_torsion_data[1]).is_zero()
assert (B_torsion * B_torsion_data[0]).is_zero()
assert (B_torsion * B_torsion_data[1]).is_zero()

# ========================================== #
#        BSIDH: https://ia.cr/2019/1145      #
# ========================================== #

def keygen(torsion_data, other_torsion_data, order):
    """
    BSIDH KeyGen following https://ia.cr/2019/1145
    """
    # Unpack torsion
    xP, xQ, xPQ = torsion_data
    xP_other, xQ_other, xPQ_other = other_torsion_data

    # Secret values
    sk = randint(0, order)
    xG = xQ.ladder_3_pt(xP, xPQ, sk)
    phi = KummerLineIsogeny(xG.parent(), xG, order)

    # Public Key
    E = phi.codomain()
    imxP, imxQ, imxPQ = [phi(X) for X in (xP_other, xQ_other, xPQ_other)]
    pk = (E, imxP, imxQ, imxPQ)

    return sk, pk

def shared_secret(sk, pk_other, order):
    """
    BSIDH Shared Secret generation following https://ia.cr/2019/1145
    """
    # Unpack pk
    K, xP, xQ, xPQ = pk_other

    # Compute shared curve
    xK = xQ.ladder_3_pt(xP, xPQ, sk)
    psi = KummerLineIsogeny(K, xK, order)
    return psi.codomain().j_invariant()


# ============================== #
#        Time the Protocol       #
# ============================== #

t0 = time.time()
skA, pkA = keygen(A_torsion_data, B_torsion_data, A_torsion)
print(f"Keygen A took: {time.time() - t0}")

t0 = time.time()
skB, pkB = keygen(B_torsion_data, A_torsion_data, B_torsion)
print(f"Keygen B took: {time.time() - t0}")

t0 = time.time()
sA = shared_secret(skA, pkB, A_torsion)
print(f"Secret A took: {time.time() - t0}")

t0 = time.time()
sB = shared_secret(skB, pkA, B_torsion)
print(f"Secret B took: {time.time() - t0}")

assert sA == sB

# Compare isogeny timings for p+1 torsion and p-1 torsion
xPA, xQA, _ = A_torsion_data
xPB, xQB, _ = B_torsion_data
compare_isogeny(PA, QA, xPA, xQA, A_torsion)
compare_isogeny(PB, QB, xPB, xQB, B_torsion)

