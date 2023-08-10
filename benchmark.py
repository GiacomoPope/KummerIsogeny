import time
import cProfile
import pstats

from sage.all import *

from kummer_line import KummerLine
from kummer_isogeny import KummerLineIsogeny_Velu, KummerLineIsogeny_VeluSqrt
from utilities import print_info

proof.all(False)


def print_comparison():
    for ell in [311, 571, 1321, 6011, 76667]:
        if ell < 10:
            continue
        cofactor = (p + 1) // ell
        R = cofactor * P
        ker = K(R)
        img = K(Q)

        print_info(f"ell = {ell}")

        t0 = time.process_time()
        psi = KummerLineIsogeny_Velu(K, ker, ell)
        print(f"Velu isogeny took: {time.process_time() - t0}")

        t0 = time.process_time()
        _ = psi(img)
        print(f"Velu image took: {time.process_time() - t0}")
        print()

        t0 = time.process_time()
        phi = KummerLineIsogeny_VeluSqrt(K, ker, ell)
        print(f"Sqrt isogeny took: {time.process_time() - t0}")

        t0 = time.process_time()
        _ = phi(img)
        print(f"Sqrt image took: {time.process_time() - t0}")
        print()

        j1 = psi.codomain().j_invariant()
        j2 = phi.codomain().j_invariant()

        assert j1 == j2


def profile_codomain():
    p = K.base_ring().characteristic()

    cofactor = (p + 1) // ell
    R = cofactor * P
    ker = K(R)

    p_codomain = cProfile.Profile()
    p_codomain.enable()

    _ = KummerLineIsogeny_VeluSqrt(K, ker, ell)

    p_codomain.disable()
    p_codomain.dump_stats("p_codomain.cProfile")
    p = pstats.Stats("p_codomain.cProfile")
    p.strip_dirs().sort_stats("cumtime").print_stats(30)


def profile_image(ell, K, P, Q):
    p = K.base_ring().characteristic()

    cofactor = (p + 1) // ell
    R = cofactor * P
    ker = K(R)
    img = K(Q)

    phi = KummerLineIsogeny_VeluSqrt(K, ker, ell)

    p_image = cProfile.Profile()
    p_image.enable()

    phi(img)

    p_image.disable()
    p_image.dump_stats("p_image.cProfile")
    p = pstats.Stats("p_image.cProfile")
    p.strip_dirs().sort_stats("cumtime").print_stats(30)


if __name__ == "__main__":
    # BSIDH Prime
    p = 0x1935BECE108DC6C0AAD0712181BB1A414E6A8AAA6B510FC29826190FE7EDA80F

    # Compute curve
    F = GF(p**2, name="i", modulus=[1, 0, 1])
    i = F.gens()[0]
    E = EllipticCurve(F, [0, 6, 0, 1, 0])
    E.set_order((p + 1) ** 2)
    K = KummerLine(E)

    # Compute A-torsion
    A_torsion = ZZ(
        2**4
        * 3
        * 7**16
        * 17**9
        * 31**8
        * 311
        * 571
        * 1321
        * 5119
        * 6011
        * 14207
        * 28477
        * 76667
    )
    A_cofactor = (p + 1) // A_torsion
    P, Q = E.gens()
    P, Q = A_cofactor * P, A_cofactor * Q

    # print_comparison()
    ell = ZZ(28477)
    profile_codomain()
