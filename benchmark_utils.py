import time

# local imports
from kummer_isogeny import (
    KummerLineIsogeny_Velu,
    KummerLineIsogeny_VeluSqrt,
    KummerLineIsogeny,
)

from utilities import EllipticCurveIsogenyFactored, print_info


def compare_isogeny(P, Q, xP, xQ, order):
    """
    A function which compares the times of isogeny computation
    and evaluation using:

    - Sage Naive isogeny E.isogeny(K, algotithm="factored)
    - Optimised Sage Isogeny written for
      https://github.com/LearningToSQI/SQISign-SageMath
    - New KummerLine x-only isogenies using KummerLineIsogeny
    """
    E = P.curve()
    L = xP.parent()

    print_info(f"Computing isogenies of degree:\n{order.factor()}")

    # Time naive SageMath isogeny evaluation
    # t0 = time.time()
    # psi = E.isogeny(P, algorithm="factored")
    # print(f"Naive SageMath codomain computation: {time.time() - t0:.5f}")
    # t0 = time.time()
    # psi(Q)
    # print(f"Naive SageMath point evaluation: {time.time() - t0:.5f}\n")

    # Time optimisation with sparse strategy and velusqrt
    t0 = time.time()
    sigma = EllipticCurveIsogenyFactored(E, P, order=order)
    print(f"Optimised SageMath codomain computation: {time.time() - t0:.5f}")
    t0 = time.time()
    sigma(Q)
    print(f"Optimised SageMath point evaluation: {time.time() - t0:.5f}\n")

    # Time x-only formula using KummerLine and KummerIsogeny classes
    t0 = time.time()
    phi = KummerLineIsogeny(L, xP, order)
    print(f"KummerLine codomain computation: {time.time() - t0:.5f}")
    t0 = time.time()
    phi(xQ)
    print(f"KummerLine point evaluation: {time.time() - t0:.5f}\n")


def compare_isogeny_factors(P, Q, xP, xQ, order):
    """
    A function which compares the times of prime-degree isogeny
    computation and evaluation using:

    - Sage velu isogeny E.isogeny(K)
    - Sage velusqrt isogeny E.isogeny(K, algorithm="velusqrt")
    - KummerLine x-only velu using KummerLineIsogeny_Velu
    - KummerLine x-only velusqrt using KummerLineIsogeny_VeluSqrt
    """
    E = P.curve()
    L = xP.parent()

    for l, _ in factor(order):
        # Compute cofactor
        k = order // l

        # Kernels of order l
        K = k * P
        xK = k * xP

        print_info(f"{l = }")
        print_info("SageMath timings", banner="-")
        # Codomain computation
        t0 = time.time()
        psi = E.isogeny(K)
        print(f"Velu codomain took: {time.time() - t0:.5f}")

        # Point evaluation
        t0 = time.time()
        psi(Q)
        print(f"Velu evaluation took: {time.time() - t0:.5f}")

        # for large l, do velusqrt
        if l > 10:
            t0 = time.time()
            psi_prime = E.isogeny(K, algorithm="velusqrt")
            print(f"Sqrt codomain took: {time.time() - t0:.5f}")
            t0 = time.time()
            psi_prime(Q)
            print(f"Sqrt evaluation took: {time.time() - t0:.5f}")

        print_info("KummerLine timings", banner="-")

        # Codomain computation
        t0 = time.time()
        phi = KummerLineIsogeny_Velu(L, xK, l)
        print(f"x-only Velu codomain took: {time.time() - t0:.5f}")

        # Point evaluation
        t0 = time.time()
        phi(xQ)
        print(f"Velu evaluation took: {time.time() - t0:.5f}")

        # for large l, do the same but with velusqrt
        if l > 10:
            t0 = time.time()
            phi_prime = KummerLineIsogeny_VeluSqrt(L, xK, l)
            print(f"x-only Sqrt codomain took: {time.time() - t0:.5f}")
            t0 = time.time()
            phi_prime(xQ)
            print(f"x-only Sqrt evaluation took: {time.time() - t0:.5f}")

        print()
