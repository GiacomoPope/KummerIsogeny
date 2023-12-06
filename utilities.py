# Python imports
import textwrap
import sys

# Sage imports
from sage.all import EllipticCurve, EllipticCurveIsogeny, ZZ
from sage.schemes.elliptic_curves.hom_velusqrt import EllipticCurveHom_velusqrt
from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite


def compute_quadratic_twist(E):
    """
    Compute the quadratic twist of E and return
    Et, together with the explicit twisting
    parameter D
    """
    K = E.base_field()
    p = K.characteristic()
    # Compute a good D for the quadratic
    # twist
    while True:
        D_twist = K.random_element()
        if not D_twist.is_square():
            break

    # Make the twisted curve, we need it to
    # be in this form so we can simply map
    # elements from Et to E later on
    _, _, _, a, b = E.a_invariants()
    Et = EllipticCurve(K, [a / D_twist**2, b / D_twist**3])

    # Make sure the twist is actually a twist
    assert Et.is_isomorphic(E.quadratic_twist())
    assert Et.order() == (p - 1) ** 2
    assert Et.is_supersingular()

    return D_twist, Et


def fix_even_torsion(P, Q, twist=False):
    """
    We have to set the torsion basis P,Q such that
    the point (0,0) lies at the bottom of Q to avoid
    making a kernel K which has (0,0) as it's order
    two point
    """
    p = P.curve().base_ring().characteristic()
    if twist:
        oo = (p - 1) // 2
    else:
        oo = (p + 1) // 2

    # Even fix
    Pa, Pb = oo * P, oo * Q

    if Pa[0] == 0:
        return Q, P
    elif Pb[0] == 0:
        return P, Q
    else:
        return P, P + Q


def EllipticCurveIsogenyFactored(E, P, order=None, velu_bound=400):
    """
    Works similarly to EllipticCurveHom_composite
    but with two main additions:

    Introduces a sparse strategy for prime power
    isogenies, taken from
    https://trac.sagemath.org/ticket/34239
    This should be default soon (9.8 maybe)

    For primes l > 400, we use velusqrt as
    the algorithm. This bound was found by testing
    in tests/test_isogenies.sage

    Additionally, we allow `order` as an optional parameter
    and `velu_bound` controls when sqrtvelu kicks in
    """

    def EllipticCurveHom_velusqrt_setorder(P):
        """
        To speed things up, we manually set the order
        assuming all curves have order (p^2 - 1)^2

        I think this is fixed for 9.8, but not everyone
        will be running the latest SageMath version.
        """
        E = P.curve()
        p = E.base().characteristic()
        E._order = ZZ((p**2 - 1) ** 2)
        return EllipticCurveHom_velusqrt(E, P)

    def evaluate_factored_isogeny(phi_list, P):
        """
        Given a list of isogenies, evaluates the
        point for each isogeny in the list
        """
        for phi in phi_list:
            P = phi(P)
        return P

    def sparse_isogeny_prime_power(P, l, e, split=0.8, velu_bound=2000):
        """
        Compute chain of isogenies quotienting
        out a point P of order l**e
        https://trac.sagemath.org/ticket/34239
        """
        if l > velu_bound:
            isogeny_algorithm = lambda Q, l: EllipticCurveHom_velusqrt_setorder(Q)
        else:
            isogeny_algorithm = lambda Q, l: EllipticCurveIsogeny(
                Q.curve(), Q, degree=l, check=False
            )

        def recursive_sparse_isogeny(Q, k):
            assert k
            if k == 1:  # base case
                return [isogeny_algorithm(Q, l)]

            k1 = int(k * split + 0.5)
            k1 = max(1, min(k - 1, k1))  # clamp to [1, k-1]

            Q1 = l**k1 * Q
            L = recursive_sparse_isogeny(Q1, k - k1)

            Q2 = evaluate_factored_isogeny(L, Q)
            R = recursive_sparse_isogeny(Q2, k1)

            return L + R

        return recursive_sparse_isogeny(P, e)

    # Ensure P is a point on E
    if P.curve() != E:
        raise ValueError(f"The supplied kernel must be a point on the curve E")

    if order:
        P._order = ZZ(order)
    cofactor = P.order()

    # Deal with isomorphisms
    if cofactor == 1:
        return EllipticCurveIsogeny(P.curve(), P)

    ϕ_list = []
    for l, e in cofactor.factor():
        # Compute point Q of order l^e
        D = ZZ(l**e)
        cofactor //= D
        Q = cofactor * P

        # Manually setting the order means
        # Sage won't try and do it for each
        # l-isogeny in the iteration
        Q._order = D

        # Use Q as kernel of degree l^e isogeny
        ψ_list = sparse_isogeny_prime_power(Q, l, e, velu_bound=velu_bound)

        # Map P through chain length e of l-isogenies
        P = evaluate_factored_isogeny(ψ_list, P)
        ϕ_list += ψ_list

    return EllipticCurveHom_composite.from_factors(ϕ_list)


def print_info(str, banner="="):
    """
    Print information with a banner to help
    with visibility during debug printing
    """
    print(banner * 80)
    info = textwrap.fill(str, 80)
    print(info.center(80))
    print(banner * 80)
