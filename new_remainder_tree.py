from sage.all import GF, prod, PolynomialRing, random_prime, proof
from sage.rings.generic import ProductTree

proof.all(False)

def product_tree_resultant(hI_tree, poly):
    r"""
    Helper function to evaluate a resultant with `h_I` quickly,
    using the product tree, taken from FastEllipticPolynomial
    sage/src/sage/schemes/elliptic_curves/hom_velusqrt.py

    Original author: Lorenz Panny (2022)
    """
    rems = hI_tree.remainders(poly)
    r = prod(rems)
    s = -1 if len(hI_tree) % 2 == 1 == poly.degree() else 1
    assert r.is_constant()
    return s * r[0]


def janky_product_tree(f, n):
    """
    Product tree of given list of polynomials
    """

    if n == 0:
        return {'left': None, 'right': None, 'poly': [1], 'deg': 0}

    if n == 1:

        # No multiplication is required
        return {
            'left': None,
            'right': None,
            'poly': f[0],
            'deg': f[0].degree()
        }

    else:

        m = n - (n // 2)
        left = janky_product_tree(f[:m], m)
        right = janky_product_tree(f[m:], n - m)
        return {
            'left': left,
            'right': right,
            'poly': left['poly'] * right['poly'],
            'deg': left['deg'] + right['deg'],

        }

def multieval_unscaled(g, ptree_f, n):
    """
    Next function computes g(x) mod f_1(x), ..., g(x) mod f_n(x)
    """

    if n == 0:
        return [[1]]

    g_mod = g % ptree_f['poly']

    if n == 1:
        # Now, we have g corresponds with the initial G but now it is modulus 
        # a leaf of the product tree of f
        return [g_mod]

    else:

        m = n - (n // 2)
        # Reducing g(x) modulo the current node polynomial
        left = multieval_unscaled(
            g_mod, ptree_f['left'], m
        )
        right = multieval_unscaled(
            g_mod, ptree_f['right'], n - m
        )
        return left + right


def poly_mul_middle(f, g):
    fg = f*g


def multieval_scaled(g, glen, f, flen, ptree_f, n):
    """
    Next functions computes the scaled remainder tree
    """

    if n == 0:
        return [[1]]

    if flen == n and glen == n and n > 1:
        fg = list(g)
    else:
        fg = poly_mul_middle(f, flen, g, glen)

    if n == 1:
        # The last coefficient should be the desire modular reduction with linear modulus
        if fg != []:
            return [[fg[-1]]]
        else:
            return [[1]]

    m = n - (n // 2)
    left = multieval_scaled(
        fg,
        flen,
        ptree_f['right']['poly'],
        ptree_f['right']['deg'] + 1,
        ptree_f['left'],
        m,
    )
    right = multieval_scaled(
        fg,
        flen,
        ptree_f['left']['poly'],
        ptree_f['left']['deg'] + 1,
        ptree_f['right'],
        n - m,
    )
    return left + right

while True:
    p = random_prime(2**1024)
    if p % 4 == 3:
        break

F = GF(p**2, name="i", modulus=[1,0,1])
R = PolynomialRing(F, name="X")
X = R.gens()[0]

d = 123
roots = [F.random_element() for _ in range(d)]
tree = ProductTree([X - root for root in roots])
