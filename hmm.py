
from sage.all import GF, prod, PolynomialRing, random_prime, proof
from sage.all import load
attach("multi_point_eval.pyx")

proof.all(False)

p_bit = 256
d = 128

while True:
    p = random_prime(2**p_bit)
    if p % 4 == 3:
        break

F = GF(p**2, name="i", modulus=[1,0,1])
R = PolynomialRing(F, name="X")
x = R.gens()[0]

f = R.random_element(degree=2*d)
g = R.random_element(degree=d)

def test_1():
    v = (f * g) % (x**d - 1)
    v = v.reverse()
    v = v % X**d
    v = v.reverse()
    return v

def test_2():
    v = (f * g) % (x**d - 1)
    v = ntl_reverse(v)
    v = v.truncate(d)
    v = ntl_reverse(v)
    return v

def test_3():
    v = (f * g) % (x**d - 1)
    return ntl_all(v, d)