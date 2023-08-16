import time
from kummer_line import KummerLine

proof.all(False)

p = 1

while p % 4 != 3:
    p = random_prime(2**1024)

F.<i> = GF(p^2, modulus=[1,0,1])
E = EllipticCurve(F, [0,6,0,1,0])

P, Q = E.random_point(),  E.random_point()
PQ = P-Q

K = KummerLine(E)
R, S, RS = map(K, [P,Q,PQ])

t0 = time.process_time_ns()
for _ in range(500):
    _ = R._double()
print((time.process_time_ns() - t0) // (500*1000))

t0 = time.process_time_ns()
for _ in range(500):
    _ = R._add(S, RS)
print((time.process_time_ns() - t0) // (500*1000))

print("-"*100)

t0 = time.process_time_ns()
X, Z = R.XZ()
A, C = R._parent.extract_constants()
for _ in range(500):
    _ = R.xDBL(X, Z, A, C)
print((time.process_time_ns() - t0) // (500*1000))

t0 = time.process_time_ns()
XP, ZP = R.XZ()
XQ, ZQ = S.XZ()
XPQ, ZPQ = RS.XZ()
for _ in range(500):
    _ = R.xADD(XP, ZP, XQ, ZQ, XPQ, ZPQ)
print((time.process_time_ns() - t0) // (500*1000))