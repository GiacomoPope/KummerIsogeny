# Kummer Isogeny

This repository contains SageMath classes for `KummerLine`, `KummerPoint` and `KummerLineIsogeny`.

- `KummerLine` is the Kummer Line associated to a Montgomery curve $C : y^2 = x(x^2 + Ax + x)$.
- `KummerPoint` is the point on this line, and is the $x$-coordinate of a point on a Montgomery curve 
  in projective coordinates $x(P) = (X : Z)$.
- `KummerLineIsogeny` computes a composite degree isogeny between two Kummer Lines

There are two main reasons to use $x$-only point arithmetic and isogeny computations:

1. **Fast**: Using $x$-only formula is more efficient. For example, we gain something like a 5-25x speed up 
   against SageMath isogenies
2. **More Torsion**: Working with x-only isogenies allows us to access rational torsion on both the elliptic curve 
   and its quadratic twist. This is used, for example, in [SQISign](https://eprint.iacr.org/2020/1240).

However, $x$-only formula are not always suitable. The biggest issue is we lose our group structure. On an elliptic 
curve, we can always compute $P + Q$ given points $P, Q$. However, for points on the Kummer Line, we can
only compute $x(P) + x(Q)$ if we have $x(P), x(Q), x(P-Q)$. Using a square-root, a $y$-coordinate can always
be recovered to get $P$ *or* $-P$, but more information is needed if this sign needs to be set. Similarly, when
we compute $\phi(xP)$ using `KummerLineIsogeny`, we can only evaluate the image of the $x$-coordinate.

The `KummerLine` and `KummerPoint` classes are simply $x$-only arithmetic on Montgomery curves packaged somewhat 
nicely. Further detail of these classes is included as comments in the code.

The `KummerLineIsogeny` classes implement:
- [A simple and compact algorithm for SIDH with arbitrary degree isogenies](https://ia.cr/2017/1198)
  by Craig Costello and Huseyin Hisil for the small, odd $\ell$ degree isogenies
- [Computing Isogenies between Montgomery Curves Using the Action of (0, 0)](https://ia.cr/2017/1198)
  by Joost Renes for the $2$-isogenies 
- A trick from [A faster way to the CSIDH](https://ia.cr/2018/782) by Michael Meyer and Steffen Reith 
  for faster codomain computations using twisted Edwards curves
- [Faster computation of isogenies of large prime degree](ttps://velusqrt.isogeny.org/) the VeluSqrt formula 
  by Daniel J. Bernstein, Luca De Feo, Antonin Leroux, Benjamin Smith for large $\ell$ degree isogenies

## Examples

### Construction of Kummer Line and Kummer Points
A KummerLine can be constructed straight from a Montgomery curve:

```py
E = EllipticCurve(F, [0,A,0,1,0])
K = KummerLine(E)
```

Or, it can be constructed from the Montgomery coefficient

```py
K = KummerLine(F, A)
```

Additionally, we allow A = (A : C) to be stored projectively and
we can construct this in the following way:

```py
K = KummerLine(F, [A, C])
```

A KummerPoint can be constructed from coordinates

```py
xP = K(X, Z)
```

Where $x(P) = (X : Z)$ is the $x$-coordinate in projective $XZ$-coordinates.

A KummerPoint can also be made straight from an elliptic curve point

```py
E = EllipticCurve(F, [0,A,0,1,0])
P = E.random_point()

K = KummerLine(E)
xP = K(P)
```

### Arithmetic of Kummer Points

We can perform arithmetic in the following way:

```py
# Take an elliptic curve and two points
F = GF(163)
A = 6
E = EllipticCurve(F, [0,A,0,1,0])
P, Q = E.random_point(), E.random_point()

# Map these to the Kummer Line and three Kummer points
K = KummerLine(E)
xP, xQ, xPQ = K(P), K(Q), K(P - Q)

# We can double points in two ways
assert 2*xP == xP.double()
assert 2*xP == K(2*P)

# We can perform differential addition
assert xP.add(xQ, xPQ) == K(P + Q)

# We can perform scalar multiplication
assert 11*xP == K(11*P)

# We can also do a three point ladder
assert xQ.ladder_3_pt(xP, xPQ, 11) == K(P + 11*Q)
```

### Kummer Line Isogenies

The code should feel familiar to those who use the Elliptic Curve isogeny classes.

The main thing to understand is that for these $x$-only formula, we cannot compute a
degree two isogeny with the point $P = (0, 0)$ as the kernel. Most of the time this
can just be avoided when constructing the torsion basis, but it's something to be aware of.

```py
# Supersingular curve
p = 163
F = GF(p^2)
A = 6
E = EllipticCurve(F, [0,A,0,1,0])
P, Q = E.gens()

# We can't use the point (0,0) in E[2] as a kernel
check = (p+1)//2 * P
if check[0] == 0:
    P, Q = Q, P

# Compute an isogeny of degree p+1
assert P.order() == 164
phi = E.isogeny(P, algorithm="factored")
EA = phi.codomain()

# Compute the same isogeny using x-only formula
K = KummerLine(E)
xP, xQ = K(P), K(Q)
psi = KummerLineIsogeny(K, xP, p+1)
EB = psi.codomain().curve()

# Isogenies should be the same up to isomorphism
assert EA.is_isomorphic(EB)

# We can also evaluate points
imQ = phi(Q)
imxQ = psi(xQ)

# Up to isomorphism and an overall sign, the images
# should be the same!
iso = EA.isomorphism_to(EB)

# We can check by comparing x-coordinates
assert iso(imQ)[0] == imxQ.x()
```

## Future Work

There's a lot that could be improved, but the main things I'm thinking about:

- Implement isomorphisms between Kummer Lines, this is simply a case of writing explicit
  isomorphisms between Montgomery curves
- If we have these isomorphisms, we can avoid avoiding the point $(0,0)$ in the kernel by
  mapping to some other isomorphic Kummer Line
- Improve the performance of the `velusqrt` formula. They seem to be underperforming by a factor
  of 5-10x!!
- Implement the ability to compose two isogenies $\phi \circ \psi$ by `phi * psi`.


## Rough Benchmarking

### SQISign

The starting motivation for this code was to write $x$-only isogenies which could be integrated 
into [https://github.com/LearningToSQI/SQISign-SageMath](https://github.com/LearningToSQI/SQISign-SageMath)
to allow our implementation SQISign to more closely match the proof of concept.

Looking at the rough benchmarks, the inclusion of this new code could have something between a 5-10x performance
improvement of the isogeny computations, which by far are the most expensive part of our implementation.

#### Original SQISign Parameters

Using the original prime $p_{6983}$ from the paper 
[SQISign: compact post-quantum signatures from quaternions and isogenies](https://eprint.iacr.org/2020/1240) 
by Luca De Feo, David Kohel, Antonin Leroux, Christophe Petit, and Benjamin Wesolowski,
we get the following results:

```
================================================================================
Computing isogenies of degree: 2^33 * 5^21 * 7^2 * 11 * 31 * 83 * 107 * 137 *
751 * 827 * 3691 * 4019 * 6983
================================================================================
Naive SageMath codomain computation: 4.08866
Naive SageMath point evaluation: 0.34477

Optimised SageMath codomain computation: 1.11831
Optimised SageMath point evaluation: 0.20448

KummerLine codomain computation: 0.18936
KummerLine point evaluation: 0.04760

================================================================================
Computing isogenies of degree: 2 * 3^53 * 43 * 103^2 * 109 * 199 * 227 * 419 *
491 * 569 * 631 * 677 * 857 * 859 * 883 * 1019 * 1171 * 1879 * 2713 * 4283
================================================================================
Naive SageMath codomain computation: 4.46888
Naive SageMath point evaluation: 0.35535

Optimised SageMath codomain computation: 1.81987
Optimised SageMath point evaluation: 0.31383

KummerLine codomain computation: 0.27545
KummerLine point evaluation: 0.06021
```

#### Updated SQISign Parameters

Using the prime $p_{3923}$ from the paper 
[New algorithms for the Deuring correspondence: Towards practical and secure SQISign signatures](https://ia.cr/2022/234)
by Luca De Feo, Antonin Leroux, Patrick Longa and Benjamin Wesolowski
we get the following results:

```
================================================================================
Computing isogenies of degree: 2^65 * 5^2 * 7 * 11 * 19 * 29^2 * 37^2 * 47 * 197
* 263 * 281 * 461 * 521 * 3923
================================================================================
Naive SageMath codomain computation: 1.60352
Naive SageMath point evaluation: 0.12905

Optimised SageMath codomain computation: 0.69290
Optimised SageMath point evaluation: 0.10337

KummerLine codomain computation: 0.16735
KummerLine point evaluation: 0.02182

================================================================================
Computing isogenies of degree: 2 * 3^65 * 13 * 17 * 43 * 79 * 157 * 239 * 271 *
283 * 307 * 563 * 599 * 607 * 619 * 743 * 827 * 941 * 2357
================================================================================
Naive SageMath codomain computation: 2.36184
Naive SageMath point evaluation: 0.22675

Optimised SageMath codomain computation: 1.29857
Optimised SageMath point evaluation: 0.19723

KummerLine codomain computation: 0.18210
KummerLine point evaluation: 0.03099
```

### BSIDH

Another example of where $x$-only isogenies are required is 
[B-SIDH: supersingular isogeny Diffie-Hellman using twisted torsion](https://eprint.iacr.org/2019/1145.pdf)
by Craig Costello. 

The file `example_BSIDH.sage` has two parameter sets which are Example 2 and Example 3 from the paper. 

```
# Example 2
# p + 1 = 2^4 * 3 * 7^16 * 17^9 * 31^8 * 311 * 571 * 1321 * 5119 * 6011 * 14207 * 28477 * 76667 * 315668179
# p - 1 = 2 * 11^18 * 19 * 23^13 * 47 * 79 * 83 * 89 * 151 * 3347 * 17449 * 33461 * 51193 * 258434945441
p = 0x1935BECE108DC6C0AAD0712181BB1A414E6A8AAA6B510FC29826190FE7EDA80F

# Example 3
# p + 1 = 2^110 * 5 * 7^2 * 67 * 223 * 4229 * 9787 * 13399 * 21521 * 32257 * 47353 * 616228535059
# p - 1 = 2 * 3^34 * 11 * 17 * 19^2 * 29 * 37 * 53^2 * 97 * 107 * 109 * 131 * 137 * 197 * 199 * 227 * 251 * 5519 * 9091 * 33997 * 38201 * 460409 * 5781011
p = 0x76042798BBFB78AEBD02490BD2635DEC131ABFFFFFFFFFFFFFFFFFFFFFFFFFFF
```

#### BSIDH: Example 2 parameters

For Alice and Bob's keygen and secret computation times, we find:

```
Alice
Keygen took: 0.74583
Secret took: 0.37905

Bob
Keygen took: 0.63818
Secret took: 0.29245
```

Comparing the performance of KummerLine isogenies with the in-built SageMath isogenies
we find:

```
================================================================================
Computing isogenies of degree: 2^4 * 3 * 7^16 * 17^9 * 31^8 * 311 * 571 * 1321 *
5119 * 6011 * 14207 * 28477 * 76667
================================================================================
Naive SageMath codomain computation: 30.01183
Naive SageMath point evaluation: 3.04108

Optimised SageMath codomain computation: 2.46512
Optimised SageMath point evaluation: 0.54163

KummerLine codomain computation: 0.37883
KummerLine point evaluation: 0.12098

================================================================================
Computing isogenies of degree: 11^18 * 19 * 23^13 * 47 * 79 * 83 * 89 * 151 *
3347 * 17449 * 33461 * 51193
================================================================================
Naive SageMath codomain computation: 24.16028
Naive SageMath point evaluation: 2.46105

Optimised SageMath codomain computation: 2.02521
Optimised SageMath point evaluation: 0.43988

KummerLine codomain computation: 0.80380
KummerLine point evaluation: 0.10258
```

### BSIDH: Example 3 parameters

For Alice and Bob's keygen and secret computation times, we find:

```
Alice
Keygen took: 0.85299
Secret took: 0.47513

Bob
Keygen took: 0.57196
Secret took: 0.29563
```

Comparing the performance of KummerLine isogenies with the in-built SageMath isogenies
we find:

```
================================================================================
Computing isogenies of degree: 2^110 * 5 * 7^2 * 67 * 223 * 4229 * 9787 * 13399
* 21521 * 32257 * 47353
================================================================================
Naive SageMath codomain computation: 29.69447
Naive SageMath point evaluation: 2.99138

Optimised SageMath codomain computation: 2.66121
Optimised SageMath point evaluation: 0.60815

KummerLine codomain computation: 0.43675
KummerLine point evaluation: 0.13993

================================================================================
Computing isogenies of degree: 3^34 * 11 * 17 * 19^2 * 29 * 37 * 53^2 * 97 * 107
* 109 * 131 * 137 * 197 * 199 * 227 * 251 * 5519 * 9091 * 33997 * 38201
================================================================================
Naive SageMath codomain computation: 20.22029
Naive SageMath point evaluation: 2.01559

Optimised SageMath codomain computation: 2.08677
Optimised SageMath point evaluation: 0.44141

KummerLine codomain computation: 0.29352
KummerLine point evaluation: 0.09691
```
