"""Exact checks for the computations in 'The whisker threshold for WARM graphs with equal rates'."""
from sympy import Rational, binomial, diff, expand, factor, log, Poly, simplify, sqrt, symbols, nsimplify

t, s, p, a, x = symbols("t s p alpha x", positive=True)

# Lemma (branch): P(t) = (1-t)(1+t)^6 - t(3/2-t)^6, s = 2t, Bernstein coefficients on [0,1].
P = (1 - t) * (1 + t) ** 6 - t * (Rational(3, 2) - t) ** 6
Ps = Poly(expand(P.subs(t, s / 2)), s)   # s here is sigma = 2t
n = Ps.degree()
c = [Ps.coeff_monomial(s**k) for k in range(n + 1)]
bern = [sum(binomial(j, k) / binomial(n, k) * c[k] for k in range(j + 1)) for j in range(n + 1)]
print("degree", n)
print("Bernstein coefficients", bern)
assert bern == [1, Rational(487, 896), Rational(165, 224), Rational(297, 224), Rational(243, 112),
                Rational(711, 224), Rational(1899, 448), Rational(665, 128)]
assert all(b > 0 for b in bern)

# Remark: at alpha = 9 and t = 1/12 the lower bound -11/12 + (13/12) 11^(1/9) is below 1/2, because 11 * 13^9 < 17^9.
assert 11 * 13**9 < 17**9
tt = Rational(1, 12)
assert (tt - 1 + (1 + tt) * ((1 - tt) / tt) ** Rational(1, 9) - Rational(1, 2)).evalf(30) < 0

# Lemma B: (2-p)^2 p - (1+p)^2 (1-p) = (2p-1)(p^2-p+1).
assert expand((2 - p) ** 2 * p - (1 + p) ** 2 * (1 - p) - (2 * p - 1) * (p**2 - p + 1)) == 0

# Lemma (depth): phi'(t) numerator.
phi = log((1 - t) / t) - a * log(1 + 1 / a - t)
dphi = diff(phi, t)
target = -(1 + (1 - a * t) * (a - a * t)) / (a * t * (1 - t) * (1 + 1 / a - t))
assert simplify(dphi - target) == 0
assert simplify(phi.subs(t, 1 / a) - log(a - 1)) == 0

# Example: alpha = 3/2, path with four edges, weights (A, B, B, A).
al = Rational(3, 2)
A = (5 + sqrt(5)) / 4
B = (5 - sqrt(5)) / 4
w = [A, B, B, A]
verts = [[0], [0, 1], [1, 2], [2, 3], [3]]          # edges at each vertex of the path y0..y4
S = [sum(w[e] ** al for e in vs) for vs in verts]
q = {(e, v): w[e] ** al / S[v] for v, vs in enumerate(verts) for e in vs}
for e in range(4):
    assert simplify(q[(e, e)] + q[(e, e + 1)] - w[e]) == 0   # equilibrium
# test vector on the two middle edges (k = 1): path y1 y2 y3, g1 = edge 1, g2 = edge 2
u = {1: 1, 2: -q[(1, 2)] / q[(2, 2)]}
Q = 0
for v, vs in enumerate(verts):
    Q += (al - 1) * sum(q[(e, v)] * u.get(e, 0) ** 2 for e in vs) - al * sum(q[(e, v)] * u.get(e, 0) for e in vs) ** 2
print("Q - (7 sqrt5 - 11)/8 =", simplify(Q - (7 * sqrt(5) - 11) / 8))
assert simplify(Q - (7 * sqrt(5) - 11) / 8) == 0

# Symmetric k-edge star: weights (k+1)/k is an equilibrium and Q is negative definite for alpha <= 2.
k, al2 = symbols("k alpha2", positive=True)
for kk in range(1, 8):
    uu = symbols(f"u0:{kk}")
    qc = Rational(1, kk)
    Qs = (al2 - 1) * qc * sum(v**2 for v in uu) - al2 * (qc * sum(uu)) ** 2 + sum(-(v**2) for v in uu)
    target_s = ((al2 - 1) / kk - 1) * sum(v**2 for v in uu) - al2 / kk**2 * sum(uu) ** 2
    assert expand(Qs - target_s) == 0
    for av in [Rational(1, 1) + Rational(j, 10) for j in range(1, 11)]:
        M = Poly(Qs.subs(al2, av), *uu)
        from sympy import Matrix
        H2 = Matrix(kk, kk, lambda i, j: diff(Qs.subs(al2, av), uu[i], uu[j]) / 2)
        assert all(ev < 0 for ev in H2.eigenvals())
    assert simplify(qc + 1 - Rational(kk + 1, kk)) == 0
print("all checks pass")
