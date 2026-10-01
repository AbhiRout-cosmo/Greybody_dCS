#!/usr/bin/env python3
"""
Supplemental Material -- Script 05: Phase-sector integral D_2 = zeta(3) - 9/8
============================================================================
Fully SYMBOLIC evaluation (no PSLQ) of the weight-three phase integral that
sets the transmission-phase coefficient c_1 = 2 pi^2 D_2.

Same rationalizing substitution u = 2/(x+1) as kappa_2, but the phase integral
needs only FOUR polygamma-closed master integrals (it avoids the hard
Euler-sum master), so it evaluates in closed form directly.

Paper cross-reference
---------------------
  D_2 = int_1^inf [ Q2(x)/(x+1) a(x) + P2(x)/(x+1) b(x) ] dx = zeta(3) - 9/8   -> Eq. (37)
  c_1 = 2 pi^2 ( zeta(3) - 9/8 )   [prefactor 2 pi^2 identified numerically]   -> Eq. (37), Table I

Master integrals used (exponent a, over [0,1]):
  int u^a du = 1/(a+1);  int u^a ln u du = -1/(a+1)^2
  int u^a Li2(u) du = zeta2/(a+1) - H_{a+1}/(a+1)^2
  int u^a ln u ln(1-u) du = -psi'(a+2)/(a+1) + H_{a+1}/(a+1)^2
with all polygammas reduced to argument 1+eps by the recurrence, then eps->0.

Requires: sympy.
"""
import sympy as sy

u, eps = sy.symbols('u epsilon', positive=True); g = sy.EulerGamma
l2 = sy.log(2); xp1 = 2/u; x = xp1 - 1
lnu = sy.log(u); Lm = sy.log(1-u); Li2 = sy.polylog(2, u)
P2 = (3*x**2 - 1)/2
Q2 = sy.Rational(1, 2)*P2*(-Lm) - sy.Rational(3, 2)*x            # ln((x+1)/(x-1)) = -ln(1-u)
a_x = sy.Rational(3, 4)*(x-1)**2 + (-lnu)                        # a(x) = 3(x-1)^2/4 + ln((x+1)/2)
lnxp1 = l2 - lnu; lnxm1 = l2 + Lm - lnu
Hx = (sy.Rational(3, 4)*((x**2/2 - x)*lnxp1 - x**2/4 + sy.Rational(3, 2)*x - sy.Rational(3, 2)*lnxp1
        - ((x-1)**2/2)*lnxm1 + (x-1)**2/4) - sy.Rational(1, 2)*Li2 - sy.Rational(3, 2)*x + sy.Rational(3, 2)*lnxp1)
b_x = -sy.Rational(21, 16) - Hx
integ = sy.expand((Q2/xp1*a_x + P2/xp1*b_x)*2/u**2)             # D_2 = int_0^1 integ du

Lu, Lw, D = sy.symbols('Lu Lw D')
poly = sy.Poly(sy.expand(integ.subs({lnu: Lu, Lm: Lw, Li2: D})), Lu, Lw, D)

# polygamma at m+eps reduced to argument 1+eps via the recurrence
DG = sy.digamma(1+eps); PG = sy.polygamma(1, 1+eps)
def dig(m):
    v = DG
    if m >= 1:
        for j in range(1, m): v += 1/(j+eps)
    else:
        for j in range(m, 1): v -= 1/(j+eps)
    return v
def trig(m):
    v = PG
    if m >= 1:
        for j in range(1, m): v -= 1/(j+eps)**2
    else:
        for j in range(m, 1): v += 1/(j+eps)**2
    return v
Hh = lambda k: dig(k+2) + g                # H_{a+1}, a=k+eps
z2 = sy.zeta(2)
def master(mono, k):
    s = k + 1 + eps
    if mono == (0, 0, 0): return 1/s
    if mono == (1, 0, 0): return -1/s**2
    if mono == (0, 0, 1): return z2/s - Hh(k)/s**2
    if mono == (1, 1, 0): return -trig(k+2)/s + Hh(k)/s**2
    raise ValueError(mono)

tot = 0
for mono, coeff in poly.terms():
    cp = sy.Poly(sy.expand(sy.expand(coeff)*u**8), u)   # clear negative powers (shift +8)
    for (pw,), cf in zip(cp.monoms(), cp.coeffs()):
        tot += cf*master(mono, pw-8)
lim = sy.limit(sy.expand(tot), eps, 0)                  # eps->0 with poles cancelling
print("Transcendental monomials (all ln2-free):",
      [tuple(m) for m, _ in poly.terms()])
print("\nD_2 (fully symbolic) =", sy.simplify(lim))
print("             expected =", sy.zeta(3) - sy.Rational(9, 8))
print("       numeric        =", sy.N(lim, 40))
print("\nc_1 = 2 pi^2 D_2 = 2 pi^2 (zeta(3)-9/8) =", sy.N(2*sy.pi**2*lim, 12),
      " [prefactor 2 pi^2 is a numerical identification]")
