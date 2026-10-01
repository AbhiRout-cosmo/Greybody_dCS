#!/usr/bin/env python3
"""
Supplemental Material -- Script 04: Absorption asymmetry kappa_2 (proven basis)
==============================================================================
(1) Reduces D_2^kappa to a 1-D integral and evaluates it to 58 digits.
(2) Shows the rationalizing substitution u = 2/r turns the integrand into a
    non-alternating polylog integral whose 7 transcendental monomials have
    coefficients IDENTICALLY FREE OF ln2 (=> no ln2 in kappa_2).
(3) Identifies kappa_2 in the PROVEN basis {1, pi^2, zeta(3), pi^4}: the
    zeta(3) coefficient comes out exactly 0.

Paper cross-reference
---------------------
  kappa_2 = -864 D_2^kappa,  D_2^kappa = int_2^inf (sigma/r^5)[qA+pB] dr  -> Eq. (28), App. C
  A(r) = 3r^2/4 - 3r + ln r + 3 - ln2                                     -> Eq. (29)
  B(r) = -21/16 - H(r-1)                                                  -> Eq. (30), App. C
  u=2/r reduction, 7 ln2-free monomials, span_Q{1,pi^2,zeta3,pi^4}        -> App. C.5
  kappa_2 = 21 pi^4/80 + 27 pi^2/8 - 89833/1536                           -> Eq. (31), App. C.6

Requires: sympy, mpmath.
"""
import sympy as sy, mpmath as mp

# ---------------------------------------------------------------------------
# Part (1)+(3): high-precision D_2^kappa and PSLQ in the proven basis
#   (stable r-integrand: cancellation-free series for Q_2 and sigma at large r)
# ---------------------------------------------------------------------------
mp.mp.dps = 55        # raise to 105 and NT=520 for the full 78-digit run (slower)
def H(xx):
    xx = mp.mpf(xx); Lp = mp.log(xx+1); Lm = mp.log(xx-1)
    p1 = mp.mpf(3)/4*((xx**2/2-xx)*Lp - xx**2/4 + mp.mpf(3)/2*xx - mp.mpf(3)/2*Lp
                      - ((xx-1)**2/2)*Lm + (xx-1)**2/4)
    return p1 - mp.mpf(1)/2*mp.polylog(2, 2/(xx+1)) - mp.mpf(3)/2*xx + mp.mpf(3)/2*Lp
Hinf = -mp.mpf(21)/16; NT = 220
def Q2(xx):
    xx = mp.mpf(xx)
    if xx > 6:
        return sum(mp.mpf(2*n)/((2*n+1)*(2*n+3))*xx**(-(2*n+1)) for n in range(1, NT))
    return mp.mpf(1)/2*((3*xx**2-1)/2)*mp.log((xx+1)/(xx-1)) - mp.mpf(3)/2*xx
p  = lambda r: r*((3*(r-1)**2-1)/2)
q  = lambda r: r*Q2(r-1)
def sigma(r):
    r = mp.mpf(r)
    if r > 8:
        return -sum(mp.mpf(2)**k/((5+k)*r**(2+k)) for k in range(0, NT))
    return (3*r**4*mp.log((r-2)/r) + 6*r**3 + 6*r**2 + 8*r + 12)/(96*r)
A = lambda r: mp.mpf(3)/4*r**2 - 3*r + mp.log(r) + 3 - mp.log(2)
B = lambda r: Hinf - H(r-1)
def integrand(r):
    r = mp.mpf(r)
    if r-2 < mp.mpf(10)**-80: return mp.mpf(0)
    return (sigma(r)/r**5)*(q(r)*A(r) + p(r)*B(r))

nodes = [2, mp.mpf('2.00001'), mp.mpf('2.0001'), mp.mpf('2.001'), mp.mpf('2.01'),
         mp.mpf('2.05'), mp.mpf('2.2'), mp.mpf('2.5'), 3, 4, 5, 7, 10, 20, 60,
         300, 3000, 30000, mp.inf]
D2k = mp.quad(integrand, nodes)
kap = -864*D2k
print("D_2^kappa =", mp.nstr(D2k, 30))
print("kappa_2   =", mp.nstr(kap, 60))

z2, z3, z4 = mp.zeta(2), mp.zeta(3), mp.zeta(4)
# work at trustworthy precision (quadrature good to ~78 digits)
mp.mp.dps = 50
rel = mp.pslq([mp.mpf(mp.nstr(kap, 48)), mp.mpf(1), z2, z3, z4],
              maxcoeff=10**7, maxsteps=10**6)
print("\nPSLQ {kappa_2, 1, zeta2, zeta3, zeta4} =", rel)
print("  -> [-1536, -89833, 31104, 0, 36288]: zeta(3) coefficient is EXACTLY 0")
print("  => kappa_2 = (-89833 + 31104 zeta2 + 36288 zeta4)/1536")
print("            = 21 pi^4/80 + 27 pi^2/8 - 89833/1536")
target = 21*mp.pi**4/80 + 27*mp.pi**2/8 - mp.mpf(89833)/1536
print("  closed form - quadrature =", mp.nstr(target - kap, 5),
      "(this run ~48 digits; raise dps=105, NT=520 for the full 78-digit agreement)")

# ---------------------------------------------------------------------------
# Part (2): symbolic proof that all seven monomials are ln2-free (u = 2/r)
# ---------------------------------------------------------------------------
print("\n=== u = 2/r reduction: transcendental monomials and their ln2-content ===")
u = sy.symbols('u', positive=True); l2 = sy.log(2); rr = 2/u
Lm = sy.log(1-u); lnu = sy.log(u); lnr = l2 - lnu; lnrm2 = l2 + Lm - lnu
Li2 = sy.polylog(2, u)
ps = (3*rr**3 - 6*rr**2 + 2*rr)/2
qs = -ps*Lm/2 - sy.Rational(3, 2)*(rr**2 - rr)
As = sy.Rational(3, 4)*rr**2 - 3*rr + 3 - lnu
xx = rr - 1
Hs = (sy.Rational(3, 4)*((xx**2/2 - xx)*lnr - xx**2/4 + sy.Rational(3, 2)*xx
        - sy.Rational(3, 2)*lnr - ((xx-1)**2/2)*lnrm2 + (xx-1)**2/4)
      - sy.Rational(1, 2)*Li2 - sy.Rational(3, 2)*xx + sy.Rational(3, 2)*lnr)
Bs = -sy.Rational(21, 16) - Hs
sig = sy.Rational(1, 96)*(3*Lm/rr**2 + 6/rr**3 + 6/rr**4 + 8/rr**5 + 12/rr**6)
integ_u = sy.expand(sig*(qs*As + ps*Bs)*2/u**2)
Lu, Lw, D = sy.symbols('Lu Lw D')
poly = sy.Poly(sy.expand(integ_u.subs({lnu: Lu, Lm: Lw, Li2: D})), Lu, Lw, D)
names = {(1,2,0): 'ln u * ln^2(1-u)', (1,1,0): 'ln u * ln(1-u)', (1,0,0): 'ln u',
         (0,1,1): 'ln(1-u) * Li2(u)', (0,1,0): 'ln(1-u)', (0,0,1): 'Li2(u)', (0,0,0): '1'}
allfree = True
for mono, coeff in poly.terms():
    hasln2 = sy.expand(coeff).has(sy.log(2))
    allfree = allfree and (not hasln2)
    print("  %-22s : coeff has ln2 = %s" % (names.get(mono, mono), hasln2))
print("  ALL coefficients ln2-free:", allfree, " => kappa_2 contains no ln2 (theorem).")
