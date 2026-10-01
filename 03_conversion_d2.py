#!/usr/bin/env python3
"""
Supplemental Material -- Script 03: Conversion coefficient d_2 (closed form)
===========================================================================
Evaluates the single conversion integral J and assembles d_2 in closed form.

Paper cross-reference
---------------------
  J = int_1^inf Q_2(x)/(x+1) dx = pi^2/12 - 3/4     -> Eqs. (24),(25); App. B (Eq. B11)
  H(x)  antiderivative of Q_2(x)/(x+1)              -> App. B.2 (Eq. B8)
  H(1) = -9/16 - pi^2/12 ,  H(inf) = -21/16         -> Eqs. (B9),(B10)
  d_2 = 1944 * (1/36) * J^2 = 54 J^2 = 3(pi^2-9)^2/8 -> Eqs. (26),(27); App. B.5

Requires: sympy, mpmath.
"""
import sympy as sy, mpmath as mp

x = sy.symbols('x', positive=True)
P2 = (3*x**2 - 1)/2
Q2 = sy.Rational(1, 2)*P2*sy.log((x+1)/(x-1)) - sy.Rational(3, 2)*x

# Hand antiderivative H(x) of Q_2(x)/(x+1)  [App. B.2, Eq. (B8)]
Lp, Lm = sy.log(x+1), sy.log(x-1)
H = (sy.Rational(3, 4)*((x**2/2 - x - sy.Rational(3, 2))*Lp - (x-1)**2/2*Lm
        - x**2/4 + sy.Rational(3, 2)*x + (x-1)**2/4)
     - sy.Rational(1, 2)*sy.polylog(2, 2/(x+1)) - sy.Rational(3, 2)*x + sy.Rational(3, 2)*Lp)

# check H'(x) = Q_2/(x+1) symbolically (collapse polylog(1,z) = -log(1-z) first)
diff_expr = sy.diff(H, x) - Q2/(x+1)
diff_expr = diff_expr.replace(lambda e: e.func == sy.polylog and e.args[0] == 1,
                              lambda e: -sy.log(1 - e.args[1]))
resid = sy.simplify(sy.expand_log(diff_expr, force=True))
print("H'(x) - Q2/(x+1) =", resid, " (expected 0)")

# endpoint values (Li2(1)=pi^2/6 at x->1; Li2(0)=0 at x->inf)
H1   = sy.Rational(-9, 16) - sy.pi**2/12          # Eq. (B9)
Hinf = sy.Rational(-21, 16)                        # Eq. (B10)
J    = Hinf - H1
print("H(1)   =", H1)
print("H(inf) =", Hinf)
print("J = H(inf) - H(1) =", sy.simplify(J), " = (pi^2-9)/12")

# numeric cross-check of J with a cancellation-free Q_2 at large argument
mp.mp.dps = 30
def Q2num(t):
    t = mp.mpf(t)
    if t > 6:                                   # series: Q_2(t) = sum_{n>=1} 2n/((2n+1)(2n+3)) t^-(2n+1)
        return sum(mp.mpf(2*n)/((2*n+1)*(2*n+3))*t**(-(2*n+1)) for n in range(1, 200))
    return mp.mpf(1)/2*((3*t**2-1)/2)*mp.log((t+1)/(t-1)) - mp.mpf(3)/2*t
Jnum = mp.quad(lambda t: Q2num(t)/(t+1), [1, 1.5, 3, 6, 30, 300, mp.inf])
print("J (quadrature) =", mp.nstr(Jnum, 20), " vs pi^2/12-3/4 =",
      mp.nstr(mp.pi**2/12 - mp.mpf(3)/4, 20))

# assembly:  d_2 = (18 sqrt6)^2 * (1/36) * J^2 = 1944/36 * J^2 = 54 J^2
d2 = sy.nsimplify(54*J**2)
print("\nd_2 = 1944 * (1/36) * J^2 = 54 J^2 =", sy.simplify(54*J**2))
print("    = 3(pi^2-9)^2/8 =", sy.simplify(sy.Rational(3, 8)*(sy.pi**2-9)**2))
print("    =", mp.nstr(mp.mpf(str(sy.N(54*J**2, 40))), 15))
print("  vertex (18 sqrt6)^2 = 1944 ; greybody ratio 1/36 (Script 02); J exact above.")
