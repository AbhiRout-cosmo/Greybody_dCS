#!/usr/bin/env python3
"""
Supplemental Material -- Script 01: Static (near-zone) solutions and Wronskians
==============================================================================
Parity-violating greybody factors of a dynamical Chern-Simons black hole.

Verifies SYMBOLICALLY (M=1, l=2) that the four static solutions annihilate
their static operators, and computes the two tortoise-coordinate Wronskians.

Paper cross-reference
---------------------
  rho = r^3                                   -> Eq. (15)
  sigma (irregular Regge-Wheeler)             -> Eq. (16)
  W_{rho,sigma} = f(rho sigma' - sigma rho')  -> Eq. (17)  [ = +1 ]
  p = r P_2(r-1),  q = r Q_2(r-1)             -> Eq. (18)
  W_{pq} = f(p q' - q p')                     -> Eq. (19)  [ = -1 ]
  static operators                            -> App. A.2 (V_RW, V_sc)

Requires: sympy.
"""
import sympy as sy

r = sy.symbols('r', positive=True)
M = 1
f = 1 - 2*M/r
L = 6                                   # l(l+1) for l=2

# --- Regge-Wheeler (axial graviton) static operator:  d/dr(f psi') = (L/r^2 - 6M/r^3) psi
def RW(psi):
    return sy.diff(f*sy.diff(psi, r), r) - (L/r**2 - 6*M/r**3)*psi

# --- massless-scalar static operator (psi = r R):    d/dr(f psi') = (L/r^2 + 2M/r^3) psi
def SC(psi):
    return sy.diff(f*sy.diff(psi, r), r) - (L/r**2 + 2*M/r**3)*psi

# Regge-Wheeler solutions [Eqs. (15),(16)]
rho = r**3
sigma = (3*r**4*sy.log((r-2)/r) + 6*r**3 + 6*r**2 + 8*r + 12)/(96*r)

# Scalar solutions [Eq. (18)]: p = r P_2(r-1), q = r Q_2(r-1)
P2 = (3*(r-1)**2 - 1)/2
Q2 = sy.Rational(1, 2)*P2*sy.log(r/(r-2)) - sy.Rational(3, 2)*(r-1)
p = r*P2
q = r*Q2

print("=== Static-operator residuals (should all be 0) ===")
print("  RW[rho]   =", sy.simplify(RW(rho)))
print("  RW[sigma] =", sy.simplify(RW(sigma)))
print("  SC[p]     =", sy.simplify(SC(p)))
print("  SC[q]     =", sy.simplify(SC(q)))

W_rs = sy.simplify(f*(rho*sy.diff(sigma, r) - sigma*sy.diff(rho, r)))
W_pq = sy.simplify(f*(p*sy.diff(q, r) - q*sy.diff(p, r)))
print("\n=== Tortoise Wronskians ===")
print("  W_{rho,sigma} =", W_rs, "  (Eq. 17, expected +1)")
print("  W_{pq}        =", W_pq, "  (Eq. 19, expected -1)")

print("\n=== Asymptotics used in the matching (Script 02) ===")
print("  sigma large-r :", sy.series(sigma, r, sy.oo, 6))     # -> -1/(5 r^2) leading
print("  p     large-r :", sy.series(p, r, sy.oo, 2))         # -> (3/2) r^3
print("  q     large-r :", sy.series(q, r, sy.oo, 4))         # -> (2/15) r^-2
print("  rho(2)=%s  p(2)=%s   (horizon values rho_h=8M^3, p_h=2M)" % (rho.subs(r, 2), p.subs(r, 2)))
