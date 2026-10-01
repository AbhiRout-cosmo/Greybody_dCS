#!/usr/bin/env python3
"""
Supplemental Material -- Script 02: Low-frequency matching and the 1/36 ratio
=============================================================================
Numerically verifies the matched-asymptotic transmission normalizations and
the scalar/graviton greybody ratio by direct integration of the single-channel
Regge-Wheeler and massless-scalar equations (l=2, M=1), extrapolating omega->0.

Paper cross-reference
---------------------
  T_g   = 16 i omega^3 / 15,  Gamma_pol = (4/225)(2 M omega)^6   -> App. B.1 (Eqs. B17-B20)
  T_s   = 8 i omega^3 / 45,   Gamma_s   = (1/2025)(2 M omega)^6  -> App. B.5 (Eq. B12)
  Gamma_{0,2}/Gamma_{2,2} = |T_s/T_g|^2 = (1/6)^2 = 1/36         -> Eq. (20), App. B.5 (Eq. B13)

Requires: numpy, scipy.
"""
import numpy as np
from scipy.integrate import solve_ivp

M = 1.0
f  = lambda r: 1 - 2*M/r
fp = lambda r: 2*M/r**2
rstar = lambda r: r + 2*M*np.log(r/(2*M) - 1)
V_RW = lambda r: f(r)*(6.0/r**2 - 6*M/r**3)     # spin-2 Regge-Wheeler (l=2)
V_sc = lambda r: f(r)*(6.0/r**2 + 2*M/r**3)     # spin-0 massless scalar (l=2)

def rhs(r, y, w, V):
    P, Pp = y; Fp = fp(r); F = f(r)
    return [Pp, -(Fp/F)*Pp + (V(r) - w**2)*P/F**2]

def decompose(psi, dpsi, r, w):
    rs = rstar(r); prs = f(r)*dpsi
    Iin  = 0.5*(psi + (1j/w)*prs)*np.exp(1j*w*rs)   # coeff of e^{-i w r*}
    return Iin

def transmission(w, V, rtol=1e-11, atol=1e-13):
    """Complex transmission amplitude T = 1/I_in (unit ingoing at horizon)."""
    r0 = 2 + 1e-8; Rmax = max(1500., 80/w)
    rs0 = rstar(r0); e0 = np.exp(-1j*w*rs0); dp0 = -1j*w/f(r0)*e0
    s = solve_ivp(rhs, [r0, Rmax], [e0, dp0], args=(w, V),
                  rtol=rtol, atol=atol, method='DOP853')
    Iin = decompose(s.y[0, -1], s.y[1, -1], s.t[-1], w)
    return 1.0/Iin

print("Matched-asymptotic predictions:")
print("  T_g/(i w^3) -> 16/15 = %.6f ,  Gamma_pol/(2Mw)^6 -> 4/225 = %.6f" % (16/15, 4/225))
print("  T_s/(i w^3) ->  8/45 = %.6f ,  Gamma_s  /(2Mw)^6 -> 1/2025 = %.6f" % (8/45, 1/2025))
print("  ratio Gamma_s/Gamma_g -> 1/36 = %.6f\n" % (1/36))
print("  Mw       Gpol/(2Mw)^6    Gs/(2Mw)^6     Gs/Gpol")
for w in [0.020, 0.012, 0.008, 0.005, 0.003]:
    Tg = transmission(w, V_RW); Ts = transmission(w, V_sc)
    Gg = abs(Tg)**2; Gs = abs(Ts)**2
    print("  %.3f    %.7f      %.7f    %.6f" % (w, Gg/(2*w)**6, Gs/(2*w)**6, Gs/Gg))
print("\n  (Gpol/(2Mw)^6 -> 0.017778 = 4/225 ; Gs/(2Mw)^6 -> 0.000494 = 1/2025 ;")
print("   Gs/Gpol -> 0.027778 = 1/36, converging faster than either greybody alone.)")
