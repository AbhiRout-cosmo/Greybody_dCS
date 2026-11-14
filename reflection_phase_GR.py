#!/usr/bin/env python3
"""
reflection_phase_GR.py
----------------------
Reproducibility script for:
  "Analytic parity-violating greybody factors of a dynamical Chern-Simons
   black hole" (arXiv:2610.10280).

Backs: Sec. V.C (Phase sector) -- the parity-EVEN general-relativistic
reflection-phase difference between the Regge-Wheeler (axial) and Zerilli
(polar) channels, i.e. the Chandrasekhar-Detweiler relation

    arg R_RW - arg R_Z = -2 arctan[ 12 M w / ((l-1)l(l+1)(l+2)) ]
                       -> -2 arctan(M w / 2)   for l = 2.

This is O(M w), zeta-independent, and present already in GR (NOT a signature of
parity violation). The script confirms, for l = 2, that the Regge-Wheeler and
Zerilli potentials share identical |T| and |R| moduli (isospectrality) while
their reflection PHASES differ by exactly the closed form above.

Requires: numpy, scipy.
"""
import numpy as np
from scipy.integrate import solve_ivp

M = 1.0
f  = lambda r: 1 - 2*M/r
fp = lambda r: 2*M/r**2
rstar = lambda r: r + 2*M*np.log(r/(2*M) - 1)

def VRW(r, l=2):
    """Regge-Wheeler (axial, spin-2) potential."""
    L = l*(l+1)
    return f(r)*(L/r**2 - 6*M/r**3)

def VZ(r, l=2):
    """Zerilli (polar) potential."""
    lam = (l-1)*(l+2)/2.0
    num = 2*lam**2*(lam+1)*r**3 + 6*lam**2*M*r**2 + 18*lam*M**2*r + 18*M**3
    den = r**3*(lam*r + 3*M)**2
    return f(r)*num/den

def rhs(r, y, w, V):
    P, Pp = y; F = f(r); Fp = fp(r)
    return [Pp, -(Fp/F)*Pp + (V(r) - w**2)*P/F**2]

def scatter(V, w, r0=2+1e-7, Rfac=60., Rfloor=800.):
    """Single-channel reflection/transmission from purely ingoing horizon data."""
    Rmax = max(Rfloor, Rfac/w)
    rs0 = rstar(r0); e0 = np.exp(-1j*w*rs0); dp0 = -1j*w/f(r0)*e0
    s = solve_ivp(rhs, [r0, Rmax], [e0, dp0], args=(w, V),
                  rtol=1e-12, atol=1e-14, method='DOP853')
    assert s.success
    rb = s.t[-1]; Psi = s.y[0, -1]; Psip = s.y[1, -1]
    rs = rstar(rb); d = f(rb)*Psip
    A = 0.5*(Psi - d/(1j*w))*np.exp(1j*w*rs)   # outgoing coeff at infinity
    B = 0.5*(Psi + d/(1j*w))*np.exp(-1j*w*rs)  # incoming coeff at infinity
    return B/A, 1.0/A                           # (reflection, transmission)

if __name__ == "__main__":
    print("l=2 GR:  Regge-Wheeler (axial) vs Zerilli (polar)")
    print(" Mw    |T_RW|     |T_Z|      ||T|eq  |R_RW|  |R_Z|   "
          "argR_RW-argR_Z   -2*atan(Mw/2)   ratio")
    for w in [0.30, 0.20, 0.12, 0.08, 0.05, 0.03, 0.02]:
        Rr, Tr = scatter(VRW, w); Rz, Tz = scatter(VZ, w)
        dphi = np.angle(Rr) - np.angle(Rz)
        dphi = (dphi + np.pi) % (2*np.pi) - np.pi
        ref = -2*np.arctan(M*w/2)               # Chandrasekhar-Detweiler, l=2
        print("%.3f  %.3e %.3e  %.1e  %.4f %.4f   %+.6f     %+.6f   %.4f" % (
            w, abs(Tr), abs(Tz), abs(abs(Tr)-abs(Tz)), abs(Rr), abs(Rz),
            dphi, ref, dphi/ref))
