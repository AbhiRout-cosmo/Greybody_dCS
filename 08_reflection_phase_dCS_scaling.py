#!/usr/bin/env python3
"""
reflection_phase_dCS_scaling.py
-------------------------------
Reproducibility script for:
  "Analytic parity-violating greybody factors of a dynamical Chern-Simons
   black hole" (arXiv:2610.10280).

Backs: Sec. V.C (Phase sector) -- the genuinely parity-VIOLATING O(zeta)
reflection-phase shift of the AXIAL channel relative to its general-relativistic
(Regge-Wheeler) baseline,

    Dphi_R = arg R_11^ax - arg R_11^RW = O( zeta (M w)^5 ) ,

i.e. strongly suppressed at low frequency (vanishing as M w -> 0), far below the
parity-even GR reflection phase [O(M w), see reflection_phase_GR.py] and the
transmission-phase coefficient c_1 [O(zeta w)].

The script (i) confirms the shift is linear in zeta (genuinely O(zeta)); and
(ii) extracts the leading frequency power, whose local exponent drifts downward
and extrapolates toward ~5 as M w -> 0 (consistent with a (M w)^5 law up to a
possible logarithm). The coupled two-channel potentials are those of Eq. (10):
  V11 = f[ L/r^2 - 6M/r^3 ]                       (axial Regge-Wheeler)
  V22 = f[ (L/r^2)(1 + 36 zeta (M/r)^6) + 2M/r^3] (pseudoscalar)
  Voff = 6 sqrt(Lambda) sqrt(zeta) M^3 f / r^5,  Lambda = (l+2)!/(l-2)!

Requires: numpy, scipy.
"""
import numpy as np
from scipy.integrate import solve_ivp

M = 1.0
f  = lambda r: 1 - 2*M/r
fp = lambda r: 2*M/r**2
rstar = lambda r: r + 2*M*np.log(r/(2*M) - 1)
TOL = dict(rtol=1e-11, atol=1e-13, method='DOP853')

def pots(l, zeta):
    L = l*(l+1); Lam = (l+2)*(l+1)*l*(l-1); Vc = 6*np.sqrt(Lam*zeta)
    return (lambda r: f(r)*(L/r**2 - 6*M/r**3),
            lambda r: f(r)*((L/r**2)*(1 + 36*zeta*(M/r)**6) + 2*M/r**3),
            lambda r: Vc*M**3*f(r)/r**5)

def dec(ps, psr, r, w):
    rs = rstar(r); prs = f(r)*psr
    return (0.5*(ps - (1j/w)*prs)*np.exp(1j*w*rs),
            0.5*(ps + (1j/w)*prs)*np.exp(-1j*w*rs))

def RW_reflect(l, w):
    """Exact zeta=0 axial (single-channel Regge-Wheeler) reflection amplitude."""
    V = lambda r: f(r)*(l*(l+1)/r**2 - 6*M/r**3)
    def rhs(r, y):
        P, Pp = y; F = f(r); Fp = fp(r)
        return [Pp, -(Fp/F)*Pp + (V(r) - w**2)*P/F**2]
    Rmax = max(600., 50./w); r0 = 2 + 1e-8
    rs0 = rstar(r0); e0 = np.exp(-1j*w*rs0); dp0 = -1j*w/f(r0)*e0
    s = solve_ivp(rhs, [r0, Rmax], [e0, dp0], **TOL); assert s.success
    A, B = dec(s.y[0, -1], s.y[1, -1], s.t[-1], w)
    return B/A

def Rax(l, w, zeta):
    """Coupled two-channel axial reflection amplitude R_11^ax at coupling zeta."""
    V11, V22, Vo = pots(l, zeta); Rmax = max(600., 50./w); r0 = 2 + 1e-8
    rs0 = rstar(r0); e0 = np.exp(-1j*w*rs0); dp0 = -1j*w/f(r0)*e0; S = {}
    def rhs2(r, y):
        P, Pp, T, Tp = y; F = f(r); Fp = fp(r)
        return [Pp, -(Fp/F)*Pp + ((V11(r)-w**2)*P + Vo(r)*T)/F**2,
                Tp, -(Fp/F)*Tp + ((V22(r)-w**2)*T + Vo(r)*P)/F**2]
    for tag, ic in [('A', [e0, dp0, 0, 0]), ('B', [0, 0, e0, dp0])]:
        s = solve_ivp(rhs2, [r0, Rmax], ic, **TOL); assert s.success
        yb = s.y[:, -1]; rb = s.t[-1]
        IA, OA = dec(yb[0], yb[1], rb, w); IB, OB = dec(yb[2], yb[3], rb, w)
        S[tag] = (np.array([IA, IB]), np.array([OA, OB]))
    IAA, IBB = S['A'][0], S['B'][0]; OAA, OBB = S['A'][1], S['B'][1]
    Mx = np.array([[IAA[0], IBB[0]], [IAA[1], IBB[1]]], dtype=complex)
    a, b = np.linalg.solve(Mx, np.array([1., 0.], dtype=complex))
    return a*OAA[0] + b*OBB[0]

def dphi(w, zeta, l=2):
    d = np.angle(Rax(l, w, zeta)) - np.angle(RW_reflect(l, w))
    return (d + np.pi) % (2*np.pi) - np.pi

if __name__ == "__main__":
    print("(i) zeta-linearity of Dphi_R at Mw=0.10  (d/zeta should be constant):")
    for z in [0.1, 0.05, 0.02]:
        print("    zeta=%.3f : Dphi_R/zeta = %+.5e" % (z, dphi(0.10, z)/z))

    print("\n(ii) leading frequency power (local exponent -> ~5 as Mw->0):")
    print("    Mw      Dphi_R/zeta       local power")
    ws = [0.20, 0.16, 0.13, 0.10, 0.08, 0.06, 0.05, 0.045, 0.04]
    zeta = 0.1; prev = None
    for w in ws:
        d = dphi(w, zeta)/zeta
        if prev is None:
            print("    %.3f   %+.5e" % (w, d))
        else:
            w0, d0 = prev
            p = np.log(abs(d/d0))/np.log(w/w0)
            print("    %.3f   %+.5e    p=%.2f" % (w, d, p))
        prev = (w, d)
