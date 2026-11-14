#!/usr/bin/env python3
"""
Smatrix_two_channel_ell2.py
---------------------------
Reproducibility script for:
  "Analytic parity-violating greybody factors of a dynamical Chern-Simons
   black hole" (arXiv:2610.10280).

Backs: Appendix E ("The complete two-channel S-matrix at l=2") and Table VII,
with the Sec. VI pointer. Builds the full infinity-incidence 2x2 reflection (R)
and transmission (T) matrices of the coupled axial-graviton (channel 1) /
pseudoscalar (channel 2) system at l=2. Both channels are massless, so the
asymptotic wavenumbers coincide and flux conservation takes the channel-symmetric
form  R^dag R + T^dag T = 1.

Verifies:
  * UNITARITY  R^dag R + T^dag T = 1  to ~1e-12 (each incidence column sums to 1);
  * RECIPROCITY  S = S^T  (R symmetric, R12 = R21 to ~1e-13);
  * the graviton-incidence column reproduces the headline result
    |T21|^2/(zeta Gpol) -> d2 = 3(pi^2-9)^2/8 = 0.283579... as Mw -> 0;
  * the pseudoscalar-incidence column: |T22|^2 = GR massless-scalar greybody at
    leading order, and the scalar->graviton conversion |T12|^2 is a distinct
    O(zeta) coefficient (|T12|/|T21| ~ 0.548, frequency-independent).

The parity-violating observables kappa2 (axial diagonal) and d2 (off diagonal)
are confined to the graviton-incidence column and are UNAFFECTED by the
scalar-incidence sector; the full S-matrix completes the scattering problem
without altering them.

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
    """Return (outgoing coeff of e^{+i w r*}, incoming coeff of e^{-i w r*})."""
    rs = rstar(r); prs = f(r)*psr
    return (0.5*(ps - (1j/w)*prs)*np.exp(1j*w*rs),
            0.5*(ps + (1j/w)*prs)*np.exp(-1j*w*rs))

def smatrix(l, w, zeta):
    """Full 2x2 infinity-incidence S-matrix. Channel 1 = axial graviton,
    channel 2 = pseudoscalar. Returns (R, T) with T the transmission into the
    horizon. Flux-normalized (equal wavenumbers)."""
    V11, V22, Vo = pots(l, zeta); Rmax = max(800., 60./w); r0 = 2 + 1e-8
    rs0 = rstar(r0); e0 = np.exp(-1j*w*rs0); dp0 = -1j*w/f(r0)*e0
    def rhs2(r, y):
        P, Pp, T, Tp = y; F = f(r); Fp = fp(r)
        return [Pp, -(Fp/F)*Pp + ((V11(r)-w**2)*P + Vo(r)*T)/F**2,
                Tp, -(Fp/F)*Tp + ((V22(r)-w**2)*T + Vo(r)*P)/F**2]
    C = np.zeros((2, 2), complex); D = np.zeros((2, 2), complex)
    # basis k: purely ingoing-at-horizon in channel k only -> horizon vector e_k
    for k, ic in enumerate([[e0, dp0, 0, 0], [0, 0, e0, dp0]]):
        s = solve_ivp(rhs2, [r0, Rmax], ic, **TOL); assert s.success
        yb = s.y[:, -1]; rb = s.t[-1]
        o1, i1 = dec(yb[0], yb[1], rb, w); o2, i2 = dec(yb[2], yb[3], rb, w)
        C[:, k] = [i1, i2]; D[:, k] = [o1, o2]
    Cinv = np.linalg.inv(C)
    T = Cinv            # horizon transmission matrix
    R = D @ Cinv        # reflection matrix
    return R, T

if __name__ == "__main__":
    d2 = 3*(np.pi**2 - 9)**2/8
    print("Two-channel S-matrix probabilities at l=2, zeta=0.1  (Table VII).")
    print(" Mw    |T11|^2    |T21|^2    |T22|^2    |T12|^2    unit.resid   |R12-R21|")
    for w in [0.20, 0.15, 0.10, 0.07, 0.05, 0.03]:
        R, T = smatrix(2, w, 0.1)
        U = R.conj().T @ R + T.conj().T @ T - np.eye(2)
        print(" %.2f  %.3e  %.3e  %.3e  %.3e   %.1e      %.1e" % (
            w, abs(T[0,0])**2, abs(T[1,0])**2, abs(T[1,1])**2, abs(T[0,1])**2,
            np.max(np.abs(U)), abs(R[0,1]-R[1,0])))

    print("\nGraviton->scalar conversion recovers d2 (ratio -> %.6f as Mw->0):" % d2)
    for w in [0.10, 0.07, 0.05, 0.03]:
        R, T = smatrix(2, w, 0.05)
        R0, T0 = smatrix(2, w, 1e-12)          # zeta->0 : polar(=axial GR) greybody
        gpol = abs(T0[0, 0])**2
        print("  Mw=%.3f : |T21|^2/(zeta Gpol) = %.5f" % (w, abs(T[1,0])**2/(gpol*0.05)))
