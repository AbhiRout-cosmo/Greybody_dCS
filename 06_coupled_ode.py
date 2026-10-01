#!/usr/bin/env python3
"""
Supplemental Material -- Script 06: Two-channel numerical integration
====================================================================
Integrates the coupled axial-graviton / pseudoscalar system with purely
ingoing horizon data, extracts the greybody factors and the O(zeta) asymmetry
coefficients, and runs the paper's numerical checks:
  (a) convergence of kappa(w),d(w) to the closed forms + flux conservation
  (b) linearity in zeta
  (c) drop-delta-V22 test  (change is O(zeta) normalized = O(zeta^2) physical)
  (d) l=3 coefficients
  (e) high-frequency coefficient c_infinity and the geometric-optics identity

Boundary conditions / extraction
--------------------------------
  horizon (r0 = 2M + eps): purely ingoing, [Psi,Theta] ~ e^{-i w r*}
  infinity (Rmax >> 1/w): decompose into e^{+-i w r*}; T_11=a (axial), T_21=b (conversion)
  kappa(w) = 1/2 (Gax-Gpol)/(zeta Gpol),  d(w) = Ggs/(zeta Gpol)
  flux residual = |T11|^2+|T21|^2+|R11|^2+|R21|^2 - 1

Paper cross-reference: Sec. VI, Fig. 2, App. C.8; kappa_2=0.394775, d_2=0.283579.
Requires: numpy, scipy.  (Runtime a few minutes; reduce grids for a quick check.)
"""
import numpy as np
from scipy.integrate import solve_ivp

M = 1.0
f  = lambda r: 1 - 2*M/r
fp = lambda r: 2*M/r**2
rstar = lambda r: r + 2*M*np.log(r/(2*M) - 1)

def pots(l, zeta, dropdV=False):
    L = l*(l+1); Lam = (l+2)*(l+1)*l*(l-1); Vc = 6*np.sqrt(Lam*zeta)
    dv = 0.0 if dropdV else 1.0
    return (lambda r: f(r)*(L/r**2 - 6*M/r**3),
            lambda r: f(r)*((L/r**2)*(1 + dv*36*zeta*(M/r)**6) + 2*M/r**3),
            lambda r: Vc*M**3*f(r)/r**5)

def rhs2(r, y, w, V11, V22, Vo):
    P, Pp, T, Tp = y; F = f(r); Fp = fp(r)
    return [Pp, -(Fp/F)*Pp + ((V11(r)-w**2)*P + Vo(r)*T)/F**2,
            Tp, -(Fp/F)*Tp + ((V22(r)-w**2)*T + Vo(r)*P)/F**2]

def rhs1(r, y, w, V):
    P, Pp = y; F = f(r); Fp = fp(r)
    return [Pp, -(Fp/F)*Pp + (V(r)-w**2)*P/F**2]

def dec(ps, psr, r, w):
    rs = rstar(r); prs = f(r)*psr
    return (0.5*(ps+(1j/w)*prs)*np.exp(1j*w*rs),
            0.5*(ps-(1j/w)*prs)*np.exp(-1j*w*rs))

def solve(l, w, zeta, rtol=1e-11, atol=1e-13, Rfac=45., Rfloor=400., dropdV=False):
    V11, V22, Vo = pots(l, zeta, dropdV); Rmax = max(Rfloor, Rfac/w); r0 = 2 + 1e-8
    rs0 = rstar(r0); e0 = np.exp(-1j*w*rs0); dp0 = -1j*w/f(r0)*e0; S = {}
    for tag, ic in [('A', [e0, dp0, 0, 0]), ('B', [0, 0, e0, dp0])]:
        s = solve_ivp(rhs2, [r0, Rmax], ic, args=(w, V11, V22, Vo), rtol=rtol, atol=atol, method='DOP853')
        yb = s.y[:, -1]; rb = s.t[-1]
        I1, O1 = dec(yb[0], yb[1], rb, w); I2, O2 = dec(yb[2], yb[3], rb, w)
        S[tag] = (np.array([I1, I2]), np.array([O1, O2]))
    IA, IB = S['A'][0], S['B'][0]; OA, OB = S['A'][1], S['B'][1]
    a, b = np.linalg.solve(np.array([[IA[0], IB[0]], [IA[1], IB[1]]]), [1., 0.])
    Rax, Rsc = a*OA[0]+b*OB[0], a*OA[1]+b*OB[1]
    Gax, Ggs = abs(a)**2, abs(b)**2
    s1 = solve_ivp(rhs1, [r0, Rmax], [e0, dp0], args=(w, V11), rtol=rtol, atol=atol, method='DOP853')
    I1, _ = dec(s1.y[0, -1], s1.y[1, -1], s1.t[-1], w); Gpol = 1/abs(I1)**2
    flux = Gax + Ggs + abs(Rax)**2 + abs(Rsc)**2 - 1.0
    return dict(kappa=0.5*(Gax-Gpol)/(zeta*Gpol), d=Ggs/(zeta*Gpol),
                Gax=Gax, Ggs=Ggs, Gpol=Gpol, flux=flux)

KAP, D2 = 0.3947752084355584, 0.2835794304

print("(a) Convergence + flux conservation  [l=2, zeta=0.01;  kappa_2=%.6f, d_2=%.6f]" % (KAP, D2))
print("  Mw      kappa(w)     d(w)        flux residual")
for w in [0.20, 0.12, 0.08, 0.05, 0.03]:
    r = solve(2, w, 0.01)
    print("  %.3f   %.6f    %.6f    %+.2e" % (w, r['kappa'], r['d'], r['flux']))

print("\n(b) Linearity in zeta  (kappa,d are zeta-independent at leading order)")
for zeta in [0.01, 0.005]:
    r = solve(2, 0.05, zeta)
    print("  zeta=%.3f : kappa=%.6f  d=%.6f" % (zeta, r['kappa'], r['d']))

print("\n(c) Drop-delta-V22 test (remove 36 zeta (M/r)^6 from V22; Mw=0.05)")
print("  zeta     dkappa/zeta   dd/zeta   (constant => O(zeta) normalized = O(zeta^2) physical)")
for zeta in [0.02, 0.01, 0.005]:
    full = solve(2, 0.05, zeta); drop = solve(2, 0.05, zeta, dropdV=True)
    print("  %.3f    %+.4f      %+.4f" % (zeta, (full['kappa']-drop['kappa'])/zeta,
                                          (full['d']-drop['d'])/zeta))

print("\n(d) l=3 coefficients (ratio method 2*val(0.01)-val(0.02), Mw=0.05)")
vals = []
for zz in [0.01, 0.02]:
    r = solve(3, 0.05, zz, Rfloor=300.)
    vals.append((r['kappa'], r['d']))
print("  kappa_3 ~ %.3f ,  d_3 ~ %.3f   (paper: 0.85(1), 0.87(1))"
      % (2*vals[0][0]-vals[1][0], 2*vals[0][1]-vals[1][1]))

print("\n(e) High-frequency c_infinity and geometric-optics identity (l=2, zeta=0.05)")
print("  Mw    asym*(Mw)^2    -Ggs/(zGpol)*(Mw)^2   [asym = (Gax-Gpol)/(z Gpol)]")
for w in [2.0, 4.0]:
    r = solve(2, w, 0.05, Rfac=40., Rfloor=80.)
    asym = (r['Gax']-r['Gpol'])/(0.05*r['Gpol'])
    print("  %.1f   %+.5f       %+.5f" % (w, asym*w**2, -r['Ggs']/(0.05*r['Gpol'])*w**2))
print("  (asym*(Mw)^2 -> -c_infinity ~ -0.053; asym = -Ggs/(zGpol) exactly => pure conversion loss)")
