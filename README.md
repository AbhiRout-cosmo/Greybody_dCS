# Supplemental Material — Reproducibility package

**Parity-violating greybody factors of a dynamical Chern–Simons black hole:
closed-form absorption asymmetry and graviton-to-pseudoscalar conversion**
A. Rout

This package reproduces every quantitative result of the paper: the two closed-form
quadrupole coefficients, the matched-asymptotic normalizations they rest on, the
phase-sector transcendental, and the two-channel numerical validation. Each script is
standalone, prints results labeled by the paper's equation/table numbers, and requires
only open-source Python packages.

## Dependencies
- Python ≥ 3.9
- `sympy` (symbolic algebra), `mpmath` (arbitrary precision), `numpy`, `scipy` (ODE integration)

Install: `pip install -r requirements.txt`

## Contents and paper cross-reference

| Script | Reproduces | Paper items | Method | Runtime |
|---|---|---|---|---|
| `01_static_solutions.py` | ρ, σ, p, q solve the static operators; W_{ρσ}=+1, W_{pq}=−1 | Eqs. (15)–(19); App. A.2, B.1, C.1 | symbolic (sympy) | seconds |
| `02_matching_and_ratio.py` | T_g=16iω³/15, Γ_pol=(4/225)(2Mω)⁶; T_s=8iω³/45; ratio → 1/36 | Eq. (20); App. B.1, B.5 | numeric (scipy), ω→0 | ~1 min |
| `03_conversion_d2.py` | J=π²/12−3/4; H(x); d₂=3(π²−9)²/8 | Eqs. (24)–(27); App. B | symbolic + quadrature | seconds |
| `04_kappa2_reduction.py` | D₂^κ; 7 ln2-free monomials (u=2/r); κ₂ in proven basis, ζ(3) coeff = 0 | Eqs. (28)–(31); App. C | mpmath (58-dig) + PSLQ + sympy | ~1–2 min |
| `05_phase_D2.py` | D₂ = ζ(3)−9/8 (fully symbolic, no PSLQ); c₁=2π²D₂ | Eq. (37); Table I | symbolic (polygamma masters) | seconds |
| `06_coupled_ode.py` | convergence + flux; linearity; drop-δV₂₂; ℓ=3; c_∞ | Sec. VI; Fig. 2; App. C.8 | numeric (scipy) | ~2–4 min |

## Notes on precision and honest status
- **d₂ [Eq. (27)]** is an exact closed form: every factor (J, the vertex 1944, the ratio 1/36)
  is derived symbolically (scripts 01–03).
- **κ₂ [Eq. (31)]** is determined analytically **in a proven basis**: script 04 shows the
  u=2/r reduction makes the integrand a non-alternating polylogarithm whose seven monomials are
  *identically* free of ln2, placing κ₂ in span_Q{1, π², ζ(3), π⁴}; PSLQ then returns the
  ζ(3) coefficient as exactly 0. Raise `dps=105, NT=520` in script 04 for the full 78-digit
  quadrature↔closed-form agreement (default settings give ~48 digits, ample for the identification).
- **D₂ = ζ(3)−9/8 [Eq. (37)]** is fully symbolic (script 05). The overall phase prefactor
  2π² is a numerical identification (four digits; script 06 part e and the paper's Table I footnote).
- **Numerical vs analytic precision:** the coupled ODE (script 06) reproduces κ₂, d₂ to ~3–4
  digits, limited by O((Mω)²ln Mω) corrections; flux is conserved to ~1e−12–1e−15. The
  high-precision certification is the analytic reduction (script 04), not the ODE.

## Units and conventions
Geometric units G=c=1; M=1 in scripts (restored by dimensional analysis). Coupling
ζ=(ℓ_CS/M)⁴ with α=1. Tortoise coordinate r*=r+2M ln(r/2M−1). Regge–Wheeler master variable
Ψ=(f/r)h₁; scalar ψ=rR. Boundary conditions and extraction are documented in script 06's header.
