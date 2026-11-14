# Supplemental Material — Reproducibility package

**Analytic parity-violating greybody factors of a dynamical Chern–Simons black hole**
A. Rout

*Package version v1.1.0 (revised manuscript). The core results — the two closed-form
quadrupole coefficients — are unchanged from v1.0.0; this version adds reproducibility
scripts for the revised phase-sector analysis and the full two-channel S-matrix.*

This package reproduces every quantitative result of the paper: the two closed-form
quadrupole coefficients, the matched-asymptotic normalizations they rest on, the
phase-sector transcendental, the general-relativistic and parity-violating reflection-phase
relations, the complete two-channel S-matrix at ℓ=2, and the two-channel numerical
validation. Each script is standalone, prints results labeled by the paper's equation/table
numbers, and requires only open-source Python packages.

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
| `07_reflection_phase_GR.py` | GR Regge–Wheeler/Zerilli reflection-phase relation arg R_RW − arg R_Z = −2 arctan(Mω/2) (equal moduli, differing phases) | Sec. V.C; Eq. (35) | numeric (scipy) | seconds |
| `08_reflection_phase_dCS.py` | O(ζ) parity-violating axial reflection-phase shift Δφ_R: ζ-linearity and ~(Mω)⁵ low-frequency suppression | Sec. V.C; Eq. (36) | numeric (scipy) | ~1 min |
| `09_Smatrix_ell2.py` | full 2×2 S-matrix at ℓ=2: unitarity, reciprocity, Table VII probabilities; recovers d₂ from graviton→scalar conversion | App. E; Table VII | numeric (scipy) | ~1 min |

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
- **Phase sector [Sec. V.C].** The reflection-phase difference between the axial (Regge–Wheeler)
  and polar (Zerilli) channels is, at ζ=0, the exact parity-even Chandrasekhar–Detweiler relation
  −2 arctan(Mω/2) [O(Mω), ζ-independent; script 07, matched to all digits]. The genuinely
  parity-violating O(ζ) shift Δφ_R is strongly suppressed, numerically ~(Mω)⁵ up to a possible
  logarithm (script 08), confirming it is negligible at low frequency — far below both the GR
  baseline and the transmission phase c₁=O(ζω).
- **Full S-matrix [App. E].** The complete two-channel S-matrix at ℓ=2 (script 09) is unitary
  and reciprocal to ~1e−12; the pseudoscalar-incidence channel (scalar greybody and
  scalar→graviton conversion) does not modify d₂ or κ₂, which reside in the graviton-incidence column.
- **Numerical vs analytic precision:** the coupled ODE (script 06) reproduces κ₂, d₂ to ~3–4
  digits, limited by O((Mω)²ln Mω) corrections; flux is conserved to ~1e−12–1e−15. The
  high-precision certification is the analytic reduction (script 04), not the ODE.

## Units and conventions
Geometric units G=c=1; M=1 in scripts (restored by dimensional analysis). Coupling
ζ=(ℓ_CS/M)⁴ with α=1. Tortoise coordinate r*=r+2M ln(r/2M−1). Regge–Wheeler master variable
Ψ=(f/r)h₁; scalar ψ=rR. Boundary conditions and extraction are documented in script 06's header.
