# ew_ymi2l51r_2

Inputs for `uvbuilder fill`. Conventions: `docs/conventions.md` (uvbuilder repository).

## Mass eigenstates

| state | PDG | kind | charge | mass (protomodel) |
| --- | --- | --- | --- | --- |
| N1 | 1000022 | real scalar | 0 | 2.57E+02 [GeV] |
| N2 | 1000023 | real scalar | 0 | 6.17E+02 [GeV] |
| C1 | 1000024 | complex scalar | 1 | 6.18E+02 [GeV] |
| N3 | 1000025 | real scalar | 0 | 1.64E+03 [GeV] |
| N4 | 1000035 | real scalar | 0 | 1.64E+03 [GeV] |
| C2 | 1000037 | complex scalar | 1 | 1.64E+03 [GeV] |

## Input parameters

| name | card | type | enters |
| --- | --- | --- | --- |
| MPhiI1Y0 | BSMMASS 1 | M² | RN1N2N3N4 |
| MPhiI2Y3 | BSMMASS 2 | M² | CC1C2, RN1N2N3N4 |
| MPhiI3Y0 | BSMMASS 3 | M² | CC1C2, RN1N2N3N4 |
| RmuPhiI1Y0PhiI2Y3 | BSMTRILINEAR 1 | Re | RN1N2N3N4 |
| ImuPhiI1Y0PhiI2Y3 | IMBSMTRILINEAR 1 | Im | RN1N2N3N4 |
| RmuPhiI2Y3PhiI3Y0 | BSMTRILINEAR 2 | Re | CC1C2, RN1N2N3N4 |
| ImuPhiI2Y3PhiI3Y0 | IMBSMTRILINEAR 2 | Im | CC1C2, RN1N2N3N4 |
| lambdaPhiI1Y0PhiI1Y0 | BSMQUARTIC 1 | coupling | RN1N2N3N4 |
| lambdaPhiI1Y0PhiI3Y0 | BSMQUARTIC 2 | coupling | RN1N2N3N4 |
| RlambdaPhiI2Y3PhiI2Y3 | BSMQUARTIC 3 | Re | RN1N2N3N4 |
| IlambdaPhiI2Y3PhiI2Y3 | IMBSMQUARTIC 3 | Im | RN1N2N3N4 |
| lambdaPhiI2Y3PhiI2Y3barP1 | BSMQUARTIC 4 | coupling | CC1C2, RN1N2N3N4 |
| lambdaPhiI2Y3PhiI2Y3barP2 | BSMQUARTIC 5 | coupling | CC1C2, RN1N2N3N4 |
| lambdaPhiI3Y0PhiI3Y0 | BSMQUARTIC 6 | coupling | CC1C2, RN1N2N3N4 |

## Mass matrices

### UMixCC1C2: complex scalar, unitary → C1, C2

```math
\begin{array}{c|cc}
 & \mathrm{PhiI2Y3q6} & \mathrm{PhiI3Y0q6} \\ \hline
\mathrm{PhiI2Y3q6bar} & M_{\mathrm{PhiI2Y3}}^{2} + \frac{\lambda_{\mathrm{PhiI2Y3PhiI2Y3barP1}} v^{2}}{4} + \frac{\sqrt{3} \lambda_{\mathrm{PhiI2Y3PhiI2Y3barP2}} v^{2}}{12} & \frac{\sqrt{6} v \left(i \mathrm{Im}\,\mu_{\mathrm{PhiI2Y3PhiI3Y0}} - \mathrm{Re}\,\mu_{\mathrm{PhiI2Y3PhiI3Y0}}\right)}{6} \\
\mathrm{PhiI3Y0q6bar} & - \frac{\sqrt{6} v \left(i \mathrm{Im}\,\mu_{\mathrm{PhiI2Y3PhiI3Y0}} + \mathrm{Re}\,\mu_{\mathrm{PhiI2Y3PhiI3Y0}}\right)}{6} & M_{\mathrm{PhiI3Y0}}^{2} - \frac{\sqrt{6} \lambda_{\mathrm{PhiI3Y0PhiI3Y0}} v^{2}}{6} \\
\end{array}
```

### UMixRN1N2N3N4: real scalar, orthogonal → N1, N2, N3, N4

```math
\begin{array}{c|cccc}
 & \mathrm{PhiI1Y0q0} & \mathrm{PhiI2Y3q0\_R} & \mathrm{PhiI2Y3q0\_I} & \mathrm{PhiI3Y0q0} \\ \hline
\mathrm{PhiI1Y0q0} & M_{\mathrm{PhiI1Y0}}^{2} + \lambda_{\mathrm{PhiI1Y0PhiI1Y0}} v^{2} & - \mathrm{Re}\,\mu_{\mathrm{PhiI1Y0PhiI2Y3}} v & \mathrm{Im}\,\mu_{\mathrm{PhiI1Y0PhiI2Y3}} v & - \frac{\sqrt{6} \lambda_{\mathrm{PhiI1Y0PhiI3Y0}} v^{2}}{12} \\
\mathrm{PhiI2Y3q0\_R} & - \mathrm{Re}\,\mu_{\mathrm{PhiI1Y0PhiI2Y3}} v & M_{\mathrm{PhiI2Y3}}^{2} + \frac{\sqrt{3} \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3PhiI2Y3}} v^{2}}{3} + \frac{\lambda_{\mathrm{PhiI2Y3PhiI2Y3barP1}} v^{2}}{4} - \frac{\sqrt{3} \lambda_{\mathrm{PhiI2Y3PhiI2Y3barP2}} v^{2}}{12} & - \frac{\sqrt{3} \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3PhiI2Y3}} v^{2}}{3} & - \frac{\sqrt{6} \mathrm{Re}\,\mu_{\mathrm{PhiI2Y3PhiI3Y0}} v}{6} \\
\mathrm{PhiI2Y3q0\_I} & \mathrm{Im}\,\mu_{\mathrm{PhiI1Y0PhiI2Y3}} v & - \frac{\sqrt{3} \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3PhiI2Y3}} v^{2}}{3} & M_{\mathrm{PhiI2Y3}}^{2} - \frac{\sqrt{3} \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3PhiI2Y3}} v^{2}}{3} + \frac{\lambda_{\mathrm{PhiI2Y3PhiI2Y3barP1}} v^{2}}{4} - \frac{\sqrt{3} \lambda_{\mathrm{PhiI2Y3PhiI2Y3barP2}} v^{2}}{12} & \frac{\sqrt{6} \mathrm{Im}\,\mu_{\mathrm{PhiI2Y3PhiI3Y0}} v}{6} \\
\mathrm{PhiI3Y0q0} & - \frac{\sqrt{6} \lambda_{\mathrm{PhiI1Y0PhiI3Y0}} v^{2}}{12} & - \frac{\sqrt{6} \mathrm{Re}\,\mu_{\mathrm{PhiI2Y3PhiI3Y0}} v}{6} & \frac{\sqrt{6} \mathrm{Im}\,\mu_{\mathrm{PhiI2Y3PhiI3Y0}} v}{6} & M_{\mathrm{PhiI3Y0}}^{2} - \frac{\sqrt{6} \lambda_{\mathrm{PhiI3Y0PhiI3Y0}} v^{2}}{6} \\
\end{array}
```

## Multiplets

| multiplet | spin | SU(2) | Y | color | global charges |
| --- | --- | --- | --- | --- | --- |
| PhiI1Y0 | Scalar | 1 | 0 | 1 | () |
| PhiI2Y3 | Scalar | 2 | 1/2 | 1 | () |
| PhiI3Y0 | Scalar | 3 | 0 | 1 | () |

## Notes

- v = 246.220569073486 GeV (G_F).
- Complex couplings: y = Re y + i Im y, card entries `Ry...`, `Iy...`.
- Real scalar blocks are in the real basis Phi = (Phi_R + i Phi_I)/sqrt(2).
