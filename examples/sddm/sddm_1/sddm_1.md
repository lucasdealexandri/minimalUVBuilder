# sddm_1

Inputs for `uvbuilder fill`. Conventions: `docs/conventions.md` (uvbuilder repository).

## Mass eigenstates

| state | PDG | kind | charge | mass (protomodel) |
| --- | --- | --- | --- | --- |
| Z1 | 1000022 | real scalar | 0 | 1.09E+02 [GeV] |
| Z2 | 1000023 | real scalar | 0 | 2.90E+02 [GeV] |
| W1 | 1000024 | complex scalar | 1 | 2.00E+02 [GeV] |
| Z3 | 1000025 | real scalar | 0 | 4.31E+02 [GeV] |

## Input parameters

| name | card | type | enters |
| --- | --- | --- | --- |
| MPhiI1Y0 | BSMMASS 1 | M² | RZ1Z2Z3 |
| MPhiI2Y3 | BSMMASS 2 | M² | CW1, RZ1Z2Z3 |
| RmuPhiI1Y0PhiI2Y3 | BSMTRILINEAR 1 | Re | RZ1Z2Z3 |
| ImuPhiI1Y0PhiI2Y3 | IMBSMTRILINEAR 1 | Im | RZ1Z2Z3 |
| lambdaPhiI1Y0PhiI1Y0 | BSMQUARTIC 1 | coupling | RZ1Z2Z3 |
| RlambdaPhiI2Y3PhiI2Y3 | BSMQUARTIC 2 | Re | RZ1Z2Z3 |
| IlambdaPhiI2Y3PhiI2Y3 | IMBSMQUARTIC 2 | Im | RZ1Z2Z3 |
| lambdaPhiI2Y3PhiI2Y3barP1 | BSMQUARTIC 3 | coupling | CW1, RZ1Z2Z3 |
| lambdaPhiI2Y3PhiI2Y3barP2 | BSMQUARTIC 4 | coupling | CW1, RZ1Z2Z3 |

## Mass matrices

### UMixCW1: complex scalar, single state → W1

$$
\begin{array}{c|c}
 & \mathrm{PhiI2Y3q6} \\ \hline
\mathrm{PhiI2Y3q6bar} & M_{\mathrm{PhiI2Y3}}^{2} + \frac{\lambda_{\mathrm{PhiI2Y3PhiI2Y3barP1}} v^{2}}{4} + \frac{\sqrt{3} \lambda_{\mathrm{PhiI2Y3PhiI2Y3barP2}} v^{2}}{12} \\
\end{array}
$$

### UMixRZ1Z2Z3: real scalar, orthogonal → Z1, Z2, Z3

$$
\begin{array}{c|ccc}
 & \mathrm{PhiI1Y0q0} & \mathrm{PhiI2Y3q0\_R} & \mathrm{PhiI2Y3q0\_I} \\ \hline
\mathrm{PhiI1Y0q0} & M_{\mathrm{PhiI1Y0}}^{2} + \lambda_{\mathrm{PhiI1Y0PhiI1Y0}} v^{2} & - \mathrm{Re}\,\mu_{\mathrm{PhiI1Y0PhiI2Y3}} v & \mathrm{Im}\,\mu_{\mathrm{PhiI1Y0PhiI2Y3}} v \\
\mathrm{PhiI2Y3q0\_R} & - \mathrm{Re}\,\mu_{\mathrm{PhiI1Y0PhiI2Y3}} v & M_{\mathrm{PhiI2Y3}}^{2} + \frac{\sqrt{3} \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3PhiI2Y3}} v^{2}}{3} + \frac{\lambda_{\mathrm{PhiI2Y3PhiI2Y3barP1}} v^{2}}{4} - \frac{\sqrt{3} \lambda_{\mathrm{PhiI2Y3PhiI2Y3barP2}} v^{2}}{12} & - \frac{\sqrt{3} \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3PhiI2Y3}} v^{2}}{3} \\
\mathrm{PhiI2Y3q0\_I} & \mathrm{Im}\,\mu_{\mathrm{PhiI1Y0PhiI2Y3}} v & - \frac{\sqrt{3} \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3PhiI2Y3}} v^{2}}{3} & M_{\mathrm{PhiI2Y3}}^{2} - \frac{\sqrt{3} \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3PhiI2Y3}} v^{2}}{3} + \frac{\lambda_{\mathrm{PhiI2Y3PhiI2Y3barP1}} v^{2}}{4} - \frac{\sqrt{3} \lambda_{\mathrm{PhiI2Y3PhiI2Y3barP2}} v^{2}}{12} \\
\end{array}
$$

## Multiplets

| multiplet | spin | SU(2) | Y | color | global charges |
| --- | --- | --- | --- | --- | --- |
| PhiI1Y0 | Scalar | 1 | 0 | 1 | () |
| PhiI2Y3 | Scalar | 2 | 1/2 | 1 | () |

## Notes

- v = 246.220569073486 GeV (G_F).
- Complex couplings: y = Re y + i Im y, card entries `Ry...`, `Iy...`.
- Real scalar blocks are in the real basis Phi = (Phi_R + i Phi_I)/sqrt(2).
