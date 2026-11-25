# ew_ymi2l51r_1

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
| MPhiI1Y0F1 | BSMMASS 1 | M² | RN1N2N3N4 |
| MPhiI1Y0F2 | BSMMASS 2 | M² | RN1N2N3N4 |
| MPhiI1Y6 | BSMMASS 3 | M² | CC1C2 |
| MPhiI2Y3 | BSMMASS 4 | M² | CC1C2, RN1N2N3N4 |
| RmuPhiI1Y0F1PhiI2Y3 | BSMTRILINEAR 1 | Re | RN1N2N3N4 |
| ImuPhiI1Y0F1PhiI2Y3 | IMBSMTRILINEAR 1 | Im | RN1N2N3N4 |
| RmuPhiI1Y0F2PhiI2Y3 | BSMTRILINEAR 2 | Re | RN1N2N3N4 |
| ImuPhiI1Y0F2PhiI2Y3 | IMBSMTRILINEAR 2 | Im | RN1N2N3N4 |
| RmuPhiI2Y3PhiI1Y6bar | BSMTRILINEAR 3 | Re | CC1C2 |
| ImuPhiI2Y3PhiI1Y6bar | IMBSMTRILINEAR 3 | Im | CC1C2 |
| lambdaPhiI1Y0F1PhiI1Y0F1 | BSMQUARTIC 1 | coupling | RN1N2N3N4 |
| lambdaPhiI1Y0F1PhiI1Y0F2 | BSMQUARTIC 2 | coupling | RN1N2N3N4 |
| lambdaPhiI1Y0F2PhiI1Y0F2 | BSMQUARTIC 3 | coupling | RN1N2N3N4 |
| lambdaPhiI1Y6PhiI1Y6bar | BSMQUARTIC 4 | coupling | CC1C2 |
| RlambdaPhiI2Y3PhiI2Y3 | BSMQUARTIC 5 | Re | RN1N2N3N4 |
| IlambdaPhiI2Y3PhiI2Y3 | IMBSMQUARTIC 5 | Im | RN1N2N3N4 |
| lambdaPhiI2Y3PhiI2Y3barP1 | BSMQUARTIC 6 | coupling | CC1C2, RN1N2N3N4 |
| lambdaPhiI2Y3PhiI2Y3barP2 | BSMQUARTIC 7 | coupling | CC1C2, RN1N2N3N4 |

## Mass matrices

### UMixCC1C2: complex scalar, unitary → C1, C2

$$
\begin{array}{c|cc}
 & \mathrm{PhiI1Y6q6} & \mathrm{PhiI2Y3q6} \\ \hline
\mathrm{PhiI1Y6q6bar} & M_{\mathrm{PhiI1Y6}}^{2} + \frac{\lambda_{\mathrm{PhiI1Y6PhiI1Y6bar}} v^{2}}{2} & \frac{\sqrt{2} v \left(i \mathrm{Im}\,\mu_{\mathrm{PhiI2Y3PhiI1Y6bar}} + \mathrm{Re}\,\mu_{\mathrm{PhiI2Y3PhiI1Y6bar}}\right)}{2} \\
\mathrm{PhiI2Y3q6bar} & \frac{\sqrt{2} v \left(- i \mathrm{Im}\,\mu_{\mathrm{PhiI2Y3PhiI1Y6bar}} + \mathrm{Re}\,\mu_{\mathrm{PhiI2Y3PhiI1Y6bar}}\right)}{2} & M_{\mathrm{PhiI2Y3}}^{2} + \frac{\lambda_{\mathrm{PhiI2Y3PhiI2Y3barP1}} v^{2}}{4} + \frac{\sqrt{3} \lambda_{\mathrm{PhiI2Y3PhiI2Y3barP2}} v^{2}}{12} \\
\end{array}
$$

### UMixRN1N2N3N4: real scalar, orthogonal → N1, N2, N3, N4

$$
\begin{array}{c|cccc}
 & \mathrm{PhiI1Y0q0F1} & \mathrm{PhiI1Y0q0F2} & \mathrm{PhiI2Y3q0\_R} & \mathrm{PhiI2Y3q0\_I} \\ \hline
\mathrm{PhiI1Y0q0F1} & M_{\mathrm{PhiI1Y0F1}}^{2} + \lambda_{\mathrm{PhiI1Y0F1PhiI1Y0F1}} v^{2} & \frac{\lambda_{\mathrm{PhiI1Y0F1PhiI1Y0F2}} v^{2}}{2} & - \mathrm{Re}\,\mu_{\mathrm{PhiI1Y0F1PhiI2Y3}} v & \mathrm{Im}\,\mu_{\mathrm{PhiI1Y0F1PhiI2Y3}} v \\
\mathrm{PhiI1Y0q0F2} & \frac{\lambda_{\mathrm{PhiI1Y0F1PhiI1Y0F2}} v^{2}}{2} & M_{\mathrm{PhiI1Y0F2}}^{2} + \lambda_{\mathrm{PhiI1Y0F2PhiI1Y0F2}} v^{2} & - \mathrm{Re}\,\mu_{\mathrm{PhiI1Y0F2PhiI2Y3}} v & \mathrm{Im}\,\mu_{\mathrm{PhiI1Y0F2PhiI2Y3}} v \\
\mathrm{PhiI2Y3q0\_R} & - \mathrm{Re}\,\mu_{\mathrm{PhiI1Y0F1PhiI2Y3}} v & - \mathrm{Re}\,\mu_{\mathrm{PhiI1Y0F2PhiI2Y3}} v & M_{\mathrm{PhiI2Y3}}^{2} + \frac{\sqrt{3} \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3PhiI2Y3}} v^{2}}{3} + \frac{\lambda_{\mathrm{PhiI2Y3PhiI2Y3barP1}} v^{2}}{4} - \frac{\sqrt{3} \lambda_{\mathrm{PhiI2Y3PhiI2Y3barP2}} v^{2}}{12} & - \frac{\sqrt{3} \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3PhiI2Y3}} v^{2}}{3} \\
\mathrm{PhiI2Y3q0\_I} & \mathrm{Im}\,\mu_{\mathrm{PhiI1Y0F1PhiI2Y3}} v & \mathrm{Im}\,\mu_{\mathrm{PhiI1Y0F2PhiI2Y3}} v & - \frac{\sqrt{3} \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3PhiI2Y3}} v^{2}}{3} & M_{\mathrm{PhiI2Y3}}^{2} - \frac{\sqrt{3} \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3PhiI2Y3}} v^{2}}{3} + \frac{\lambda_{\mathrm{PhiI2Y3PhiI2Y3barP1}} v^{2}}{4} - \frac{\sqrt{3} \lambda_{\mathrm{PhiI2Y3PhiI2Y3barP2}} v^{2}}{12} \\
\end{array}
$$

## Multiplets

| multiplet | spin | SU(2) | Y | colour | global charges |
| --- | --- | --- | --- | --- | --- |
| PhiI1Y0F1 | Scalar | 1 | 0 | 1 | () |
| PhiI1Y0F2 | Scalar | 1 | 0 | 1 | () |
| PhiI1Y6 | Scalar | 1 | 1 | 1 | () |
| PhiI2Y3 | Scalar | 2 | 1/2 | 1 | () |

## Notes

- v = 246.220569073486 GeV (G_F).
- Complex couplings: y = Re y + i Im y, card entries `Ry...`, `Iy...`.
- Real scalar blocks are in the real basis Phi = (Phi_R + i Phi_I)/sqrt(2).
