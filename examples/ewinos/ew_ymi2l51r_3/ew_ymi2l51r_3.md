# ew_ymi2l51r_3

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
| MPhiI2Y3F1 | BSMMASS 1 | M² | CC1C2, RN1N2N3N4 |
| MPhiI2Y3F2 | BSMMASS 2 | M² | CC1C2, RN1N2N3N4 |
| RlambdaPhiI2Y3F1PhiI2Y3F1 | BSMQUARTIC 1 | Re | RN1N2N3N4 |
| IlambdaPhiI2Y3F1PhiI2Y3F1 | IMBSMQUARTIC 1 | Im | RN1N2N3N4 |
| lambdaPhiI2Y3F1PhiI2Y3F1barP1 | BSMQUARTIC 2 | coupling | CC1C2, RN1N2N3N4 |
| lambdaPhiI2Y3F1PhiI2Y3F1barP2 | BSMQUARTIC 3 | coupling | CC1C2, RN1N2N3N4 |
| RlambdaPhiI2Y3F1PhiI2Y3F2 | BSMQUARTIC 4 | Re | RN1N2N3N4 |
| IlambdaPhiI2Y3F1PhiI2Y3F2 | IMBSMQUARTIC 4 | Im | RN1N2N3N4 |
| RlambdaPhiI2Y3F1PhiI2Y3F2barP1 | BSMQUARTIC 5 | Re | CC1C2, RN1N2N3N4 |
| IlambdaPhiI2Y3F1PhiI2Y3F2barP1 | IMBSMQUARTIC 5 | Im | CC1C2, RN1N2N3N4 |
| RlambdaPhiI2Y3F1PhiI2Y3F2barP2 | BSMQUARTIC 6 | Re | CC1C2, RN1N2N3N4 |
| IlambdaPhiI2Y3F1PhiI2Y3F2barP2 | IMBSMQUARTIC 6 | Im | CC1C2, RN1N2N3N4 |
| RlambdaPhiI2Y3F2PhiI2Y3F2 | BSMQUARTIC 7 | Re | RN1N2N3N4 |
| IlambdaPhiI2Y3F2PhiI2Y3F2 | IMBSMQUARTIC 7 | Im | RN1N2N3N4 |
| lambdaPhiI2Y3F2PhiI2Y3F2barP1 | BSMQUARTIC 8 | coupling | CC1C2, RN1N2N3N4 |
| lambdaPhiI2Y3F2PhiI2Y3F2barP2 | BSMQUARTIC 9 | coupling | CC1C2, RN1N2N3N4 |

## Mass matrices

### UMixCC1C2: complex scalar, unitary → C1, C2

```math
\begin{array}{c|cc}
 & \mathrm{PhiI2Y3q6F1} & \mathrm{PhiI2Y3q6F2} \\ \hline
\mathrm{PhiI2Y3q6F1bar} & M_{\mathrm{PhiI2Y3F1}}^{2} + \frac{\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F1barP1}} v^{2}}{4} + \frac{\sqrt{3} \lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F1barP2}} v^{2}}{12} & \frac{v^{2} \left(- 3 i \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2barP1}} + 3 \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2barP1}} - \sqrt{3} \left(i \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2barP2}} - \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2barP2}}\right)\right)}{12} \\
\mathrm{PhiI2Y3q6F2bar} & \frac{v^{2} \left(3 i \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2barP1}} + 3 \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2barP1}} + \sqrt{3} \left(i \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2barP2}} + \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2barP2}}\right)\right)}{12} & M_{\mathrm{PhiI2Y3F2}}^{2} + \frac{\lambda_{\mathrm{PhiI2Y3F2PhiI2Y3F2barP1}} v^{2}}{4} + \frac{\sqrt{3} \lambda_{\mathrm{PhiI2Y3F2PhiI2Y3F2barP2}} v^{2}}{12} \\
\end{array}
```

### UMixRN1N2N3N4: real scalar, orthogonal → N1, N2, N3, N4

```math
\begin{array}{c|cccc}
 & \mathrm{PhiI2Y3q0F1\_R} & \mathrm{PhiI2Y3q0F1\_I} & \mathrm{PhiI2Y3q0F2\_R} & \mathrm{PhiI2Y3q0F2\_I} \\ \hline
\mathrm{PhiI2Y3q0F1\_R} & M_{\mathrm{PhiI2Y3F1}}^{2} + \frac{\sqrt{3} \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F1}} v^{2}}{3} + \frac{\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F1barP1}} v^{2}}{4} - \frac{\sqrt{3} \lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F1barP2}} v^{2}}{12} & - \frac{\sqrt{3} \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F1}} v^{2}}{3} & \frac{v^{2} \left(2 \sqrt{3} \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2}} + 3 \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2barP1}} - \sqrt{3} \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2barP2}}\right)}{12} & \frac{v^{2} \left(- 2 \sqrt{3} \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2}} + 3 \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2barP1}} - \sqrt{3} \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2barP2}}\right)}{12} \\
\mathrm{PhiI2Y3q0F1\_I} & - \frac{\sqrt{3} \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F1}} v^{2}}{3} & M_{\mathrm{PhiI2Y3F1}}^{2} - \frac{\sqrt{3} \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F1}} v^{2}}{3} + \frac{\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F1barP1}} v^{2}}{4} - \frac{\sqrt{3} \lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F1barP2}} v^{2}}{12} & \frac{v^{2} \left(- 2 \sqrt{3} \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2}} - 3 \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2barP1}} + \sqrt{3} \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2barP2}}\right)}{12} & \frac{v^{2} \left(- 2 \sqrt{3} \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2}} + 3 \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2barP1}} - \sqrt{3} \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2barP2}}\right)}{12} \\
\mathrm{PhiI2Y3q0F2\_R} & \frac{v^{2} \left(2 \sqrt{3} \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2}} + 3 \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2barP1}} - \sqrt{3} \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2barP2}}\right)}{12} & \frac{v^{2} \left(- 2 \sqrt{3} \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2}} - 3 \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2barP1}} + \sqrt{3} \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2barP2}}\right)}{12} & M_{\mathrm{PhiI2Y3F2}}^{2} + \frac{\sqrt{3} \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3F2PhiI2Y3F2}} v^{2}}{3} + \frac{\lambda_{\mathrm{PhiI2Y3F2PhiI2Y3F2barP1}} v^{2}}{4} - \frac{\sqrt{3} \lambda_{\mathrm{PhiI2Y3F2PhiI2Y3F2barP2}} v^{2}}{12} & - \frac{\sqrt{3} \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3F2PhiI2Y3F2}} v^{2}}{3} \\
\mathrm{PhiI2Y3q0F2\_I} & \frac{v^{2} \left(- 2 \sqrt{3} \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2}} + 3 \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2barP1}} - \sqrt{3} \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2barP2}}\right)}{12} & \frac{v^{2} \left(- 2 \sqrt{3} \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2}} + 3 \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2barP1}} - \sqrt{3} \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3F1PhiI2Y3F2barP2}}\right)}{12} & - \frac{\sqrt{3} \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3F2PhiI2Y3F2}} v^{2}}{3} & M_{\mathrm{PhiI2Y3F2}}^{2} - \frac{\sqrt{3} \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3F2PhiI2Y3F2}} v^{2}}{3} + \frac{\lambda_{\mathrm{PhiI2Y3F2PhiI2Y3F2barP1}} v^{2}}{4} - \frac{\sqrt{3} \lambda_{\mathrm{PhiI2Y3F2PhiI2Y3F2barP2}} v^{2}}{12} \\
\end{array}
```

## Multiplets

| multiplet | spin | SU(2) | Y | color | global charges |
| --- | --- | --- | --- | --- | --- |
| PhiI2Y3F1 | Scalar | 2 | 1/2 | 1 | () |
| PhiI2Y3F2 | Scalar | 2 | 1/2 | 1 | () |

## Notes

- v = 246.220569073486 GeV (G_F).
- Complex couplings: y = Re y + i Im y, card entries `Ry...`, `Iy...`.
- Real scalar blocks are in the real basis Phi = (Phi_R + i Phi_I)/sqrt(2).
