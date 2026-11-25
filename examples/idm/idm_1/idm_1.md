# idm_1

Inputs for `uvbuilder fill`. Conventions: `docs/conventions.md` (uvbuilder repository).

## Mass eigenstates

| state | PDG | kind | charge | mass (protomodel) |
| --- | --- | --- | --- | --- |
| H0 | 1000022 | real scalar | 0 | 1.00E+02 [GeV] |
| A0 | 1000023 | real scalar | 0 | 2.00E+02 [GeV] |
| Hp | 1000024 | complex scalar | 1 | 2.50E+02 [GeV] |

## Input parameters

| name | card | type | enters |
| --- | --- | --- | --- |
| MPhiI2Y3 | BSMMASS 1 | M² | CHp, RH0A0 |
| RlambdaPhiI2Y3PhiI2Y3 | BSMQUARTIC 1 | Re | RH0A0 |
| IlambdaPhiI2Y3PhiI2Y3 | IMBSMQUARTIC 1 | Im | RH0A0 |
| lambdaPhiI2Y3PhiI2Y3barP1 | BSMQUARTIC 2 | coupling | CHp, RH0A0 |
| lambdaPhiI2Y3PhiI2Y3barP2 | BSMQUARTIC 3 | coupling | CHp, RH0A0 |

## Mass matrices

### UMixCHp: complex scalar, single state → Hp

$$
\begin{array}{c|c}
 & \mathrm{PhiI2Y3q6} \\ \hline
\mathrm{PhiI2Y3q6bar} & M_{\mathrm{PhiI2Y3}}^{2} + \frac{\lambda_{\mathrm{PhiI2Y3PhiI2Y3barP1}} v^{2}}{4} + \frac{\sqrt{3} \lambda_{\mathrm{PhiI2Y3PhiI2Y3barP2}} v^{2}}{12} \\
\end{array}
$$

### UMixRH0A0: real scalar, orthogonal → H0, A0

$$
\begin{array}{c|cc}
 & \mathrm{PhiI2Y3q0\_R} & \mathrm{PhiI2Y3q0\_I} \\ \hline
\mathrm{PhiI2Y3q0\_R} & M_{\mathrm{PhiI2Y3}}^{2} + \frac{\sqrt{3} \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3PhiI2Y3}} v^{2}}{3} + \frac{\lambda_{\mathrm{PhiI2Y3PhiI2Y3barP1}} v^{2}}{4} - \frac{\sqrt{3} \lambda_{\mathrm{PhiI2Y3PhiI2Y3barP2}} v^{2}}{12} & - \frac{\sqrt{3} \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3PhiI2Y3}} v^{2}}{3} \\
\mathrm{PhiI2Y3q0\_I} & - \frac{\sqrt{3} \mathrm{Im}\,\lambda_{\mathrm{PhiI2Y3PhiI2Y3}} v^{2}}{3} & M_{\mathrm{PhiI2Y3}}^{2} - \frac{\sqrt{3} \mathrm{Re}\,\lambda_{\mathrm{PhiI2Y3PhiI2Y3}} v^{2}}{3} + \frac{\lambda_{\mathrm{PhiI2Y3PhiI2Y3barP1}} v^{2}}{4} - \frac{\sqrt{3} \lambda_{\mathrm{PhiI2Y3PhiI2Y3barP2}} v^{2}}{12} \\
\end{array}
$$

## Multiplets

| multiplet | spin | SU(2) | Y | color | global charges |
| --- | --- | --- | --- | --- | --- |
| PhiI2Y3 | Scalar | 2 | 1/2 | 1 | () |

## Notes

- v = 246.220569073486 GeV (G_F).
- Complex couplings: y = Re y + i Im y, card entries `Ry...`, `Iy...`.
- Real scalar blocks are in the real basis Phi = (Phi_R + i Phi_I)/sqrt(2).
