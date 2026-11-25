# stop_2

Inputs for `uvbuilder fill`. Conventions: `docs/conventions.md` (uvbuilder repository).

## Mass eigenstates

| state | PDG | kind | charge | mass (protomodel) |
| --- | --- | --- | --- | --- |
| T1 | 1000006 | complex scalar | 2/3 | 1.00E+03 [GeV] |
| X1 | 1000022 | Majorana fermion | 0 | 2.00E+02 [GeV] |

## Input parameters

| name | card | type | enters |
| --- | --- | --- | --- |
| MChiI1Y0 | BSMMASS 1 | mass | MX1 |
| MPhiC3I1Y4 | BSMMASS 2 | M² | CT1 |
| RzPhiC3I1Y4ChiI1Y0cRbar | BSMYUKAWA 1 | Re | - |
| IzPhiC3I1Y4ChiI1Y0cRbar | IMBSMYUKAWA 1 | Im | - |
| RzPhiC3I1Y4ChiI1Y0tRbar | BSMYUKAWA 2 | Re | - |
| IzPhiC3I1Y4ChiI1Y0tRbar | IMBSMYUKAWA 2 | Im | - |
| RzPhiC3I1Y4ChiI1Y0uRbar | BSMYUKAWA 3 | Re | - |
| IzPhiC3I1Y4ChiI1Y0uRbar | IMBSMYUKAWA 3 | Im | - |
| lambdaPhiC3I1Y4PhiC3I1Y4bar | BSMQUARTIC 1 | coupling | CT1 |

## Mass matrices

### UMixCT1: complex scalar, single state → T1

$$
\begin{array}{c|c}
 & \mathrm{PhiC3I1Y4q4} \\ \hline
\mathrm{PhiC3I1Y4q4bar} & M_{\mathrm{PhiC3I1Y4}}^{2} + \frac{\lambda_{\mathrm{PhiC3I1Y4PhiC3I1Y4bar}} v^{2}}{2} \\
\end{array}
$$

### UMixMX1: Majorana fermion, single state → X1

$$
\begin{array}{c|c}
 & \mathrm{ChiI1Y0q0} \\ \hline
\mathrm{ChiI1Y0q0} & M_{\mathrm{ChiI1Y0}} \\
\end{array}
$$

## Multiplets

| multiplet | spin | SU(2) | Y | colour | global charges |
| --- | --- | --- | --- | --- | --- |
| PhiC3I1Y4 | Scalar | 1 | 2/3 | 3 | () |
| ChiI1Y0Lw | LH WeylSpinor | 1 | 0 | 1 | () |

## Notes

- v = 246.220569073486 GeV (G_F).
- Complex couplings: y = Re y + i Im y, card entries `Ry...`, `Iy...`.
