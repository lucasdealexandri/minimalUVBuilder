# stop_1

Inputs for `uvbuilder fill`. Conventions: `docs/conventions.md` (uvbuilder repository).

## Mass eigenstates

| state | PDG | kind | charge | mass (protomodel) |
| --- | --- | --- | --- | --- |
| T1 | 1000006 | Dirac fermion | 2/3 | 1.00E+03 [GeV] |
| X1 | 1000022 | real scalar | 0 | 2.00E+02 [GeV] |

## Input parameters

| name | card | type | enters |
| --- | --- | --- | --- |
| MChiC3I1Y4ChiC3bI1Ym4 | BSMMASS 1 | mass | DT1 |
| MPhiI1Y0 | BSMMASS 2 | M² | RX1 |
| RzChiC3I1Y4PhiI1Y0cRbar | BSMYUKAWA 1 | Re | - |
| IzChiC3I1Y4PhiI1Y0cRbar | IMBSMYUKAWA 1 | Im | - |
| RzChiC3I1Y4PhiI1Y0tRbar | BSMYUKAWA 2 | Re | - |
| IzChiC3I1Y4PhiI1Y0tRbar | IMBSMYUKAWA 2 | Im | - |
| RzChiC3I1Y4PhiI1Y0uRbar | BSMYUKAWA 3 | Re | - |
| IzChiC3I1Y4PhiI1Y0uRbar | IMBSMYUKAWA 3 | Im | - |
| lambdaPhiI1Y0PhiI1Y0 | BSMQUARTIC 1 | coupling | RX1 |

## Mass matrices

### UMixDT1: Dirac fermion, single state → T1

$$
\begin{array}{c|c}
 & \mathrm{ChiC3bI1Ym4qm4} \\ \hline
\mathrm{ChiC3I1Y4q4} & M_{\mathrm{ChiC3I1Y4ChiC3bI1Ym4}} \\
\end{array}
$$

### UMixRX1: real scalar, single state → X1

$$
\begin{array}{c|c}
 & \mathrm{PhiI1Y0q0} \\ \hline
\mathrm{PhiI1Y0q0} & M_{\mathrm{PhiI1Y0}}^{2} + \lambda_{\mathrm{PhiI1Y0PhiI1Y0}} v^{2} \\
\end{array}
$$

## Multiplets

| multiplet | spin | SU(2) | Y | color | global charges |
| --- | --- | --- | --- | --- | --- |
| PhiI1Y0 | Scalar | 1 | 0 | 1 | () |
| ChiC3bI1Ym4Lw | LH WeylSpinor | 1 | -2/3 | 3* | () |
| ChiC3I1Y4Lw | LH WeylSpinor | 1 | 2/3 | 3 | () |

## Notes

- v = 246.220569073486 GeV (G_F).
- Complex couplings: y = Re y + i Im y, card entries `Ry...`, `Iy...`.
- Real scalar blocks are in the real basis Phi = (Phi_R + i Phi_I)/sqrt(2).
