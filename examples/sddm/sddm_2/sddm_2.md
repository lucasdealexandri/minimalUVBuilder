# sddm_2

Inputs for `uvbuilder fill`. Conventions: `docs/conventions.md` (uvbuilder repository).

## Mass eigenstates

| state | PDG | kind | charge | mass (protomodel) |
| --- | --- | --- | --- | --- |
| Z1 | 1000022 | Majorana fermion | 0 | 1.09E+02 [GeV] |
| Z2 | 1000023 | Majorana fermion | 0 | 2.90E+02 [GeV] |
| W1 | 1000024 | Dirac fermion | 1 | 2.00E+02 [GeV] |
| Z3 | 1000025 | Majorana fermion | 0 | 4.31E+02 [GeV] |

## Input parameters

| name | card | type | enters |
| --- | --- | --- | --- |
| MChiI1Y0 | BSMMASS 1 | mass | MZ1Z2Z3 |
| MChiI2Y3ChiI2Ym3 | BSMMASS 2 | mass | DW1, MZ1Z2Z3 |
| RyChiI1Y0ChiI2Y3 | BSMYUKAWA 1 | Re | MZ1Z2Z3 |
| IyChiI1Y0ChiI2Y3 | IMBSMYUKAWA 1 | Im | MZ1Z2Z3 |
| RyChiI1Y0ChiI2Ym3 | BSMYUKAWA 2 | Re | MZ1Z2Z3 |
| IyChiI1Y0ChiI2Ym3 | IMBSMYUKAWA 2 | Im | MZ1Z2Z3 |

## Mass matrices

### UMixDW1: Dirac fermion, single state → W1

$$
\begin{array}{c|c}
 & \mathrm{ChiI2Ym3qm6} \\ \hline
\mathrm{ChiI2Y3q6} & M_{\mathrm{ChiI2Y3ChiI2Ym3}} \\
\end{array}
$$

### UMixMZ1Z2Z3: Majorana fermion, Takagi → Z1, Z2, Z3

$$
\begin{array}{c|ccc}
 & \mathrm{ChiI1Y0q0} & \mathrm{ChiI2Y3q0} & \mathrm{ChiI2Ym3q0} \\ \hline
\mathrm{ChiI1Y0q0} & M_{\mathrm{ChiI1Y0}} & - \frac{\sqrt{2} v y_{\mathrm{ChiI1Y0ChiI2Y3}}}{2} & \frac{\sqrt{2} v y_{\mathrm{ChiI1Y0ChiI2Ym3}}}{2} \\
\mathrm{ChiI2Y3q0} & - \frac{\sqrt{2} v y_{\mathrm{ChiI1Y0ChiI2Y3}}}{2} & 0 & - M_{\mathrm{ChiI2Y3ChiI2Ym3}} \\
\mathrm{ChiI2Ym3q0} & \frac{\sqrt{2} v y_{\mathrm{ChiI1Y0ChiI2Ym3}}}{2} & - M_{\mathrm{ChiI2Y3ChiI2Ym3}} & 0 \\
\end{array}
$$

## Multiplets

| multiplet | spin | SU(2) | Y | color | global charges |
| --- | --- | --- | --- | --- | --- |
| ChiI1Y0Lw | LH WeylSpinor | 1 | 0 | 1 | () |
| ChiI2Ym3Lw | LH WeylSpinor | 2 | -1/2 | 1 | () |
| ChiI2Y3Lw | LH WeylSpinor | 2 | 1/2 | 1 | () |

## Notes

- v = 246.220569073486 GeV (G_F).
- Complex couplings: y = Re y + i Im y, card entries `Ry...`, `Iy...`.
