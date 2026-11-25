# ew_ymi2l51r_4

Inputs for `uvbuilder fill`. Conventions: `docs/conventions.md` (uvbuilder repository).

## Mass eigenstates

| state | PDG | kind | charge | mass (protomodel) |
| --- | --- | --- | --- | --- |
| N1 | 1000022 | Majorana fermion | 0 | 2.57E+02 [GeV] |
| N2 | 1000023 | Majorana fermion | 0 | 6.17E+02 [GeV] |
| C1 | 1000024 | Dirac fermion | 1 | 6.18E+02 [GeV] |
| N3 | 1000025 | Majorana fermion | 0 | 1.64E+03 [GeV] |
| N4 | 1000035 | Majorana fermion | 0 | 1.64E+03 [GeV] |
| C2 | 1000037 | Dirac fermion | 1 | 1.64E+03 [GeV] |

## Input parameters

| name | card | type | enters |
| --- | --- | --- | --- |
| MChiI1Y0F1 | BSMMASS 1 | mass | MN1N2N3N4 |
| MChiI1Y0F2 | BSMMASS 2 | mass | MN1N2N3N4 |
| MChiI1Y6ChiI1Ym6 | BSMMASS 3 | mass | DC1C2 |
| MChiI2Y3ChiI2Ym3 | BSMMASS 4 | mass | DC1C2, MN1N2N3N4 |
| RyChiI1Y0F1ChiI2Y3 | BSMYUKAWA 1 | Re | MN1N2N3N4 |
| IyChiI1Y0F1ChiI2Y3 | IMBSMYUKAWA 1 | Im | MN1N2N3N4 |
| RyChiI1Y0F1ChiI2Ym3 | BSMYUKAWA 2 | Re | MN1N2N3N4 |
| IyChiI1Y0F1ChiI2Ym3 | IMBSMYUKAWA 2 | Im | MN1N2N3N4 |
| RyChiI1Y0F2ChiI2Y3 | BSMYUKAWA 3 | Re | MN1N2N3N4 |
| IyChiI1Y0F2ChiI2Y3 | IMBSMYUKAWA 3 | Im | MN1N2N3N4 |
| RyChiI1Y0F2ChiI2Ym3 | BSMYUKAWA 4 | Re | MN1N2N3N4 |
| IyChiI1Y0F2ChiI2Ym3 | IMBSMYUKAWA 4 | Im | MN1N2N3N4 |
| RyChiI2Y3ChiI1Ym6 | BSMYUKAWA 5 | Re | DC1C2 |
| IyChiI2Y3ChiI1Ym6 | IMBSMYUKAWA 5 | Im | DC1C2 |
| RyChiI2Ym3ChiI1Y6 | BSMYUKAWA 6 | Re | DC1C2 |
| IyChiI2Ym3ChiI1Y6 | IMBSMYUKAWA 6 | Im | DC1C2 |

## Mass matrices

### UMixDC1C2: Dirac fermion, SVD → C1, C2

$$
\begin{array}{c|cc}
 & \mathrm{ChiI1Ym6qm6} & \mathrm{ChiI2Ym3qm6} \\ \hline
\mathrm{ChiI1Y6q6} & M_{\mathrm{ChiI1Y6ChiI1Ym6}} & - \frac{\sqrt{2} v y_{\mathrm{ChiI2Ym3ChiI1Y6}}}{2} \\
\mathrm{ChiI2Y3q6} & \frac{\sqrt{2} v y_{\mathrm{ChiI2Y3ChiI1Ym6}}}{2} & M_{\mathrm{ChiI2Y3ChiI2Ym3}} \\
\end{array}
$$

### UMixMN1N2N3N4: Majorana fermion, Takagi → N1, N2, N3, N4

$$
\begin{array}{c|cccc}
 & \mathrm{ChiI1Y0q0F1} & \mathrm{ChiI1Y0q0F2} & \mathrm{ChiI2Y3q0} & \mathrm{ChiI2Ym3q0} \\ \hline
\mathrm{ChiI1Y0q0F1} & M_{\mathrm{ChiI1Y0F1}} & 0 & - \frac{\sqrt{2} v y_{\mathrm{ChiI1Y0F1ChiI2Y3}}}{2} & \frac{\sqrt{2} v y_{\mathrm{ChiI1Y0F1ChiI2Ym3}}}{2} \\
\mathrm{ChiI1Y0q0F2} & 0 & M_{\mathrm{ChiI1Y0F2}} & - \frac{\sqrt{2} v y_{\mathrm{ChiI1Y0F2ChiI2Y3}}}{2} & \frac{\sqrt{2} v y_{\mathrm{ChiI1Y0F2ChiI2Ym3}}}{2} \\
\mathrm{ChiI2Y3q0} & - \frac{\sqrt{2} v y_{\mathrm{ChiI1Y0F1ChiI2Y3}}}{2} & - \frac{\sqrt{2} v y_{\mathrm{ChiI1Y0F2ChiI2Y3}}}{2} & 0 & - M_{\mathrm{ChiI2Y3ChiI2Ym3}} \\
\mathrm{ChiI2Ym3q0} & \frac{\sqrt{2} v y_{\mathrm{ChiI1Y0F1ChiI2Ym3}}}{2} & \frac{\sqrt{2} v y_{\mathrm{ChiI1Y0F2ChiI2Ym3}}}{2} & - M_{\mathrm{ChiI2Y3ChiI2Ym3}} & 0 \\
\end{array}
$$

## Multiplets

| multiplet | spin | SU(2) | Y | colour | global charges |
| --- | --- | --- | --- | --- | --- |
| ChiI1Ym6Lw | LH WeylSpinor | 1 | -1 | 1 | () |
| ChiI1Y0F1Lw | LH WeylSpinor | 1 | 0 | 1 | () |
| ChiI1Y0F2Lw | LH WeylSpinor | 1 | 0 | 1 | () |
| ChiI1Y6Lw | LH WeylSpinor | 1 | 1 | 1 | () |
| ChiI2Ym3Lw | LH WeylSpinor | 2 | -1/2 | 1 | () |
| ChiI2Y3Lw | LH WeylSpinor | 2 | 1/2 | 1 | () |

## Notes

- v = 246.220569073486 GeV (G_F).
- Complex couplings: y = Re y + i Im y, card entries `Ry...`, `Iy...`.
