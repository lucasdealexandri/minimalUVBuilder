# idm_2

Inputs for `uvbuilder fill`. Conventions: `docs/conventions.md` (uvbuilder repository).

## Mass eigenstates

| state | PDG | kind | charge | mass (protomodel) |
| --- | --- | --- | --- | --- |
| H0 | 1000022 | Majorana fermion | 0 | 1.00E+02 [GeV] |
| A0 | 1000023 | Majorana fermion | 0 | 2.00E+02 [GeV] |
| Hp | 1000024 | Dirac fermion | 1 | 2.50E+02 [GeV] |

## Input parameters

| name | card | type | enters |
| --- | --- | --- | --- |
| MChiI2Y3ChiI2Ym3 | BSMMASS 1 | mass | DHp, MH0A0 |

## Mass matrices

### UMixDHp: Dirac fermion, single state → Hp

$$
\begin{array}{c|c}
 & \mathrm{ChiI2Ym3qm6} \\ \hline
\mathrm{ChiI2Y3q6} & M_{\mathrm{ChiI2Y3ChiI2Ym3}} \\
\end{array}
$$

### UMixMH0A0: Majorana fermion, Takagi → H0, A0

$$
\begin{array}{c|cc}
 & \mathrm{ChiI2Y3q0} & \mathrm{ChiI2Ym3q0} \\ \hline
\mathrm{ChiI2Y3q0} & 0 & - M_{\mathrm{ChiI2Y3ChiI2Ym3}} \\
\mathrm{ChiI2Ym3q0} & - M_{\mathrm{ChiI2Y3ChiI2Ym3}} & 0 \\
\end{array}
$$

## Multiplets

| multiplet | spin | SU(2) | Y | colour | global charges |
| --- | --- | --- | --- | --- | --- |
| ChiI2Ym3Lw | LH WeylSpinor | 2 | -1/2 | 1 | () |
| ChiI2Y3Lw | LH WeylSpinor | 2 | 1/2 | 1 | () |

## Notes

- v = 246.220569073486 GeV (G_F).
