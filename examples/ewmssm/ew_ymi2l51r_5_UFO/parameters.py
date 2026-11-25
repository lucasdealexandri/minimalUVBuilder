# This file was automatically created by FeynRules 2.3.49
# Mathematica version: 14.3.0 for Mac OS X ARM (64-bit) (July 8, 2025)
# Date: Wed 25 Nov 2026 09:43:14



from object_library import all_parameters, Parameter


from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot

# This is a default parameter object representing 0.
ZERO = Parameter(name = 'ZERO',
                 nature = 'internal',
                 type = 'real',
                 value = '0.0',
                 texname = '0')

# User-defined parameters.
MChiI1Y0 = Parameter(name = 'MChiI1Y0',
                     nature = 'external',
                     type = 'real',
                     value = 1,
                     texname = '\\text{MChiI1Y0}',
                     lhablock = 'BSMMASS',
                     lhacode = [ 1 ])

MChiI2Y3ChiI2Ym3 = Parameter(name = 'MChiI2Y3ChiI2Ym3',
                             nature = 'external',
                             type = 'real',
                             value = 1,
                             texname = '\\text{MChiI2Y3ChiI2Ym3}',
                             lhablock = 'BSMMASS',
                             lhacode = [ 2 ])

MChiI3Y0 = Parameter(name = 'MChiI3Y0',
                     nature = 'external',
                     type = 'real',
                     value = 1,
                     texname = '\\text{MChiI3Y0}',
                     lhablock = 'BSMMASS',
                     lhacode = [ 3 ])

RyChiI1Y0ChiI2Y3 = Parameter(name = 'RyChiI1Y0ChiI2Y3',
                             nature = 'external',
                             type = 'real',
                             value = 1,
                             texname = '\\text{RyChiI1Y0ChiI2Y3}',
                             lhablock = 'BSMYUKAWA',
                             lhacode = [ 1 ])

RyChiI1Y0ChiI2Ym3 = Parameter(name = 'RyChiI1Y0ChiI2Ym3',
                              nature = 'external',
                              type = 'real',
                              value = 1,
                              texname = '\\text{RyChiI1Y0ChiI2Ym3}',
                              lhablock = 'BSMYUKAWA',
                              lhacode = [ 2 ])

RyChiI2Y3ChiI3Y0 = Parameter(name = 'RyChiI2Y3ChiI3Y0',
                             nature = 'external',
                             type = 'real',
                             value = 1,
                             texname = '\\text{RyChiI2Y3ChiI3Y0}',
                             lhablock = 'BSMYUKAWA',
                             lhacode = [ 3 ])

RyChiI2Ym3ChiI3Y0 = Parameter(name = 'RyChiI2Ym3ChiI3Y0',
                              nature = 'external',
                              type = 'real',
                              value = 1,
                              texname = '\\text{RyChiI2Ym3ChiI3Y0}',
                              lhablock = 'BSMYUKAWA',
                              lhacode = [ 4 ])

cabi = Parameter(name = 'cabi',
                 nature = 'external',
                 type = 'real',
                 value = 0.227736,
                 texname = '\\theta _c',
                 lhablock = 'CKMBLOCK',
                 lhacode = [ 1 ])

IyChiI1Y0ChiI2Y3 = Parameter(name = 'IyChiI1Y0ChiI2Y3',
                             nature = 'external',
                             type = 'real',
                             value = 0,
                             texname = '\\text{IyChiI1Y0ChiI2Y3}',
                             lhablock = 'IMBSMYUKAWA',
                             lhacode = [ 1 ])

IyChiI1Y0ChiI2Ym3 = Parameter(name = 'IyChiI1Y0ChiI2Ym3',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{IyChiI1Y0ChiI2Ym3}',
                              lhablock = 'IMBSMYUKAWA',
                              lhacode = [ 2 ])

IyChiI2Y3ChiI3Y0 = Parameter(name = 'IyChiI2Y3ChiI3Y0',
                             nature = 'external',
                             type = 'real',
                             value = 0,
                             texname = '\\text{IyChiI2Y3ChiI3Y0}',
                             lhablock = 'IMBSMYUKAWA',
                             lhacode = [ 3 ])

IyChiI2Ym3ChiI3Y0 = Parameter(name = 'IyChiI2Ym3ChiI3Y0',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{IyChiI2Ym3ChiI3Y0}',
                              lhablock = 'IMBSMYUKAWA',
                              lhacode = [ 4 ])

IUMixDC1C21x1 = Parameter(name = 'IUMixDC1C21x1',
                          nature = 'external',
                          type = 'real',
                          value = 0,
                          texname = '\\text{IUMixDC1C21x1}',
                          lhablock = 'IMUMIXDC1C2',
                          lhacode = [ 1, 1 ])

IUMixDC1C21x2 = Parameter(name = 'IUMixDC1C21x2',
                          nature = 'external',
                          type = 'real',
                          value = 0,
                          texname = '\\text{IUMixDC1C21x2}',
                          lhablock = 'IMUMIXDC1C2',
                          lhacode = [ 1, 2 ])

IUMixDC1C22x1 = Parameter(name = 'IUMixDC1C22x1',
                          nature = 'external',
                          type = 'real',
                          value = 0,
                          texname = '\\text{IUMixDC1C22x1}',
                          lhablock = 'IMUMIXDC1C2',
                          lhacode = [ 2, 1 ])

IUMixDC1C22x2 = Parameter(name = 'IUMixDC1C22x2',
                          nature = 'external',
                          type = 'real',
                          value = 0,
                          texname = '\\text{IUMixDC1C22x2}',
                          lhablock = 'IMUMIXDC1C2',
                          lhacode = [ 2, 2 ])

IUMixMN1N2N3N41x1 = Parameter(name = 'IUMixMN1N2N3N41x1',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{IUMixMN1N2N3N41x1}',
                              lhablock = 'IMUMIXMN1N2N3N4',
                              lhacode = [ 1, 1 ])

IUMixMN1N2N3N41x2 = Parameter(name = 'IUMixMN1N2N3N41x2',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{IUMixMN1N2N3N41x2}',
                              lhablock = 'IMUMIXMN1N2N3N4',
                              lhacode = [ 1, 2 ])

IUMixMN1N2N3N41x3 = Parameter(name = 'IUMixMN1N2N3N41x3',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{IUMixMN1N2N3N41x3}',
                              lhablock = 'IMUMIXMN1N2N3N4',
                              lhacode = [ 1, 3 ])

IUMixMN1N2N3N41x4 = Parameter(name = 'IUMixMN1N2N3N41x4',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{IUMixMN1N2N3N41x4}',
                              lhablock = 'IMUMIXMN1N2N3N4',
                              lhacode = [ 1, 4 ])

IUMixMN1N2N3N42x1 = Parameter(name = 'IUMixMN1N2N3N42x1',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{IUMixMN1N2N3N42x1}',
                              lhablock = 'IMUMIXMN1N2N3N4',
                              lhacode = [ 2, 1 ])

IUMixMN1N2N3N42x2 = Parameter(name = 'IUMixMN1N2N3N42x2',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{IUMixMN1N2N3N42x2}',
                              lhablock = 'IMUMIXMN1N2N3N4',
                              lhacode = [ 2, 2 ])

IUMixMN1N2N3N42x3 = Parameter(name = 'IUMixMN1N2N3N42x3',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{IUMixMN1N2N3N42x3}',
                              lhablock = 'IMUMIXMN1N2N3N4',
                              lhacode = [ 2, 3 ])

IUMixMN1N2N3N42x4 = Parameter(name = 'IUMixMN1N2N3N42x4',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{IUMixMN1N2N3N42x4}',
                              lhablock = 'IMUMIXMN1N2N3N4',
                              lhacode = [ 2, 4 ])

IUMixMN1N2N3N43x1 = Parameter(name = 'IUMixMN1N2N3N43x1',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{IUMixMN1N2N3N43x1}',
                              lhablock = 'IMUMIXMN1N2N3N4',
                              lhacode = [ 3, 1 ])

IUMixMN1N2N3N43x2 = Parameter(name = 'IUMixMN1N2N3N43x2',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{IUMixMN1N2N3N43x2}',
                              lhablock = 'IMUMIXMN1N2N3N4',
                              lhacode = [ 3, 2 ])

IUMixMN1N2N3N43x3 = Parameter(name = 'IUMixMN1N2N3N43x3',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{IUMixMN1N2N3N43x3}',
                              lhablock = 'IMUMIXMN1N2N3N4',
                              lhacode = [ 3, 3 ])

IUMixMN1N2N3N43x4 = Parameter(name = 'IUMixMN1N2N3N43x4',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{IUMixMN1N2N3N43x4}',
                              lhablock = 'IMUMIXMN1N2N3N4',
                              lhacode = [ 3, 4 ])

IUMixMN1N2N3N44x1 = Parameter(name = 'IUMixMN1N2N3N44x1',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{IUMixMN1N2N3N44x1}',
                              lhablock = 'IMUMIXMN1N2N3N4',
                              lhacode = [ 4, 1 ])

IUMixMN1N2N3N44x2 = Parameter(name = 'IUMixMN1N2N3N44x2',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{IUMixMN1N2N3N44x2}',
                              lhablock = 'IMUMIXMN1N2N3N4',
                              lhacode = [ 4, 2 ])

IUMixMN1N2N3N44x3 = Parameter(name = 'IUMixMN1N2N3N44x3',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{IUMixMN1N2N3N44x3}',
                              lhablock = 'IMUMIXMN1N2N3N4',
                              lhacode = [ 4, 3 ])

IUMixMN1N2N3N44x4 = Parameter(name = 'IUMixMN1N2N3N44x4',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{IUMixMN1N2N3N44x4}',
                              lhablock = 'IMUMIXMN1N2N3N4',
                              lhacode = [ 4, 4 ])

IVMixDC1C21x1 = Parameter(name = 'IVMixDC1C21x1',
                          nature = 'external',
                          type = 'real',
                          value = 0,
                          texname = '\\text{IVMixDC1C21x1}',
                          lhablock = 'IMVMIXDC1C2',
                          lhacode = [ 1, 1 ])

IVMixDC1C21x2 = Parameter(name = 'IVMixDC1C21x2',
                          nature = 'external',
                          type = 'real',
                          value = 0,
                          texname = '\\text{IVMixDC1C21x2}',
                          lhablock = 'IMVMIXDC1C2',
                          lhacode = [ 1, 2 ])

IVMixDC1C22x1 = Parameter(name = 'IVMixDC1C22x1',
                          nature = 'external',
                          type = 'real',
                          value = 0,
                          texname = '\\text{IVMixDC1C22x1}',
                          lhablock = 'IMVMIXDC1C2',
                          lhacode = [ 2, 1 ])

IVMixDC1C22x2 = Parameter(name = 'IVMixDC1C22x2',
                          nature = 'external',
                          type = 'real',
                          value = 0,
                          texname = '\\text{IVMixDC1C22x2}',
                          lhablock = 'IMVMIXDC1C2',
                          lhacode = [ 2, 2 ])

aEWM1 = Parameter(name = 'aEWM1',
                  nature = 'external',
                  type = 'real',
                  value = 127.9,
                  texname = '\\text{aEWM1}',
                  lhablock = 'SMINPUTS',
                  lhacode = [ 1 ])

Gf = Parameter(name = 'Gf',
               nature = 'external',
               type = 'real',
               value = 0.0000116637,
               texname = 'G_f',
               lhablock = 'SMINPUTS',
               lhacode = [ 2 ])

aS = Parameter(name = 'aS',
               nature = 'external',
               type = 'real',
               value = 0.1184,
               texname = '\\alpha _s',
               lhablock = 'SMINPUTS',
               lhacode = [ 3 ])

RUMixDC1C21x1 = Parameter(name = 'RUMixDC1C21x1',
                          nature = 'external',
                          type = 'real',
                          value = 1,
                          texname = '\\text{RUMixDC1C21x1}',
                          lhablock = 'UMIXDC1C2',
                          lhacode = [ 1, 1 ])

RUMixDC1C21x2 = Parameter(name = 'RUMixDC1C21x2',
                          nature = 'external',
                          type = 'real',
                          value = 0,
                          texname = '\\text{RUMixDC1C21x2}',
                          lhablock = 'UMIXDC1C2',
                          lhacode = [ 1, 2 ])

RUMixDC1C22x1 = Parameter(name = 'RUMixDC1C22x1',
                          nature = 'external',
                          type = 'real',
                          value = 0,
                          texname = '\\text{RUMixDC1C22x1}',
                          lhablock = 'UMIXDC1C2',
                          lhacode = [ 2, 1 ])

RUMixDC1C22x2 = Parameter(name = 'RUMixDC1C22x2',
                          nature = 'external',
                          type = 'real',
                          value = 1,
                          texname = '\\text{RUMixDC1C22x2}',
                          lhablock = 'UMIXDC1C2',
                          lhacode = [ 2, 2 ])

RUMixMN1N2N3N41x1 = Parameter(name = 'RUMixMN1N2N3N41x1',
                              nature = 'external',
                              type = 'real',
                              value = 1,
                              texname = '\\text{RUMixMN1N2N3N41x1}',
                              lhablock = 'UMIXMN1N2N3N4',
                              lhacode = [ 1, 1 ])

RUMixMN1N2N3N41x2 = Parameter(name = 'RUMixMN1N2N3N41x2',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{RUMixMN1N2N3N41x2}',
                              lhablock = 'UMIXMN1N2N3N4',
                              lhacode = [ 1, 2 ])

RUMixMN1N2N3N41x3 = Parameter(name = 'RUMixMN1N2N3N41x3',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{RUMixMN1N2N3N41x3}',
                              lhablock = 'UMIXMN1N2N3N4',
                              lhacode = [ 1, 3 ])

RUMixMN1N2N3N41x4 = Parameter(name = 'RUMixMN1N2N3N41x4',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{RUMixMN1N2N3N41x4}',
                              lhablock = 'UMIXMN1N2N3N4',
                              lhacode = [ 1, 4 ])

RUMixMN1N2N3N42x1 = Parameter(name = 'RUMixMN1N2N3N42x1',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{RUMixMN1N2N3N42x1}',
                              lhablock = 'UMIXMN1N2N3N4',
                              lhacode = [ 2, 1 ])

RUMixMN1N2N3N42x2 = Parameter(name = 'RUMixMN1N2N3N42x2',
                              nature = 'external',
                              type = 'real',
                              value = 1,
                              texname = '\\text{RUMixMN1N2N3N42x2}',
                              lhablock = 'UMIXMN1N2N3N4',
                              lhacode = [ 2, 2 ])

RUMixMN1N2N3N42x3 = Parameter(name = 'RUMixMN1N2N3N42x3',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{RUMixMN1N2N3N42x3}',
                              lhablock = 'UMIXMN1N2N3N4',
                              lhacode = [ 2, 3 ])

RUMixMN1N2N3N42x4 = Parameter(name = 'RUMixMN1N2N3N42x4',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{RUMixMN1N2N3N42x4}',
                              lhablock = 'UMIXMN1N2N3N4',
                              lhacode = [ 2, 4 ])

RUMixMN1N2N3N43x1 = Parameter(name = 'RUMixMN1N2N3N43x1',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{RUMixMN1N2N3N43x1}',
                              lhablock = 'UMIXMN1N2N3N4',
                              lhacode = [ 3, 1 ])

RUMixMN1N2N3N43x2 = Parameter(name = 'RUMixMN1N2N3N43x2',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{RUMixMN1N2N3N43x2}',
                              lhablock = 'UMIXMN1N2N3N4',
                              lhacode = [ 3, 2 ])

RUMixMN1N2N3N43x3 = Parameter(name = 'RUMixMN1N2N3N43x3',
                              nature = 'external',
                              type = 'real',
                              value = 1,
                              texname = '\\text{RUMixMN1N2N3N43x3}',
                              lhablock = 'UMIXMN1N2N3N4',
                              lhacode = [ 3, 3 ])

RUMixMN1N2N3N43x4 = Parameter(name = 'RUMixMN1N2N3N43x4',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{RUMixMN1N2N3N43x4}',
                              lhablock = 'UMIXMN1N2N3N4',
                              lhacode = [ 3, 4 ])

RUMixMN1N2N3N44x1 = Parameter(name = 'RUMixMN1N2N3N44x1',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{RUMixMN1N2N3N44x1}',
                              lhablock = 'UMIXMN1N2N3N4',
                              lhacode = [ 4, 1 ])

RUMixMN1N2N3N44x2 = Parameter(name = 'RUMixMN1N2N3N44x2',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{RUMixMN1N2N3N44x2}',
                              lhablock = 'UMIXMN1N2N3N4',
                              lhacode = [ 4, 2 ])

RUMixMN1N2N3N44x3 = Parameter(name = 'RUMixMN1N2N3N44x3',
                              nature = 'external',
                              type = 'real',
                              value = 0,
                              texname = '\\text{RUMixMN1N2N3N44x3}',
                              lhablock = 'UMIXMN1N2N3N4',
                              lhacode = [ 4, 3 ])

RUMixMN1N2N3N44x4 = Parameter(name = 'RUMixMN1N2N3N44x4',
                              nature = 'external',
                              type = 'real',
                              value = 1,
                              texname = '\\text{RUMixMN1N2N3N44x4}',
                              lhablock = 'UMIXMN1N2N3N4',
                              lhacode = [ 4, 4 ])

RVMixDC1C21x1 = Parameter(name = 'RVMixDC1C21x1',
                          nature = 'external',
                          type = 'real',
                          value = 1,
                          texname = '\\text{RVMixDC1C21x1}',
                          lhablock = 'VMIXDC1C2',
                          lhacode = [ 1, 1 ])

RVMixDC1C21x2 = Parameter(name = 'RVMixDC1C21x2',
                          nature = 'external',
                          type = 'real',
                          value = 0,
                          texname = '\\text{RVMixDC1C21x2}',
                          lhablock = 'VMIXDC1C2',
                          lhacode = [ 1, 2 ])

RVMixDC1C22x1 = Parameter(name = 'RVMixDC1C22x1',
                          nature = 'external',
                          type = 'real',
                          value = 0,
                          texname = '\\text{RVMixDC1C22x1}',
                          lhablock = 'VMIXDC1C2',
                          lhacode = [ 2, 1 ])

RVMixDC1C22x2 = Parameter(name = 'RVMixDC1C22x2',
                          nature = 'external',
                          type = 'real',
                          value = 1,
                          texname = '\\text{RVMixDC1C22x2}',
                          lhablock = 'VMIXDC1C2',
                          lhacode = [ 2, 2 ])

ymdo = Parameter(name = 'ymdo',
                 nature = 'external',
                 type = 'real',
                 value = 0.00504,
                 texname = '\\text{ymdo}',
                 lhablock = 'YUKAWA',
                 lhacode = [ 1 ])

ymup = Parameter(name = 'ymup',
                 nature = 'external',
                 type = 'real',
                 value = 0.00255,
                 texname = '\\text{ymup}',
                 lhablock = 'YUKAWA',
                 lhacode = [ 2 ])

yms = Parameter(name = 'yms',
                nature = 'external',
                type = 'real',
                value = 0.101,
                texname = '\\text{yms}',
                lhablock = 'YUKAWA',
                lhacode = [ 3 ])

ymc = Parameter(name = 'ymc',
                nature = 'external',
                type = 'real',
                value = 1.27,
                texname = '\\text{ymc}',
                lhablock = 'YUKAWA',
                lhacode = [ 4 ])

ymb = Parameter(name = 'ymb',
                nature = 'external',
                type = 'real',
                value = 4.7,
                texname = '\\text{ymb}',
                lhablock = 'YUKAWA',
                lhacode = [ 5 ])

ymt = Parameter(name = 'ymt',
                nature = 'external',
                type = 'real',
                value = 172,
                texname = '\\text{ymt}',
                lhablock = 'YUKAWA',
                lhacode = [ 6 ])

yme = Parameter(name = 'yme',
                nature = 'external',
                type = 'real',
                value = 0.000511,
                texname = '\\text{yme}',
                lhablock = 'YUKAWA',
                lhacode = [ 11 ])

ymm = Parameter(name = 'ymm',
                nature = 'external',
                type = 'real',
                value = 0.10566,
                texname = '\\text{ymm}',
                lhablock = 'YUKAWA',
                lhacode = [ 13 ])

ymtau = Parameter(name = 'ymtau',
                  nature = 'external',
                  type = 'real',
                  value = 1.777,
                  texname = '\\text{ymtau}',
                  lhablock = 'YUKAWA',
                  lhacode = [ 15 ])

MZ = Parameter(name = 'MZ',
               nature = 'external',
               type = 'real',
               value = 91.1876,
               texname = '\\text{MZ}',
               lhablock = 'MASS',
               lhacode = [ 23 ])

Me = Parameter(name = 'Me',
               nature = 'external',
               type = 'real',
               value = 0.000511,
               texname = '\\text{Me}',
               lhablock = 'MASS',
               lhacode = [ 11 ])

MMU = Parameter(name = 'MMU',
                nature = 'external',
                type = 'real',
                value = 0.10566,
                texname = '\\text{MMU}',
                lhablock = 'MASS',
                lhacode = [ 13 ])

MTA = Parameter(name = 'MTA',
                nature = 'external',
                type = 'real',
                value = 1.777,
                texname = '\\text{MTA}',
                lhablock = 'MASS',
                lhacode = [ 15 ])

MU = Parameter(name = 'MU',
               nature = 'external',
               type = 'real',
               value = 0.00255,
               texname = 'M',
               lhablock = 'MASS',
               lhacode = [ 2 ])

MC = Parameter(name = 'MC',
               nature = 'external',
               type = 'real',
               value = 1.27,
               texname = '\\text{MC}',
               lhablock = 'MASS',
               lhacode = [ 4 ])

MT = Parameter(name = 'MT',
               nature = 'external',
               type = 'real',
               value = 172,
               texname = '\\text{MT}',
               lhablock = 'MASS',
               lhacode = [ 6 ])

MD = Parameter(name = 'MD',
               nature = 'external',
               type = 'real',
               value = 0.00504,
               texname = '\\text{MD}',
               lhablock = 'MASS',
               lhacode = [ 1 ])

MS = Parameter(name = 'MS',
               nature = 'external',
               type = 'real',
               value = 0.101,
               texname = '\\text{MS}',
               lhablock = 'MASS',
               lhacode = [ 3 ])

MB = Parameter(name = 'MB',
               nature = 'external',
               type = 'real',
               value = 4.7,
               texname = '\\text{MB}',
               lhablock = 'MASS',
               lhacode = [ 5 ])

MH = Parameter(name = 'MH',
               nature = 'external',
               type = 'real',
               value = 125,
               texname = '\\text{MH}',
               lhablock = 'MASS',
               lhacode = [ 25 ])

MN1 = Parameter(name = 'MN1',
                nature = 'external',
                type = 'real',
                value = 256.924777,
                texname = '\\text{MN1}',
                lhablock = 'MASS',
                lhacode = [ 1000022 ])

MN2 = Parameter(name = 'MN2',
                nature = 'external',
                type = 'real',
                value = 617.426187,
                texname = '\\text{MN2}',
                lhablock = 'MASS',
                lhacode = [ 1000023 ])

MN3 = Parameter(name = 'MN3',
                nature = 'external',
                type = 'real',
                value = 1638.2255,
                texname = '\\text{MN3}',
                lhablock = 'MASS',
                lhacode = [ 1000025 ])

MN4 = Parameter(name = 'MN4',
                nature = 'external',
                type = 'real',
                value = 1641.63702,
                texname = '\\text{MN4}',
                lhablock = 'MASS',
                lhacode = [ 1000035 ])

MC1 = Parameter(name = 'MC1',
                nature = 'external',
                type = 'real',
                value = 617.58599,
                texname = '\\text{MC1}',
                lhablock = 'MASS',
                lhacode = [ 1000024 ])

MC2 = Parameter(name = 'MC2',
                nature = 'external',
                type = 'real',
                value = 1642.09149,
                texname = '\\text{MC2}',
                lhablock = 'MASS',
                lhacode = [ 1000037 ])

WZ = Parameter(name = 'WZ',
               nature = 'external',
               type = 'real',
               value = 2.4952,
               texname = '\\text{WZ}',
               lhablock = 'DECAY',
               lhacode = [ 23 ])

WW = Parameter(name = 'WW',
               nature = 'external',
               type = 'real',
               value = 2.085,
               texname = '\\text{WW}',
               lhablock = 'DECAY',
               lhacode = [ 24 ])

WT = Parameter(name = 'WT',
               nature = 'external',
               type = 'real',
               value = 1.50833649,
               texname = '\\text{WT}',
               lhablock = 'DECAY',
               lhacode = [ 6 ])

WH = Parameter(name = 'WH',
               nature = 'external',
               type = 'real',
               value = 0.00407,
               texname = '\\text{WH}',
               lhablock = 'DECAY',
               lhacode = [ 25 ])

WN1 = Parameter(name = 'WN1',
                nature = 'external',
                type = 'real',
                value = 0.,
                texname = '\\text{WN1}',
                lhablock = 'DECAY',
                lhacode = [ 1000022 ])

WN2 = Parameter(name = 'WN2',
                nature = 'external',
                type = 'real',
                value = 6.1742618700000005,
                texname = '\\text{WN2}',
                lhablock = 'DECAY',
                lhacode = [ 1000023 ])

WN3 = Parameter(name = 'WN3',
                nature = 'external',
                type = 'real',
                value = 16.382255,
                texname = '\\text{WN3}',
                lhablock = 'DECAY',
                lhacode = [ 1000025 ])

WN4 = Parameter(name = 'WN4',
                nature = 'external',
                type = 'real',
                value = 16.4163702,
                texname = '\\text{WN4}',
                lhablock = 'DECAY',
                lhacode = [ 1000035 ])

WC1 = Parameter(name = 'WC1',
                nature = 'external',
                type = 'real',
                value = 6.175859900000001,
                texname = '\\text{WC1}',
                lhablock = 'DECAY',
                lhacode = [ 1000024 ])

WC2 = Parameter(name = 'WC2',
                nature = 'external',
                type = 'real',
                value = 16.4209149,
                texname = '\\text{WC2}',
                lhablock = 'DECAY',
                lhacode = [ 1000037 ])

aEW = Parameter(name = 'aEW',
                nature = 'internal',
                type = 'real',
                value = '1/aEWM1',
                texname = '\\alpha _{\\text{EW}}')

G = Parameter(name = 'G',
              nature = 'internal',
              type = 'real',
              value = '2*cmath.sqrt(aS)*cmath.sqrt(cmath.pi)',
              texname = 'G')

CKM1x1 = Parameter(name = 'CKM1x1',
                   nature = 'internal',
                   type = 'complex',
                   value = 'cmath.cos(cabi)',
                   texname = '\\text{CKM1x1}')

CKM1x2 = Parameter(name = 'CKM1x2',
                   nature = 'internal',
                   type = 'complex',
                   value = 'cmath.sin(cabi)',
                   texname = '\\text{CKM1x2}')

CKM1x3 = Parameter(name = 'CKM1x3',
                   nature = 'internal',
                   type = 'complex',
                   value = '0',
                   texname = '\\text{CKM1x3}')

CKM2x1 = Parameter(name = 'CKM2x1',
                   nature = 'internal',
                   type = 'complex',
                   value = '-cmath.sin(cabi)',
                   texname = '\\text{CKM2x1}')

CKM2x2 = Parameter(name = 'CKM2x2',
                   nature = 'internal',
                   type = 'complex',
                   value = 'cmath.cos(cabi)',
                   texname = '\\text{CKM2x2}')

CKM2x3 = Parameter(name = 'CKM2x3',
                   nature = 'internal',
                   type = 'complex',
                   value = '0',
                   texname = '\\text{CKM2x3}')

CKM3x1 = Parameter(name = 'CKM3x1',
                   nature = 'internal',
                   type = 'complex',
                   value = '0',
                   texname = '\\text{CKM3x1}')

CKM3x2 = Parameter(name = 'CKM3x2',
                   nature = 'internal',
                   type = 'complex',
                   value = '0',
                   texname = '\\text{CKM3x2}')

CKM3x3 = Parameter(name = 'CKM3x3',
                   nature = 'internal',
                   type = 'complex',
                   value = '1',
                   texname = '\\text{CKM3x3}')

yChiI1Y0ChiI2Y3 = Parameter(name = 'yChiI1Y0ChiI2Y3',
                            nature = 'internal',
                            type = 'complex',
                            value = 'complex(0,1)*IyChiI1Y0ChiI2Y3 + RyChiI1Y0ChiI2Y3',
                            texname = '\\text{yChiI1Y0ChiI2Y3}')

yChiI1Y0ChiI2Ym3 = Parameter(name = 'yChiI1Y0ChiI2Ym3',
                             nature = 'internal',
                             type = 'complex',
                             value = 'complex(0,1)*IyChiI1Y0ChiI2Ym3 + RyChiI1Y0ChiI2Ym3',
                             texname = '\\text{yChiI1Y0ChiI2Ym3}')

yChiI2Y3ChiI3Y0 = Parameter(name = 'yChiI2Y3ChiI3Y0',
                            nature = 'internal',
                            type = 'complex',
                            value = 'complex(0,1)*IyChiI2Y3ChiI3Y0 + RyChiI2Y3ChiI3Y0',
                            texname = '\\text{yChiI2Y3ChiI3Y0}')

yChiI2Ym3ChiI3Y0 = Parameter(name = 'yChiI2Ym3ChiI3Y0',
                             nature = 'internal',
                             type = 'complex',
                             value = 'complex(0,1)*IyChiI2Ym3ChiI3Y0 + RyChiI2Ym3ChiI3Y0',
                             texname = '\\text{yChiI2Ym3ChiI3Y0}')

UMixDC1C21x1 = Parameter(name = 'UMixDC1C21x1',
                         nature = 'internal',
                         type = 'complex',
                         value = 'complex(0,1)*IUMixDC1C21x1 + RUMixDC1C21x1',
                         texname = '\\text{UMixDC1C21x1}')

UMixDC1C21x2 = Parameter(name = 'UMixDC1C21x2',
                         nature = 'internal',
                         type = 'complex',
                         value = 'complex(0,1)*IUMixDC1C21x2 + RUMixDC1C21x2',
                         texname = '\\text{UMixDC1C21x2}')

UMixDC1C22x1 = Parameter(name = 'UMixDC1C22x1',
                         nature = 'internal',
                         type = 'complex',
                         value = 'complex(0,1)*IUMixDC1C22x1 + RUMixDC1C22x1',
                         texname = '\\text{UMixDC1C22x1}')

UMixDC1C22x2 = Parameter(name = 'UMixDC1C22x2',
                         nature = 'internal',
                         type = 'complex',
                         value = 'complex(0,1)*IUMixDC1C22x2 + RUMixDC1C22x2',
                         texname = '\\text{UMixDC1C22x2}')

UMixMN1N2N3N41x1 = Parameter(name = 'UMixMN1N2N3N41x1',
                             nature = 'internal',
                             type = 'complex',
                             value = 'complex(0,1)*IUMixMN1N2N3N41x1 + RUMixMN1N2N3N41x1',
                             texname = '\\text{UMixMN1N2N3N41x1}')

UMixMN1N2N3N41x2 = Parameter(name = 'UMixMN1N2N3N41x2',
                             nature = 'internal',
                             type = 'complex',
                             value = 'complex(0,1)*IUMixMN1N2N3N41x2 + RUMixMN1N2N3N41x2',
                             texname = '\\text{UMixMN1N2N3N41x2}')

UMixMN1N2N3N41x3 = Parameter(name = 'UMixMN1N2N3N41x3',
                             nature = 'internal',
                             type = 'complex',
                             value = 'complex(0,1)*IUMixMN1N2N3N41x3 + RUMixMN1N2N3N41x3',
                             texname = '\\text{UMixMN1N2N3N41x3}')

UMixMN1N2N3N41x4 = Parameter(name = 'UMixMN1N2N3N41x4',
                             nature = 'internal',
                             type = 'complex',
                             value = 'complex(0,1)*IUMixMN1N2N3N41x4 + RUMixMN1N2N3N41x4',
                             texname = '\\text{UMixMN1N2N3N41x4}')

UMixMN1N2N3N42x1 = Parameter(name = 'UMixMN1N2N3N42x1',
                             nature = 'internal',
                             type = 'complex',
                             value = 'complex(0,1)*IUMixMN1N2N3N42x1 + RUMixMN1N2N3N42x1',
                             texname = '\\text{UMixMN1N2N3N42x1}')

UMixMN1N2N3N42x2 = Parameter(name = 'UMixMN1N2N3N42x2',
                             nature = 'internal',
                             type = 'complex',
                             value = 'complex(0,1)*IUMixMN1N2N3N42x2 + RUMixMN1N2N3N42x2',
                             texname = '\\text{UMixMN1N2N3N42x2}')

UMixMN1N2N3N42x3 = Parameter(name = 'UMixMN1N2N3N42x3',
                             nature = 'internal',
                             type = 'complex',
                             value = 'complex(0,1)*IUMixMN1N2N3N42x3 + RUMixMN1N2N3N42x3',
                             texname = '\\text{UMixMN1N2N3N42x3}')

UMixMN1N2N3N42x4 = Parameter(name = 'UMixMN1N2N3N42x4',
                             nature = 'internal',
                             type = 'complex',
                             value = 'complex(0,1)*IUMixMN1N2N3N42x4 + RUMixMN1N2N3N42x4',
                             texname = '\\text{UMixMN1N2N3N42x4}')

UMixMN1N2N3N43x1 = Parameter(name = 'UMixMN1N2N3N43x1',
                             nature = 'internal',
                             type = 'complex',
                             value = 'complex(0,1)*IUMixMN1N2N3N43x1 + RUMixMN1N2N3N43x1',
                             texname = '\\text{UMixMN1N2N3N43x1}')

UMixMN1N2N3N43x2 = Parameter(name = 'UMixMN1N2N3N43x2',
                             nature = 'internal',
                             type = 'complex',
                             value = 'complex(0,1)*IUMixMN1N2N3N43x2 + RUMixMN1N2N3N43x2',
                             texname = '\\text{UMixMN1N2N3N43x2}')

UMixMN1N2N3N43x3 = Parameter(name = 'UMixMN1N2N3N43x3',
                             nature = 'internal',
                             type = 'complex',
                             value = 'complex(0,1)*IUMixMN1N2N3N43x3 + RUMixMN1N2N3N43x3',
                             texname = '\\text{UMixMN1N2N3N43x3}')

UMixMN1N2N3N43x4 = Parameter(name = 'UMixMN1N2N3N43x4',
                             nature = 'internal',
                             type = 'complex',
                             value = 'complex(0,1)*IUMixMN1N2N3N43x4 + RUMixMN1N2N3N43x4',
                             texname = '\\text{UMixMN1N2N3N43x4}')

UMixMN1N2N3N44x1 = Parameter(name = 'UMixMN1N2N3N44x1',
                             nature = 'internal',
                             type = 'complex',
                             value = 'complex(0,1)*IUMixMN1N2N3N44x1 + RUMixMN1N2N3N44x1',
                             texname = '\\text{UMixMN1N2N3N44x1}')

UMixMN1N2N3N44x2 = Parameter(name = 'UMixMN1N2N3N44x2',
                             nature = 'internal',
                             type = 'complex',
                             value = 'complex(0,1)*IUMixMN1N2N3N44x2 + RUMixMN1N2N3N44x2',
                             texname = '\\text{UMixMN1N2N3N44x2}')

UMixMN1N2N3N44x3 = Parameter(name = 'UMixMN1N2N3N44x3',
                             nature = 'internal',
                             type = 'complex',
                             value = 'complex(0,1)*IUMixMN1N2N3N44x3 + RUMixMN1N2N3N44x3',
                             texname = '\\text{UMixMN1N2N3N44x3}')

UMixMN1N2N3N44x4 = Parameter(name = 'UMixMN1N2N3N44x4',
                             nature = 'internal',
                             type = 'complex',
                             value = 'complex(0,1)*IUMixMN1N2N3N44x4 + RUMixMN1N2N3N44x4',
                             texname = '\\text{UMixMN1N2N3N44x4}')

VMixDC1C21x1 = Parameter(name = 'VMixDC1C21x1',
                         nature = 'internal',
                         type = 'complex',
                         value = 'complex(0,1)*IVMixDC1C21x1 + RVMixDC1C21x1',
                         texname = '\\text{VMixDC1C21x1}')

VMixDC1C21x2 = Parameter(name = 'VMixDC1C21x2',
                         nature = 'internal',
                         type = 'complex',
                         value = 'complex(0,1)*IVMixDC1C21x2 + RVMixDC1C21x2',
                         texname = '\\text{VMixDC1C21x2}')

VMixDC1C22x1 = Parameter(name = 'VMixDC1C22x1',
                         nature = 'internal',
                         type = 'complex',
                         value = 'complex(0,1)*IVMixDC1C22x1 + RVMixDC1C22x1',
                         texname = '\\text{VMixDC1C22x1}')

VMixDC1C22x2 = Parameter(name = 'VMixDC1C22x2',
                         nature = 'internal',
                         type = 'complex',
                         value = 'complex(0,1)*IVMixDC1C22x2 + RVMixDC1C22x2',
                         texname = '\\text{VMixDC1C22x2}')

MW = Parameter(name = 'MW',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(MZ**2/2. + cmath.sqrt(MZ**4/4. - (aEW*cmath.pi*MZ**2)/(Gf*cmath.sqrt(2))))',
               texname = 'M_W')

ee = Parameter(name = 'ee',
               nature = 'internal',
               type = 'real',
               value = '2*cmath.sqrt(aEW)*cmath.sqrt(cmath.pi)',
               texname = 'e')

sw2 = Parameter(name = 'sw2',
                nature = 'internal',
                type = 'real',
                value = '1 - MW**2/MZ**2',
                texname = '\\text{sw2}')

cw = Parameter(name = 'cw',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(1 - sw2)',
               texname = 'c_w')

sw = Parameter(name = 'sw',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(sw2)',
               texname = 's_w')

g1 = Parameter(name = 'g1',
               nature = 'internal',
               type = 'real',
               value = 'ee/cw',
               texname = 'g_1')

gw = Parameter(name = 'gw',
               nature = 'internal',
               type = 'real',
               value = 'ee/sw',
               texname = 'g_w')

vev = Parameter(name = 'vev',
                nature = 'internal',
                type = 'real',
                value = '(2*MW*sw)/ee',
                texname = '\\text{vev}')

lam = Parameter(name = 'lam',
                nature = 'internal',
                type = 'real',
                value = 'MH**2/(2.*vev**2)',
                texname = '\\text{lam}')

yb = Parameter(name = 'yb',
               nature = 'internal',
               type = 'real',
               value = '(ymb*cmath.sqrt(2))/vev',
               texname = '\\text{yb}')

yc = Parameter(name = 'yc',
               nature = 'internal',
               type = 'real',
               value = '(ymc*cmath.sqrt(2))/vev',
               texname = '\\text{yc}')

ydo = Parameter(name = 'ydo',
                nature = 'internal',
                type = 'real',
                value = '(ymdo*cmath.sqrt(2))/vev',
                texname = '\\text{ydo}')

ye = Parameter(name = 'ye',
               nature = 'internal',
               type = 'real',
               value = '(yme*cmath.sqrt(2))/vev',
               texname = '\\text{ye}')

ym = Parameter(name = 'ym',
               nature = 'internal',
               type = 'real',
               value = '(ymm*cmath.sqrt(2))/vev',
               texname = '\\text{ym}')

ys = Parameter(name = 'ys',
               nature = 'internal',
               type = 'real',
               value = '(yms*cmath.sqrt(2))/vev',
               texname = '\\text{ys}')

yt = Parameter(name = 'yt',
               nature = 'internal',
               type = 'real',
               value = '(ymt*cmath.sqrt(2))/vev',
               texname = '\\text{yt}')

ytau = Parameter(name = 'ytau',
                 nature = 'internal',
                 type = 'real',
                 value = '(ymtau*cmath.sqrt(2))/vev',
                 texname = '\\text{ytau}')

yup = Parameter(name = 'yup',
                nature = 'internal',
                type = 'real',
                value = '(ymup*cmath.sqrt(2))/vev',
                texname = '\\text{yup}')

muH = Parameter(name = 'muH',
                nature = 'internal',
                type = 'real',
                value = 'cmath.sqrt(lam*vev**2)',
                texname = '\\mu')

I1a11 = Parameter(name = 'I1a11',
                  nature = 'internal',
                  type = 'complex',
                  value = 'ydo*complexconjugate(CKM1x1)',
                  texname = '\\text{I1a11}')

I1a12 = Parameter(name = 'I1a12',
                  nature = 'internal',
                  type = 'complex',
                  value = 'ydo*complexconjugate(CKM2x1)',
                  texname = '\\text{I1a12}')

I1a13 = Parameter(name = 'I1a13',
                  nature = 'internal',
                  type = 'complex',
                  value = 'ydo*complexconjugate(CKM3x1)',
                  texname = '\\text{I1a13}')

I1a21 = Parameter(name = 'I1a21',
                  nature = 'internal',
                  type = 'complex',
                  value = 'ys*complexconjugate(CKM1x2)',
                  texname = '\\text{I1a21}')

I1a22 = Parameter(name = 'I1a22',
                  nature = 'internal',
                  type = 'complex',
                  value = 'ys*complexconjugate(CKM2x2)',
                  texname = '\\text{I1a22}')

I1a23 = Parameter(name = 'I1a23',
                  nature = 'internal',
                  type = 'complex',
                  value = 'ys*complexconjugate(CKM3x2)',
                  texname = '\\text{I1a23}')

I1a31 = Parameter(name = 'I1a31',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yb*complexconjugate(CKM1x3)',
                  texname = '\\text{I1a31}')

I1a32 = Parameter(name = 'I1a32',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yb*complexconjugate(CKM2x3)',
                  texname = '\\text{I1a32}')

I1a33 = Parameter(name = 'I1a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yb*complexconjugate(CKM3x3)',
                  texname = '\\text{I1a33}')

I2a11 = Parameter(name = 'I2a11',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yup*complexconjugate(CKM1x1)',
                  texname = '\\text{I2a11}')

I2a12 = Parameter(name = 'I2a12',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yc*complexconjugate(CKM2x1)',
                  texname = '\\text{I2a12}')

I2a13 = Parameter(name = 'I2a13',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yt*complexconjugate(CKM3x1)',
                  texname = '\\text{I2a13}')

I2a21 = Parameter(name = 'I2a21',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yup*complexconjugate(CKM1x2)',
                  texname = '\\text{I2a21}')

I2a22 = Parameter(name = 'I2a22',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yc*complexconjugate(CKM2x2)',
                  texname = '\\text{I2a22}')

I2a23 = Parameter(name = 'I2a23',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yt*complexconjugate(CKM3x2)',
                  texname = '\\text{I2a23}')

I2a31 = Parameter(name = 'I2a31',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yup*complexconjugate(CKM1x3)',
                  texname = '\\text{I2a31}')

I2a32 = Parameter(name = 'I2a32',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yc*complexconjugate(CKM2x3)',
                  texname = '\\text{I2a32}')

I2a33 = Parameter(name = 'I2a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yt*complexconjugate(CKM3x3)',
                  texname = '\\text{I2a33}')

I3a11 = Parameter(name = 'I3a11',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM1x1*yup',
                  texname = '\\text{I3a11}')

I3a12 = Parameter(name = 'I3a12',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM1x2*yup',
                  texname = '\\text{I3a12}')

I3a13 = Parameter(name = 'I3a13',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM1x3*yup',
                  texname = '\\text{I3a13}')

I3a21 = Parameter(name = 'I3a21',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM2x1*yc',
                  texname = '\\text{I3a21}')

I3a22 = Parameter(name = 'I3a22',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM2x2*yc',
                  texname = '\\text{I3a22}')

I3a23 = Parameter(name = 'I3a23',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM2x3*yc',
                  texname = '\\text{I3a23}')

I3a31 = Parameter(name = 'I3a31',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM3x1*yt',
                  texname = '\\text{I3a31}')

I3a32 = Parameter(name = 'I3a32',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM3x2*yt',
                  texname = '\\text{I3a32}')

I3a33 = Parameter(name = 'I3a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM3x3*yt',
                  texname = '\\text{I3a33}')

I4a11 = Parameter(name = 'I4a11',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM1x1*ydo',
                  texname = '\\text{I4a11}')

I4a12 = Parameter(name = 'I4a12',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM1x2*ys',
                  texname = '\\text{I4a12}')

I4a13 = Parameter(name = 'I4a13',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM1x3*yb',
                  texname = '\\text{I4a13}')

I4a21 = Parameter(name = 'I4a21',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM2x1*ydo',
                  texname = '\\text{I4a21}')

I4a22 = Parameter(name = 'I4a22',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM2x2*ys',
                  texname = '\\text{I4a22}')

I4a23 = Parameter(name = 'I4a23',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM2x3*yb',
                  texname = '\\text{I4a23}')

I4a31 = Parameter(name = 'I4a31',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM3x1*ydo',
                  texname = '\\text{I4a31}')

I4a32 = Parameter(name = 'I4a32',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM3x2*ys',
                  texname = '\\text{I4a32}')

I4a33 = Parameter(name = 'I4a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM3x3*yb',
                  texname = '\\text{I4a33}')

