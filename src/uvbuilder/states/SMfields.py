"""
.. module:: SMparticleDefinitions
   :synopsis: Defines the SM particles.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from uvbuilder.states.field import  Field
from uvbuilder.states.multiplet import Multiplet
from uvbuilder.tools.reps import WEYL_LEFT, WEYL_RIGHT, SCALAR, VECTOR, SU3Rep, SU2Rep, U1Rep, BrokenRep
from uvbuilder.tools.reps import GlobalU1s

#  SM particles
# We only consider one generation, which is sufficient for determining
# the exotics representations

#  Lepton doublets:
lL1 = Multiplet(
    rep=BrokenRep(
        colordim=SU3Rep(0, 0),
        spin=WEYL_LEFT,
        isospin=SU2Rep(1,),
        hypercharge=U1Rep(-1,2),
        globalU1s=GlobalU1s(L=1),  # Lepton number L=1
    ),
    label='lL1',
    isSM=True,
    index=1,
    labels={U1Rep(-1,2) : "eL", U1Rep(1,2) : "nueL"},
    fillMultiplet=True,
)
lL2 = Multiplet(
    rep=BrokenRep(
        colordim=SU3Rep(0, 0),
        spin=WEYL_LEFT,
        isospin=SU2Rep(1,),
        hypercharge=U1Rep(-1,2),
        globalU1s=GlobalU1s(L=1),  # Lepton number L=1
    ),
    label='lL2',
    isSM=True,
    index=2,
    labels={U1Rep(-1,2) : "muL", U1Rep(1,2) : "numL"},
    fillMultiplet=True,
)

lL3 = Multiplet(
    rep=BrokenRep(
        colordim=SU3Rep(0, 0),
        spin=WEYL_LEFT,
        isospin=SU2Rep(1,),
        hypercharge=U1Rep(-1,2),
        globalU1s=GlobalU1s(L=1),  # Lepton number L=1
    ),
    label='lL3',
    isSM=True,
    index=3,
    labels={U1Rep(-1,2) : "tauL", U1Rep(1,2) : "nutauL"},
    fillMultiplet=True,
)

# Lepton singlets:
eR = Field(
    rep=BrokenRep(
        colordim=SU3Rep(0, 0),
        spin=WEYL_RIGHT,
        isospin=SU2Rep(0,),
        hypercharge=U1Rep(-1,1),
        i3=U1Rep(0),
        globalU1s=GlobalU1s(L=1),  # Lepton number L=1
    ),
    index=1,
    isSM=True,
    label="eR",
)


muR = Field(
    rep=BrokenRep(
        colordim=SU3Rep(0, 0),
        spin=WEYL_RIGHT,
        isospin=SU2Rep(0,),
        hypercharge=U1Rep(-1,1),
        i3=U1Rep(0),
        globalU1s=GlobalU1s(L=1),  # Lepton number L=1
    ),
    index=2,
    isSM=True,
    label="muR",
)


tauR = Field(
    rep=BrokenRep(
        colordim=SU3Rep(0, 0),
        spin=WEYL_RIGHT,
        isospin=SU2Rep(0,),
        hypercharge=U1Rep(-1,1),
        i3=U1Rep(0),
        globalU1s=GlobalU1s(L=1),  # Lepton number L=1
    ),
    index=3,
    isSM=True,
    label="tauR",
)

# -----------------------------------------------------------------------------------------------------------------------------

# Quark doublets:
QL1 = Multiplet(
    rep=BrokenRep(
        colordim=SU3Rep(1, 0),
        globalU1s=GlobalU1s(B=U1Rep(1,3)),  # Baryon number B=1/3
        spin=WEYL_LEFT,
        isospin=SU2Rep(1,),
        hypercharge=U1Rep(1,6),
    ),
    label='QL1',
    isSM=True,
    index=1,
    labels={U1Rep(1,2) : "uL", U1Rep(-1,2) : "dL"},
    fillMultiplet=True,
)

QL2 = Multiplet(
    rep=BrokenRep(
        colordim=SU3Rep(1, 0),
        globalU1s=GlobalU1s(B=U1Rep(1,3)),  # Baryon number B=1/3
        spin=WEYL_LEFT,
        isospin=SU2Rep(1,),
        hypercharge=U1Rep(1,6),
    ),
    label='QL2',
    isSM=True,
    index=2,
    labels={U1Rep(1,2) : "cL", U1Rep(-1,2) : "sL"},
    fillMultiplet=True,
)

QL3 = Multiplet(
    rep=BrokenRep(
        colordim=SU3Rep(1, 0),
        globalU1s=GlobalU1s(B=U1Rep(1,3)),  # Baryon number B=1/3
        spin=WEYL_LEFT,
        isospin=SU2Rep(1,),
        hypercharge=U1Rep(1,6),
    ),
    label='QL3',
    isSM=True,
    index=3,
    labels={U1Rep(1,2) : "tL", U1Rep(-1,2) : "bL"},
    fillMultiplet=True,
)

# Quark singlets
uR = Field(
    rep=BrokenRep(
        colordim=SU3Rep(1, 0),
        globalU1s=GlobalU1s(B=U1Rep(1,3)),  # Baryon number B=1/3
        spin=WEYL_RIGHT,
        isospin=SU2Rep(0,),
        hypercharge=U1Rep(2,3),
        i3=U1Rep(0),
    ),
    isSM=True,
    index=1,
    label="uR",
)

dR = Field(
    rep=BrokenRep(
        colordim=SU3Rep(1, 0),
        globalU1s=GlobalU1s(B=U1Rep(1,3)),  # Baryon number B=1/3
        spin=WEYL_RIGHT,
        isospin=SU2Rep(0,),
        hypercharge=U1Rep(-1,3),
        i3=U1Rep(0),
    ),
    isSM=True,
    index=1,
    label="dR",
)

cR = Field(
    rep=BrokenRep(
        colordim=SU3Rep(1, 0),
        globalU1s=GlobalU1s(B=U1Rep(1,3)),  # Baryon number B=1/3
        spin=WEYL_RIGHT,
        isospin=SU2Rep(0,),
        hypercharge=U1Rep(2,3),
        i3=U1Rep(0),
    ),
    isSM=True,
    index=2,
    label="cR",
)

sR = Field(
    rep=BrokenRep(
        colordim=SU3Rep(1, 0),
        globalU1s=GlobalU1s(B=U1Rep(1,3)),  # Baryon number B=1/3
        spin=WEYL_RIGHT,
        isospin=SU2Rep(0,),
        hypercharge=U1Rep(-1,3),
        i3=U1Rep(0),
    ),
    isSM=True,
    index=2,
    label="sR",
)

tR = Field(
    rep=BrokenRep(
        colordim=SU3Rep(1, 0),
        globalU1s=GlobalU1s(B=U1Rep(1,3)),  # Baryon number B=1/3
        spin=WEYL_RIGHT,
        isospin=SU2Rep(0,),
        hypercharge=U1Rep(2,3),
        i3=U1Rep(0),
    ),
    isSM=True,
    index=3,
    label="tR",
)

bR = Field(
    rep=BrokenRep(
        colordim=SU3Rep(1, 0),
        globalU1s=GlobalU1s(B=U1Rep(1,3)),  # Baryon number B=1/3
        spin=WEYL_RIGHT,
        isospin=SU2Rep(0,),
        hypercharge=U1Rep(-1,3),
        i3=U1Rep(0),
    ),
    isSM=True,
    index=3,
    label="bR",
)

# -----------------------------------------------------------------------------------------------------------------------------

# Gauge bosons:
G = Field(
    rep=BrokenRep(
        colordim=SU3Rep(1, 1),
        isospin=SU2Rep(0,),
        spin=VECTOR,
        hypercharge=U1Rep(0,1),
        i3=U1Rep(0),
    ),
    isSM=True,
    label="g",
)
B = Field(
    rep=BrokenRep(
        colordim=SU3Rep(0, 0),
        spin=VECTOR,
        isospin=SU2Rep(0,),
        hypercharge=U1Rep(0,1),
        i3=U1Rep(0),
    ),
    isSM=True,
    label="B",
)

Wplus = Field(
    rep=BrokenRep(
        colordim=SU3Rep(0, 0),
        spin=VECTOR,
        isospin=SU2Rep(2,),
        hypercharge=U1Rep(0,1),
        i3=U1Rep(1),
    ),
    isSM=True,
    label="W+",
)
Wminus = Wplus.conjugateField  # i3=-1, auto-label "W-"
W3 = Field(
    rep=BrokenRep(
        colordim=SU3Rep(0, 0),
        spin=VECTOR,
        isospin=SU2Rep(2,),
        hypercharge=U1Rep(0,1),
        i3=U1Rep(0),
    ),
    isSM=True,
    label="W3",
)
Wv = Multiplet(
    rep=BrokenRep(
        colordim=SU3Rep(0, 0),
        spin=VECTOR,
        isospin=SU2Rep(2,),
        hypercharge=U1Rep(0,1),
    ),
    label='W',
    isSM=True,
    fields=[Wplus, W3, Wminus],
    fillMultiplet=True
)

# -----------------------------------------------------------------------------------------------------------------------------

# Higgs:

HiggsDoublet = Multiplet(
    rep=BrokenRep(
        colordim=SU3Rep(0,0),
        spin=SCALAR,
        isospin=SU2Rep(1,),
        hypercharge=U1Rep(1,2)
    ),
    label='Phi',
    labels={U1Rep(-1,2): 'Phi[2]', U1Rep(1,2): 'Phi[1]'},
    isSM=True,
    fillMultiplet=True
)

Higgs = Field(
    rep=BrokenRep(
        colordim=SU3Rep(0, 0),
        spin=SCALAR,
        isospin=SU2Rep(1,),
        hypercharge=U1Rep(1,2),
        i3=-U1Rep(1/2),
    ),
    isSM=True,
    label="higgs",
)

LH_SM_DOUBLETS = [lL1, lL2, lL3, QL1, QL2, QL3]
LH_SM_SINGLETS = [f.conjugateField for f in (eR, muR, tauR, uR, dR, cR, sR, tR, bR)]
LH_SM_SINGLETMULTIPLETS = [
    Multiplet(
        rep=f.rep.multiplet_rep, 
        label=f.label, 
        labels={U1Rep(0): f.label}, 
        index=f.index,
        isSM=True,
        fillMultiplet=True) 
        for f in LH_SM_SINGLETS]

SM_MULTIPLETS = LH_SM_DOUBLETS + LH_SM_SINGLETMULTIPLETS + [HiggsDoublet, HiggsDoublet.conjugateMultiplet]
