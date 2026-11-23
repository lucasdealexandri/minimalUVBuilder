"""
.. module:: SMparticleDefinitions
   :synopsis: Defines the SM particles.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from uvbuilder.states.physicsUnits import MeV, GeV
from uvbuilder.states.particle import Particle
from uvbuilder.states.stateDict import StateDict
from uvbuilder.tools.reps import DIRAC_SPINOR, SCALAR, VECTOR, SU3Rep, U1Rep, UnbrokenRep
from uvbuilder.tools.reps import GlobalU1s
from uvbuilder.states.SMfields import lL1, lL2, lL3, eR, muR, tauR, QL1, QL2, QL3, uR, cR, tR, dR, sR, bR, B, Wv, Higgs, G

# SM particles
#--------------------------------
#  Charged leptons:
em = Particle(
    isSM=True,
    label="e-",
    pdg=11,
    mass=0.5 * MeV,
    rep=UnbrokenRep(eCharge=U1Rep(-1), colordim=SU3Rep(0, 0), spin=DIRAC_SPINOR,
                    globalU1s=GlobalU1s(L=1)),  # Lepton number L=1
    fields=[lL1[U1Rep(-1,2)],eR],
)
muon = Particle(
    isSM=True,
    label="mu-",
    pdg=13,
    mass=106 * MeV,
    rep=UnbrokenRep(eCharge=U1Rep(-1), colordim=SU3Rep(0, 0), spin=DIRAC_SPINOR,
                    globalU1s=GlobalU1s(L=1)),
    fields=[lL2[U1Rep(-1,2)],muR],
) 

tau = Particle(
    isSM=True,
    label="tau-",
    pdg=15,
    mass=1.77 * GeV,
    rep=UnbrokenRep(eCharge=U1Rep(-1), colordim=SU3Rep(0, 0), spin=DIRAC_SPINOR,
                    globalU1s=GlobalU1s(L=1)),
    fields=[lL3[U1Rep(-1,2)],tauR],
) 

# Neutrinos
nue = Particle(
    isSM=True,
    label="nue",
    pdg=12,
    mass=0.0 * MeV,
    rep=UnbrokenRep(eCharge=U1Rep(0), colordim=SU3Rep(0, 0), spin=DIRAC_SPINOR,
                    globalU1s=GlobalU1s(L=1)), # L=1 keeps nu distinct from its anti-particle
    fields=[lL1[U1Rep(1,2)]],
)
num = Particle(
    isSM=True,
    label="numu",
    pdg=14,
    mass=0.0 * MeV,
    rep=UnbrokenRep(eCharge=U1Rep(0), colordim=SU3Rep(0, 0), spin=DIRAC_SPINOR,
                    globalU1s=GlobalU1s(L=1)),
    fields=[lL2[U1Rep(1,2)]],
)
nutau = Particle(
    isSM=True,
    label="nutau",
    pdg=16,
    mass=0.0 * MeV,
    rep=UnbrokenRep(eCharge=U1Rep(0), colordim=SU3Rep(0, 0), spin=DIRAC_SPINOR,
                    globalU1s=GlobalU1s(L=1)),
    fields=[lL3[U1Rep(1,2)]],
)


# -----------------------------------------------------------------------------------------------------------------------------

# Light quarks:
d = Particle(
    isSM=True,
    label="d",
    pdg=1,
    mass=0.0 * MeV,
    rep=UnbrokenRep(eCharge=U1Rep(-1,3), colordim=SU3Rep(1, 0), spin=DIRAC_SPINOR,
                    globalU1s=GlobalU1s(B=U1Rep(1,3))),  # Baryon number B=1/3
    fields=[QL1[U1Rep(-1,2)],dR],
)

s = Particle(
    isSM=True,
    label="s",
    pdg=3,
    mass=0.0 * MeV,
    rep=UnbrokenRep(eCharge=U1Rep(-1,3), colordim=SU3Rep(1, 0), spin=DIRAC_SPINOR,
                    globalU1s=GlobalU1s(B=U1Rep(1,3))),  # Baryon number B=1/3
    fields=[QL2[U1Rep(-1,2)],sR],
)

b = Particle(
    isSM=True,
    label="b",
    pdg=5,
    mass=4.7 * GeV,
    rep=UnbrokenRep(eCharge=U1Rep(-1,3), colordim=SU3Rep(1, 0), spin=DIRAC_SPINOR,
                    globalU1s=GlobalU1s(B=U1Rep(1,3))),  # Baryon number B=1/3
    fields=[QL3[U1Rep(-1,2)],bR],
)


u = Particle(
    isSM=True,
    label="u",
    pdg=2,
    mass=0.0 * MeV,
    rep=UnbrokenRep(eCharge=U1Rep(2,3), colordim=SU3Rep(1, 0), spin=DIRAC_SPINOR,
                    globalU1s=GlobalU1s(B=U1Rep(1,3))),  # Baryon number B=1/3
    fields=[QL1[U1Rep(1,2)],uR],
)

c = Particle(
    isSM=True,
    label="c",
    pdg=4,
    mass=0.0 * MeV,
    rep=UnbrokenRep(eCharge=U1Rep(2,3), colordim=SU3Rep(1, 0), spin=DIRAC_SPINOR,
                    globalU1s=GlobalU1s(B=U1Rep(1,3))),  # Baryon number B=1/3
    fields=[QL2[U1Rep(1,2)],cR],
)

t = Particle(
    isSM=True,
    label="t",
    pdg=6,
    mass=172 * GeV,
    rep=UnbrokenRep(eCharge=U1Rep(2,3), colordim=SU3Rep(1, 0), spin=DIRAC_SPINOR,
                    globalU1s=GlobalU1s(B=U1Rep(1,3))),  # Baryon number B=1/3
    fields=[QL3[U1Rep(1,2)],tR],
)


# -----------------------------------------------------------------------------------------------------------------------------

# Gauge bosons:
g = Particle(
    isSM=True,
    label="g",
    pdg=21,
    mass=0.0 * MeV,
    rep=UnbrokenRep(eCharge=U1Rep(0), colordim=SU3Rep(1, 1), spin=VECTOR),
    totalwidth=0.0 * GeV,
    fields=[G],
)
photon = Particle(
    isSM=True,
    label="photon",
    pdg=22,
    mass=0.0 * MeV,
    rep=UnbrokenRep(eCharge=U1Rep(0), colordim=SU3Rep(0, 0), spin=VECTOR),
    totalwidth=0.0 * GeV,
    fields=[B, Wv[0]],
)
Z = Particle(
    isSM=True,
    label="Z",
    pdg=23,
    mass=91.0 * GeV,
    rep=UnbrokenRep(eCharge=U1Rep(0), colordim=SU3Rep(0, 0), spin=VECTOR),
    totalwidth=2.5 * GeV,
    fields=[Wv[0],B],
)

W = Particle(
    isSM=True,
    label="W+",
    pdg=24,
    mass=80.0 * GeV,
    rep=UnbrokenRep(eCharge=U1Rep(1), colordim=SU3Rep(0, 0), spin=VECTOR),
    totalwidth=2.0 * GeV,
    fields=[Wv[U1Rep(1)]], # W+ is simply the W+ field
)
# -----------------------------------------------------------------------------------------------------------------------------

# Higgs:
higgs = Particle(
    isSM=True,
    label="higgs",
    pdg=25,
    mass=125.0 * GeV,
    rep=UnbrokenRep(eCharge=U1Rep(0), colordim=SU3Rep(0, 0), spin=SCALAR),
    totalwidth=0.0 * GeV,
    fields = [Higgs],
)

leptons = [em, muon, tau]
gauge = [g, photon, Z, W]

jetList = [u, c, t,
           d, s, b]

nuList = [nue, num, nutau]


SMparticles = leptons + nuList + gauge + jetList + [higgs]
for ptc in SMparticles[:]:
    if ptc.antiParticle is not None:
        SMparticles.append(ptc.antiParticle)
# SMList = sorted(list(set(SMList)), key = lambda ptc: ptc.pdg)

for ptc in SMparticles:
    assert isinstance(ptc.eCharge, U1Rep)
    assert isinstance(ptc.colordim, SU3Rep)

# Store particles in StateDict, which only stores particles, but
# allows to retrieve the conjugate of a particle if needed.
SMparticlesDict = StateDict({ptc.pdg : ptc for ptc in SMparticles})