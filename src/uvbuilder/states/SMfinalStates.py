"""
.. module:: SMfinalStates
   :synopsis: Defines SM final states (inclusive particles).

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from uvbuilder.states.SMparticles import SMparticlesDict
from uvbuilder.states.finalState import FinalState
from uvbuilder.states.stateDict import StateDict

# SM particles
#--------------------------------
#  Charged leptons are not grouped: e, mu and tau lead to distinct final states in the detector


# Jet:
pdgs = [1,2,3,4,21,-1,-2,-3,-4]
jet = FinalState(label='jet',particleList=[SMparticlesDict[pdg] for pdg in pdgs])

# Neutrino
pdgs = [12,14,16,-12,-14,-16]
nu = FinalState(label='nu',particleList=[SMparticlesDict[pdg] for pdg in pdgs])

finalStates = [jet,nu]

# Use a StateDict to store the final states, which only stores particles/final states, 
# but allows to #retrieve the conjugate of a particle/final state if needed.
SMfinalStatesDict = StateDict()
for fs in finalStates:
    pdgList = fs.getPDGs()
    if any(pdg in SMfinalStatesDict for pdg in pdgList):
        raise ValueError(f"Some PDG in {fs} is already contained in another final state!")
    for pdg in pdgList:
        SMfinalStatesDict[pdg] = fs

# If a particle is not contained in a final state,
# create a dummy final state holding the particle (one-to-one mapping)
for pdg,ptc in SMparticlesDict.items():
    if pdg in SMfinalStatesDict:
        continue
    fs = FinalState(label=ptc.label,particleList=[ptc])
    SMfinalStatesDict[pdg] = fs