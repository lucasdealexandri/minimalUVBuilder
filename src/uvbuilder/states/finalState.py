"""
.. module:: finalState
   :synopsis: Defines the FinalState class, which represents a signature/inclusive particle (i.e. jet, lepton, neutrino, ...)

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

from typing import List,Union
from uvbuilder.states.particle import Particle
from uvbuilder.tools.reps import AllReps


class FinalState(object):
    """
    An instance of this class represents a final state signature.

    The properties are: label, particleList
    """
    __slots__ = ["label", "particleList", "selfConjugate", "antiFinalState", "pdg", "isSM"]
    


    def __init__(self, label : str, particleList : List[Particle], 
                 selfConjugate : Union[None,bool] = None):
        """
        Creates a final state objects, which is simply a wrapper for a list of particles.

        Possible properties for arguments.
        label: str, e.g. 'jet'
        particleList: list of particle objects corresponding to the final state, 
                      e.g. if the final state is a jet particleList = [u, d, s, c, b, g].
        selfConjugate: If None, it will be defined based on the particles in the list. ]
                       If True, the anti-final state is defined to be its own anti-finalstate.
                       If False, the anti-final state is created to be a new final state with all particles 
                       replaced by their anti-particles.
        """

        # Initialize attributes
        self.label : str = label
        self.particleList : List[Particle] = sorted(particleList)

        # Define the final state to be its own anti-particle if all particles in the list 
        # are self-conjugate or if the anti-particle of each particle is also in the list
        if selfConjugate is None:
            self.selfConjugate = True
            for ptc in self.particleList:
                if ptc.selfConjugate:
                    continue
                elif ptc.antiParticle is None:
                    continue
                if ptc.antiParticle not in self.particleList:
                    self.selfConjugate = False
                    break
        else:
            self.selfConjugate = selfConjugate

        self.isSM : bool = all(p.isSM for p in self.particleList)
        self.pdg : List[int] = self.getPDGs()
        if self.selfConjugate:
            self.antiFinalState : FinalState = self
        else:
            particleListC = [p.antiParticle for p in self.particleList 
                             if p.antiParticle is not None]
            # Set selfConjugate to True for the anti-final state to avoid a recursive loop.
            fsC = FinalState(label = self.conjugate_label(),
                             particleList=particleListC,
                             selfConjugate=True)
            self.antiFinalState = fsC
            fsC.antiFinalState = self
            # Reset to correct value
            fsC.selfConjugate = False

    
    def __str__(self) -> str:
        return self.label

    def __repr__(self) -> str:
        return self.__str__()
    
    def __hash__(self) -> int:
        return id(self)
    
    def __eq__(self, other) -> bool:
        """
        Compare particles according to their particles.
        """

        if type(self) != type(other):
            return False
        elif len(self.particleList) != len(other.particleList):
            return False
        else:
            return self.particleList == other.particleList
    
    def __neq__(self, other) -> bool:
        return not (self == other)
    
    def __lt__(self,other) -> bool:
        if isinstance(other,Particle):
            return max(self.pdg) < other.pdg
        else:
            return self.pdg < other.pdg

    def __gt__(self,other) -> bool:
        if isinstance(other,Particle):
            return max(self.pdg) > other.pdg
        else:
            return self.pdg > other.pdg
        
    def pretty_print(self) -> str:
        """
        Converts the final state to a more readable string, 
        listing the particles it contains.
        """

        pStr = [str(p) for p in self.particleList]
        fsStr = f'{self.label}({",".join(sorted(pStr))})'

        return fsStr
    
    def conjugate_label(self) -> str:
        """
        Returns the label of the conjugate final state.
        """

        if self.selfConjugate:
            return self.label
        if self.label.endswith("+"):
            return self.label[:-1] + "-"
        elif self.label.endswith("-"):
            return self.label[:-1] + "+"
        elif self.label.endswith("bar"):
            return self.label[:-3]
        elif self.label.endswith("~"):
            return self.label[:-1]
        else:            
            return self.label + "~"
    
    def copy(self) -> "FinalState":
        """
        Create a copy of the particle.
        """

        particleList_copy = [ptc.copy() for ptc in self.particleList]

        self_copy = FinalState(label=self.label,
                                particleList=particleList_copy)
        
        
        return self_copy

    def getPDGs(self) -> List[int]:

        pdgList = sorted([ptc.pdg for ptc in self.particleList])

        return pdgList

    def getParticles(self, uniqueReps=False):
        """
        Returns a list of particles.
        If uniqueReps = True, it will return only one particle 
        for each distinct representation (when all quantum numbers are defined).
        """

        if not uniqueReps:
            return self.particleList[:]
        
        uniqueReps = []
        particleList = []
        for ptc in self.particleList:
            rep = ptc.representation
            # Always include particles for which representations have not been completely defined yet
            if any(not isinstance(r,AllReps) for r in rep.values()):
                particleList.append(ptc)
            else:
                rep = frozenset(rep.items())
                if rep in uniqueReps:
                    continue
                uniqueReps.append(rep)
                particleList.append(ptc)

        return particleList


    

