'''ParticleVertex: a decay vertex of particles, the building block of a ParticleModel.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
'''

from uvbuilder.states.particle import Particle
from typing import List,Dict
from itertools import product
from uvbuilder.tools.reps import AllReps


class ParticleVertex(object):
    """
    Class describing a Particle vertex, which contains Particle objects
    """
    
    def __init__(
        self, 
        mother: Particle, 
        daughters: List[Particle],
        br: float = 0.0
        ):
        
        self.mother : Particle = mother
        self.daughters : List[Particle] = sorted(list(daughters[:]))
        self.br : float = br

        if not self.hasParticles():
            raise ValueError(f"Error creating ParticleVertex. Not all particles are Particle objects")
        
    def __str__(self) -> str:
        vStr = f"{self.mother} -> {','.join([str(d) for d in self.daughters])}"
        return vStr
    
    def __repr__(self) -> str:
        return str(self)
    
    def __eq__(self, other) -> bool:
        """
        Define a comparison between Particle vertices based
        on the particles appearing in the vertex.
        The daughter particles are assumed to be sorted.
        """

        if type(self) != type(other):
            return False
        elif len(self) != len(other):
            return False
        elif self.knownQuantumNumbers() != other.knownQuantumNumbers():
            return False
        elif self.mother != other.mother:
            return False
        elif self.daughters != other.daughters:
            return False
        elif self.br != other.br:
            return False
        return True
    
    def __ne__(self, other) -> bool:
        return not (self == other)

    def __len__(self):
        return len(self.daughters) + 1

    def hasParticles(self):

        onlyParticles = all(isinstance(ptc,Particle)
                            for ptc in [self.mother] + self.daughters)
        return onlyParticles
    
    @property
    def representation(self) -> Dict[str,AllReps]:
        """
        Returns a dictionary with the representations appearing in the vertex
        """

        dList = sorted(self.daughters, key = lambda d : 
                       list(d.representation.items()))
        repDict = {'mother' : self.mother.representation,
                   'daughters' : [d.representation 
                                         for d in dList]}
        
        return repDict

    def conjugate(self) -> "ParticleVertex":
        """
        Creates a new vertex with its particle conjugates.
        Replace the particles by their anti-particles.
        """
        
        if self.hasParticles():
            assert isinstance(self.mother,Particle)
            if self.mother.antiParticle is None:
                newMother = self.mother
            else:
                newMother = self.mother.antiParticle
            newDaughters = []
            for d in self.daughters:
                assert isinstance(d,Particle)
                if d.antiParticle is None:
                    newDaughter = d
                else:
                    newDaughter = d.antiParticle
                newDaughters.append(newDaughter)
        else:
            raise ValueError("Can not conjugate a vertex not containing Particle objects")
        
        newVertex = ParticleVertex(mother=newMother,
                                daughters=newDaughters,
                                br=self.br)
        
        return newVertex
    
    def getParticles(self) -> List[Particle]:
        """
        Returns all particles appearing in the vertex.
        """

        particleList = [self.mother] + self.daughters[:]

        return particleList
    
    def knownQuantumNumbers(self,qnumber="eCharge") -> int:
        """
        Function that maps a vertex to an integer representing how many 'known'
        numbers there are.

        """
        tally = 0
        for particle in self.getParticles():
            if hasattr(particle, qnumber):
                if getattr(particle,qnumber) is not None:
                    tally += 1

        return tally

    def copy(self) -> "ParticleVertex":
        """
        Creates a copy of the vertex with the same particles and branching ratio.
        """

        newMother = self.mother.copy()
        newDaughters = [d.copy() for d in self.daughters]
        newVertex = ParticleVertex(mother=newMother,
                                   daughters=newDaughters,
                                   br=self.br)
        
        return newVertex

    