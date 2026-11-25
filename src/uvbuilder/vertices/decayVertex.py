"""
.. module:: decayVertex
   :synopsis: Defines the DecayVertex class, which represents a vertex of final states (i.e. jet, l-,...)
              or particles (photon,W+,...)

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""


from uvbuilder.states.finalState import FinalState
from uvbuilder.states.particle import Particle
from uvbuilder.states.stateDict import StateDict, States
from uvbuilder.vertices.particleVertex import ParticleVertex
from itertools import product
from typing import List,Dict



class DecayVertex(object):
    """
    Class describing a decay vertex, which contains Particle or FinalState objects
    """
    
    def __init__(
        self, 
        mother: int | States, 
        daughters: List[int | States],
        pdgDict: StateDict = StateDict(),
        br: float = 0.0,
        ):
        
        if isinstance(mother,int):
            if mother not in pdgDict:
                raise ValueError(f"Error creating DecayVertex. {mother} not found in {pdgDict}")
            ptc = pdgDict[mother]
            if not isinstance(ptc,FinalState):
                raise ValueError(f"Trying to create DecayVertex, but {ptc} with type {type(ptc)} is not allowed.")                
            mother = ptc
        
        daughterList = []
        for d in daughters:
            if isinstance(d,int):
                if d not in pdgDict:
                    raise ValueError(f"Error creating DecayVertex. {d} not found in {pdgDict}")
                ptc = pdgDict[d]
                if not isinstance(ptc,FinalState):
                    raise ValueError(f"Trying to create DecayVertex, but {ptc} with type {type(ptc)} is not allowed.")                 
                daughterList.append(ptc)
            else:
                daughterList.append(d)

        assert isinstance(mother,FinalState)
        self.mother : FinalState = mother
        self.daughters : List[FinalState] = sorted(daughterList)
        self.br : float = br

        if not self.hasFinalStates():
            raise ValueError(f"Error creating DecayVertex. Not all particles are Particle or FinalState objects")

    def __str__(self) -> str:
        vStr = f"{self.mother} -> {','.join([str(d) for d in self.daughters])}"
        return vStr
    
    def __repr__(self) -> str:
        return str(self)
    
    def __eq__(self, other) -> bool:
        
        if type(self) != type(other):
            return False
        if self.mother != other.mother:
            return False
        if self.daughters != other.daughters:
            return False
        if self.br != other.br:
            return False
        return True

    def __ne__(self, other) -> bool:
        return not (self == other)
    
    def __gt__(self, other) -> bool:
        
        if self.mother != other.mother:
            return (self.mother > other.mother)
        if self.daughters != other.daughters:
            return (self.daughters > other.daughters)
        return self.br > other.br

    def __lt__(self, other) -> bool:
        return not (self > other)
    
    def __len__(self):
        return len(self.daughters) + 1
    
    @classmethod
    def fromDict(cls,vertexDict : Dict, pdgDict : StateDict = StateDict()) -> List["DecayVertex"]:
        '''One vertex per decay of {mother: {(daughters): BR}}, the PDGs resolved with
        pdgDict.'''
        vertexList = []
        for mom, decays in vertexDict.items():
            for daughters,br in decays.items():
                vertex = DecayVertex(mother=mom,
                                     daughters=daughters,
                                     br=br,
                                     pdgDict=pdgDict)
                vertexList.append(vertex)
        
        return vertexList
   
    def copy(self) -> "DecayVertex":
        """
        Creates a copy of the vertex with the same particles and branching ratio.
        """

        if isinstance(self.mother,int):
            newMother = self.mother
        else:
            newMother = self.mother.copy()
        newDaughters = []
        for d in self.daughters:
            if isinstance(d,int):
                newDaughters.append(d)
            else:
                newDaughters.append(d.copy())
        newVertex = DecayVertex(mother=newMother,
                                daughters=newDaughters,
                                br=self.br)
        
        return newVertex

    def hasFinalStates(self):
        '''True if the mother and the daughters are all FinalState objects.'''
        onlyParticles = all(isinstance(ptc,FinalState) 
                            for ptc in [self.mother] + self.daughters)
        return onlyParticles
    
    def getFinalStates(self) -> List[FinalState]:
        """
        Get a list of all final states appearing in the vertex.
        Does not expand final states into their particle content.
        """

        finalStates = [self.mother] 
        for fs in self.daughters:
            if any(fs is f for f in finalStates):
                continue
            finalStates.append(fs)

        return sorted(finalStates)
    
    def getParticles(self) -> List[Particle]:
        """
        Returns all particles appearing in the vertex. If the vertex contains
        FinalStates (i.e. inclusive particles such as 'jet'), include the particles
        belonging to the FinalState.
        """

        particleList = []
        if isinstance(self.mother,FinalState):
            particleList += self.mother.particleList[:]
        else:
            particleList.append(self.mother)
        for daughter in self.daughters:
            if isinstance(daughter,FinalState):
                particleList += daughter.particleList[:]
            else:
                particleList.append(daughter)

        return particleList
    
    def conjugate(self) -> "DecayVertex":
        """
        Creates a new vertex with its final state conjugates.
        Replaces the particles by their anti-particles.
        """
        
        if self.hasFinalStates():
            assert isinstance(self.mother,(FinalState,Particle))
            newMother = self.mother.getConjugate()
            newDaughters = []
            for d in self.daughters:
                assert isinstance(d,(Particle,FinalState))
                if d.selfConjugate:
                    newDaughter = d
                else:
                    newDaughter = d.getConjugate()
                newDaughters.append(newDaughter)
        else:
            raise ValueError("Can not conjugate a vertex not containing Particle objects")

        if newMother is None:
            raise ValueError(f"Trying to conjugate vertex, but conjugate for {self.mother} is not defined.")
        if any(d is None for d in newDaughters):
            raise ValueError(f"Trying to conjugate vertex, but conjugate for {self.daughters} is not defined.")

        
        newVertex = DecayVertex(mother=newMother,
                                daughters=newDaughters,
                                br=self.br)
        
        return newVertex

    def expandToParticleVertices(self, uniqueReps : bool = False) -> List[ParticleVertex]:
        """
        Expands the vertex into a list of ParticleVertex objects containing only Particle objects.
        All possible particle combinations from the FinalStates are considered.
        If uniqueReps is True, only particles with distinct representations are considered
        when expanding the vertex.
        """

        vertexList = []
        for mother in self.mother.getParticles(uniqueReps=uniqueReps):
            for dauList in product(*[d.getParticles(uniqueReps=uniqueReps) for d in self.daughters]):
                vertexList.append(ParticleVertex(mother=mother, daughters=list(dauList), br=self.br))

        return vertexList