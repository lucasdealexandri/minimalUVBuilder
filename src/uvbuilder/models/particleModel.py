#!/usr/bin/env python3


from pyexpat import model
from typing import List, Union, Dict
from dataclasses import replace
from uvbuilder.states.particle import Particle, family_key
from uvbuilder.states.stateDict import StateDict
from uvbuilder.vertices.particleVertex import ParticleVertex
from collections import Counter
import itertools
from uvbuilder.tools.reps import (WEYL_LEFT, WEYL_RIGHT, DiracSpinor, WeylSpinor,
                        Scalar, UnbrokenRep, U1Rep)
from uvbuilder.tools.components import WeylComponent, ScalarComponent
from uvbuilder.tools.fr_utils import (LEFT_SUFFIX, RIGHT_SUFFIX, chirality_suffix,
                            formatVariable, formatClassBlock, fr_quantum_numbers,
                            parameter_block, merge_chunks)
from uvbuilder.states.physicsUnits import GeV



class ParticleModel:

    __slots__ = ["vertices","pdgDict","_parent"]

    def __init__(self,vertices: List[ParticleVertex] | None = None):
        
        if vertices is None:
            vertices = []
        
        # Sort vertices according to PDGs
        sorted_vertices = sorted(vertices, key = lambda v: v.getParticles())
        self.vertices : List[ParticleVertex] = sorted_vertices
        self._parent : Union[None,"ParticleModel"] = None

        # Create dictionary with all particles appearing in vertices
        self.pdgDict : StateDict = StateDict()
        for v in vertices:
            particles = v.getParticles()
            for p in particles:
                self.pdgDict.setdefault(p.pdg,p)

    @property
    def bsm_particles(self) -> List[Particle]:
        return [p for p in self.getParticles() if not p.isSM]

    def to_feynrules(self, tracker) -> Dict[str, List[str]]:
        '''FeynRules chunk for the physical (mass eigenstate) sector.

        Groups the BSM particles by representation (same rep = one family);
        each group emits the Weyl component W classes and the F/S class with
        ``ClassMembers``, masses/widths and a generation index (``FGen``).
        The index ``tracker`` is duck-typed (class numbers, gauge/family
        indices).
        '''

        chunks = []
        particles = sorted(self.bsm_particles, key=lambda p: (p.rep, family_key(p)))
        for _, group in itertools.groupby(particles, key=lambda p: p.rep):
            chunks.append(self._physical_group_chunk(list(group), tracker))
        return merge_chunks(*chunks)

    def _physical_group_chunk(self, members: List[Particle],
                              tracker) -> Dict[str, List[str]]:

        first = members[0]
        majorana = first.selfConjugate
        rep = first.rep
        single = len(members) == 1
        container = str(members[0]) if single else self._container_name(members)
        family = None
        if not single:
            family = tracker.family(self._family_name(members), len(members))

        chunk: Dict[str, List[str]] = {'Classes': []}

        # Weyl component classes
        if isinstance(rep.spin, DiracSpinor):
            chiralities = [(WEYL_LEFT, LEFT_SUFFIX)]
            if not majorana: chiralities += [(WEYL_RIGHT, RIGHT_SUFFIX)]
        elif isinstance(rep.spin, WeylSpinor):
            chiralities = [(rep.spin, chirality_suffix(rep.spin))]
        else:
            chiralities = []
        for spin, suffix in chiralities:
            component = WeylComponent(f'{container}{suffix}',
                                      replace(rep, spin=spin))
            chunk['Classes'].extend(
                component.to_feynrules(tracker, family)['Classes'])

        # The 4-component / scalar class
        spin_label = formatVariable(rep.spin)
        header = f"{spin_label}[{tracker.class_number(spin_label)}]"
        blockDict = {
            'ClassName': container,
            'SelfConjugate': 'True' if majorana else 'False',
        }
        # Self-conjugate classes carry no quantum numbers
        if not majorana:
            blockDict['QuantumNumbers'] = fr_quantum_numbers({'Q': rep.eCharge}, rep.globalU1s)
        if isinstance(rep.spin, DiracSpinor):
            if majorana: blockDict['WeylComponents'] = f"{container}{LEFT_SUFFIX}"
            else:
                blockDict['WeylComponents'] = f"{{{container}{LEFT_SUFFIX}, {container}{RIGHT_SUFFIX}}}"
        if not single:
            blockDict['ClassMembers'] = '{' + ','.join(str(p) for p in members) + '}'
            blockDict['FlavorIndex'] = family.name
        indices = [str(index) for index in tracker.gauge(rep)]
        if family is not None:
            indices.append(str(family))
        if indices:
            blockDict['Indices'] = f"{{{','.join(indices)}}}"
        masses = [f"{{M{p.label}, {p.mass.asNumber(GeV)}}}"
                  for p in members if p.mass is not None]
        # W = 1% of mass for initial guess / prompt decays
        # Except for DM, where W = 0.
        dm_candidate = min(self.bsm_particles, key=lambda p: p.mass)
        widths = [f"{{W{p.label}, {0 * p.mass.asNumber(GeV)}}}" 
                  if p == dm_candidate else f"{{W{p.label}, {0.01 * p.mass.asNumber(GeV)}}}"
                  for p in members if p.mass is not None]
        # Single-member class: the particle's own {M, value}, no class-level symbol
        if masses:
            blockDict['Mass'] = masses[0] if single else '{' + f'M{container}, ' + ', '.join(masses) + '}'
        if widths:
            blockDict['Width'] = widths[0] if single else '{' + f'W{container}, ' + ', '.join(widths) + '}'
        blockDict['PDG'] = '{' + ', '.join(str(p.pdg) for p in members) + '}'
        chunk['Classes'].append(formatClassBlock(header, blockDict))

        return chunk

    @staticmethod
    def _container_name(members: List[Particle]) -> str:
        spin = members[0].rep.spin
        if isinstance(spin, DiracSpinor):
            name = 'Chi'
        elif isinstance(spin, Scalar):
            name = 'Phi'
        else:
            name = 'ChiL'
        return name + ParticleModel._gauge_name(members)

    @staticmethod
    def _gauge_name(members: List[Particle]) -> str:
        '''Tag of a physical class: charge x 3, color (b = conjugate), non-zero global
        charges, e.g. Q2C3 (Q = 2/3, color 3), Qm1C3b (Q = -1/3, color 3bar),
        Q0C1D1 (neutral, DarkCharge 1).'''
        from uvbuilder.models.multipletModel import _global_u1_tag   # local: avoids circular import
        rep = members[0].rep
        q3 = 3 * U1Rep(rep.eCharge)
        if q3.denominator != 1:
            raise ValueError(f'Electric charge {rep.eCharge} of {members[0]} is not a multiple of 1/3')
        n = q3.numerator
        charge = f'Q{n}' if n >= 0 else f'Qm{-n}'
        color = f'C{rep.colordim.dim}' + ('b' if rep.colordim.lower > rep.colordim.upper else '')
        return charge + color + _global_u1_tag(rep.globalU1s)

    @staticmethod
    def _family_name(members: List[Particle]) -> str:
        if isinstance(members[0].rep.spin, (DiracSpinor, WeylSpinor)):
            name = 'F'
        else:
            name = 'S'
        return name + 'Gen' + ParticleModel._gauge_name(members)

    @property
    def bsm_components(self) -> List[WeylComponent | ScalarComponent]:
        components = []

        for particle in self.bsm_particles:
            for component in particle.components:

                # Avoid adding duplicates
                if component in components: continue

                # CONVENTION: Only use LH Weyl components
                if component.rep.spin == WEYL_LEFT:
                    components.append(component)

                elif component.rep.spin == WEYL_RIGHT:
                    # Majorana case: just skip

                    # Dirac case: add if not already in.
                    if not component.rep.is_realInternalRep: 
                        if component.dagger not in components:
                            components.append(component.dagger)

                elif isinstance(component.rep.spin, Scalar):
                    components.append(component)

        return components

    @property
    def unbroken_reps(self) -> List[UnbrokenRep]:

        return [component.rep for component in self.bsm_components]

    def __str__(self):
        """
        Returns the particles and their quantum numbers (if assigned)
        """

        modelStr = "ParticleModel"
        labels = [ptc.label for ptc in self.pdgDict.values() if not ptc.isSM]

        return f"{modelStr}({','.join(labels)})"

    def __repr__(self):
        return str(self)

    def describe(self):
        """
        Returns a more detailed description of the model, 
        including the quantum numbers of the particles.
        """

        modelStr = "ParticleModel"
        ptc_description = [ptc.describe() for ptc in self.pdgDict.values() if not ptc.isSM]

        return f"{modelStr} ({', '.join(ptc_description)})"     

    def __eq__(self, other) -> bool:
        """
        Two models are defined as equal if they have the same vertices
        and the particles appearing in the vertices are equal.
        """
        
        if type(self) != type(other):
            return False
        if len(self.vertices) != len(other.vertices):
            return False
        
        for pdg,ptc in self.pdgDict.items():
            if pdg not in other.pdgDict:
                return False
            if ptc != other.pdgDict[pdg]:
                return False
        
        for v1,v2 in zip(self.vertices,other.vertices):
            if v1 != v2:
                return False
        
        return True
        

    def pretty_print(self) -> str:
        """
        Converts the model to a more readable string
        listing the particles and their representations under the SM
        gauge group.
        """

        modelStr = "ParticleModel"
        ptc_description = [ptc.describe() for ptc in self.pdgDict.values() if not ptc.isSM]
        
        return f"{modelStr}:\n  {'\n  '.join(ptc_description)}"   

    def copy(self) -> "ParticleModel":
        
        # Make a copy of all particles
        newDict = self.pdgDict.copy()
        newVertices = []
        # Create new vertices with the original particles replaced by their copies
        for vertex in self.vertices:
            new_mother = newDict[vertex.mother.pdg]
            new_daughters = [newDict[daughter.pdg] for daughter in vertex.daughters]
            assert isinstance(new_mother,Particle)
            newVertex = ParticleVertex(mother=new_mother,
                                       daughters=new_daughters,
                                       br=vertex.br)
            newVertices.append(newVertex)
        
        # Create a new ParticleModel, keeping the vertex order (explodeQuantityHelper
        # walks the vertices by index)
        newModel = ParticleModel()
        newModel.vertices = newVertices
        newModel.pdgDict = newDict
        if self._parent is None:  
            newModel._parent = self
        else:
            newModel._parent = self._parent

        return newModel


    def getModelSize(self) -> int:
        """
        Compute the size of a ParticleModel based on the number
        of defined attributes.
        
        """

        n = 0
        for ptc in self.pdgDict.values():
            n += len(ptc) # the length of a particle corresponds to how many qnumbers have been set

        return n
    
    def hasQNumbers(self, ignoreSM=True) -> List[str]:
        """
        Returns a list of particle qnumbers which have been defined for ALL particles.
        If ignoreSM = True, it will ignore the SM particles
        """

        qnumbers = None
        for ptc in self.pdgDict.values():
            if ptc.isSM and ignoreSM:
                continue
            if qnumbers is None:
                qnumbers = set(ptc.hasQNumbers())
            else:
                qnumbers = qnumbers & set(ptc.hasQNumbers())

        if qnumbers is None:
            return []
        else:
            return sorted(list(qnumbers))

    def diffAttrs(self, other : "ParticleModel") ->List[str]:
        """
        Return a list with the names of attributes which were
        defined in self and not in other and vice-versa.
        """
        
        m1_qnumbers = set(self.hasQNumbers())
        m2_qnumbers = set(other.hasQNumbers())
        diff_atts = (m1_qnumbers.symmetric_difference(m2_qnumbers))
        return sorted(list(diff_atts))
    
    def getParticles(self) -> List[Particle]:
        """
        Returns a list of all particles appearing in the model.
        Their anti-particles are not included in the list.
        """
        
        return sorted(list(self.pdgDict.values()))
        
    def sortVerticesByKnownQNumber(self, qnumber: str):
        """
        Sort the vertices according to the number of particles for which
        the quantum number have already been set.
        """
        
        sorted_vertices : List[ParticleVertex] = sorted(self.vertices[:],
                            key=lambda v: v.knownQuantumNumbers(qnumber=qnumber), 
                            reverse=True)
        self.vertices = sorted_vertices[:]


    def signature_order(self, rep):
        return (rep.eCharge, rep.globalU1s)

    @property
    def decay_signature(self) -> Counter:
        '''Decay signature is a Counter that has (p1.rep, p2.rep): count'''

        seen_vertices = set()
        signature = Counter()

        for vertex in self.vertices:

            bsm_particles = [p for p in vertex.getParticles() if not p.isSM]

            # Maybe throw an error? But works for 2HDM
            if len(bsm_particles) != 2:
                continue
            
            p1, p2 = bsm_particles

            if (p1, p2) in seen_vertices:
                continue

            seen_vertices.add((p1, p2))

            rep1, rep2 = sorted([p.components[0].rep for p in (p1, p2)], key=lambda rep: self.signature_order(rep))
            signature[(rep1, rep2)] += 1

            if p1.selfConjugate and p2.selfConjugate:
                continue
            else:
                signature[(rep2.conjInternalRep, rep1.conjInternalRep)] += 1

        return signature
