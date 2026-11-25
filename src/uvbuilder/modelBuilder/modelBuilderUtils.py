#!/usr/bin/env python3
'''Filling of the quantum numbers of a ParticleModel from its decay vertices: at each
vertex an unknown value is one of the products of the other two particles', and each
possibility gives a new model.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
'''

from typing import List, Sequence, Iterator, Union

from uvbuilder.states.particle import Particle
from uvbuilder.tools.reps import spin_prod, LorentzRep, SU3Rep, SU2Rep, U1Rep, AllReps, GlobalU1s, VECTOR
from uvbuilder.models.particleModel import ParticleModel
from uvbuilder.vertices.particleVertex import ParticleVertex


def productDispatch(quantity: str, particle_1: Particle, particle_2: Particle, inv: bool = False
) -> Sequence[AllReps | int]:
    """
    Aux. function called in the `explodeQuantity` method of `ParticleModel`.
    """
    if not hasattr(particle_1,quantity):
        raise ValueError(f"Trying to compute the product of {quantity}, but it has not been defined for {particle_1}")
    if not hasattr(particle_2,quantity):
        raise ValueError(f"Trying to compute the product of {quantity}, but it has not been defined for {particle_2}")
    x1 = getattr(particle_1,quantity)
    x2 = getattr(particle_2,quantity)
    if (x1 is None) or (x2 is None):
        raise ValueError(f"Trying to compute the product of {quantity} with {x1} and {x2} for particles {particle_1} and {particle_2}.")

    if isinstance(x1,int):
        result: int = x1+x2
        if inv:
            result = result*(-1)
        result : List[int] = [result]
    elif isinstance(x1,U1Rep):
        result : U1Rep = x1+x2
        if inv:
            result = result*(-1)
        result : List[U1Rep] = [result]
    elif isinstance(x1, GlobalU1s):
        result: GlobalU1s = x1 + x2
        if inv:
            result = -result
        result: List[GlobalU1s] = [result]
    elif isinstance(x1,SU2Rep):
        # No inv needed for SU(2)
        result: List[SU2Rep] = x1*x2
    elif isinstance(x1,SU3Rep):
        result: List[SU3Rep] = x1*x2
        if inv:
            result = [c.inv() for c in result[:]]
    elif isinstance(x1,LorentzRep):
        result: List[LorentzRep] = spin_prod(x1,x2)
    else:
        raise Exception(f"Can not compute a product of type {type(x1)}")
    
    return result


def explodeQuantity(model : ParticleModel, 
                    qnumber: str) -> Iterator[ParticleModel]:
    """
    Compute all possible values for qnumber and generates all ParticleModels
    for each possibility.
    """
    
    # Sort vertices according to the number of particles for which
    # the quantum number have already been set.
    model.sortVerticesByKnownQNumber(qnumber=qnumber)

    for newModel in explodeQuantityHelper(model,qnumber=qnumber):
        if newModel is not model:
            newModel._parent = model
            yield newModel
    
def explodeQuantityHelper(
    model : ParticleModel, qnumber: str, vertex_number: int = 0
) -> Iterator[ParticleModel]:
    """Returns a generator of ParticleModel objects with allowed quantum numbers for
    `quantity` filled in at each vertex.

    Iterates through the vertices one by one recursively, increments
    `vertex_number` at each step. Should be called only through
    `explodeQuantity` since this sorts the list (in reverse order) first.
    The remaining vertices are re-sorted at each step, so the next vertex is
    always one with the most known values.
    """

    # Recursive base case
    if vertex_number == len(model.vertices):
        yield model
    else:
        # Re-sort the remaining vertices: the next one has the most known values
        model.vertices[vertex_number:] = sorted(model.vertices[vertex_number:],
            key=lambda v: v.knownQuantumNumbers(qnumber=qnumber), reverse=True)
        # The vertex being exploded (or checked)
        vertex : ParticleVertex = model.vertices[vertex_number]

        # Are we setting a quantum number value for the first time or just
        # checking for consistency?
        is_check = vertex.knownQuantumNumbers(qnumber=qnumber) == len(vertex)

        # Get vertex particles so they all appear in the same side of the vertex
        # ( i.e. X -> Y + a --> [X,anti-Y,anti-a])
        particles = [vertex.mother]
        for d in vertex.daughters:
            if isinstance(d,Particle) and d.antiParticle is not None:
                particles.append(d.antiParticle)
            else:
                raise ValueError(f"Trying to compute the product of {qnumber} for vertex {vertex},"
                                 f" but particle {d} does not have a defined conjugate field/particle,"
                                 f" which is needed to put all particles on the same side of the vertex.")
        # Extract particle with unknown quantum number as p1
        if not is_check:
            particles = sorted(particles, key=lambda x: getattr(x, qnumber) is not None)
        if len(particles) != 3:
                raise ValueError(f"Can only deal with two body decays and not with {vertex}")
        p1,p2,p3 = particles
        
        # Make sure all particles have been properly defined
        assert p1 is not None
        assert p2 is not None
        assert p3 is not None
        assert hasattr(p2, qnumber)
        assert hasattr(p3, qnumber)

        # Possible assignments for p1 quantum number
        possibleQuantumNumbers = productDispatch(quantity=qnumber, 
                                                    particle_1=p2, particle_2=p3, 
                                                    inv=True)
        for quantumNumber in possibleQuantumNumbers:
            if not is_check:
                copy_ = model.copy()
            else:
                copy_  = model
            try:
                # Set attribute to p1 in new model:
                assert isinstance(p1,Particle)
                setattr(copy_.pdgDict[p1.pdg], qnumber, quantumNumber)
                
                for newModel in explodeQuantityHelper(copy_,
                                                       qnumber=qnumber, 
                                                       vertex_number=vertex_number + 1):
                    yield newModel
            except AttributeError:
                # Skip the possibility if it is inconsistent
                pass

def fillECharge(model : ParticleModel) -> Iterator[ParticleModel]:
    '''Models with the electric charges filled (explodeQuantity).'''
    return explodeQuantity(model, qnumber="eCharge")

def fillGlobalU1s(model: ParticleModel) -> Iterator[ParticleModel]:
    '''Models with the global charges filled (explodeQuantity).'''
    return explodeQuantity(model, qnumber="globalU1s")
    
def explodeColordim(model : ParticleModel) -> Iterator[ParticleModel]:
    '''Models with the color representations filled, one per possibility (explodeQuantity).'''
    return explodeQuantity(model,qnumber="colordim")

def explodeSpin(model : ParticleModel) -> Iterator[ParticleModel]:
    '''Models with the spins filled, one per possibility (explodeQuantity).'''
    return explodeQuantity(model, qnumber="spin")


