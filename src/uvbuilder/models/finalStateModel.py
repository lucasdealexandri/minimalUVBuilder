#!/usr/bin/env python3
"""
.. module:: finalStateModel
   :synopsis: Defines the FinalStateModel class, which holds the decay vertices of a protomodel.

.. moduleauthor:: John Gargalionis <johngargalionis@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

from typing import List, Dict, Union

from uvbuilder.states.particle import Particle
from uvbuilder.states.finalState import FinalState
from uvbuilder.states.stateDict import StateDict
from uvbuilder.states.SMfinalStates import SMfinalStatesDict as SMdict
from uvbuilder.states.physicsUnits import GeV


from uvbuilder.vertices.decayVertex import DecayVertex



class FinalStateModel(object):
    """
    Decay vertices of a protomodel and the states appearing in them (pdgDict: PDG ->
    Particle for the BSM states, FinalState for the SM ones).
    """

    __slots__ = ["vertices","pdgDict"]

    def __init__(self,vertices: List[DecayVertex] | None = None):
        
        if vertices is None:
            vertices = []
        
        self.vertices = vertices

        # Create dictionary with all particles/final states appearing in the vertices
        self.pdgDict : StateDict = StateDict()
        for v in vertices:
            for fs in v.getFinalStates():
                for pdg in fs.getPDGs():
                    self.pdgDict.setdefault(pdg, fs)
            

    @classmethod
    def fromDict(cls, modelDict : Dict, pdgDict :  StateDict = SMdict) -> Union[None,"FinalStateModel"]:

        """
        Use the dictionary to create a FinalStateModel.

        :param modelDict: Model dictionary with decays ({mother: {(daughters): BR}}) and,
                          optionally, masses ({pdg: GeV}) and labels ({pdg: label}; a BSM
                          state without a label is named X<pdg>). Decays of SM particles
                          are skipped; decays must be two-body (an off-shell decay is
                          written with its virtual boson, e.g. (1000022, 23)).
        :param pdgDict: Dictionary to convert PDGs to Particle/FinalState objects (SM states). If a PDG is not found, assumes it is a BSM particle.

        :return: FinalStateModel dictionary
        """

        # modelDict should at least contain decays
        if 'decays' not in modelDict or len(modelDict['decays']) == 0:
            raise ValueError(f"Error creating FinalStateModel from {modelDict}. No decays found")
        
        # SM decays carry no information on the BSM states
        decays = {mom: d for mom, d in modelDict['decays'].items() if abs(mom) not in pdgDict}
        many = [f'{mom} -> {" ".join(map(str, d))}' for mom, ds in decays.items()
                for d in ds if len(d) != 2]
        if many:
            raise ValueError('only two-body decays are supported; write an off-shell decay '
                             'with its virtual boson, e.g. (1000022, 23) instead of '
                             '(1000022, 12, -12): ' + ', '.join(many))
        masses = modelDict.get('masses',{})
        xsections = modelDict.get('xsections',{})
        labels = modelDict.get('labels',{})

        # Compile a list of PDGs appearing in the modelDict:
        pdgList = set([])
        for pdg,decay in decays.items():
            pdgList.add(pdg)
            for daughter_pdgs in decay.keys():
                pdgList.update(daughter_pdgs)
        pdgList.update(masses.keys())        
        for prod_pdgs in xsections:
            pdgList.update(prod_pdgs)

        # Consider only particles (positive pdg) since 
        # they already contain the information about their anti-particles
        pdgList = set([abs(pdg) for pdg in list(pdgList)[:]])

        # Create a dictionary with all model final states;
        # If the pdg is not found in pdgDict, create a new final state containing
        # the BSM particle
        model_pdgDict : StateDict = StateDict()
        for  pdg in pdgList:
            if pdg in model_pdgDict:
                continue
            elif pdg in pdgDict:
                finalState = pdgDict[pdg].copy()
            else:
                bsm_ptc = Particle(pdg=pdg,
                                label=labels.get(pdg,f'X{pdg}'),
                                isSM=False)
                mass = masses.get(pdg,None)
                if mass is not None:
                    bsm_ptc.mass = mass*GeV
                finalState = FinalState(label=bsm_ptc.label, 
                                        particleList=[bsm_ptc])
            # Make sure all the final state pdgs point to the same object
            if not isinstance(finalState,FinalState):
                raise ValueError(f"Error creating FinalStateModel from {modelDict}. Object for pdg {pdg} is not a FinalState.")
            
            for ptc in finalState.getParticles():
                if ptc.isSM:
                    continue
                if ptc.pdg in masses:
                    ptc.mass = masses[ptc.pdg]*GeV
                if ptc.pdg in labels:
                    ptc.label = labels[ptc.pdg]
            for pdg in finalState.getPDGs():
                model_pdgDict[pdg] = finalState
        
        # Finally create vertices
        decayVertices = DecayVertex.fromDict(decays,pdgDict=model_pdgDict)
        
        # Create FinalStateModel:
        fsModel = FinalStateModel(vertices=decayVertices)

        # Early, warn-only name check (the FeynRules writer does the renaming)
        from uvbuilder.tools.fr_names import warn_labels, particle_symbols
        bsm_labels = sorted({p.label for fs in model_pdgDict.values()
                             for p in fs.getParticles() if not p.isSM and p.pdg > 0})
        warn_labels(bsm_labels, particle_symbols, 'Particle')

        return fsModel

    
    def __str__(self):
        """
        Simple string representation listing the final states.
        """

        # Get List of final states (ignore conjugate final states and SM final states)
        finalStates = []
        for fs in self.pdgDict.values():
            if fs.isSM:
                continue
            elif fs.selfConjugate:
                finalStates.append(fs)
            elif fs.antiFinalState not in finalStates:
                finalStates.append(fs)

        BSMstates = [str(fs) for fs in finalStates]
        bsmStr = ','.join(sorted(BSMstates))
        mstr = f'FinalStateModel: ({bsmStr})'

        return mstr
    
    def __repr__(self) -> str:
        return str(self)

    def pretty_print(self) -> str:
        """
        Converts the model to a more readable string
        listing the particles, their masses and their decays.
        """
        
        # Get list of final states (ignore conjugate final states and SM final states)        
        finalStates = []
        for fs in self.pdgDict.values():
            if fs.isSM:
                continue
            finalStates.append(fs)
        
        pStr = [str(fs) for fs in finalStates]
        # Construct particle labels for pretty printing
        vStr = [str(v) for v in self.vertices]
        
        modelStr = f'FinalStateModel: {','.join(sorted(pStr))}'
        mStr = [f'{fs}: {fs.particleList[0].mass}' for fs in sorted(finalStates, key=str)
                if fs.particleList[0].mass is not None]
        if mStr:
            modelStr += '\n  Masses: ' + ', '.join(mStr)
        if vStr:
            modelStr += '\n  Vertices: ' + ', '.join(sorted(vStr))

        return modelStr    
    
    def getFinalStates(self) -> List[FinalState | Particle]:
        """
        Get a list of all final states appearing in the model.
        Does not expand final states into their particle content and does not
        include anti-particles.
        """

        finalStates = []
        for fs in self.pdgDict.values():
            if any(fs is f for f in finalStates):
                continue
            finalStates.append(fs)

        return sorted(finalStates)
    
    def copy(self) -> "FinalStateModel":
        """
        Create a copy of the model.
        The copy contains the same vertices, but with copies of the original final states.
        """

        pdgDict = self.pdgDict.copy()
        new_vertices = []
        for v in self.vertices:
            new_vertex = DecayVertex(mother=pdgDict[v.mother.getPDGs()[0]],
                                     daughters=[pdgDict[d.getPDGs()[0]] for d in v.daughters],
                                     br=v.br,
                                     pdgDict=pdgDict)
            new_vertices.append(new_vertex)
        return FinalStateModel(vertices=new_vertices)


