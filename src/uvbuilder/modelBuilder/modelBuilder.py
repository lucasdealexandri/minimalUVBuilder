#!/usr/bin/env python3

import os
from dataclasses import replace
from typing import Dict,Union,List,Iterator,Dict,Optional
from itertools import product, chain

from uvbuilder.states.finalState import FinalState
from uvbuilder.states.particle import Particle
from uvbuilder.states.stateDict import StateDict
from uvbuilder.tools.reps import SCALAR,DIRAC_SPINOR,U1Rep,SU3Rep, AllReps, UnbrokenRep, GlobalU1s

from uvbuilder.vertices.particleVertex import ParticleVertex

from uvbuilder.models.finalStateModel import FinalStateModel
from uvbuilder.models.particleModel import ParticleModel
from uvbuilder.models.multipletModel import (MultipletModel, default_component_label,
                                            multiplet_model_from_partition)
from uvbuilder.models.feynrulesmodel import FeynRulesModel
from uvbuilder.partition.partitions_builder import PartitionsBuilder

from uvbuilder.modelBuilder.modelBuilderUtils import fillECharge,explodeColordim,explodeSpin,fillGlobalU1s

from uvbuilder.states.SMfinalStates import SMfinalStatesDict

import logging
logger = logging.getLogger(__name__)


# TODO: WHEN SETTING SPIN, SELECT ONLY SCALAR or DIRAC_SPINOR for BSM particles. Erro pra setar rep do modelo padrão


class ModelBuilder(object):
    """
    Class holding all the methods required for building models from the basic input data.
    """
   
    
    def __init__(self,
                 modelData: Union[Dict,FinalStateModel,None] = None,
                 dm_pdg : int = 1000022,
                 possible_dm_reps: Dict[str,List[AllReps]] = {}, 
                 SMpdgDict: StateDict = SMfinalStatesDict,
                 ) -> None:
        
        
        self.pdgDict : StateDict = StateDict(dict(SMpdgDict.items()))
        self.finalStateModel : Union[FinalStateModel,None] = None
        dm_default_reps = {'spin': [SCALAR, DIRAC_SPINOR],
                            'eCharge' : U1Rep(0),
                            'colordim' : SU3Rep(0,0),
                            'globalU1s': [GlobalU1s(), GlobalU1s(D=1)]
                           }
        # Define possible DM representations
        # (use default values if not defined)
        self.dm_pdg : int = dm_pdg
        self.possible_dm_reps = dict(possible_dm_reps.items())
        for key,value in dm_default_reps.items():            
            self.possible_dm_reps.setdefault(key,value)

        logger.debug(f'Created Model Builder:\n {self.pretty_print()}')

        if modelData is not None:
            self.load(modelData)

        self.check_types()

    def load(self, modelData : Union[Dict,FinalStateModel]):
        """
        Load the input data.
        It can be:
           * a dictionary with masses, decays and labels (FinalStateModel.fromDict)
           * a FinalStateModel object
        """

        if isinstance(modelData,dict):
            self.finalStateModel = FinalStateModel.fromDict(modelData,
                                                            pdgDict=self.pdgDict)
        elif isinstance(modelData,FinalStateModel):
            self.finalStateModel = modelData
        else:
            raise ValueError(f'modelData ({type(modelData)}) should be a dictionary or FinalStateModel object')

    def __repr__(self):
        return f'ModelBuilder for {self.finalStateModel}'
    
    def __str__(self):
        return self.__repr__()
        
    def pretty_print(self,addDecays=False):

        modelStr = "ModelBuilder"
        modelStr += '\nPossible DM representations:\n'
        modelStr += ', '.join([f'{key}: {self.possible_dm_reps[key]}' for key in self.possible_dm_reps])

        if self.finalStateModel is not None:
            modelStr += f'\n{self.finalStateModel.pretty_print(addDecays=addDecays)}'
        
        return modelStr

    def check_types(self):

        modelData = self.finalStateModel
        possible_dmn_reps = self.possible_dm_reps


        # Enforce type checking
        if not isinstance(modelData,Union[None,FinalStateModel]):
            raise ValueError(f'finalStateModel ({type(modelData)}) should be a FinalStateModel object')
        if not isinstance(possible_dmn_reps,dict):
            raise ValueError(f'possible_spin_reps ({type(possible_dmn_reps)}) should be a dictionary')
        elif len(possible_dmn_reps) == 0:
            raise ValueError(f'possible_dm_reps can not be an empty dictionary.')
        
    def getDMParticles(self, finalStateModel : FinalStateModel) -> List[Particle]:
        """
        Creates a list of Particle objects containing all possible DM particles
        (with all possible combinations of unbroken quantum numbers (spin, electric charge and color) 
        defined in self.possible_dm_reps)
        """

        dm_label = finalStateModel.pdgDict[self.dm_pdg].label
        
       # Get possible particle quantum numbers
        qnumbers = {q : self.possible_dm_reps.get(q,None) 
                    for q in Particle.QUANTUM_NUMBERS}

        values = []
        for q,qnumber in qnumbers.items():
            if qnumber is None:
                raise ValueError(f'No value defined for DM {q}.'
                                 'Check the definition of possible_dm_reps.')
            if isinstance(qnumber, AllReps):
                qnumber = [qnumber]
            elif not isinstance(qnumber, list):
                raise ValueError(
                    f'Value for DM {q} should be a list of representations or a single representation.'
                    'Check the definition of possible_dm_reps.'
                )
            values.append(qnumber)

        dm_particles = []
        for comb in product(*values):
            qDict = dict(zip(qnumbers.keys(), comb))
            rep = UnbrokenRep(**qDict)
            dm_ptc = Particle(rep=rep,label=dm_label, pdg=self.dm_pdg)
            src = finalStateModel.pdgDict[self.dm_pdg].getParticles()[0]
            dm_ptc.mass = src.mass
            dm_particles.append(dm_ptc)

        return dm_particles

    def explodeDMParticles(self, finalStateModel : FinalStateModel) -> Iterator[FinalStateModel]:
        """
        Expand the DM particle in the final states contained in self.finalStateModel 
        to all possible DM particles defined in self.possible_dm_reps, generating a list of FinalStateModels.
        """

        if finalStateModel is None:
            raise ValueError('No FinalStateModel defined. Please load a model before generating ParticleModels.')

        dm_particles = self.getDMParticles(finalStateModel)
        
        for dm_particle in dm_particles:
            newFinalStateModel = finalStateModel.copy()
            dm_finalState = newFinalStateModel.pdgDict[self.dm_pdg]
            # Set the particle in the DM final state to dm_particle
            assert isinstance(dm_finalState,FinalState)
            dm_finalState.particleList = [dm_particle]
            
            yield newFinalStateModel

    def finalStatesToParticles(self,finalStateModel: FinalStateModel, 
                               uniqueReps=False) -> Iterator[ParticleModel]:
        """
        Expand the final states contained in self.finalStateModel to ParticleModels
        (i.e. X -> jet + Y is expanded to X -> u + Y, X -> d + Y,...).
        If uniqueReps = True the expansion will be only to particles with distinct quantum numbers (when all quantum numbers are defined).
        Note, however, that particles with the same quantum numbers may be mapped to distinct fields when generating
        FieldModels, so this option should be used with caution.
        """


        if not isinstance(finalStateModel, FinalStateModel):
            raise ValueError(f'FinalStateModel should be a FinalStateModel object and not {type(finalStateModel)}.')
        
        expandedVertices = []
        for vertex in finalStateModel.vertices:
            expandedVertices.append(vertex.expandToParticleVertices(uniqueReps=uniqueReps))

        generatedParticleModels = []
        for vertices in product(*expandedVertices):
            particleModel = ParticleModel(vertices=list(vertices))
            if particleModel not in generatedParticleModels:
                generatedParticleModels.append(particleModel)
                yield particleModel
            else:
                continue

    def setParticleQNumbers(self, particleModel: ParticleModel) -> Iterator[ParticleModel]:
        """
        Fill the quantum numbers of the particles (eCharge, colordim and spin) in the ParticleModel based on the decay vertices and the quantum numbers of the final states.
        Since there can be multiple possibilities for the quantum numbers of a given particle, 
        this method generates a list of ParticleModels.
        """

        eChargeResults = list(fillECharge(particleModel))
        assert len(eChargeResults) == 1, f'Inconsistent model: {len(eChargeResults)} different eCharge assignments found for the particles in the model.'

        for eChargeRes in eChargeResults:
            model, conserved = self._project_global_u1s(eChargeRes)
            filled = fillGlobalU1s(model) if conserved else [model]
            for globalU1Res in filled:
                for colourRes in explodeColordim(globalU1Res):
                    for model in explodeSpin(colourRes):
                        assert isinstance(model, ParticleModel)
                        yield model

    def _project_global_u1s(self, particleModel: ParticleModel):
        """Fix the global U(1)s tracked in this model and prepare their filling.

        The conserved charges are the ones tracked by the DM candidate (a charge
        set explicitly, even to 0, is conserved). SM particles (the model's own
        copies) are projected onto them; the other BSM particles are marked
        unresolved (``None``) for fillGlobalU1s, or set to nothing-tracked when
        no charge is conserved.
        Returns (model copy, frozenset of conserved charge names).
        """
        model = particleModel.copy()
        dm = model.pdgDict[self.dm_pdg]
        conserved = dm.globalU1s.tracked
        for particle in model.getParticles():
            if particle is dm:
                continue
            if particle.isSM:
                particle.rep = replace(particle.rep, globalU1s=particle.globalU1s.project(conserved))
            else:
                particle.rep = replace(particle.rep, globalU1s=None if conserved else GlobalU1s())
        return model, conserved

    def generateParticleModels(self, finalStateModel: Union[FinalStateModel, None] = None, 
                               uniqueReps: bool = False) -> Iterator[ParticleModel]:
        """
        Converts the FinalStateModel to a list of ParticleModels, each one containing a set of particles 
        with well defined spin, electric charge and color consistent with the decay vertices.
        If multiple possibilities for the DM particle are defined in self.possible_dm_reps, 
        it will generate a list of ParticleModels for each of these possibilities.
        """

        if finalStateModel is None:
            finalStateModel = self.finalStateModel
        
        if finalStateModel is None:
            raise ValueError('No FinalStateModel defined. Please load a model before generating ParticleModels.')
        
        # Go over all possibple DM particles and generate a list of ParticleModels for each of them.
        for fsModel in self.explodeDMParticles(finalStateModel):
            # Go over all possible expansions of the final states to particles and generate a list of ParticleModels for each of them.
            for particleModel in self.finalStatesToParticles(fsModel,uniqueReps=uniqueReps):
                # Explode the quantum numbers of the particles in the model and generate a list of ParticleModels for each of them.
                # (if ParticleModel is inconsistent, skip it)
                try:
                    for particleModelFilled in self.setParticleQNumbers(particleModel):
                        yield particleModelFilled
                except Exception as error:
                    logger.debug(f'{particleModel} rejected: {error}')

    def generateMultipletModels(self, finalStateModel: Union[FinalStateModel, None] = None,
                                uniqueReps: bool = False,
                                labeler=default_component_label,
                                validOnly: bool = True) -> Iterator[MultipletModel]:
        """
        Generates the MultipletModels of all ParticleModels (see generateParticleModels),
        one for each partition of the BSM particles into gauge multiplets.
        The ParticleModel of each result is stored in its particle_model attribute.

        :param labeler: function naming the gauge fields of each multiplet, applied to
                        the valid models only (the others keep the default labels)
        :param validOnly: if True, yields only the models that reproduce all decay vertices
        """

        for particleModel in self.generateParticleModels(finalStateModel, uniqueReps=uniqueReps):
            for partition in PartitionsBuilder(particleModel).get_partitions():
                multipletModel = multiplet_model_from_partition(partition, particle_model=particleModel)
                valid = multipletModel.is_valid
                if validOnly and not valid:
                    continue
                if valid and labeler is not default_component_label:
                    multipletModel = multiplet_model_from_partition(partition, labeler=labeler,
                                                                    particle_model=particleModel)
                yield multipletModel

    def generateFeynRulesModels(self, finalStateModel: Union[FinalStateModel, None] = None,
                                uniqueReps: bool = False,
                                labeler=default_component_label,
                                outputFile: Optional[str] = None,
                                name: Optional[str] = None,
                                **kwargs) -> Iterator[FeynRulesModel]:
        """
        Generates one FeynRulesModel for each possible mixing of each valid MultipletModel
        (see generateMultipletModels).

        :param outputFile: file name of the models (default <name>.fr, or model.fr without
                           a name); with more than one model a running number is added
                           before the extension (model.fr -> model_1.fr, model_2.fr, ...)
        :param name: model name (M$ModelName), numbered like the files (name_1, name_2,
                     ...); without it every model is named as FeynRulesModel's default
        :param kwargs: additional arguments of FeynRulesModel (sm_file, strict, author,
                       version)
        """

        root, extension = os.path.splitext(outputFile or f'{name or "model"}.fr')
        pairs = ((multipletModel, mixings)
                 for multipletModel in self.generateMultipletModels(
                     finalStateModel, uniqueReps=uniqueReps, labeler=labeler)
                 for mixings in multipletModel.possible_mixings)
        # a single model is unnumbered: look one model ahead before yielding the first
        first, second = next(pairs, None), next(pairs, None)
        if first is None:
            return
        if second is None:
            numbered = [(None, first)]
        else:
            numbered = enumerate(chain([first, second], pairs), 1)
        for number, (multipletModel, mixings) in numbered:
            suffix = '' if number is None else f'_{number}'
            names = {'name': f'{name}{suffix}'} if name else {}
            yield FeynRulesModel(multipletModel, multipletModel.particle_model, mixings,
                                 f'{root}{suffix}{extension}', **names, **kwargs)
