import itertools
from typing import List
from uvbuilder.tools.reps import UnbrokenRep, BrokenRep, WEYL_RIGHT
import networkx as nx
from networkx import DiGraph
from collections import Counter
from uvbuilder.models.particleModel import ParticleModel
from uvbuilder.models.multipletModel import (MultipletModel, default_component_label,
                                       multiplet_model_from_partition)


class PartitionsBuilder:

    def __init__(self, particle_model: ParticleModel):
        self.particle_model = particle_model

        self.bsm_particles = self.particle_model.bsm_particles
        
        self.unbroken_reps = self.particle_model.unbroken_reps

    def charge_graph(self, unbroken_reps: List[dict]) -> DiGraph:

        unique_reps = []
        for rep in unbroken_reps:
            if rep not in unique_reps:
                unique_reps.append(rep)

        graph = DiGraph()
        graph.add_nodes_from(unique_reps)

        for rep1, rep2 in itertools.combinations(unique_reps,2):

            if rep1.can_connect(rep2):

                if rep1.eCharge < rep2.eCharge:
                    graph.add_edge(rep1, rep2)
                else:
                    graph.add_edge(rep2, rep1)

        return graph


    def building_blocks(self, charge_graph: DiGraph) -> List[List[UnbrokenRep]]:

        blocks = []
        for rep1, rep2 in itertools.product(charge_graph.nodes, charge_graph.nodes):

            simple_paths = list(nx.all_simple_paths(charge_graph, rep1, rep2))
            if simple_paths:
                blocks.append(*simple_paths)

        return blocks


    def valid_partition(self, partition: List[BrokenRep]) -> bool:
        '''
        A Given partition is valid if every broken representation has a
        conjugate representation in it (vectorlike condition).
        This must strictly account for multiplicity to allow mass terms.
        '''
        
        # A valid vector-like partition means the count of all reps 
        # is exactly equal to the count of all their conjugates.
        return Counter(partition) == Counter(rep.conjInternalRep for rep in partition)


    def get_partitions(self, unbroken_reps: List[UnbrokenRep] = None) -> List[BrokenRep]:

        '''
        Returns all valid broken reps a list of unbroken reps can have come from.
        '''

        unbroken_reps = self.unbroken_reps if unbroken_reps is None else unbroken_reps
        valid_partitions = []
        seen_signatures = set()

        # Recursive step
        def _search(
                current_unbroken_reps: List[UnbrokenRep], 
                current_partition: List[BrokenRep]):

            if not current_unbroken_reps:
                signature = frozenset(Counter(current_partition).items())

                if signature not in seen_signatures and self.valid_partition(current_partition):
                    seen_signatures.add(signature)
                    valid_partitions.append(current_partition)
                return 

            charge_graph = self.charge_graph(current_unbroken_reps)
            anchor = next(iter(charge_graph.nodes))

            for block in self.building_blocks(charge_graph):

                if anchor in block:

                    new_rep = BrokenRep.from_block(block)
                    remaining = list(
                        (Counter(current_unbroken_reps) - Counter(block)).elements())

                    _search(remaining, current_partition + [new_rep])

        _search(unbroken_reps, [])

        return valid_partitions

    def get_multipletModels(self, labeler=default_component_label):
        '''Each partition becomes a
        MultipletModel whose components are named gauge fields (labeler).'''
        return [multiplet_model_from_partition(partition, labeler=labeler,
                                               particle_model=self.particle_model)
                for partition in self.get_partitions()]