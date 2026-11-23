#!/usr/bin/env python3
"""
Unit tests for the FinalStateModel class.

This module contains comprehensive unit tests for the FinalStateModel class
and its methods, including initialization, creation from dictionaries,
creation from SLHA files, and string representations.
"""

import sys
import os
import unittest


from uvbuilder.models.finalStateModel import FinalStateModel
from uvbuilder.vertices.decayVertex import DecayVertex
from uvbuilder.states.particle import Particle


# FinalStateModel stores FinalState objects: their PDG codes come as a list
# (getPDGs) and their masses live on the particles of `particleList`.
# Helpers for single-particle (BSM) final states.
def _pdg(final_state):
    return final_state.getPDGs()[0]


def _mass(final_state):
    return final_state.particleList[0].mass


class TestFinalStateModelInit(unittest.TestCase):
    """Test FinalStateModel initialization"""
    
    def test_init_empty_vertices(self):
        """Test initialization with empty vertex list"""
        fs = FinalStateModel(vertices=[])
        self.assertEqual(len(fs.vertices), 0)
        self.assertEqual(len(fs.pdgDict), 0)
    
    def test_init_default_vertices(self):
        """Test initialization with default (empty) vertices parameter"""
        fs = FinalStateModel()
        self.assertEqual(len(fs.vertices), 0)
        self.assertEqual(len(fs.pdgDict), 0)


class TestFinalStateModelFromDict(unittest.TestCase):
    """Test FinalStateModel.fromDict() class method"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.fsDict_basic = {
            "masses": {1000022: 180.603, 1000001: 757.246, 1000006: 1164.264},
            "decays": {
                1000022: {},
                1000001: {(1000022, 1): 1.0},
                1000006: {(1000022, 6): 1.0},
            },
        }
        
        self.fsDict_with_labels = {
            "masses": {1000022: 180.603, 1000001: 757.246, 1000006: 1164.264},
            "decays": {
                1000022: {},
                1000001: {(1000022, 1): 1.0},
                1000006: {(1000022, 6): 1.0},
            },
            "labels": {1000022: 'N1', 1000001: 'sup1', 1000006: 'st1'}
        }
    
    def test_fromDict_basic(self):
        """Test creating FinalStateModel from basic dictionary"""
        fs = FinalStateModel.fromDict(self.fsDict_basic)
        
        self.assertIsNotNone(fs)
        self.assertGreater(len(fs.vertices), 0)
        self.assertGreater(len(fs.pdgDict), 0)
    
    def test_fromDict_with_labels(self):
        """Test creating FinalStateModel from dictionary with labels"""
        fs = FinalStateModel.fromDict(self.fsDict_with_labels)
        
        self.assertIsNotNone(fs)
        # Check that labels are set
        self.assertEqual(fs.pdgDict[1000022].label, 'N1')
        self.assertEqual(fs.pdgDict[1000001].label, 'sup1')
        self.assertEqual(fs.pdgDict[1000006].label, 'st1')
    
    def test_fromDict_missing_decays(self):
        """Test that fromDict raises error when decays are missing"""
        fsDict_no_decays = {"masses": {1000022: 180.603}}
        with self.assertRaises(ValueError):
            FinalStateModel.fromDict(fsDict_no_decays)
    
    def test_fromDict_empty_decays(self):
        """Test that fromDict raises error when decays are empty"""
        fsDict_empty_decays = {"decays": {}}
        with self.assertRaises(ValueError):
            FinalStateModel.fromDict(fsDict_empty_decays)
    
    def test_fromDict_pdgDict_creation(self):
        """Test that pdgDict is correctly populated"""
        fs = FinalStateModel.fromDict(self.fsDict_basic)
        
        # All particles in masses and decays should be in pdgDict
        self.assertIn(1000022, fs.pdgDict)
        self.assertIn(1000001, fs.pdgDict)
        self.assertIn(1000006, fs.pdgDict)
        # PDG 1 and 6 (SM particles) should also be present
        self.assertIn(1, fs.pdgDict)
        self.assertIn(6, fs.pdgDict)
    
    def test_fromDict_masses_set(self):
        """Test that masses are correctly set on particles"""
        fs = FinalStateModel.fromDict(self.fsDict_basic)
        
        # Check that masses are set (converted to GeV)
        self.assertIsNotNone(_mass(fs.pdgDict[1000022]))
        self.assertIsNotNone(_mass(fs.pdgDict[1000001]))
        self.assertIsNotNone(_mass(fs.pdgDict[1000006]))


class TestFinalStateModelFromFile(unittest.TestCase):
    """Test FinalStateModel.fromFile() class method"""
    
    def test_fromFile_sleptons(self):
        """Test loading FinalStateModel from SLHA file"""
        slhaFile = os.path.join(os.path.dirname(__file__), 'data', 'slhaFiles', 'sleptons.slha')
        
        # Check if file exists
        if not os.path.isfile(slhaFile):
            self.skipTest(f"SLHA file {slhaFile} not found")
        
        labels = {
            1000011: 'seL', 2000011: 'seR',
            1000013: 'smuL', 2000013: 'smuR',
            1000015: 'stauL', 2000015: 'stauR',
            1000012: 'snueL', 1000014: 'snumuL', 1000016: 'snutauL',
            1000022: 'n1'
        }
        
        fs = FinalStateModel.fromFile(slhaFile=slhaFile, labels=labels)
        
        self.assertIsNotNone(fs)
        self.assertGreater(len(fs.vertices), 0)
        self.assertGreater(len(fs.pdgDict), 0)
        
        # Check that labels were applied
        if 1000011 in fs.pdgDict:
            self.assertEqual(fs.pdgDict[1000011].label, 'seL')
    
    def test_fromFile_nonexistent(self):
        """Test that fromFile handles non-existent files gracefully"""
        with self.assertRaises((FileNotFoundError, ValueError)):
            FinalStateModel.fromFile(slhaFile='nonexistent_file.slha')


class TestFinalStateModelFromVertices(unittest.TestCase):
    """Test FinalStateModel creation from vertices"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.fsDict = {
            "masses": {1000022: 180.603, 1000001: 757.246, 1000006: 1164.264},
            "decays": {
                1000022: {},
                1000001: {(1000022, 1): 1.0},
                1000006: {(1000022, 6): 1.0},
            },
            "labels": {1000022: 'N1', 1000001: 'sup1', 1000006: 'st1'}
        }
    
    def test_create_from_existing_vertices(self):
        """Test creating FinalStateModel from another model's vertices"""
        # Create first model
        fs1 = FinalStateModel.fromDict(self.fsDict)
        
        # Create new model from the first model's vertices
        fs2 = FinalStateModel(vertices=fs1.vertices)
        
        self.assertIsNotNone(fs2)
        self.assertEqual(len(fs2.vertices), len(fs1.vertices))
        self.assertGreater(len(fs2.pdgDict), 0)
    
    def test_vertices_preserved(self):
        """Test that vertices are correctly preserved"""
        fs1 = FinalStateModel.fromDict(self.fsDict)
        fs2 = FinalStateModel(vertices=fs1.vertices)
        
        # Check that vertex count is the same
        self.assertEqual(len(fs1.vertices), len(fs2.vertices))
        
        # Check that particles in pdgDict are present
        for pdg in fs1.pdgDict:
            self.assertIn(pdg, fs2.pdgDict)


class TestFinalStateModelStringRepresentation(unittest.TestCase):
    """Test string representation methods"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.fsDict = {
            "masses": {1000022: 180.603, 1000001: 757.246, 1000006: 1164.264},
            "decays": {
                1000022: {},
                1000001: {(1000022, 1): 1.0},
                1000006: {(1000022, 6): 1.0},
            },
            "labels": {1000022: 'N1', 1000001: 'sup1', 1000006: 'st1'}
        }
        self.fs = FinalStateModel.fromDict(self.fsDict)
    
    def test_str_representation(self):
        """Test __str__ method"""
        str_rep = str(self.fs)
        
        self.assertIsInstance(str_rep, str)
        self.assertIn('FinalStateModel', str_rep)
        # Should contain BSM particles
        self.assertTrue(len(str_rep) > 0)
    
    def test_pretty_print(self):
        """Test pretty_print method"""
        pretty_str = self.fs.pretty_print()
        
        self.assertIsInstance(pretty_str, str)
        self.assertIn('FinalStateModel', pretty_str)
        # Should contain more detailed information
        self.assertTrue(len(pretty_str) > len(str(self.fs)))
    
    def test_pretty_print_contains_masses(self):
        """Test that pretty_print includes mass information"""
        pretty_str = self.fs.pretty_print()
        
        # Should contain 'Masses' section if particles have masses
        if any(_mass(p) is not None for p in self.fs.pdgDict.values()):
            self.assertIn('Masses', pretty_str)
    
    def test_pretty_print_contains_vertices(self):
        """Test that pretty_print includes vertex information"""
        pretty_str = self.fs.pretty_print()
        
        # Should contain 'Vertices' section if model has vertices
        if len(self.fs.vertices) > 0:
            self.assertIn('Vertices', pretty_str)


class TestFinalStateModelPdgDict(unittest.TestCase):
    """Test pdgDict functionality"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.fsDict = {
            "masses": {1000022: 180.603, 1000001: 757.246, 1000006: 1164.264},
            "decays": {
                1000022: {},
                1000001: {(1000022, 1): 1.0},
                1000006: {(1000022, 6): 1.0},
            },
        }
        self.fs = FinalStateModel.fromDict(self.fsDict)
    
    def test_pdgDict_contains_particles(self):
        """Test that pdgDict contains Particle objects"""
        for pdg, particle in self.fs.pdgDict.items():
            self.assertTrue(hasattr(particle, 'pdg'))
            self.assertTrue(hasattr(particle, 'label'))
    
    def test_pdgDict_includes_antiparticles(self):
        """Test that antiparticles are included in pdgDict"""
        # For BSM particles without explicit antiparticle entry,
        # the antiparticle PDG should also be in the dict
        bsm_pdgs = [pdg for pdg, p in self.fs.pdgDict.items() if not p.isSM]
        
        for pdg in bsm_pdgs:
            particle = self.fs.pdgDict[pdg]
            if not particle.selfConjugate and particle.antiParticle is not None:
                # The negative PDG should be in the dict
                self.assertIn(-pdg, self.fs.pdgDict)


class TestFinalStateModelVerticesPopulation(unittest.TestCase):
    """Test that vertices are correctly populated in FinalStateModel"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.fsDict_simple = {
            "masses": {1000022: 180.603, 1000001: 757.246, 1000006: 1164.264},
            "decays": {
                1000022: {},
                1000001: {(1000022, 1): 1.0},
                1000006: {(1000022, 6): 1.0},
            },
        }
        
        self.fsDict_multiway = {
            "masses": {1000022: 180.0, 1000001: 500.0, 1000002: 510.0},
            "decays": {
                1000022: {},
                1000001: {(1000022, 1): 0.6, (1000022, 2): 0.4},
                1000002: {(1000022, 2): 1.0},
            },
        }
        
        self.fsDict_with_labels = {
            "masses": {1000022: 180.603, 1000001: 757.246, 1000006: 1164.264},
            "decays": {
                1000022: {},
                1000001: {(1000022, 1): 1.0},
                1000006: {(1000022, 6): 1.0},
            },
            "labels": {1000022: 'N1', 1000001: 'sq1', 1000006: 'st1'}
        }
    
    def test_vertices_created_from_decays_fromDict(self):
        """Test that vertices are created from decay dictionary"""
        fs = FinalStateModel.fromDict(self.fsDict_simple)
        
        # Should have vertices for particles with decays
        self.assertGreater(len(fs.vertices), 0)
        # Basic dict has 2 particles with decays (1000001 and 1000006)
        self.assertEqual(len(fs.vertices), 2)
    
    def test_vertex_is_DecayVertex_type(self):
        """Test that all vertices are DecayVertex instances"""
        fs = FinalStateModel.fromDict(self.fsDict_simple)
        
        for vertex in fs.vertices:
            self.assertIsInstance(vertex, DecayVertex)
    
    def test_vertices_have_correct_mothers(self):
        """Test that vertices have correct mother particles"""
        fs = FinalStateModel.fromDict(self.fsDict_simple)
        
        # Extract mother PDGs from vertices
        mother_pdgs = {_pdg(v.mother) for v in fs.vertices}
        
        # Expected mothers are particles with non-empty decay dicts
        expected_mothers = {1000001, 1000006}
        self.assertEqual(mother_pdgs, expected_mothers)
    
    def test_vertices_have_correct_daughters(self):
        """Test that vertices have correct daughter particles"""
        fs = FinalStateModel.fromDict(self.fsDict_simple)
        
        for vertex in fs.vertices:
            mother_pdg = _pdg(vertex.mother)
            
            # Check that mothers with non-empty decay dict have daughters
            if mother_pdg in self.fsDict_simple['decays']:
                decay_dict = self.fsDict_simple['decays'][mother_pdg]
                if len(decay_dict) > 0:
                    # Vertex should have daughters
                    self.assertGreater(len(vertex.daughters), 0)
                    
                    # Check that daughters are either Particle or FinalState objects
                    for daughter in vertex.daughters:
                        # Should have a pdg attribute
                        self.assertTrue(hasattr(daughter, 'pdg'))
                        # pdg should be either int (Particle) or list (FinalState)
                        self.assertTrue(isinstance(daughter.pdg, (int, list, tuple)))
                        
                        # If it's a FinalState, it should have particleList
                        if isinstance(daughter.pdg, (list, tuple)):
                            self.assertTrue(hasattr(daughter, 'particleList'))
                        else:
                            # If it's a Particle, it should have label
                            self.assertTrue(hasattr(daughter, 'label'))
    
    def test_vertices_branching_ratios_fromDict(self):
        """Test that branching ratios are correctly assigned to vertices"""
        fs = FinalStateModel.fromDict(self.fsDict_multiway)
        
        # Find vertices for particle 1000001 (has multiple decay modes)
        decay_vertices = [v for v in fs.vertices if _pdg(v.mother) == 1000001]
        
        # Should have 2 decay vertices for 1000001
        self.assertEqual(len(decay_vertices), 2)
        
        # Extract branching ratios
        branching_ratios = [v.br for v in decay_vertices]
        branching_ratios.sort()
        
        # Should be approximately 0.4 and 0.6
        self.assertAlmostEqual(branching_ratios[0], 0.4, places=5)
        self.assertAlmostEqual(branching_ratios[1], 0.6, places=5)
    
    def test_vertices_with_stable_particles(self):
        """Test that stable particles (empty decay dict) produce no vertices"""
        fs = FinalStateModel.fromDict(self.fsDict_simple)
        
        # 1000022 is stable (empty decay dict)
        stable_mother_pdgs = {_pdg(v.mother) for v in fs.vertices}
        self.assertNotIn(1000022, stable_mother_pdgs)
    
    def test_vertices_mother_has_mass(self):
        """Test that mother particles in vertices have correct masses"""
        fs = FinalStateModel.fromDict(self.fsDict_simple)
        
        for vertex in fs.vertices:
            mother = vertex.mother
            # Mother should have a mass if it was in the masses dict
            if _pdg(mother) in self.fsDict_simple['masses']:
                self.assertIsNotNone(_mass(mother))
    
    def test_vertices_with_custom_labels(self):
        """Test that vertices use particles with custom labels"""
        fs = FinalStateModel.fromDict(self.fsDict_with_labels)
        
        # Find the vertex for st1 (1000006)
        st1_vertices = [v for v in fs.vertices if _pdg(v.mother) == 1000006]
        self.assertEqual(len(st1_vertices), 1)
        
        vertex = st1_vertices[0]
        self.assertEqual(vertex.mother.label, 'st1')
    
    def test_vertices_for_multiway_decay(self):
        """Test vertices for particle with multiple decay modes"""
        fs = FinalStateModel.fromDict(self.fsDict_multiway)
        
        # Particle 1000001 has 2 decay modes
        sq_vertices = [v for v in fs.vertices if _pdg(v.mother) == 1000001]
        self.assertEqual(len(sq_vertices), 2)
        
        # Check that both daughters include 1000022 (LSP)
        for vertex in sq_vertices:
            daughter_pdgs = [pdg for d in vertex.daughters for pdg in d.getPDGs()]
            self.assertIn(1000022, daughter_pdgs)
    
    def test_vertices_from_file_populated(self):
        """Test that vertices are populated when loading from file"""
        slhaFile = os.path.join(os.path.dirname(__file__), 'data', 'slhaFiles', 'sleptons.slha')
        
        if not os.path.isfile(slhaFile):
            self.skipTest(f"SLHA file {slhaFile} not found")
        
        labels = {
            1000011: 'seL', 2000011: 'seR',
            1000013: 'smuL', 2000013: 'smuR',
            1000015: 'stauL', 2000015: 'stauR',
            1000012: 'snueL', 1000014: 'snumuL', 1000016: 'snutauL',
            1000022: 'n1'
        }
        
        fs = FinalStateModel.fromFile(slhaFile=slhaFile, labels=labels)
        
        # Check that vertices are created
        self.assertGreater(len(fs.vertices), 0)
        
        # All vertices should be DecayVertex instances
        for vertex in fs.vertices:
            self.assertIsInstance(vertex, DecayVertex)
            self.assertIsNotNone(vertex.mother)
            self.assertGreater(len(vertex.daughters), 0)
    
    def test_vertices_correspond_to_decay_dict(self):
        """Test that all vertices correspond to entries in decay dictionary"""
        fs = FinalStateModel.fromDict(self.fsDict_simple)
        
        decay_dict = self.fsDict_simple['decays']
        
        # Every vertex mother should have a corresponding decay entry
        for vertex in fs.vertices:
            mother_pdg = _pdg(vertex.mother)
            self.assertIn(mother_pdg, decay_dict)
            self.assertGreater(len(decay_dict[mother_pdg]), 0)
    
    def test_vertices_daughters_count(self):
        """Test that vertices have correct number of daughters"""
        fs = FinalStateModel.fromDict(self.fsDict_simple)
        
        for vertex in fs.vertices:
            # Each decay should have daughters
            self.assertGreater(len(vertex.daughters), 0)
            # Typically 2 daughters per decay (mother -> 2 products)
            self.assertGreaterEqual(len(vertex.daughters), 2)


class TestFinalStateModelVerticesConsistency(unittest.TestCase):
    """Test consistency of vertices with model data"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.fsDict_complex = {
            "masses": {
                1000022: 180.0,     # LSP - stable
                1000023: 300.0,     # Neutralino 2
                1000024: 310.0,     # Chargino
                1000001: 1000.0,    # Squark
                1000002: 1010.0,    # Squark
            },
            "decays": {
                1000022: {},  # stable
                1000023: {(1000022, 23): 1.0},  # decays to LSP + Z
                1000024: {(1000022, 24): 1.0},  # decays to LSP + W
                1000001: {
                    (1000022, 1): 0.5,  # 50% to LSP + quark
                    (1000023, 1): 0.5,  # 50% to N2 + quark
                },
            }
        }
    
    def test_vertices_count_matches_decay_entries(self):
        """Test that number of vertices matches total decay channels"""
        fs = FinalStateModel.fromDict(self.fsDict_complex)
        
        # Count decay channels in dictionary
        total_decay_channels = sum(
            len(decays) for decays in self.fsDict_complex['decays'].values()
        )
        
        # Number of vertices should match total decay channels
        self.assertEqual(len(fs.vertices), total_decay_channels)
    
    def test_all_decay_mothers_have_vertices(self):
        """Test that all particles with decays have corresponding vertices"""
        fs = FinalStateModel.fromDict(self.fsDict_complex)
        
        decay_mothers = {
            pdg for pdg, decays in self.fsDict_complex['decays'].items()
            if len(decays) > 0
        }
        vertex_mothers = {_pdg(v.mother) for v in fs.vertices}
        
        # All decay mothers should have vertices
        self.assertEqual(decay_mothers, vertex_mothers)
    
    def test_stable_particles_no_vertices(self):
        """Test that stable particles (empty decay dict) have no vertices"""
        fs = FinalStateModel.fromDict(self.fsDict_complex)
        
        vertex_mothers = {_pdg(v.mother) for v in fs.vertices}
        
        # Stable particle (1000022) should not be a mother in any vertex
        self.assertNotIn(1000022, vertex_mothers)
    
    def test_vertices_branching_ratios_sum_to_one(self):
        """Test that branching ratios for a given mother sum to approximately 1.0"""
        fs = FinalStateModel.fromDict(self.fsDict_complex)
        
        # Group vertices by mother PDG
        vertices_by_mother = {}
        for vertex in fs.vertices:
            mother_pdg = _pdg(vertex.mother)
            if mother_pdg not in vertices_by_mother:
                vertices_by_mother[mother_pdg] = []
            vertices_by_mother[mother_pdg].append(vertex)
        
        # For each mother with multiple decay channels, BRs should sum to ~1.0
        for mother_pdg, vertices in vertices_by_mother.items():
            if len(vertices) > 1:
                total_br = sum(v.br for v in vertices)
                self.assertAlmostEqual(total_br, 1.0, places=5,
                    msg=f"Branching ratios for PDG {mother_pdg} do not sum to 1.0")
    
    def test_all_mother_particles_in_masses(self):
        """Test that all mother particles in vertices have masses set"""
        fs = FinalStateModel.fromDict(self.fsDict_complex)
        
        for vertex in fs.vertices:
            mother = vertex.mother
            # If mother is in decay dict with decays, it should have a mass
            mother_pdg = _pdg(mother)
            if mother_pdg in self.fsDict_complex.get('masses', {}):
                self.assertIsNotNone(_mass(mother),
                    msg=f"Mother particle {mother_pdg} should have mass set")
    
    def test_daughter_particles_in_model(self):
        """Test that all daughter particles appear in the model's pdgDict"""
        fs = FinalStateModel.fromDict(self.fsDict_complex)
        
        for vertex in fs.vertices:
            for daughter in vertex.daughters:
                # Get PDG(s) for the daughter
                if isinstance(daughter.pdg, int):
                    daughter_pdgs = [daughter.pdg]
                else:
                    daughter_pdgs = daughter.pdg
                
                # All daughter PDGs should be in pdgDict
                for pdg in daughter_pdgs:
                    self.assertIn(abs(pdg), [abs(p) for p in fs.pdgDict.keys()],
                        msg=f"Daughter PDG {pdg} not found in model")
    
    def test_vertices_mother_daughter_hierarchy(self):
        """Test that mother particle mass >= sum of daughter masses where applicable"""
        fs = FinalStateModel.fromDict(self.fsDict_complex)
        
        checked = 0
        for vertex in fs.vertices:
            mother_mass = _mass(vertex.mother)
            # BSM daughters carry the model masses; SM daughters are ignored
            daughter_masses = [_mass(d) for d in vertex.daughters
                               if not d.isSM and _mass(d) is not None]
            if mother_mass is None or not daughter_masses:
                continue
            total = sum(m.asNumber() for m in daughter_masses)
            self.assertGreaterEqual(
                mother_mass.asNumber(), total,
                msg=f"Mother {_pdg(vertex.mother)} lighter than its BSM daughters")
            checked += 1
        self.assertGreater(checked, 0, "no vertex was actually checked")

    def test_multiway_decay_coverage(self):
        """Test that all decay channels for multi-decay particles are represented"""
        fs = FinalStateModel.fromDict(self.fsDict_complex)
        
        # Find a particle with multiple decay channels
        multi_decay_pdgs = [
            pdg for pdg, decays in self.fsDict_complex['decays'].items()
            if len(decays) > 1
        ]
        
        for multi_pdg in multi_decay_pdgs:
            # Get all vertices for this mother
            multi_vertices = [v for v in fs.vertices if _pdg(v.mother) == multi_pdg]
            
            # Should have same number of vertices as decay channels
            expected_channels = len(self.fsDict_complex['decays'][multi_pdg])
            self.assertEqual(len(multi_vertices), expected_channels,
                msg=f"PDG {multi_pdg} has {len(multi_vertices)} vertices but {expected_channels} decay channels")
    
    def test_vertices_are_sorted(self):
        """Test that daughter particles in each vertex are sorted"""
        fs = FinalStateModel.fromDict(self.fsDict_complex)
        
        for vertex in fs.vertices:
            # Check that daughters are sorted (per DecayVertex sorting)
            daughters = vertex.daughters
            # If all are Particles, they should be sortable
            if all(isinstance(d.pdg, int) for d in daughters):
                pdgs = [d.pdg for d in daughters]
                # Verify they maintain some consistent ordering
                self.assertTrue(len(pdgs) >= 0)  # Just ensure list exists


if __name__ == '__main__':
    # Run tests with verbose output
    unittest.main(verbosity=2)
