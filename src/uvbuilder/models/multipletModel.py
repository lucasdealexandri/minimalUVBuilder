"""
.. module:: multipletModel
   :synopsis: Defines the MultipletModel class, a model of SU(2) multiplets and their couplings.

.. moduleauthor:: John Gargalionis <johngargalionis@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Lucas Heck <lucasdealexandri.heck@gmail.com>
"""

from collections import Counter
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Sequence, Tuple, Optional, Set

import itertools
import warnings
import networkx as nx
import sympy as sp
from networkx import Graph

from sympy.printing.mathematica import mathematica_code
from uvbuilder.tools.operators import Operator
from uvbuilder.tools.cg_tensors import SU2CGTensor, su2_singlet_paths
from uvbuilder.models.particleModel import ParticleModel
from uvbuilder.states.field import Field
from uvbuilder.states.multiplet import Multiplet
from uvbuilder.states.particle import Particle
from uvbuilder.states.SMfields import SM_MULTIPLETS, HiggsDoublet
from uvbuilder.tools.reps import su2_conj_phase
from uvbuilder.tools.reps import (BrokenRep, Scalar, U1Rep, UnbrokenRep, WEYL_LEFT,
                        WEYL_RIGHT, SU3Rep, Vector, WeylSpinor, DiracSpinor)
from uvbuilder.tools.fr_utils import (formatVariable, formatClassBlock, merge_chunks,
                            getSU2Index, HIGGS_TILDE, toggle_dagger, DAGGER_SUFFIX,
                            fr_quantum_numbers, RIGHT_SUFFIX)
from uvbuilder.partition.particle_mixing import ParticleMixing, MixingMatrix
from uvbuilder.vertices.particleVertex import ParticleVertex


# ----------------------------------------------------------------------
# Default component naming (overridable at partition -> MultipletModel build)
# ----------------------------------------------------------------------


def _spin_letter(spin) -> str:
    '''Spin letter of a default component label: S (scalar), V (vector), LW (Weyl).'''
    name = type(spin).__name__.lower()
    if 'scalar' in name:
        return 'S'
    if 'vector' in name:
        return 'V'
    return 'LW'


def _spin_prefix(spin) -> str:
    '''Label prefix, as in the physical class names: Chi (fermions; the chirality
    is the Lw/Rw suffix of Multiplet.fr_symbol), Phi (scalars), V (vectors).'''
    if isinstance(spin, (WeylSpinor, DiracSpinor)):
        return 'Chi'
    if isinstance(spin, Scalar):
        return 'Phi'
    if isinstance(spin, Vector):
        return 'V'
    return 'Chi'


def _default_charge_suffix(q: U1Rep) -> str:
    '''Charge superscript of a default component label: ^0, ^+, ^{-2}, ^{+1/3}, ...'''
    if q.numerator == 0:
        return '^0'
    sign = '+' if q.numerator > 0 else '-'
    if q.denominator == 1:
        if abs(q.numerator) == 1: return f'^{sign}'
        return f'^{{{sign}{abs(q.numerator)}}}'
    return f'^{{{sign}{abs(q.numerator)}/{q.denominator}}}'


def sympy_component_label(
        rep: BrokenRep, 
        i3: Optional[U1Rep] = None, 
        family_index: Optional[int] = None) -> str:
    '''Human-readable gauge-field label for one component, e.g. LW2(-1/2)0.
    Not a valid bare Mathematica symbol (contains '(', ')', '^', ...); use
    fr_component_label for labels that port directly to FeynRules.'''
    colorname = f'[{rep.colordim}]' if rep.colordim != SU3Rep(0,0) else ''
    base = f'{_spin_letter(rep.spin)}{colorname}{rep.isospin.dim}({rep.hypercharge})'
    suffix = f'_{family_index}' if family_index is not None else ''
    if not rep.globalU1s.is_zero: base += f'{rep.globalU1s}'
    if i3 is None: return base + suffix
    else:
        q = U1Rep(i3) + rep.hypercharge
        return f'{base}{_default_charge_suffix(q)}{suffix}'


def _int_charge(q: U1Rep) -> int:
    '''Electric charge / hypercharge scaled by 3, so values like +/-1/2, +/-, 1/3
    all map onto distinct integers (e.g. 6q).'''
    assert 6 * q == (6*q).numerator, f'Not integer charge {6*q}'
    return 6 * q


def _charge_tag(q: U1Rep) -> str:
    '''Mathematica-valid scalar charge tag: y0, y3, ym3, ... (hypercharge) or
    q0, q6, qm3, ... (electric charge). Values scaled by 6 to stay integer.'''
    n = _int_charge(q)
    if n == 0:
        return '0'
    if n < 0:
        return f'm{abs(n)}'
    return f'{n}'


_UQN_PREFIX = {'L': 'L', 'B': 'B', 'D': 'D'}
_UQN_SCALE = {'B': 3}


def _global_u1_tag(uqn) -> str:
    '''Mathematica-valid tag for the extra quantum numbers, one per charge in
    unit-adjusted to stay integral: L=1 -> L1, D=1 -> D1, B=1/3 -> B1.
    Only non-zero tracked charges are tagged.'''
    tag = ''
    for name, charge in uqn.items():
        if charge == 0:
            continue
        prefix = _UQN_PREFIX.get(name, name.lower())
        val = charge * _UQN_SCALE.get(name, 1)
        if val < 0:
            tag += f'{prefix}m{abs(val)}'
        else:
            tag += f'{prefix}{val}'
    return tag


def fr_component_label(
        rep: BrokenRep, 
        i3: Optional[U1Rep] = None, 
        family_index: Optional[int] = None) -> str:
    '''Mathematica-portable gauge-field label, {Phi|Chi}[C{dim}[b]]I{dim}Y{6Y}{globals}[F{n}]
    for a multiplet, plus q{6q} for one component (e.g. ChiC3I2Y1q4).

    Encodes the spin (Phi scalar, Chi fermion; Weyl chirality is the Lw/Rw suffix
    added by Multiplet.fr_symbol), color (dimension, conjugate marked with 'b',
    omitted for singlets), isospin dimension, hypercharge, non-zero tracked
    global charges (B scaled by 3), the family index and the component electric
    charge q = i3 + Y (charges scaled by 6, 'm' = minus). The result is a valid
    bare Mathematica symbol: only ASCII letters and digits. Example labels::

        PhiI2Y3       # scalar doublet, Y = 1/2 (component PhiI2Y3q6: q = 1)
        ChiC3I1Y4B1   # fermion, color 3, singlet, Y = 2/3, B = 1/3
        ChiC3bI1Ym4Bm1  # its conjugate rep (color 3bar, Y = -2/3, B = -1/3)
        PhiI1Y0F2     # second family of a real singlet
    '''
    prefix = _spin_prefix(rep.spin)

    if rep.colordim == SU3Rep(0, 0):
        color = ''
    else:
        bar = 'b' if rep.colordim.lower > rep.colordim.upper else ''
        color = f'C{rep.colordim.dim}{bar}'

    isospin = f'I{rep.isospin.dim}'
    hyper = f'Y{_charge_tag(rep.hypercharge)}'

    family = f'F{family_index}' if family_index is not None else ''

    base = f'{prefix}{color}{isospin}{hyper}'
    if not rep.globalU1s.is_zero:
        base += _global_u1_tag(rep.globalU1s)
    if i3 is None:
        return base + family
    q = U1Rep(i3) + rep.hypercharge
    return f'{base}q{_charge_tag(q)}{family}'


default_component_label = fr_component_label

def default_operator_parameter(
        multiplets: List[Multiplet]):
    '''Default parameter of the operator of these multiplets (docs/conventions.md):
    M^2 (M^2/2 for a real scalar) or M (M/2 for a Majorana mass), mu, y, z or lambda,
    followed by the multiplet labels.'''

    anchor = multiplets[0]

    if len(multiplets) == 2:
        if isinstance(anchor.spin, Scalar):
            kind = 'scalar bilinear'
        else: kind = 'fermionic bilinear'
    elif len(multiplets) == 3: 
        if any(H in multiplets for H in (HiggsDoublet, HiggsDoublet.conjugateMultiplet)):
            if isinstance(anchor.spin, Scalar):
                kind = 'scalar cubic'
            else: kind = 'yukawa'
        else: kind = 'dark yukawa'
    elif len(multiplets) == 4:
        kind = 'quartic'
    else: raise NotImplementedError(f'no default parameter name for {multiplets}')

    if kind == 'scalar bilinear':
        # mass term of a multiplet with its conjugate: the mass, squared
        declared = max(multiplets, key=lambda m: m.rep)
        return sp.Symbol(f'M{declared.label}')**2 / (2 if declared.selfConjugate else 1)

    complex_bilinear_label = multiplets[1].label if multiplets[1].label != anchor.label else ""
    default_symbols = {
        'fermionic bilinear': f'M{anchor.label}{complex_bilinear_label}', 
        'scalar cubic': f'mu{"".join(m.label for m in multiplets[:2])}',
        'yukawa': f'y{"".join(m.label for m in multiplets[:2])}',
        'dark yukawa': f'z{"".join(m.label for m in multiplets)}',
        'quartic': f'lambda{"".join(m.label for m in multiplets[:2])}'}

    self_conj_factor = 1
    if kind == 'fermionic bilinear' and anchor.selfConjugate:
        self_conj_factor = 2

    return sp.Symbol(default_symbols[kind]) / self_conj_factor


# ----------------------------------------------------------------------
# MultipletModel
# ----------------------------------------------------------------------

class MultipletModel:
    '''A BSM sector: the multiplets (each one family of a multiplet representation, with
    labeled fields) and their couplings. interaction_graph has the labeled fields as nodes
    and the vertices of the BSM sector as edges ({source, smfields, coeff, param}); the
    mass sector (edges whose SM fields are all scalars, hence carry the vev) and the
    family- and representation-level multiplet graphs are derived from it.'''

    def __init__(self,
                 multiplets: Sequence[Multiplet],
                 particle_model: Optional[ParticleModel] = None):

        multiplets = list(multiplets)
        for m in multiplets:
            if not m.is_filled():
                raise ValueError(
                    f'MultipletModel requires filled multiplets (no missing '
                    f'fields); {m.rep} has an empty component.')

        self.particle_model = particle_model
        self.multiplets = multiplets
        self.fields = [f for m in multiplets for f in m.fields]
        self.counter = Counter(m.rep for m in multiplets)
        self._families = {rep: [m for m in self.multiplets if m.rep == rep]
                          for rep in self.counter}
        # Field labels must be unique within each multiplet's families, so that
        # field_symbol is well-defined (the builder appends the family index).
        for rep, fams in self._families.items():
            if len(fams) > 1:
                for i3 in rep.isospin.i3_range(reverse=True):
                    labels = [m.components[i3].label for m in fams]
                    if len(set(labels)) != len(labels):
                        raise ValueError(
                            f'Families of {rep} share a component label '
                            f'({labels[0]!r}); disambiguate with the family '
                            f'index (see multiplet_model_from_partition).')

        # Complex scalars always come with their conjugate (guaranteed by the
        # partition builder); only one of the two is declared in FeynRules.
        reps = set(self.counter)
        for m in multiplets:
            if (isinstance(m.rep.spin, Scalar) and not m.selfConjugate
                    and m.rep.conjRep not in reps):
                raise ValueError(f'Complex scalar multiplet {m.label} ({m.rep}) '
                                 f'without its conjugate in the model')

        self.sm_multiplets = list(SM_MULTIPLETS)

        # Reverse mapping: field → its parent multiplet
        self._field_to_multiplet = {f: m for m in multiplets for f in m.fields}

    # -- field access -----------------------------------------------------

    @property
    def operators(self) -> List[Operator]:
        '''The renormalizable operators: two BSM multiplets and up to two SM ones, one
        operator per independent SU(2) invariant (suffix P<k>). Operators whose singlet
        vanishes or whose color contraction is not implemented are dropped.'''
        # Cache operators to avoid recalculating
        if not hasattr(self, '_operators'):

            self._operators = []

            for bsm_multiplets in itertools.combinations_with_replacement(
                self.multiplets, 2):

                m1, m2 = bsm_multiplets

                for sm_amount in range(3):

                    if (sm_amount == 0 
                        and m1.rep == m2.rep.conjInternalRep 
                        and m1.index != m2.index):
                        continue

                    for sm_multiplets in itertools.combinations_with_replacement(
                        self.sm_multiplets, sm_amount):

                        multiplets = bsm_multiplets + sm_multiplets
                        parameter = default_operator_parameter(multiplets)

                        # one operator per independent SU(2) invariant of these fields
                        su2 = SU2CGTensor([m.isospin for m in multiplets])
                        paths = su2_singlet_paths(su2.isospins) if su2.dims else []
                        candidates = ([Operator(*multiplets, parameter=parameter)] if len(paths) < 2 else
                                      [Operator(*multiplets, parameter=sp.Symbol(f'{parameter}P{k}'), path=path)
                                       for k, path in enumerate(paths, 1)])

                        # A given operator might contain a singlet, but have its
                        # Clebsch-Gordan coefficients equaling 0. This is not a valid operator.
                        kept = [op for op in candidates if op.dim <= 4 and op.contains_singlet]
                        # An operator without an implemented color contraction cannot be
                        # written, so it is not part of the model.
                        for op in kept:
                            if op.su3_tensor.is_notimplemented:
                                warnings.warn(f'{op.pretty_print()}: color contraction '
                                              f'{tuple(m.colordim for m in op.multiplets)} '
                                              f'not implemented; operator dropped.')
                        kept = [op for op in kept if not op.su3_tensor.is_notimplemented]
                        if len(kept) == 1:
                            kept[0].parameter = parameter   # single invariant: no suffix
                        self._operators.extend(kept)

        return self._operators

    @property
    def lagrangian_operators(self) -> List[Operator]:
        '''The operators written in the Lagrangian: of each pair O, O^dagger only the
        first in self.operators is kept (the other is its h.c.).'''
        if not hasattr(self, '_lagrangian_operators'):
            def key(ms, conj=False):
                return tuple(sorted((m.isSM, (m.rep.conjRep if conj else m.rep), m.index)
                                    for m in ms))
            seen, kept = set(), []
            for op in self.operators:
                own, conj = key(op.multiplets), key(op.multiplets, conj=True)
                if conj != own and conj in seen:
                    continue
                seen.add(own)
                kept.append(op)
            self._lagrangian_operators = kept
        return self._lagrangian_operators

    @property
    def interactions_graph(self) -> Graph:
        '''Fields joined by a coupling. Each edge holds in `emits` what the coupling can
        emit at tree level: 'W' for the gauge link of neighboring components, 'h' for
        a term with a neutral Higgs component (v + h)/sqrt(2), and the other SM fields
        of the term (e.g. tR of a Yukawa). Terms with a charged Higgs component (the
        Goldstone eaten by the W) emit nothing.'''
        # Cache interactions_graph to avoid recalculating
        if not hasattr(self, '_interactions_graph'):
            graph = Graph()

            def emit(f1, f2, what):
                if not graph.has_edge(f1, f2):
                    graph.add_edge(f1, f2, emits=set())
                graph.edges[f1, f2]['emits'] |= what

            for m in self.multiplets:
                for f1, f2 in itertools.pairwise(m.fields):
                    graph.add_nodes_from((f1, f2))
                    emit(f1, f2, {'W'})

            neutral = {HiggsDoublet.components[U1Rep(-1, 2)],
                       HiggsDoublet.conjugateMultiplet.components[U1Rep(1, 2)]}
            higgs = set(HiggsDoublet.fields) | set(HiggsDoublet.conjugateMultiplet.fields)
            for op in self.operators:
                for fop in op.expand_singlets(at_vev=True):
                    f1, f2, *sm = fop.fields
                    emit(f1, f2, set())
                for fop in op.expand_singlets():
                    f1, f2, *sm = fop.fields
                    if fop.at_vev() is None:        # charged Higgs component
                        continue
                    others = {f for f in sm if f not in higgs}
                    emit(f1, f2, others or ({'h'} if any(f in neutral for f in sm) else set()))

            self._interactions_graph = graph

        return self._interactions_graph

    def _emits(self, mother, daughter, sm: Particle) -> bool:
        '''A coupling of the fields of two BSM particles (mother, daughter: sets of
        fields) emits the SM particle sm at tree level: W from a gauge link, Z if the
        shared fields have different T3 (its coupling is diagonal in fields and Q is
        common), h and SM fermions from an operator term; photon and gluon never.'''
        pdg = abs(sm.pdg)
        if pdg in (21, 22):
            return False
        if pdg == 23:
            return len({f.i3 for f in mother & daughter}) > 1
        if pdg == 24:
            wanted = {'W'}
        elif pdg == 25:
            wanted = {'h'}
        else:
            from uvbuilder.states.SMparticles import SMparticlesDict
            fields = {t.field.label for t in SMparticlesDict[pdg].fieldExpansion}
            wanted = None
        for f in mother:
            for g in daughter:
                if not self.interactions_graph.has_edge(f, g):
                    continue
                emits = self.interactions_graph.edges[f, g]['emits']
                if wanted is not None and emits & wanted:
                    return True
                if wanted is None and any(isinstance(e, Field) and e.label.removesuffix('bar')
                                          in fields for e in emits):
                    return True
        return False

    @property
    def mass_sector(self) -> Graph:
        '''Graph of the fields joined by a mass term, couplings to the Higgs vev included.'''
        # Cache mass_sector to avoid recalculating
        if not hasattr(self, '_mass_sector'):
            graph = Graph()
            for op in self.operators:
                if op.forms_bilinear:
                    for fop in op.expand_singlets():
                        f1, f2, *_ = fop.fields
                        if fop.at_vev() is not None:
                            graph.add_nodes_from([f1, f2])
                            graph.add_edge(f1, f2)
            self._mass_sector = graph

        return self._mass_sector

    @property
    def mass_clusters(self) -> List[Set[Field]]:
        '''Connected components of mass_sector: the sets of fields that can mix.'''
        # Cache mass_clusters to avoid recalculating
        if not hasattr(self, '_mass_clusters'):
            self._mass_clusters = list(
                nx.connected_components(self.mass_sector))
        return self._mass_clusters

    @property
    def interaction_clusters(self) -> List[Set[Field]]:
        '''Connected components of interactions_graph.'''
        # Cache interactions_clusters to avoid recalculating
        if not hasattr(self, '_interactions_clusters'):
            self._interactions_clusters = list(
                nx.connected_components(self.interactions_graph))
        return self._interactions_clusters

    def _match_clusters(self, clusters: List[Set[Field]], particles: List[Particle]):
        '''Recursively assigns mass eigenstates (Particles) to gauge-basis mass blocks.
        
        Matches whole particles to clusters based on their unbroken representations.
        Treating particles as rigid bundles prevents unphysical cross-mixing between
        components of different particles across disconnected mass sectors.
        '''
        if not clusters:
            yield []
            return

        cluster = clusters[0]
        cluster_reps = Counter(f.rep.unbroken for f in cluster)
        
        # A cluster of size N can be formed by at most N particles
        max_k = min(len(cluster), len(particles))
        
        for k in range(1, max_k + 1):
            for combo in itertools.combinations(particles, k):
                
                combo_components = []
                for p in combo:
                    # Extract independent gauge-basis components (folds Majorana/Real DoFs)
                    if isinstance(p.rep.spin, DiracSpinor):
                        if p.rep.is_realInternalRep:
                            indep_comps = set(c for c in p.components if c.spin == WEYL_LEFT)
                        else: indep_comps = set(c.dagger if c.spin == WEYL_RIGHT else c for c in p.components)
                    else:
                        indep_comps = set(c for c in p.components)
                    combo_components.extend(indep_comps)
                    
                # Accept the assignment if the combined particle DoFs match the cluster
                if Counter(c.rep for c in combo_components) == cluster_reps:
                    remaining_particles = [p for p in particles if p not in combo]
                    
                    for rest in self._match_clusters(clusters[1:], remaining_particles):
                        yield [(cluster, combo)] + rest


    def _to_mixings(self, 
                   matchings: List[Tuple[Set[Field], Tuple[Particle]]]) -> ParticleMixing:
        '''
        Converts a cluster assignment into a ParticleMixing BiDict.
        Maps a physical Paticle to Gauge Field components.
        '''
        mixing = ParticleMixing()

        for match in matchings:
            fields, particles = match
            for p in particles:
                mixing[p] = [
                    f for f in fields 
                    if f.rep.unbroken == p.components[0].rep
                    or f.rep.unbroken.conjInternalRep == p.components[0].rep] 

        return mixing


    def valid_decays(self, mixing: ParticleMixing) -> List[ParticleVertex]:
        '''Decay vertices realized by `mixing`. A vertex with two BSM particles
        needs a gauge coupling or an operator joining a field of each (tree level);
        one with more BSM particles needs them all in one interaction cluster.
        The emitted SM particles are not checked. Vertices with fewer than two
        BSM particles (X -> SM SM) are never realized: every operator contains a
        pair of BSM fields (implicit Z2), so they are skipped here and rejected
        explicitly in `valid_mixing`.'''

        valid = []

        for vertex in self.particle_model.vertices:
            # the model's particle for a conjugated daughter (the mixing is keyed by particles)
            bsm_particles = [self.particle_model.pdgDict[abs(p.pdg)]
                             for p in vertex.getParticles() if not p.isSM]

            if len(bsm_particles) < 2: continue

            if len(bsm_particles) == 2:
                mother, daughter = (mixing[p] for p in bsm_particles)
                sm = [p for p in vertex.getParticles() if p.isSM]
                if len(sm) == 1 and self._emits(mother, daughter, sm[0]):
                    valid.append(vertex)
                continue

            for cluster in self.interaction_clusters:
                if all(mixing[p].issubset(cluster) for p in bsm_particles):
                    valid.append(vertex)
        
        return valid
    

    def valid_mixing(self, mixing: ParticleMixing) -> bool:
        '''Valid iff (i) every BSM particle is assigned, (ii) particles of the same
        rep (one physical FeynRules class) sit in the same mixing block, so the
        block's column index is the class family index, and (iii) every decay
        vertex is realized by a tree-level coupling (which excludes decays with
        fewer than two BSM particles, see `valid_decays`).'''
        # Z2-even decays (X -> SM SM) cannot come from BSM-pair operators
        if any(sum(not p.isSM for p in v.getParticles()) < 2
               for v in self.particle_model.vertices):
            return False
        bsm = self.particle_model.bsm_particles
        one_block_per_class = all(p in mixing for p in bsm) and all(
            mixing[p] == mixing[q]
            for p, q in itertools.combinations(bsm, 2) if p.rep == q.rep)
        return (one_block_per_class
                and len(self.valid_decays(mixing)) == len(self.particle_model.vertices))

    @property
    def possible_mixings(self) -> List[List[ParticleMixing]]:
        '''The assignments of the BSM particles to the mixing blocks that pass valid_mixing.'''
        matched_clusters = self._match_clusters(
            self.mass_clusters, self.particle_model.bsm_particles)

        all_mixings = []

        for possibility in matched_clusters:
            mixings = self._to_mixings(possibility)
            if self.valid_mixing(mixings): all_mixings.append(mixings)

        return all_mixings

    @property
    def is_valid(self) -> bool:
        '''True if a mixing realizes every decay with tree-level couplings
        (valid_mixing). The masses are not checked: a valid model may not reproduce
        the protomodel spectrum, e.g. a doublet pair alone has a single mass
        parameter, so its neutral and charged states are degenerate.'''
        return any(self.possible_mixings)

    def field_symbol(self, field: Field) -> sp.Symbol:
        '''Sympy symbol for the field (its label, which carries the family index
        when the multiplet appears more than once).'''
        return sp.Symbol(field.label)

    def to_feynrules(self, tracker, mixings: ParticleMixing) -> Dict[str, List[str]]:
        '''FeynRules chunk for the whole gauge sector: the multiplet classes
        (W blocks + kinetic terms), the Dirac pairs (right-handed partners
        and F classes) and the operators (interaction terms, couplings, CG
        tensor definitions).'''

        from uvbuilder.tools.fr_utils import merge_chunks, formatClassBlock

        chunks = []
        emitted = set()

        operators = self.lagrangian_operators

        phases = self.mass_phases(operators, mixings, tracker)
        tracker.mass_phases = {p: symbol for p, (symbol, _) in phases.items()}
        chunks.append({'Parameters': [parameter for _, parameter in phases.values()]})

        for multiplet in self.non_conjugated_multiplets:
            rh = self.rh_conjugate_symbol(multiplet) if multiplet.is_color_conjugated else None
            chunks.append(multiplet.to_feynrules(tracker, mixings, rh_symbol=rh))
            emitted.add(multiplet.fr_symbol())

        # Conjugate partners with SU(2) dim > 1 enter as psi^c = C psi-bar of the
        # declared partner (a tilde class), not as the bare psi-bar.
        reps = {m.rep for m in self.multiplets}
        declared = {m.label: m for m in self.non_conjugated_multiplets}
        conjugate_symbols = {}
        for m in self.multiplets:
            if m.isospin.dim > 1 and is_conjugate_partner(m.rep, reps):
                partner = declared[toggle_dagger(m.label)]
                name = partner.fr_symbol() + TILDE_SUFFIX
                conjugate_symbols[m] = name
                chunks.append({'Classes': [tilde_class(partner, name, tracker)]})
        tracker.conjugate_symbols = conjugate_symbols

        for op in operators:
            chunks.append(op.to_feynrules(tracker))

        # H~ class, emitted once, only if some operator uses the conjugate Higgs
        # and the SM model file does not define it already (SMMod.fr does)
        sm_map = getattr(tracker, 'sm', None)
        sm_defines = sm_map is not None and sm_map.defines_class(HIGGS_TILDE)
        if not sm_defines and any(HiggsDoublet.conjugateMultiplet in op.multiplets
                                  for op in self.operators):
            chunks.append({'Classes': [tilde_class(HiggsDoublet, HIGGS_TILDE, tracker)]})

        return merge_chunks(*chunks)

    @staticmethod
    def mass_phases(operators, mixings, tracker) -> dict:
        '''{particle: (symbol, FRParameter)} for each fermion mass eigenstate alone in
        its block. Its Lagrangian mass, -m chi eta + h.c. (Dirac) or -1/2 m chi chi + h.c.
        (Majorana), is m = +-M with M a (real) mass parameter; the phase P = sign(m) (on
        the right-handed partner) or P = sqrt(sign(m)) (on the field) makes the mass
        positive. With the convention of the mass terms (Multiplet.mass_partner) P = 1
        for M > 0, except for Dirac pairs whose SU(2) singlet fixes the opposite sign.
        Massless states get no phase.'''
        from uvbuilder.tools.fr_utils import FRParameter
        blocks = {f: mm for mm in mixings.get_mixing_matrices() for f in mm.f_components}
        rename = {sp.Symbol(old): sp.Symbol(new)
                  for old, new in getattr(tracker, 'parameter_names', {}).items()}
        rename[sp.Symbol('v')] = sp.Symbol('vev')

        def owner(field):
            '''(particle, conjugate side) of a Weyl field expanding into a single state.'''
            block = blocks.get(field)
            if block is None or block.dim != 1 or not isinstance(field.rep.spin, WeylSpinor):
                return None
            (particle,) = mixings.inverse[field]
            return particle, field.rep.unbroken != particle.components[0].rep

        masses = {}
        for op in operators:
            for term in op.expand_singlets(at_vev=True) + op.hermitian_conjugate(at_vev=True):
                if len(term.fields) != 2:
                    continue
                owners = [owner(f) for f in term.fields]
                if None in owners or owners[0][0] is not owners[1][0]:
                    continue
                particle = owners[0][0]
                sides = {side for _, side in owners}
                majorana = particle.selfConjugate and sides == {False}
                if not (majorana or (not particle.selfConjugate and sides == {False, True})):
                    continue
                if op.sector != 'Mass':
                    raise ValueError(f'{particle} alone in its block gets a mass from {op}')
                coeff = term.coeff.xreplace(rename)
                masses[particle] = masses.get(particle, 0) + (2 * coeff if majorana else coeff)

        phases = {}
        for particle, mass in masses.items():
            mass = sp.simplify(mass)
            if mass == 0:
                continue
            sign, parameter = mass.as_coeff_Mul()
            if not isinstance(parameter, sp.Symbol):
                raise ValueError(f'mass of {particle} is not one mass parameter: {mass}')
            # sign(m) as m/|m|: FeynRules exports Abs and Sqrt to the UFO, not Sign
            value = f'{"-" if sign < 0 else ""}{parameter}/Abs[{parameter}]'
            value = f'Sqrt[{value}]' if particle.selfConjugate else value
            symbol = f'P{particle.label}'
            kind = 'Majorana' if particle.selfConjugate else 'Dirac'
            phases[particle] = (symbol, FRParameter(symbol, 'Phase', (
                ('ParameterType', 'Internal'), ('Value', value),
                ('ComplexParameter', True),
                ('Description', f'"Phase of {particle.label}: positive {kind} mass"')),
                (symbol, 0)))
        return phases

    def rh_conjugate_symbol(self, m: Multiplet) -> str:
        '''FeynRules name of the RH conjugate of a color-conjugated multiplet: the
        default label of the conjugate rep + Rw, e.g. ChiC3bI1Ym4Bm1 -> ChiC3I1Y4B1Rw.'''
        family = m.index if self.counter[m.rep] > 1 else None
        return fr_component_label(m.rep.conjRep, family_index=family) + RIGHT_SUFFIX

    def generate_mixing_matrices(self, tracker=None) -> list:
        '''FeynRules mixing/mass definitions between physical (mass) and
        gauge (flavor) eigenstates.

        Placeholder: the mixing machinery will live here (the index
        ``tracker`` is duck-typed, as in ``Multiplet.generate_fr_kinetic``).
        '''
        return []

    @property
    def non_conjugated_multiplets(self) -> List[Multiplet]:
        '''
        Removes Multiplets that can be achieved by Complex Conjugation of another Multiplet. 
        Keeps the Multiplet with "largest rep", defined by the canonical_order in the Reps 
        classes. Spinors are taken to be Left Weyl, and as such, are never removed from 
        the list.
        '''

        reps = {m.rep for m in self.multiplets}
        keep = [m for m in self.multiplets if not is_conjugate_partner(m.rep, reps)]
        return sorted(keep, key=lambda m: m.rep)


    def pretty_print(self) -> str:
        '''Summary of the multiplets.'''

        string = f'Multiplet Model with {len(self.counter)} Multiplets:'
        indentation = 4
        for rep, count in self.counter.items():
            string += f'\n{indentation*' '}{count} x {rep}'

        return string

    def __repr__(self):
        return self.pretty_print()



TILDE_SUFFIX = 'tilde'


def tilde_class(multiplet: Multiplet, name: str, tracker) -> str:
    '''Unphysical class for psi^c = C psi-bar, the conjugate of a declared multiplet in
    the standard basis: component -m is (-1)^(j+m) times the conjugate of component m,
    e.g. H~ = i tau^2 H*: {Phitilde[1] :> Phibar[2], Phitilde[2] :> -Phibar[1]}.'''
    conj = multiplet.conjugateMultiplet
    bar = multiplet.fr_symbol() + DAGGER_SUFFIX
    indices = tracker.gauge(multiplet.rep)
    colored = multiplet.colordim.dim > 1
    defs = []
    for field in conj.fields:
        sign = '-' if conj.phase(field.i3) == -1 else ''
        lhs = (['c_'] if colored else []) + [str(getSU2Index(field))]
        rhs = (['c'] if colored else []) + [str(getSU2Index(field.conjugateField))]
        defs.append(f'{name}[{", ".join(lhs)}] :> {sign}{bar}[{", ".join(rhs)}]')
    blockDict = {
        'ClassName': name,
        'QuantumNumbers': fr_quantum_numbers({'Y': conj.rep.hypercharge}, conj.rep.globalU1s),
        'SelfConjugate': 'False',
        'Unphysical': 'True',
        'Indices': f"{{{','.join(str(i) for i in indices)}}}",
        'FlavorIndex': indices[-1].name,
        'Definitions': '{' + ', '.join(defs) + '}',
    }
    return formatClassBlock(f"S[{tracker.class_number('S')}]", blockDict)


def is_conjugate_partner(rep: BrokenRep, reps) -> bool:
    '''Complex scalar rep whose (larger) conjugate rep is also present: it is
    named after, and declared in FeynRules through, that partner (its symbol is
    the partner's automatic FeynRules conjugate <label>bar).'''
    return isinstance(rep.spin, Scalar) and rep.conjRep in reps and rep < rep.conjRep


# ----------------------------------------------------------------------
# partition -> MultipletModel
# ----------------------------------------------------------------------

def multiplet_model_from_partition(
        partition: Sequence[BrokenRep],
        labeler: Callable[[BrokenRep, Optional[U1Rep], Optional[int]], str] = default_component_label,
        particle_model: Optional[ParticleModel] = None) -> MultipletModel:
    """Build a labeled MultipletModel from a rep-level partition.

    Each BrokenRep (multiplet rep) becomes one Multiplet; repeated reps become
    families and get a '_1', '_2', ... family index appended to every component
    label. Component labels come from labeler(rep, i3, family_index) and every
    slot is filled (fillMultiplet=True), so the result always passes the filled check.

    Labels are validated before construction: no BSM field label may clash with
    an SM field label, a particle label, or another BSM field label, and
    identical multiplets must carry a family index.
    """
    counts, reps = Counter(partition), set(partition)
    seen: Counter = Counter()

    generated: List[Optional[Tuple[str, Dict[U1Rep, str], Optional[int]]]] = [None] * len(partition)
    family_of: List[Optional[int]] = []

    # Pass 1: the labeler names everything except conjugate scalar partners
    for n, rep in enumerate(partition):
        if not rep.is_multiplet_rep:
            raise ValueError(
                f'Partition entries must be multiplet reps (i3=None), got {rep}')
        family_index = seen[rep] + 1 if counts[rep] > 1 else None
        seen[rep] += 1
        family_of.append(family_index)
        if is_conjugate_partner(rep, reps):
            continue
        generated[n] = (labeler(rep, family_index=family_index),
                        {i3: labeler(rep, i3, family_index)
                         for i3 in rep.isospin.i3_range(reverse=True)},
                        family_index)

    # Pass 2: conjugate partners take the conjugate names of the same-family
    # partner (component i3 <- partner component -i3). Names only: no phase
    # signs, the partner stays an independent multiplet on the Python side.
    kept = {(partition[n], g[2]): g for n, g in enumerate(generated) if g is not None}
    for n, rep in enumerate(partition):
        if generated[n] is None:
            k_label, k_labels, family_index = kept[(rep.conjRep, family_of[n])]
            generated[n] = (toggle_dagger(k_label),
                            {i3: toggle_dagger(k_labels[-i3])
                             for i3 in rep.isospin.i3_range(reverse=True)},
                            family_index)

    _validate_labels(partition, generated, particle_model)

    # Early, warn-only name check (the FeynRules writer does the renaming)
    from uvbuilder.tools.fr_names import warn_labels, multiplet_symbols
    from uvbuilder.tools.fr_utils import chirality_suffix
    own = [(label, rep) for (label, _, _), rep in zip(generated, partition)
           if not is_conjugate_partner(rep, reps)]
    warn_labels([label + (chirality_suffix(rep.spin) if isinstance(rep.spin, WeylSpinor) else '')
                 for label, rep in own], multiplet_symbols, 'Multiplet')

    multiplets = []
    for (label, labels, family_index), rep in zip(generated, partition):
        m = Multiplet(rep=rep, label=label, labels=labels, fillMultiplet=True,
                      index=partition[:len(multiplets)].count(rep) + 1)
        if is_conjugate_partner(rep, reps):
            # enters the operators as psi^c = C psi-bar of its partner: same phases
            m._phases = {i3: su2_conj_phase(rep.isospin, -i3) for i3 in m._phases}
        multiplets.append(m)
    return MultipletModel(multiplets, particle_model)


# ----------------------------------------------------------------------
# Label validation
# ----------------------------------------------------------------------

def _validate_labels(
    partition: Sequence[BrokenRep],
    generated: List[Tuple[str, Dict[U1Rep, str], Optional[int]]],
    particle_model: Optional[ParticleModel]
) -> None:
    """Check that generated BSM labels don't clash with SM or particle labels.

    Parameters
    ----------
    partition : the rep-level partition (one entry per multiplet).
    generated : list of (multiplet_label, {i3: component_label}, family_index)
                tuples, one per entry in *partition*, already produced by the
                labeler.
    particle_model : optional ParticleModel whose SM and BSM particle labels
                     must also be avoided.

    Raises
    ------
    ValueError  with a descriptive message on the first clash found.
    """
    # --- collect every BSM component label the labeler produced ---------------
    bsm_component_labels: Set[str] = set()
    for _, comp_labels, _ in generated:
        bsm_component_labels.update(comp_labels.values())

    # --- 1. clash with SM fields (multiplet labels + component labels) -------
    sm_labels: Set[str] = set()
    for m in SM_MULTIPLETS:
        sm_labels.add(m.label)
        for f in m.fields:
            sm_labels.add(f.label)
            sm_labels.add(f.conjugateField.label)

    sm_clashes = bsm_component_labels & sm_labels
    if sm_clashes:
        raise ValueError(
            f'BSM field labels clash with SM field labels: '
            f'{sorted(sm_clashes)}')

    # --- 2. clash with particle labels (SM + BSM, incl. anti-particles) ------
    particle_labels: Set[str] = set()
    if particle_model is not None:
        for p in particle_model.getParticles():
            particle_labels.add(p.label)
            if p.antiParticle is not None:
                particle_labels.add(p.antiParticle.label)
            if p.fieldExpansion is not None:
                for fterm in p.fieldExpansion:
                    particle_labels.add(fterm.field.label)
                    particle_labels.add(fterm.field.conjugateField.label)

    part_clashes = bsm_component_labels & particle_labels
    if part_clashes:
        raise ValueError(
            f'BSM field labels clash with particle labels: '
            f'{sorted(part_clashes)}')

    # --- 3. exact replicas must carry a family index -------------------------
    # If two entries produce identical (rep, component_labels), the labeler
    # received a family_index but ignored it.
    seen: Dict[Tuple, Optional[int]] = {}
    for i, (_, comp_labels, family_index) in enumerate(generated):
        key = (partition[i], tuple(sorted(comp_labels.items())))
        if key in seen:
            raise ValueError(
                f'Multiplet {partition[i]} appears more than once with '
                f'identical component labels; the labeler received '
                f'family_index={family_index} but did not use it to '
                f'disambiguate.')
        seen[key] = family_index

    # --- 4. internal uniqueness (no two BSM fields share a label) -------------
    all_comp_labels = [l for _, cl, _ in generated for l in cl.values()]
    dupes = {l for l, c in Counter(all_comp_labels).items() if c > 1}
    if dupes:
        raise ValueError(
            f'Duplicate BSM component labels across multiplets: '
            f'{sorted(dupes)}')
