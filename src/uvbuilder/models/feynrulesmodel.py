#!/usr/bin/env python3

from __future__ import annotations  # To avoid "ParticleModel" string
from typing import List, Dict, Tuple, Union, Optional
from itertools import groupby
from contextlib import contextmanager
from dataclasses import replace
from datetime import date
import os
import re
import warnings

from uvbuilder.tools.reps import (U1Rep, SU2Rep, SU3Rep, DiracSpinor, Scalar, Vector,
                        WeylSpinor, BrokenRep, UnbrokenRep, WEYL_LEFT, WEYL_RIGHT)
from uvbuilder.models.multipletModel import MultipletModel, Multiplet
from uvbuilder.models.particleModel import ParticleModel
from uvbuilder.tools.operators import canonical_sort
from uvbuilder.tools.fr_utils import (SECTORS, FRChunk, FRParameter, merge_chunks, formatVariable,
                            formatClassBlock, chirality_suffix, toggle_dagger)
from uvbuilder.tools.fr_names import (DEFAULT_SM_FILE, Renamer, SMClassMap, load_model,
                            reserved_names, particle_symbols, multiplet_symbols,
                            parameter_symbols)
from uvbuilder.states.particle import Particle


class FRIndex(object):
    """
    Describes a single FeynRules Index and its style.
    """

    _instances = {}

    def __new__(cls, name, **kwargs) -> FRIndex:
        if name in FRIndex._instances:
            return FRIndex._instances[name]
        newIndex = super(FRIndex, cls).__new__(cls)
        newIndex.__init__(name, **kwargs)
        FRIndex._instances[name] = newIndex
        return newIndex

    def __init__(self, name: str, range: int, index_type):

        self.name = name
        self.range = range
        self.index_type = index_type
        self.full_name = f"Index[{self.name}]"
        if index_type == SU3Rep:
            self.style = 'c'
        elif index_type == SU2Rep:
            self.style = 'i'
        elif index_type == Particle:
            self.style = 'A'
        elif index_type == Multiplet:
            self.style = 'I'
        elif index_type == DiracSpinor:
            self.style = 'sp'
        elif index_type == WeylSpinor:
            self.style = 'sp'
        elif index_type == Vector:
            self.style = 'mu'
        elif index_type == Scalar:
            self.style = None
        else:
            self.style = 'a'

    def __getnewargs_ex__(self):
        # unpickling goes through __new__, which keeps one instance per name
        return (self.name,), {'range': self.range, 'index_type': self.index_type}

    def __repr__(self) -> str:
        return self.full_name

    def __cmp__(self, other) -> int:
        if self.name == other.name:
            return 0
        elif self.name > other.name:
            return 1
        else:
            return -1

    def __eq__(self, other) -> bool:
        return self.__cmp__(other) == 0

    def __neq__(self, other) -> bool:
        return self.__cmp__(other) != 0

    def __hash__(self) -> int:
        return hash(self.name)

    @property
    def is_gauge(self) -> bool:
        return self.index_type in [SU3Rep, SU2Rep]

    @property
    def is_lorentz(self) -> bool:
        return self.index_type in [DiracSpinor, WeylSpinor, Vector, Scalar]

    @property
    def is_flavor(self) -> bool:
        return (not self.is_gauge) and (not self.is_lorentz)

    def to_feynrules(self) -> str:
        return f"IndexRange[{self.full_name}] = Range[{self.range}];"


class IndexTracker(object):
    '''Central index brain for the FeynRulesModel.

    Owns the canonical FRIndex definitions (name/range/type, memoized and
    consistency-checked) and hands out unique contraction labels and class
    numbers on demand. Every object that needs indices (Multiplets,
    operators, classes) receives this tracker as a duck-typed argument.
    Contraction labels are reset per Lagrangian term by ``reset_symbols()``
    (uniqueness is only required within one term); index definitions and
    class numbers stay global.

    ``sm`` (optional ``tools.fr_names.SMClassMap``) tells the chunk providers
    how SM multiplets are named in the SM model file; indices that file
    declares are not re-declared in ``indices_block``.
    '''

    def __init__(self, sm: Optional[SMClassMap] = None):
        self._indices: Dict[str, FRIndex] = {}
        self._symbol_counters: Dict[str, int] = {}
        self._class_counters: Dict[str, int] = {}
        self.sm = sm
        # operator parameter renames (old -> new), set by FeynRulesModel
        self.parameter_names: Dict[str, str] = {}
        # conjugate-partner multiplet -> tilde class symbol, set by MultipletModel
        self.conjugate_symbols: Dict = {}

    # -- canonical index definitions -----------------------------------

    def index(self, name: str, range: int, index_type) -> FRIndex:
        '''Return the canonical FRIndex for name/range/type.'''
        if name in self._indices:
            existing = self._indices[name]
            if (existing.range, existing.index_type) != (range, index_type):
                raise ValueError(
                    f'Conflicting index definitions for {name}: '
                    f'({existing.range}, {existing.index_type}) vs '
                    f'({range}, {index_type})')
            return existing
        index = FRIndex(name, range=range, index_type=index_type)
        self._indices[name] = index
        return index

    def gauge(self, rep) -> List[FRIndex]:
        '''Gauge FRIndex objects for a BrokenRep (colour, then isospin).'''
        indices = []
        if rep.colordim.dim > 1:
            indices.append(self.index(
                formatVariable(rep.colordim), rep.colordim.dim, SU3Rep))
        if hasattr(rep, 'isospin') and rep.isospin.dim > 1:
            indices.append(self.index(
                formatVariable(rep.isospin), rep.isospin.dim, SU2Rep))
        return indices

    def lorentz(self, spin) -> Union[FRIndex, None]:
        '''Lorentz/spin index for a spin type (None for scalars).'''
        if isinstance(spin, Scalar):
            return None
        return self.index(formatVariable(spin), 4, type(spin))

    def family(self, name: str, n: int) -> FRIndex:
        return self.index(name, n, Particle)

    def multiplet(self, n: int) -> FRIndex:
        return self.index('GenMultiplet', n, Multiplet)

    # -- unique labels and numbering -----------------------------------

    def new_symbol(self, style: str) -> str:
        '''Unique contraction label for an index style, e.g. sp1, sp2, i1...'''
        self._symbol_counters[style] = self._symbol_counters.get(style, 0) + 1
        return f'{style}{self._symbol_counters[style]}'

    def reset_symbols(self) -> None:
        '''Start a new Lagrangian term: contraction labels restart (i1, i2, ...)
        so different terms may reuse the same labels.'''
        self._symbol_counters = {}

    def class_number(self, kind: str) -> int:
        '''Next class number for a spin kind (W, F, S, V).'''
        self._class_counters[kind] = self._class_counters.get(kind, 100) + 1
        return self._class_counters[kind]

    # -- output ---------------------------------------------------------

    @property
    def all_indices(self) -> List[FRIndex]:
        return list(self._indices.values())

    def indices_block(self) -> str:
        indicesHeader = "(* ************************** *)\n"
        indicesHeader += "(* *****    Indices     ***** *)\n"
        indicesHeader += "(* ************************** *)\n\n"

        index_lines = []
        style_lines = []
        sm_indices = self.sm.index_names if self.sm is not None else frozenset()
        for index in self.all_indices:
            if index.name in sm_indices:     # declared by the SM file: never override
                continue
            if index.is_lorentz:             # mu, sp: names of built-in Lorentz/spin indices, not index types
                continue
            index_str = index.to_feynrules()
            if index_str in index_lines:
                continue
            index_lines.append(index_str)
            style_lines.append(f"IndexStyle[{index.name},{index.style}];")

        return indicesHeader + "\n".join(index_lines + style_lines)


MODEL_NAME = re.compile(r'[A-Za-z][A-Za-z0-9_]*')


class FeynRulesModel(object):
    '''Upper-level wrapper: assembles the .fr file from the chunks that each
    object produces through its ``to_feynrules(tracker)`` method. The index
    tracker is the central brain passed down to every chunk provider.'''

    def __init__(self, multiplet_model: MultipletModel,
                 particle_model: ParticleModel, mixings,
                 output_file: str, sm_file=DEFAULT_SM_FILE,
                 strict: bool = False, name: str = 'UVtest',
                 author: str = 'Andre Lessa & Lucas Heck', version: str = '1.0') -> None:
        '''name is M$ModelName and the default UFO folder (<name>_UFO). MG5 imports
        that folder as a Python package: the name must be a valid identifier (letters,
        digits, underscores, starting with a letter) and should not be the name of an
        existing Python module.'''

        if not MODEL_NAME.fullmatch(name):
            raise ValueError(f'invalid model name {name!r}: letters, digits and '
                             'underscores, starting with a letter')
        self.name = name
        self.author = author
        self.version = version
        today = date.today()
        self.date = today.strftime("%d.%m.%Y")

        self.output_file = output_file
        self.multiplet_model = multiplet_model
        self.mixings = mixings
        self.particle_model = particle_model
        # SM model file the output is loaded with (default SMMod.fr): SM
        # multiplets are written as its classes, its indices are not
        # re-declared, and its symbols are reserved.
        self.sm_file = sm_file
        self.sm_model = load_model(sm_file)
        self.reserved = reserved_names(sm_file)
        self.index_tracker = IndexTracker(sm=SMClassMap(self.sm_model))

        # Mixing matrices referenced by the emitted (non-conjugated) multiplets
        emitted = {f for m in self.multiplet_model.non_conjugated_multiplets
                   for f in m.fields}
        self._mixing_matrices = [mm for mm in mixings.get_mixing_matrices()
                                 if emitted & set(mm.f_components)]

        # Invalid or reserved labels/parameters are renamed here, for the output
        # only; the model objects keep their labels.
        renamer = Renamer(self.reserved)
        self.renames: Dict[str, Dict[str, str]] = self._plan_renames(renamer)
        messages = [f'{kind} {old} -> {new}' for kind, table in self.renames.items()
                    for old, new in table.items()]
        if messages:
            text = (f'Invalid or reserved names ({os.path.basename(str(sm_file))}) '
                    f'renamed in the FeynRules output: ' + ', '.join(messages))
            if strict:
                raise ValueError(text)
            warnings.warn(text, stacklevel=2)
        self.index_tracker.parameter_names = dict(self.renames['parameter'])

        with self._relabeled():
            self.pm_chunks: FRChunk = self.particle_model.to_feynrules(
                self.index_tracker)
            self.mm_chunks: FRChunk = self.multiplet_model.to_feynrules(
                self.index_tracker, mixings)
            self.mix_chunks: FRChunk = merge_chunks(*(
                mm.to_feynrules(self.index_tracker) for mm in self._mixing_matrices))
            self._chunks: FRChunk = merge_chunks(self.pm_chunks, self.mm_chunks,
                                                 self.mix_chunks)
            self._description = self._render_model_description()
            unfixable = self._check_fixed_names(renamer)
        if unfixable:
            text = ('Generated names that cannot be renamed clash: '
                    + '; '.join(unfixable))
            if strict:
                raise ValueError(text)
            warnings.warn(text, stacklevel=2)

    # -- names -----------------------------------------------------------

    def _partner(self, multiplet):
        '''The conjugate-partner scalar multiplet named after ``multiplet``
        (label + 'bar'), if the model has one.'''
        return next((m for m in self.multiplet_model.multiplets
                     if m is not multiplet and m.label == toggle_dagger(multiplet.label)), None)

    def _plan_renames(self, renamer: Renamer) -> Dict[str, Dict[str, str]]:
        '''Renames for the three kinds of generated roots, in a fixed order:
        particle labels (PDG order), multiplet labels (declaration order),
        operator parameters (alphabetical). Each root keeps its name when all
        derived symbols are valid and free.'''
        self._particle_renames, self._multiplet_renames = {}, {}
        table = {'particle': {}, 'multiplet': {}, 'parameter': {}}

        for p in sorted(self.particle_model.bsm_particles, key=lambda p: p.pdg):
            new = renamer.claim(p.label, particle_symbols, f'particle {p.pdg}')
            if new != p.label:
                self._particle_renames[p] = new
                table['particle'][p.label] = new

        for m in self.multiplet_model.non_conjugated_multiplets:
            suffix = (chirality_suffix(m.rep.spin)
                      if isinstance(m.rep.spin, WeylSpinor) else '')
            new = renamer.claim(m.label, lambda c: multiplet_symbols(c + suffix),
                                f'multiplet {m.label}')
            if new != m.label:
                self._multiplet_renames[m] = new
                table['multiplet'][m.label] = new

        names = sorted({str(s) for op in self.multiplet_model.operators
                        for s in op.parameter.free_symbols})
        for name in names:
            new = renamer.claim(name, parameter_symbols, 'operator parameter')
            if new != name:
                table['parameter'][name] = new
        return table

    @contextmanager
    def _relabeled(self):
        '''Temporarily give particles and multiplets their output names (and keep
        antiparticles, conjugate multiplets and conjugate-partner scalars
        consistent); everything rendered inside picks the new names up. Safe
        for dicts/sets: Particle and Multiplet hash by identity.'''
        saved = []

        def set_label(obj, new):
            saved.append((obj, obj._label))
            obj._label = new

        try:
            for p, new in self._particle_renames.items():
                set_label(p, new)
                anti = p.antiParticle
                if anti is not None and anti is not p:
                    set_label(anti, p.conjugate_label())
            for m, new in self._multiplet_renames.items():
                partner = self._partner(m)
                set_label(m, new)
                if m.conjugateMultiplet is not m:
                    set_label(m.conjugateMultiplet, toggle_dagger(new))
                if partner is not None:
                    set_label(partner, toggle_dagger(new))
                    if partner.conjugateMultiplet is not partner:
                        set_label(partner.conjugateMultiplet, new)
            yield
        finally:
            for obj, old in reversed(saved):
                obj._label = old

    def _check_fixed_names(self, renamer: Renamer) -> List[str]:
        '''Names built from reps rather than labels (multi-member class and
        family names) and mixing-matrix parameters cannot be renamed; check them.'''
        derived = []
        particles = sorted(self.particle_model.bsm_particles, key=lambda p: p.rep)
        for _, group in groupby(particles, key=lambda p: p.rep):
            members = list(group)
            if len(members) > 1:
                container = ParticleModel._container_name(members)
                derived += [(container, True), (container + 'bar', False),
                            (container + 'Lw', False), (container + 'Rw', False),
                            ('M' + container, True), ('W' + container, True),
                            (ParticleModel._family_name(members), False)]
        for m in self.multiplet_model.non_conjugated_multiplets:
            if m.is_colour_conjugated:
                rh = self.multiplet_model.rh_conjugate_symbol(m)
                derived += [(rh, False), (rh + 'bar', False)]
        for mm in self._mixing_matrices:
            if mm.dim == 1:          # single-member blocks emit no matrix
                continue
            name = mm.standard_symbol
            derived += [(name, False), ('R' + name, True), ('I' + name, True),
                        (name.upper(), True), ('IM' + name.upper(), True)]
        return renamer.check(derived, 'fixed names')

    # -- thin wrappers ---------------------------------------------------

    def get_classes_definitions(self) -> str:
        '''M$ClassesDescription block (physical then gauge eigenstates).'''

        classesHeader = "(* ************************** *)\n"
        classesHeader += "(* **** Particle classes **** *)\n"
        classesHeader += "(* ************************** *)\n\n"

        sections = []
        if self.pm_chunks.get('Classes'):
            section = "\n(***  Physical (mass eigenstates) Fields ***)\n"
            section += ',\n'.join(self.pm_chunks['Classes'])
            sections.append(section)
        if self.mm_chunks.get('Classes'):
            section = "\n(***  Multiplet (gauge eigenstates) Fields ***)\n"
            section += ',\n'.join(self.mm_chunks['Classes'])
            sections.append(section)

        frBlock = classesHeader
        frBlock += "M$ClassesDescription = {\n"
        frBlock += ',\n'.join(sections)
        frBlock += "\n};"
        return frBlock

    def parameters(self) -> List[FRParameter]:
        '''M$Parameters in output order: by sector (SECTORS order), External then
        Internal, sorted by name (real and imaginary parts adjacent). OrderBlock
        numbers the scalar External entries within each block; the two parts of a
        complex coupling share their number (BLOCK, IMBLOCK).'''
        parameters = self._chunks.get('Parameters', [])
        ordered = []
        for sector, _, _ in SECTORS:
            for kind in ('External', 'Internal'):
                numbers = {}
                for p in sorted((p for p in parameters if p.sector == sector and p.kind == kind),
                                key=lambda p: p.sort_key):
                    options = dict(p.options)
                    if kind == 'External' and 'BlockName' in options and 'Indices' not in options:
                        number = numbers.setdefault(p.base, len(numbers) + 1)
                        p = replace(p, options=p.options + (('OrderBlock', number),))
                    ordered.append(p)
        return ordered

    def get_parameters(self) -> str:
        '''M$Parameters rendered from self.parameters(), one commented group per
        sector and kind.'''
        header = "(* ************************** *)\n"
        header += "(* *****   Parameters   ***** *)\n"
        header += "(* ************************** *)\n\n"
        titles = {key: title for key, title, _ in SECTORS}
        groups = [f'  (* {titles[sector]}: {kind.lower()} parameters *)\n'
                  + ',\n'.join(p.render() for p in group)
                  for (sector, kind), group in groupby(self.parameters(),
                                                       key=lambda p: (p.sector, p.kind))]
        return header + 'M$Parameters = {\n\n' + ',\n\n'.join(groups) + '\n\n};'

    def get_indices(self) -> str:
        return self.index_tracker.indices_block()

    def _flavor_expand(self) -> str:
        expand = {index.name for index in self.index_tracker.all_indices
                  if index.index_type in (Multiplet, Particle, SU2Rep)}
        # SM SU(2) indices, as in SM.fr's FlavorExpand -> {SU2W, SU2D}: without
        # them the unphysical Wi and Phi are never expanded into W/Z/A and H/G0/GP.
        expand.update({'SU2W', 'SU2D'})
        return f"FlavorExpand -> {{{', '.join(sorted(expand))}}}"

    def get_lagrangian(self) -> str:
        '''LBSM = LBSMKin + LBSMMass + LBSMInt; the interaction terms grouped by
        sector, each group sorted.'''
        indent = 4*' '

        def piece(name, groups, block=None):
            lines = []
            for title, terms in groups:
                if len(groups) > 1:
                    lines.append(f'{indent}(* {title} *)')
                terms = sorted(t.strip() for t in terms)
                lines += [f'{indent}{t if t.startswith(("+", "-")) else "+ " + t}' for t in terms]
            if not lines:
                return f'{name} := 0;'
            body = 'ExpandIndices[\n' + '\n'.join(lines) + f'\n  , {self._flavor_expand()}]'
            return f'{name} := ' + (f'Block[{{{block}}}, {body}];' if block else f'{body};')

        chunks = self._chunks
        interactions = [(title, chunks.get(f'{key}Terms', [])) for key, title, _ in SECTORS[1:5]]
        interactions = [(title, terms) for title, terms in interactions if terms]
        out = ["(* ************************** *)",
               "(* *****   Lagrangian   ***** *)",
               "(* ************************** *)", "",
               "(* Kinetic terms *)",
               piece('LBSMKin', [('', chunks.get('KineticTerms', []))], block='mu'), "",
               "(* Mass terms *)",
               piece('LBSMMass', [('', chunks.get('MassTerms', []))]), "",
               "(* Interactions *)",
               piece('LBSMInt', interactions), "",
               "LBSM := LBSMKin + LBSMMass + LBSMInt;"]
        return '\n'.join(out) + '\n'

    # -- model description ----------------------------------------------

    def get_model_description(self) -> str:
        '''Banner listing the multiplets and physical states in the model
        (rendered once, with the output names).'''
        return self._description

    def _render_model_description(self) -> str:

        comm_left = "(*"
        comm_right = "*)\n"
        sep = comm_left + ("*" * 61) + comm_right
        header = sep + sep
        header += "(***** FeynRules model file automatically generated      *******)\n"
        header += sep
        header += comm_left + f"{" Multiplet States ":*^61}" + comm_right
        header += sep

        data = []
        data.append(['Multiplet', 'SU3C,SU2L,Y,Lorentz', 'Particles'])
        for multiplet in self.multiplet_model.non_conjugated_multiplets:
            rep = [multiplet.rep.colordim, multiplet.rep.isospin,
                   multiplet.rep.hypercharge, multiplet.rep.spin]
            rep = ','.join([str(r) for r in rep])
            data.append([multiplet.fr_symbol(), rep,
                         ', '.join(f.label for f in multiplet.fields)])

        column_widths = [max(len(str(item)) for item in col) + 1
                         for col in zip(*data)]
        header_format = " | ".join(f"{{:<{width}}}" for width in column_widths)
        header += comm_left + ' ' + header_format.format(*data[0]) + comm_right
        header += comm_left + ' ' + ("-" * (sum(column_widths) + (len(column_widths) - 1) * 3)) + comm_right

        for row in data[1:]:
            header += comm_left + ' ' + header_format.format(*row) + comm_right

        header += sep
        header += sep
        header += comm_left + f"{" Physical States ":*^61}" + comm_right
        header += sep

        data = []
        data.append(['States', 'SU3C,Q,Lorentz', 'Particles'])
        particles = sorted(self.particle_model.bsm_particles,
                           key=lambda p: p.rep)
        for _, group in groupby(particles, key=lambda p: p.rep):
            members = list(group)
            rep = [members[0].rep.colordim, members[0].rep.eCharge,
                   members[0].rep.spin]
            rep = ','.join([str(r) for r in rep])
            data.append([', '.join(p.label for p in members), rep,
                         ', '.join(p.label for p in members)])

        column_widths = [max(len(str(item)) for item in col) + 1
                         for col in zip(*data)]
        header_format = " | ".join(f"{{:<{width}}}" for width in column_widths)
        header += comm_left + ' ' + header_format.format(*data[0]) + comm_right
        header += comm_left + ' ' + ("-" * (sum(column_widths) + (len(column_widths) - 1) * 3)) + comm_right

        for row in data[1:]:
            header += comm_left + ' ' + header_format.format(*row) + comm_right

        header += sep
        return header

    def get_model_info(self) -> str:
        '''Basic model information (name, author, ...).'''
        return '\n'.join([f'M$ModelName = "{self.name}";', '',
                          'M$Information = {',
                          f'  Authors -> {{"{self.author}"}},',
                          f'  Version -> "{self.version}",',
                          f'  Date    -> "{self.date}"',
                          '};', '', 'FeynmanGauge = True;'])

    def to_feynrules(self) -> str:
        return '\n\n'.join([self.get_model_description() + '\n',
                            self.get_model_info(), self.get_indices(), self.get_parameters(),
                            self.get_classes_definitions(), self.get_lagrangian()])

    def save2file(self) -> None:

        if os.path.dirname(self.output_file):
            os.makedirs(os.path.dirname(self.output_file), exist_ok=True)

        with open(self.output_file, 'w') as f:
            f.write(self.to_feynrules())
