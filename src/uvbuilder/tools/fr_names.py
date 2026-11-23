#!/usr/bin/env python3
'''FeynRules model-file parser and reserved-name registry.

Two layers:

* ``FRModelFile.parse(path)`` reads a FeynRules ``.fr`` file *statically* (no
  Mathematica needed): comments are stripped, the file is split into top-level
  statements, and the blocks FeynRules cares about are parsed into plain data:
  particle classes (``M$ClassesDescription``), parameters (``M$Parameters``),
  gauge groups (``M$GaugeGroups``), index declarations (``IndexRange`` /
  ``IndexStyle``) and top-level definitions (``LSM := ...``, ``FeynmanGauge = ...``).
  Values that are not needed are kept as raw text.

* ``ReservedNames`` collects every symbol a set of model files defines (class
  names, automatic conjugates ``Xbar``, class members, mass/width symbols,
  parameters and their ``ParameterName`` aliases, block names, indices and index
  styles, gauge-group symbols, top-level definitions, particle-name strings) plus
  the language-level names that must not be redefined (Mathematica / FeynRules
  symbols, Python keywords and builtins, UFO module names).  Each entry remembers
  where it comes from (``file:line`` or list name) so conflicts can be reported.

Comparison rules (see ``ReservedNames.clashes``): names that end up in the UFO /
MG5 (parameters, mass/width symbols, aliases, block names, particle names) are
compared case-insensitively, because MG5 lowercases them; everything else, and
all language-level lists (Mathematica and Python are case-sensitive), is compared
exactly.

The parser is deliberately a subset parser: it understands the structure used by
SM.fr / SMMod.fr / mssm.fr style files (``K[n] == {opt -> value, ...}`` entries,
``name == {...}`` parameters and gauge groups) and ignores everything else.
'''

from __future__ import annotations

import builtins
import itertools
import keyword
import re
import warnings
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

# ----------------------------------------------------------------------
# Symbol validity
# ----------------------------------------------------------------------

VALID_SYMBOL = re.compile(r'[A-Za-z][A-Za-z0-9]*')


def is_valid_symbol(name: str) -> bool:
    '''True if ``name`` is usable as a FeynRules / UFO symbol: a letter followed by
    letters and digits only (no ``_``, ``$``, ``+``, ``-``; no leading digit).'''
    return VALID_SYMBOL.fullmatch(name) is not None


# ----------------------------------------------------------------------
# Low-level scanning (comments, strings, bracket depth)
# ----------------------------------------------------------------------

def strip_comments(text: str) -> str:
    '''Replace Mathematica comments ``(* ... *)`` (nested) by spaces, keeping
    newlines, so that character offsets and line numbers are preserved. Comment
    markers inside strings are left alone.'''
    out = list(text)
    i, n, depth, in_str = 0, len(text), 0, False
    while i < n:
        c = text[i]
        if depth == 0 and in_str:
            if c == '\\':
                i += 2
                continue
            if c == '"':
                in_str = False
            i += 1
            continue
        if text.startswith('(*', i):
            depth += 1
            out[i] = out[i + 1] = ' '
            i += 2
            continue
        if depth > 0:
            if text.startswith('*)', i):
                depth -= 1
                out[i] = out[i + 1] = ' '
                i += 2
                continue
            if c != '\n':
                out[i] = ' '
            i += 1
            continue
        if c == '"':
            in_str = True
        i += 1
    return ''.join(out)


def _depths(text: str) -> List[int]:
    '''Bracket depth of every character ((), [], {}); -1 inside strings.
    An opening bracket gets the outer depth, a closing one the outer depth too.'''
    out: List[int] = []
    depth, in_str, i = 0, False, 0
    while i < len(text):
        c = text[i]
        if in_str:
            out.append(-1)
            if c == '\\' and i + 1 < len(text):
                out.append(-1)
                i += 2
                continue
            if c == '"':
                in_str = False
        elif c == '"':
            in_str = True
            out.append(-1)
        elif c in '([{':
            out.append(depth)
            depth += 1
        elif c in ')]}':
            depth -= 1
            out.append(depth)
        else:
            out.append(depth)
        i += 1
    return out


def _split_top(text: str, sep: str, base: int = 0) -> List[Tuple[int, str]]:
    '''Split ``text`` at top-level occurrences of ``sep``. Returns
    ``(absolute offset, piece)`` pairs; empty pieces are dropped.'''
    depths = _depths(text)
    pieces, start, i = [], 0, 0
    while i < len(text):
        if depths[i] == 0 and text.startswith(sep, i):
            pieces.append((base + start, text[start:i]))
            i += len(sep)
            start = i
            continue
        i += 1
    pieces.append((base + start, text[start:]))
    return [(off, p) for off, p in pieces if p.strip()]


def _find_top(text: str, ops: Sequence[str]) -> Optional[Tuple[int, str]]:
    '''First top-level occurrence of any of ``ops``: (index, op) or None.'''
    depths = _depths(text)
    for i in range(len(text)):
        if depths[i] != 0:
            continue
        for op in ops:
            if text.startswith(op, i):
                return i, op
    return None


_CONTINUATION_START = tuple('+-*/.=><|&^@,;:')
_COMPLETE_END = re.compile(r'[\]\}\)A-Za-z0-9"$]\s*$')


def split_statements(text: str) -> List[Tuple[int, str]]:
    '''Split a comment-free file into top-level statements: at ``;`` and at
    top-level newlines that end a complete expression (Mathematica package
    semantics, approximated: the statement so far ends with a bracket, symbol,
    number or string, and the next line does not start with an operator).'''
    depths = _depths(text)
    out, start, i, n = [], 0, 0, len(text)
    while i < n:
        c = text[i]
        if depths[i] == 0 and c == ';':
            out.append((start, text[start:i]))
            start = i + 1
        elif depths[i] == 0 and c == '\n':
            current = text[start:i]
            rest = text[i + 1:].lstrip()
            if current.strip() and _COMPLETE_END.search(current) and rest \
                    and not rest.startswith(_CONTINUATION_START):
                out.append((start, current))
                start = i + 1
        i += 1
    out.append((start, text[start:]))
    return [(off + len(p) - len(p.lstrip()), p.strip()) for off, p in out if p.strip()]


_STRING = re.compile(r'"(?:[^"\\]|\\.)*"')
_SYMBOL_RE = re.compile(r'(?<![A-Za-z0-9$`\\])([A-Za-z$][A-Za-z0-9$]*)(?![A-Za-z0-9$`_])')


def symbols_in(value: str) -> List[str]:
    '''Symbols occurring in a raw value (strings, numbers, pattern variables such
    as ``sp_`` and ``\\[Alpha]`` named characters are ignored), in order, unique.'''
    blanked = _STRING.sub(lambda m: ' ' * len(m.group()), value)
    blanked = re.sub(r'\\\[[A-Za-z]+\]', ' ', blanked)
    seen: Dict[str, None] = {}
    for m in _SYMBOL_RE.finditer(blanked):
        seen.setdefault(m.group(1), None)
    return list(seen)


def strings_in(value: str) -> List[str]:
    '''Contents of the string literals in a raw value, in order.'''
    return [s[1:-1] for s in _STRING.findall(value)]


def list_items(value: str) -> List[str]:
    '''Top-level items of ``{a, b, ...}``; a non-list value is a one-item list.'''
    v = value.strip()
    if v.startswith('{') and v.endswith('}'):
        return [p.strip() for _, p in _split_top(v[1:-1], ',')]
    return [v] if v else []


def split_rule(item: str) -> Optional[Tuple[str, str, str]]:
    '''``lhs -> rhs`` / ``lhs :> rhs`` at top level: (lhs, op, rhs) or None.'''
    found = _find_top(item, ('->', ':>'))
    if found is None:
        return None
    i, op = found
    return item[:i].strip(), op, item[i + len(op):].strip()


def _head(expr: str) -> Optional[str]:
    '''Head symbol of ``f[...]`` or the symbol itself.'''
    m = re.match(r'\s*([A-Za-z$][A-Za-z0-9$]*)', expr)
    return m.group(1) if m else None


def _fraction(value: str) -> Optional[Fraction]:
    try:
        return Fraction(value.replace(' ', ''))
    except (ValueError, ZeroDivisionError):
        return None


# ----------------------------------------------------------------------
# Parsed objects
# ----------------------------------------------------------------------

@dataclass(frozen=True)
class FRIndex:
    '''An ``IndexRange`` / ``IndexStyle`` declaration.'''
    name: str
    range: Optional[str] = None      # raw right-hand side, e.g. 'Unfold[Range[2]]'
    style: Optional[str] = None
    line: Optional[int] = None       # line of the IndexRange (else of the IndexStyle)

    @property
    def unfold(self) -> Optional[bool]:
        if self.range is None:
            return None
        if self.range.startswith('NoUnfold'):
            return False
        return self.range.startswith('Unfold')

    @property
    def size(self) -> Optional[int]:
        m = re.search(r'Range\[\s*(\d+)\s*\]', self.range or '')
        return int(m.group(1)) if m else None


@dataclass(frozen=True)
class FRClass:
    '''One entry of ``M$ClassesDescription``, e.g. ``F[2] == {...}``.'''
    kind: str                         # V, F, S, U, W (or VSF/CSF in superfield models)
    number: int
    name: str
    line: int
    options: Dict[str, str] = field(default_factory=dict, compare=False, hash=False)

    # --- simple accessors -------------------------------------------------
    def _bool(self, key: str, default: bool = False) -> bool:
        v = self.options.get(key)
        return default if v is None else v.strip() == 'True'

    @property
    def unphysical(self) -> bool:
        return self._bool('Unphysical')

    @property
    def self_conjugate(self) -> bool:
        return self._bool('SelfConjugate')

    @property
    def chirality(self) -> Optional[str]:
        v = self.options.get('Chirality')
        return v.strip() if v else None

    @property
    def members(self) -> Tuple[str, ...]:
        return tuple(symbols_in(self.options.get('ClassMembers', '')))

    @property
    def indices(self) -> Tuple[str, ...]:
        '''Index types in declaration order, e.g. ('SU2D', 'Generation', 'Colour').'''
        return tuple(m.group(1) for m in
                     re.finditer(r'Index\[\s*([A-Za-z][A-Za-z0-9]*)', self.options.get('Indices', '')))

    @property
    def flavor_index(self) -> Optional[str]:
        v = self.options.get('FlavorIndex')
        return v.strip() if v else None

    @property
    def quantum_numbers(self) -> Dict[str, Fraction]:
        out: Dict[str, Fraction] = {}
        for item in list_items(self.options.get('QuantumNumbers', '')):
            rule = split_rule(item)
            if rule and _fraction(rule[2]) is not None:
                out[rule[0]] = _fraction(rule[2])
        return out

    @property
    def weyl_components(self) -> Tuple[str, ...]:
        return tuple(symbols_in(self.options.get('WeylComponents', '')))

    @property
    def pdg(self) -> Tuple[int, ...]:
        return tuple(int(x) for x in re.findall(r'-?\d+', self.options.get('PDG', '')))

    @property
    def particle_names(self) -> Tuple[str, ...]:
        return tuple(strings_in(self.options.get('ParticleName', '')))

    @property
    def antiparticle_names(self) -> Tuple[str, ...]:
        return tuple(strings_in(self.options.get('AntiParticleName', '')))

    @property
    def mass_symbols(self) -> Tuple[str, ...]:
        return tuple(s for s in symbols_in(self.options.get('Mass', ''))
                     if s not in ('Internal', 'External'))

    @property
    def width_symbols(self) -> Tuple[str, ...]:
        return tuple(s for s in symbols_in(self.options.get('Width', ''))
                     if s not in ('Internal', 'External'))

    @property
    def conjugate_names(self) -> Tuple[str, ...]:
        '''Automatic FeynRules conjugates (``Xbar``) of the class and its members.'''
        if self.self_conjugate:
            return ()
        return tuple(n + 'bar' for n in (self.name,) + self.members)


@dataclass(frozen=True)
class FRParameter:
    '''One entry of ``M$Parameters``, e.g. ``MW == {...}``.'''
    name: str
    line: int
    options: Dict[str, str] = field(default_factory=dict, compare=False, hash=False)

    @property
    def external(self) -> bool:
        return self.options.get('ParameterType', '').strip() == 'External'

    @property
    def block(self) -> Optional[str]:
        v = self.options.get('BlockName')
        return v.strip() if v else None

    @property
    def indices(self) -> Tuple[str, ...]:
        return tuple(m.group(1) for m in
                     re.finditer(r'Index\[\s*([A-Za-z][A-Za-z0-9]*)', self.options.get('Indices', '')))

    @property
    def aliases(self) -> Tuple[str, ...]:
        '''``ParameterName`` aliases: ``G`` or ``{yl[1,1] -> ye, ...}`` -> ('ye', ...).'''
        v = self.options.get('ParameterName')
        if v is None:
            return ()
        out = []
        for item in list_items(v):
            rule = split_rule(item)
            out += symbols_in(rule[2] if rule else item)
        return tuple(out)


@dataclass(frozen=True)
class FRGaugeGroup:
    '''One entry of ``M$GaugeGroups``, e.g. ``SU2L == {...}``.'''
    name: str
    line: int
    options: Dict[str, str] = field(default_factory=dict, compare=False, hash=False)

    @property
    def symbols(self) -> Tuple[str, ...]:
        '''Symbols the group introduces: option values (booleans excluded) and the
        heads of its ``Definitions`` rules (e.g. Ta, FSU2L).'''
        out: Dict[str, None] = {}
        for key, value in self.options.items():
            if key == 'Definitions':
                for item in list_items(value):
                    rule = split_rule(item)
                    head = _head(rule[0]) if rule else None
                    if head:
                        out.setdefault(head, None)
                continue
            for s in symbols_in(value):
                if s not in ('True', 'False'):
                    out.setdefault(s, None)
        return tuple(out)


def _parse_options(body: str) -> Dict[str, str]:
    '''``{key -> value, ...}`` (outer braces already removed) -> {key: raw value}.'''
    out: Dict[str, str] = {}
    for _, item in _split_top(body, ','):
        rule = split_rule(item)
        if rule is not None:
            out[rule[0]] = rule[2]
    return out


def _entries(block: str, base: int) -> Iterator[Tuple[int, str, str]]:
    '''Entries ``lhs == {body}`` of a top-level list: (offset, lhs, body).'''
    for off, item in _split_top(block, ',', base):
        found = _find_top(item, ('==',))
        if found is None:
            continue
        lhs, rhs = item[:found[0]].strip(), item[found[0] + 2:].strip()
        if rhs.startswith('{') and rhs.endswith('}'):
            yield off + len(item) - len(item.lstrip()), lhs, rhs[1:-1]


# ----------------------------------------------------------------------
# Model file
# ----------------------------------------------------------------------

@dataclass
class FRModelFile:
    '''Static content of one ``.fr`` file.'''
    source: str
    classes: Dict[str, FRClass] = field(default_factory=dict)
    parameters: Dict[str, FRParameter] = field(default_factory=dict)
    gauge_groups: Dict[str, FRGaugeGroup] = field(default_factory=dict)
    indices: Dict[str, FRIndex] = field(default_factory=dict)
    definitions: Dict[str, int] = field(default_factory=dict)   # top-level name -> line
    skipped: List[Tuple[int, str]] = field(default_factory=list)  # (line, text) not understood

    @classmethod
    def parse(cls, path) -> 'FRModelFile':
        path = Path(path)
        return cls.from_text(path.read_text(encoding='utf-8', errors='replace'), source=path.name)

    @classmethod
    def from_text(cls, text: str, source: str = '<string>') -> 'FRModelFile':
        model = cls(source=source)
        clean = strip_comments(text)
        line_starts = [0] + [m.end() for m in re.finditer('\n', clean)]

        def line(offset: int) -> int:
            lo, hi = 0, len(line_starts)
            while lo + 1 < hi:
                mid = (lo + hi) // 2
                if line_starts[mid] <= offset:
                    lo = mid
                else:
                    hi = mid
            return lo + 1

        for off, stmt in split_statements(clean):
            model._statement(off, stmt, line)
        return model

    # --- statements -------------------------------------------------------
    def _statement(self, off: int, stmt: str, line) -> None:
        m = re.match(r'IndexRange\[\s*Index\[\s*([A-Za-z][A-Za-z0-9]*)\s*\]\s*\]\s*=\s*(.*)$', stmt, re.S)
        if m:
            old = self.indices.get(m.group(1), FRIndex(m.group(1)))
            self.indices[m.group(1)] = FRIndex(m.group(1), m.group(2).strip(), old.style, line(off))
            return
        m = re.match(r'IndexStyle\[\s*([A-Za-z][A-Za-z0-9]*)\s*,\s*([A-Za-z][A-Za-z0-9]*)\s*\]$', stmt)
        if m:
            old = self.indices.get(m.group(1), FRIndex(m.group(1)))
            self.indices[m.group(1)] = FRIndex(m.group(1), old.range, m.group(2), old.line or line(off))
            return

        found = _find_top(stmt, (':=', '='))
        if found is None:            # a bare expression (Do[...], Print[...]): defines nothing we track
            return
        if stmt[found[0]:found[0] + 2] == '==':
            self.skipped.append((line(off), stmt[:80]))
            return
        i, op = found
        lhs, rhs = stmt[:i].strip(), stmt[i + len(op):].strip()
        head = _head(lhs)
        if head is None:
            self.skipped.append((line(off), stmt[:80]))
            return

        rhs_off = off + i + len(op) + (len(stmt[i + len(op):]) - len(stmt[i + len(op):].lstrip()))
        if head == 'M$ClassesDescription':
            self._classes(rhs, rhs_off, line)
        elif head == 'M$Parameters':
            self._parameters(rhs, rhs_off, line)
        elif head == 'M$GaugeGroups':
            self._gauge_groups(rhs, rhs_off, line)
        self.definitions.setdefault(head, line(off))

    @staticmethod
    def _block_body(rhs: str) -> Optional[str]:
        rhs = rhs.strip()
        return rhs[1:-1] if rhs.startswith('{') and rhs.endswith('}') else None

    def _classes(self, rhs: str, off: int, line) -> None:
        body = self._block_body(rhs)
        if body is None:
            return
        for eoff, lhs, ebody in _entries(body, off + 1):
            m = re.fullmatch(r'([A-Za-z]+)\[\s*(\d+)\s*\]', lhs)
            options = _parse_options(ebody)
            name = options.get('ClassName', '').strip()
            if not m or not name:
                self.skipped.append((line(eoff), lhs))
                continue
            self.classes[name] = FRClass(m.group(1), int(m.group(2)), name, line(eoff), options)

    def _parameters(self, rhs: str, off: int, line) -> None:
        body = self._block_body(rhs)
        if body is None:
            return
        for eoff, lhs, ebody in _entries(body, off + 1):
            if not re.fullmatch(r'[A-Za-z$][A-Za-z0-9$]*', lhs):
                self.skipped.append((line(eoff), lhs))
                continue
            self.parameters[lhs] = FRParameter(lhs, line(eoff), _parse_options(ebody))

    def _gauge_groups(self, rhs: str, off: int, line) -> None:
        body = self._block_body(rhs)
        if body is None:
            return
        for eoff, lhs, ebody in _entries(body, off + 1):
            self.gauge_groups[lhs] = FRGaugeGroup(lhs, line(eoff), _parse_options(ebody))

    # --- queries ----------------------------------------------------------
    def find_classes(self, kind: Optional[str] = None, chirality: Optional[str] = None,
                     unphysical: Optional[bool] = None, **quantum_numbers) -> List[FRClass]:
        '''Classes matching all given criteria; ``quantum_numbers`` are compared
        exactly, a missing quantum number counting as 0 (FeynRules' convention).'''
        out = []
        for c in self.classes.values():
            if kind is not None and c.kind != kind:
                continue
            if chirality is not None and c.chirality != chirality:
                continue
            if unphysical is not None and c.unphysical != unphysical:
                continue
            qn = c.quantum_numbers
            if any(qn.get(k, Fraction(0)) != Fraction(v) for k, v in quantum_numbers.items()):
                continue
            out.append(c)
        return out

    def symbols(self) -> Iterator[Tuple[str, str, int]]:
        '''Every symbol this file defines: (symbol, kind, line). ``kind`` is one of
        class, conjugate, member, mass, width, particle name, parameter, alias,
        block, index, index style, gauge, definition.'''
        for c in self.classes.values():
            yield c.name, 'class', c.line
            for m in c.members:
                yield m, 'member', c.line
            for b in c.conjugate_names:
                yield b, 'conjugate', c.line
            for s in c.mass_symbols:
                yield s, 'mass', c.line
            for s in c.width_symbols:
                yield s, 'width', c.line
            for s in c.particle_names + c.antiparticle_names:
                yield s, 'particle name', c.line
        for p in self.parameters.values():
            yield p.name, 'parameter', p.line
            for a in p.aliases:
                yield a, 'alias', p.line
            if p.block:
                yield p.block, 'block', p.line
        for i in self.indices.values():
            yield i.name, 'index', i.line or 0
            if i.style:
                yield i.style, 'index style', i.line or 0
        for g in self.gauge_groups.values():
            yield g.name, 'gauge', g.line
            for s in g.symbols:
                yield s, 'gauge', g.line
        for name, ln in self.definitions.items():
            yield name, 'definition', ln


# ----------------------------------------------------------------------
# Language-level reserved names
# ----------------------------------------------------------------------

RESERVED_DIR = Path(__file__).with_name('reserved')

# Used only when no generated Mathematica/FeynRules list is present in RESERVED_DIR
# (see scripts/dump_reserved_symbols.wl). Short, hand-picked, deliberately incomplete.
FALLBACK_MATHEMATICA = '''
C D E I K N O Pi Degree Infinity True False None Null All Automatic
Re Im Abs Arg Conjugate Sqrt Exp Log Sin Cos Tan Cot Sec Csc Sinh Cosh Tanh ArcSin ArcCos ArcTan
Plus Times Power Dot List Rule Set Module Block With If Do Table Range Sum Product Function
Expand Simplify Series Solve Part Length Map Apply First Last Head Symbol String Print
Gamma Beta Zeta Erf Max Min Mod Floor Ceiling Round Sign Delta Derivative Integrate D
'''.split()

FALLBACK_FEYNRULES = '''
LoadModel FeynmanRules WriteUFO WriteFeynArtsOutput ExpandIndices WeylToDirac OptimizeIndex
FS DC HC del Ga ProjM ProjP ProjMinus ProjPlus si sibar Eps IndexDelta Index IndexRange IndexStyle
Unfold NoUnfold PauliSigma FR Internal External anti left right CC
ClassName ClassMembers SelfConjugate Unphysical Mass Width Indices FlavorIndex QuantumNumbers PDG
ParticleName AntiParticleName PropagatorLabel PropagatorType PropagatorArrow FullName Definitions
WeylComponents Chirality Left Right Ghost Goldstone Value Values BlockName OrderBlock InteractionOrder
TeX Description ParameterName ParameterType ComplexParameter Unitary Hermitian Orthogonal
Abelian CouplingConstant GaugeBoson Charge StructureConstant Representations SymmetricTensor
Straight Sine Forward Backward GhostDash ScalarDash Spin Lorentz Spin1 Spin2 Lorentz1 Lorentz2
GhostNumber FeynmanGauge
V F S U W T R
'''.split()

UFO_NAMES = '''
Particle Parameter Vertex Coupling Lorentz Decay CouplingOrder FormFactor Function
all_particles all_parameters all_vertices all_couplings all_lorentz all_orders all_decays all_functions
complexconjugate re im csc sec acsc asec cot cmath ZERO P L C V CO Param
'''.split()


def _read_symbol_list(path: Path) -> List[str]:
    '''One symbol per line; ``#`` lines are headers/comments.'''
    out = []
    for raw in path.read_text(encoding='utf-8').splitlines():
        s = raw.strip()
        if s and not s.startswith('#'):
            out.append(s.split('`')[-1])
    return out


def language_lists(reserved_dir: Path = RESERVED_DIR) -> Dict[str, List[str]]:
    '''Language-level reserved names, by list name. Generated Mathematica /
    FeynRules lists (``mathematica*.txt``, ``feynrules*.txt``) are unioned; if a
    family has no file, the hand-picked fallback is used. Python names come from
    the running interpreter.'''
    lists: Dict[str, List[str]] = {}
    for family, fallback in (('mathematica', FALLBACK_MATHEMATICA), ('feynrules', FALLBACK_FEYNRULES)):
        files = sorted(reserved_dir.glob(f'{family}*.txt')) if reserved_dir.is_dir() else []
        if files:
            for f in files:
                lists[f.name] = _read_symbol_list(f)
        else:
            lists[f'{family} (fallback)'] = list(fallback)
    lists['python keywords'] = list(keyword.kwlist) + list(getattr(keyword, 'softkwlist', []))
    lists['python builtins'] = [n for n in dir(builtins) if not n.startswith('_')]
    lists['UFO'] = list(UFO_NAMES)
    return lists


# ----------------------------------------------------------------------
# Registry
# ----------------------------------------------------------------------

# kinds whose names reach the UFO / MG5 and are therefore compared case-insensitively
UFO_BOUND_KINDS = frozenset({'mass', 'width', 'parameter', 'alias', 'block', 'particle name'})


@dataclass(frozen=True)
class Origin:
    '''Where a reserved symbol comes from.'''
    symbol: str
    kind: str            # a FRModelFile.symbols() kind, or 'language'
    where: str           # 'SMMod.fr:289' or a list name such as 'python keywords'

    @property
    def case_insensitive(self) -> bool:
        return self.kind in UFO_BOUND_KINDS

    def __str__(self) -> str:
        return f'{self.symbol} ({self.kind}, {self.where})'


class ReservedNames:
    '''Registry of names that generated symbols must not reuse.'''

    def __init__(self):
        self._exact: Dict[str, List[Origin]] = {}
        self._lower: Dict[str, List[Origin]] = {}
        self.models: List[FRModelFile] = []

    # --- building -----------------------------------------------------------
    def add(self, symbol: str, kind: str, where: str) -> None:
        origin = Origin(symbol, kind, where)
        self._exact.setdefault(symbol, []).append(origin)
        self._lower.setdefault(symbol.lower(), []).append(origin)

    def add_model(self, model: FRModelFile) -> None:
        self.models.append(model)
        for symbol, kind, ln in model.symbols():
            self.add(symbol, kind, f'{model.source}:{ln}')

    @classmethod
    def from_fr(cls, *paths, language: bool = True,
                reserved_dir: Path = RESERVED_DIR) -> 'ReservedNames':
        '''Registry of the symbols defined in the given ``.fr`` files, plus (by
        default) the language-level lists.'''
        return cls.from_models(*(FRModelFile.parse(p) for p in paths),
                               language=language, reserved_dir=reserved_dir)

    @classmethod
    def from_models(cls, *models: FRModelFile, language: bool = True,
                    reserved_dir: Path = RESERVED_DIR) -> 'ReservedNames':
        '''As ``from_fr``, from already parsed files.'''
        reg = cls()
        for m in models:
            reg.add_model(m)
        if language:
            for list_name, names in language_lists(reserved_dir).items():
                for n in names:
                    reg.add(n, 'language', list_name)
        return reg

    # --- queries ------------------------------------------------------------
    def copy(self) -> 'ReservedNames':
        other = ReservedNames()
        other._exact = {k: list(v) for k, v in self._exact.items()}
        other._lower = {k: list(v) for k, v in self._lower.items()}
        other.models = list(self.models)
        return other

    def __contains__(self, symbol: str) -> bool:
        return symbol in self._exact

    def origins(self, symbol: str) -> List[Origin]:
        return list(self._exact.get(symbol, []))

    def clashes(self, symbol: str, ufo_bound: bool = False) -> List[Origin]:
        '''Reserved entries that ``symbol`` would collide with. Exact matches always
        count. Case-only matches count when either side reaches the UFO (the
        registry entry's kind is UFO-bound, or ``ufo_bound`` says the new symbol
        will be written to the UFO) — except for language-level names, which are
        case-sensitive.'''
        out = list(self._exact.get(symbol, []))
        for o in self._lower.get(symbol.lower(), []):
            if o.symbol == symbol or o.kind == 'language':
                continue
            if o.case_insensitive or ufo_bound:
                out.append(o)
        return out

    def is_free(self, symbol: str, ufo_bound: bool = False) -> bool:
        '''Valid symbol with no clash.'''
        return is_valid_symbol(symbol) and not self.clashes(symbol, ufo_bound)


# ----------------------------------------------------------------------
# SM model file: default path, cached loading, SM-class map
# ----------------------------------------------------------------------

DEFAULT_SM_FILE = Path(__file__).resolve().parents[1] / 'data' / 'SMMod.fr'

_MODEL_CACHE: Dict[Tuple[str, int], FRModelFile] = {}


def load_model(path=DEFAULT_SM_FILE) -> FRModelFile:
    '''``FRModelFile.parse`` cached per (resolved path, modification time).
    The returned object is shared: treat it as read-only.'''
    path = Path(path).resolve()
    key = (str(path), path.stat().st_mtime_ns)
    if key not in _MODEL_CACHE:
        _MODEL_CACHE[key] = FRModelFile.parse(path)
    return _MODEL_CACHE[key]


class SMClassMap:
    '''Maps our SM multiplets onto the unphysical classes of an SM model file,
    by quantum numbers:

    * a LH Weyl multiplet matches a ``W`` class of the same chirality, Y and
      global charges (e.g. lL2 -> ``LLw``, QL1 -> ``QLw``), or else a class of the
      opposite chirality with the conjugate charges, used through its automatic
      conjugate (e.g. our tRbar, LH with Y = -2/3 -> ``uRwbar``);
    * a scalar matches an ``S`` class the same way (Phi -> ``Phi``, our conjugate
      Higgs -> ``Phitilde``).

    Besides Y, only the global charges (LeptonNumber, BaryonNumber, ...) that the
    file actually uses are compared, and the non-generation index types must agree
    (SU2D for doublets, Colour for coloured fields). The family index of the
    multiplet (``multiplet.index``) fills the ``Generation`` index.
    '''

    def __init__(self, model: FRModelFile):
        self.model = model
        self._cache: Dict[object, Tuple[FRClass, bool]] = {}
        used = set()
        for c in model.classes.values():
            used.update(c.quantum_numbers)
        from uvbuilder.tools.fr_utils import GLOBAL_U1_QN   # local: fr_utils has no dependency on this module
        self._globals = {short: long for short, long in GLOBAL_U1_QN.items() if long in used}

    def defines_class(self, name: str) -> bool:
        return name in self.model.classes

    @property
    def index_names(self) -> frozenset:
        '''Index types declared by the SM file (IndexRange/IndexStyle).'''
        return frozenset(self.model.indices)

    # --- matching ------------------------------------------------------------
    def _signature(self, rep, conjugate: bool) -> Dict[str, Fraction]:
        def frac(x) -> Fraction:
            if x is None:
                return Fraction(0)
            return Fraction(x.numerator, x.denominator) if hasattr(x, 'denominator') else Fraction(x)

        sign = -1 if conjugate else 1
        qn = {'Y': sign * frac(rep.hypercharge)}
        for short, long in self._globals.items():
            qn[long] = sign * frac(rep.globalU1s[short] if rep.globalU1s is not None else None)
        return qn

    def lookup(self, multiplet) -> Tuple[FRClass, bool]:
        '''(class, conjugated): the class representing ``multiplet``, and whether
        the multiplet is that class' conjugate. Raises LookupError if there is no
        unique match.'''
        from uvbuilder.tools.reps import WeylSpinor, Scalar, WEYL_LEFT
        rep = multiplet.rep
        key = (rep, multiplet.index)
        if key in self._cache:
            return self._cache[key]

        if isinstance(rep.spin, WeylSpinor):
            kind = 'W'
            chir = 'Left' if rep.spin == WEYL_LEFT else 'Right'
            other = 'Right' if chir == 'Left' else 'Left'
        elif isinstance(rep.spin, Scalar):
            kind, chir, other = 'S', None, None
        else:
            raise LookupError(f'No SM-class translation for spin {rep.spin} ({multiplet.label})')

        wanted = set()
        if rep.isospin.dim > 1:
            from uvbuilder.tools.fr_utils import su2_index_type
            wanted.add(su2_index_type(rep.isospin.dim))
        if rep.colordim.dim == 3:
            wanted.add('Colour')
        elif rep.colordim.dim == 8:
            wanted.add('Gluon')
        elif rep.colordim.dim > 1:
            raise LookupError(f'No SM index type for colour {rep.colordim} ({multiplet.label})')

        def candidates(chirality, conjugate):
            out = []
            for c in self.model.find_classes(kind=kind, chirality=chirality, unphysical=True,
                                             **self._signature(rep, conjugate)):
                if set(c.indices) - {'Generation'} == wanted:
                    out.append(c)
            return out

        direct = candidates(chir, False)
        found = [(c, False) for c in direct] or [(c, True) for c in candidates(other, True)]
        if len(found) != 1:
            what = 'no' if not found else 'several (' + ', '.join(c.name for c, _ in found) + ')'
            raise LookupError(
                f'{self.model.source}: {what} class matching SM multiplet {multiplet.label} '
                f'({kind}, chirality {chir}, {self._signature(rep, False)}, indices {sorted(wanted)})')
        self._cache[key] = found[0]
        return found[0]

    # --- rendering -----------------------------------------------------------
    def symbol(self, multiplet, bar: bool = False) -> str:
        '''FeynRules symbol of the multiplet (``bar`` requests its conjugate).'''
        cls, conjugated = self.lookup(multiplet)
        return cls.name + ('bar' if conjugated != bar else '')

    def render(self, multiplet, spin: Optional[str] = None, colour: Optional[str] = None,
               isospin: Optional[str] = None, bar: bool = False) -> str:
        '''The multiplet with its indices in the class' declared order, e.g.
        ``QLw[sp1, i1, 3, c1]`` or ``uRwbar[sp1, 3, c2]``.'''
        cls, _ = self.lookup(multiplet)
        values = {'Colour': colour, 'Gluon': colour, 'Sextet': colour,
                  'Generation': str(multiplet.index)}
        idx = [spin] if cls.kind == 'W' else []
        for index_type in cls.indices:
            value = isospin if index_type.startswith('SU2') else values.get(index_type)
            if value is None:
                raise ValueError(f'No {index_type} index given for {multiplet.label} ({cls.name})')
            idx.append(value)
        if cls.kind == 'W' and spin is None:
            raise ValueError(f'No spin index given for {multiplet.label} ({cls.name})')
        name = self.symbol(multiplet, bar)
        return f'{name}[{", ".join(idx)}]' if idx else name


# ----------------------------------------------------------------------
# Renaming generated symbols: applied when the FeynRules file is
# written; earlier stages only warn.
# ----------------------------------------------------------------------

_NON_SYMBOL = re.compile(r'[^A-Za-z0-9]')


def fix_symbol(name: str) -> str:
    '''Deterministic repair of an invalid symbol: ``+`` -> ``P``, ``-`` -> ``M``,
    other invalid characters dropped, a leading digit (or nothing left) gets an
    ``X`` prefix. Valid symbols are returned unchanged.'''
    s = _NON_SYMBOL.sub('', name.replace('+', 'P').replace('-', 'M'))
    if not s or not s[0].isalpha():
        s = 'X' + s
    return s


Derived = List[Tuple[str, bool]]      # (symbol, reaches the UFO)


def particle_symbols(label: str) -> Derived:
    '''Symbols the FeynRules output derives from a BSM particle label: the class
    (also the UFO particle name) and its conjugate, the Weyl components and their
    conjugates, the mass and width parameters.'''
    return [(label, True), (label + 'bar', False),
            (label + 'Lw', False), (label + 'Rw', False),
            (label + 'Lwbar', False), (label + 'Rwbar', False),
            ('M' + label, True), ('W' + label, True)]


def multiplet_symbols(fr_symbol: str) -> Derived:
    '''A multiplet's FeynRules symbol, its automatic conjugate and its tilde class
    (see multipletModel.tilde_class).'''
    return [(fr_symbol, False), (fr_symbol + 'bar', False), (fr_symbol + 'tilde', False)]


def parameter_symbols(name: str) -> Derived:
    '''An operator parameter and the real/imaginary parts of a complex one.'''
    return [(name, True), ('R' + name, True), ('I' + name, True)]


class Renamer:
    '''Assigns generated roots (particle labels, multiplet labels, parameter
    names) the first valid, unreserved variant, and remembers every symbol it
    hands out so generated names cannot collide with each other either.'''

    def __init__(self, reserved: ReservedNames):
        self.reserved = reserved.copy()

    def problems(self, derived: Derived) -> List[str]:
        out = []
        for symbol, ufo in derived:
            if not is_valid_symbol(symbol):
                out.append(f'{symbol!r} is not a valid symbol')
            else:
                out += [f'{symbol!r} clashes with {o}' for o in self.reserved.clashes(symbol, ufo)]
        return out

    def _register(self, derived: Derived, where: str) -> None:
        for symbol, ufo in derived:
            self.reserved.add(symbol, 'parameter' if ufo else 'class', f'generated ({where})')

    def claim(self, root: str, derive, where: str) -> str:
        '''The root itself if all its derived symbols are free, else ``fix_symbol``
        of it, else that with the smallest integer suffix that works.'''
        base = fix_symbol(root)
        for candidate in itertools.chain([root, base], (f'{base}{k}' for k in itertools.count(1))):
            derived = derive(candidate)
            if not self.problems(derived):
                self._register(derived, where)
                return candidate
        raise AssertionError('unreachable')

    def check(self, derived: Derived, where: str) -> List[str]:
        '''Check (and register) symbols that cannot be renamed; returns the problems.'''
        found = self.problems(derived)
        self._register(derived, where)
        return found


_RESERVED_CACHE: Dict[Tuple[str, int], ReservedNames] = {}


def reserved_names(sm_file=DEFAULT_SM_FILE) -> ReservedNames:
    '''Registry for an SM model file plus the language lists, cached per file
    version. Shared: use ``.copy()`` (``Renamer`` does) before adding to it.'''
    model = load_model(sm_file)
    path = Path(sm_file).resolve()
    key = (str(path), path.stat().st_mtime_ns)
    if key not in _RESERVED_CACHE:
        _RESERVED_CACHE[key] = ReservedNames.from_models(model)
    return _RESERVED_CACHE[key]


def warn_labels(labels: Iterable[str], derive, what: str, sm_file=DEFAULT_SM_FILE,
                stacklevel: int = 3) -> List[str]:
    '''Early, warn-only check of labels (nothing is changed): lists the labels
    whose derived symbols are invalid or reserved. They are renamed when the
    FeynRules file is written.'''
    try:
        renamer = Renamer(reserved_names(sm_file))
    except OSError:
        return []
    bad = []
    for label in labels:
        if not is_valid_symbol(label):
            bad.append(f'{label!r}: not a valid symbol')
            continue
        problems = renamer.problems(derive(label))
        if problems:
            bad.append(f'{label!r}: ' + '; '.join(problems[:3]) + (' ...' if len(problems) > 3 else ''))
    if bad:
        warnings.warn(f'{what} labels will be renamed when the FeynRules file is written '
                      f'({Path(sm_file).name}): ' + ' | '.join(bad), stacklevel=stacklevel)
    return bad
