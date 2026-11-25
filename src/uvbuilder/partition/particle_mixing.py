"""
.. module:: particle_mixing
   :synopsis: Defines the MixingMatrix and ParticleMixing classes, the mixing of the mass eigenstates.

.. moduleauthor:: Lucas Heck <lucasdealexandri.heck@gmail.com>
"""
# We use UserDict instead of inheriting from dict for more predictable behavior
# See, e.g., https://treyhunner.com/2019/04/why-you-shouldnt-inherit-from-list-and-dict-in-python/
from collections import UserDict
from uvbuilder.tools.reps import WeylSpinor, Scalar
from typing import Sequence, Union
from uvbuilder.tools.components import ScalarComponent, WeylComponent
from uvbuilder.states.particle import family_key



class MixingMatrix:
    '''Mixing block between field components (rows) and mass eigenstates (columns).
    Blocks of real scalars are real orthogonal, in the real basis of the block
    (real_basis).'''

    def __init__(self, 
                 p_components: Union[Sequence[ScalarComponent], Sequence[WeylComponent]], 
                 f_components: Union[Sequence[ScalarComponent], Sequence[WeylComponent]],
                 particles: Sequence = ()):
        '''
        :param p_components: components of the mass eigenstates (columns)
        :param f_components: gauge components of the multiplets (rows)
        :param particles: the mass eigenstates, in mass order (their family index)
        '''
        self.p_components = p_components
        self.f_components = f_components
        self.particles = tuple(particles)   # mass eigenstates, mass order (family index)

        self.validation()

        # Once validated we know the representation
        # of all particles and fields are the same.
        self.rep = self.p_components[0].rep

    def validation(self) -> None:
        '''
        Validations:
        1. The length of the p components equals the length of the f components
        2. All the reps in the p components are equal
        3. All the unbroken reps in the f components are equal
        4. All the unbroken reps between the p and f components are equal
        '''

        p0 = self.p_components[0]
        f0 = self.f_components[0]

        same_len = len(self.p_components) == len(self.f_components)
        same_p_reps = all(x.rep == p0.rep for x in self.p_components)
        same_f_unbroken_reps = all(
            x.rep.unbroken == f0.rep.unbroken for x in self.f_components)
        same_p_f_unbroken_reps = f0.rep.unbroken == p0.rep

        if not same_len:
            raise ValueError(f"Different number of particle components ({len(self.p_components)}) and field components ({len(self.f_components)})")

        if not same_p_reps:
            raise ValueError(f"All reps of the particles should match.\n{[p.rep for p in self.p_components]}")

        if not same_f_unbroken_reps:
            raise ValueError(f"All unbroken reps of the fields should match.\n{[f.rep.unbroken for f in self.f_components]}")

        if not same_p_f_unbroken_reps:
            raise ValueError(f"Field <-> Particle rep mismatch:\n{[p.rep for p in self.p_components]}\n{[f.rep.unbroken for f in self.f_components]}")
        
        return 

    @property
    def type(self) -> str:
        '''MajoranaSpinor or DiracSpinor (Weyl components of a real or complex
        representation), RealScalar or ComplexScalar.'''
        spin = self.rep.spin
        realInternalRep = self.rep.is_realInternalRep

        if isinstance(spin, WeylSpinor):
            return {
                True: "MajoranaSpinor",
                False: "DiracSpinor"
            }[realInternalRep]

        elif isinstance(spin, Scalar):
            return {
                True: "RealScalar",
                False: "ComplexScalar"
            }[realInternalRep]

        else: raise TypeError(f"Only WeylSpinor and Scalars are allowed. Received {spin}")

    @property
    def dim(self) -> int:
        return len(self.p_components)

    @property
    def is_right_block(self) -> bool:
        '''True for the block of RH-bar components (Dirac V matrix).'''
        return self.p_components[0].rep != self.particles[0].components[0].rep

    @property
    def standard_symbol(self) -> str:
        '''Unique name {U|V}Mix{type}{particle labels}, e.g. UMixMN1N2N3.'''
        symbol_translation = {'RealScalar': 'R', 'ComplexScalar': 'C', 
                              'MajoranaSpinor': 'M', 'DiracSpinor': 'D'}
        lr = 'V' if self.is_right_block else 'U'
        labels = ''.join(str(p) for p in self.particles)
        return f'{lr}Mix{symbol_translation[self.type]}{labels}'

    def to_feynrules(self, tracker) -> dict[str, list[str]]:
        '''FeynRules chunk for this mixing matrix (MSSM convention): real and
        imaginary External parts and the Internal unitary matrix
        M[i_, j_] :> RM[i, j] + I*IM[i, j]. Real scalars mix through a real
        orthogonal matrix: a single real External M (no imaginary block).
        Indexed by the family index of the physical class; defaults to the
        identity. Single-member blocks need no matrix (empty chunk).'''
        if self.dim == 1:
            return {}
        from uvbuilder.models.particleModel import ParticleModel   # local: avoids circular import
        from uvbuilder.tools.fr_utils import FRParameter

        family = tracker.family(ParticleModel._family_name(list(self.particles)), self.dim)
        name, block = self.standard_symbol, self.standard_symbol.upper()
        indices = f'{{{family}, {family}}}'
        n = range(1, self.dim + 1)
        identity = lambda m, d: '{' + ',\n '.join(
            ', '.join(f'{m}[{i}, {j}] -> {d * (i == j)}' for j in n) for i in n) + '}'
        states = ', '.join(p.label for p in self.particles)

        def entry(symbol, sort, **options):
            return FRParameter(symbol, 'Mixing', tuple(options.items()), (name, sort))

        if self.type == 'RealScalar':
            return {'Parameters': [entry(
                name, 0, ParameterType='External', BlockName=block, Indices=indices,
                Unitary=True, Value=identity(name, 1), ComplexParameter=False,
                Description=f'"Mixing matrix of {states} (real orthogonal)"')]}
        return {'Parameters': [
            entry(f'R{name}', 1, ParameterType='External', BlockName=block, Indices=indices,
                  Value=identity(f'R{name}', 1), ComplexParameter=False,
                  Description=f'"Mixing matrix of {states} (real part)"'),
            entry(f'I{name}', 2, ParameterType='External', BlockName=f'IM{block}', Indices=indices,
                  Value=identity(f'I{name}', 0), ComplexParameter=False,
                  Description=f'"Mixing matrix of {states} (imaginary part)"'),
            entry(name, 0, ParameterType='Internal', Indices=indices, Unitary=True,
                  Value=f'{{{name}[i_, j_] :> R{name}[i, j] + I*I{name}[i, j]}}',
                  ComplexParameter=True, Description=f'"Mixing matrix of {states}"')]}

    def get_index(self, field):
        '''Row of a gauge component in the block, or None if it is not in it.'''
        if field not in self.f_components: return None

        return self.f_components.index(field)

    def real_basis(self, declared):
        '''(labels, C): the real basis rho of the block and the change of basis
        phi = C rho, phi the field components (rows). Only neutral scalars (real scalar
        blocks) have a non-trivial one, the basis of CP eigenstates if CP is conserved:
        a component Phi of a declared multiplet and its conjugate PhiBar are
            Phi = (Phi_R + i Phi_I)/sqrt(2),   PhiBar = (Phi_R - i Phi_I)/sqrt(2).
        rho keeps the order of the rows: Phi_R is listed in the row of Phi and Phi_I in
        the row of PhiBar, which fixes the rows of the orthogonal mixing matrix. Real
        fields are kept.'''
        import sympy as sp
        labels, C = [f.label for f in self.f_components], sp.eye(self.dim)
        if self.type != 'RealScalar':
            return labels, C
        for r, f in enumerate(self.f_components):
            if f in declared and not f.selfConjugate:
                if f.conjugateField not in self.f_components:
                    raise ValueError(f'{f} mixes into real states without its conjugate')
                b = self.get_index(f.conjugateField)
                labels[r], labels[b] = f'{f.label}_R', f'{f.label}_I'
                C[r, r], C[r, b] = 1 / sp.sqrt(2), sp.I / sp.sqrt(2)
                C[b, r], C[b, b] = 1 / sp.sqrt(2), -sp.I / sp.sqrt(2)
        return labels, C

    def to_real_basis(self, V, declared):
        '''V with the field components of the block written in its real basis,
        phi = C rho (real_basis).'''
        import sympy as sp
        labels, C = self.real_basis(declared)
        phi = [sp.Symbol(f.label) for f in self.f_components]
        rho = sp.Matrix([sp.Symbol(label) for label in labels])
        return V.xreplace({p: e for p, e in zip(phi, C * rho) if e != p})

    def mass_matrix(self, V, declared, right=None):
        '''(M, quadratic form) of the block from V, the two-field terms at the vev
        (L = -V - h.c.) with every real block in its real basis (to_real_basis):
          Majorana, real scalar   hessian of V in rho          (-L = 1/2 rho^T M rho [+ h.c.])
          Dirac                   d2V / dU_i dV_j, right = V block  (-L = U^T M V + h.c.)
          complex scalar          d2V / dphibar_i dphi_j       (-L = phi^dagger M^2 phi)'''
        import sympy as sp
        symbols = lambda fields: [sp.Symbol(f.label) for f in fields]
        rows = symbols(self.f_components)
        if self.type == 'DiracSpinor':
            columns = symbols(right.f_components)
            M = sp.Matrix(len(rows), len(columns), lambda i, j: sp.diff(V, rows[i], columns[j]))
            return M, (sp.Matrix([rows]) * M * sp.Matrix(columns))[0]
        if self.type == 'ComplexScalar':
            bars = symbols(f.conjugateField for f in self.f_components)
            M = sp.Matrix(len(rows), len(rows), lambda i, j: sp.diff(V, bars[i], rows[j]))
            return M, (sp.Matrix([bars]) * M * sp.Matrix(rows))[0]
        rho = [sp.Symbol(label) for label in self.real_basis(declared)[0]]
        M = sp.simplify(sp.hessian(V, rho))
        return M, (sp.Matrix([rho]) * M * sp.Matrix(rho))[0] / 2

class ParticleMixing(UserDict):
    """
    Mixing of a model: {mass eigenstate: set of gauge components it mixes with}, with its
    inverse (gauge component -> mass eigenstates) and the mixing blocks it defines
    (get_mixing_matrices).
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
    
    def __setitem__(self, key, value):
        
        if isinstance(value, set):
            value = value
            
        elif isinstance(value, (list, tuple)):
            value = set(value)
            
        else: value = {value}
        
        super().__setitem__(key, value)
        
        
    def __getitem__(self, key):
        # Retrieve the value using the parent's logic
        value = super().__getitem__(key)
        
        # Otherwise return the full set 
        return value
        
    
    def __delitem__(self, key):
        super().__delitem__(key)
    
        
    def __repr__(self):
        return super().__repr__()
    
   
    @property 
    def inverse(self):
        '''{gauge component: set of mass eigenstates}.'''
        inv_map = ParticleMixing() # Return a ParticleMapping instance for consistency
        
        for particle, fields in self.data.items():
            for f in fields:
                # We can now use standard assignment because __setitem__ handles the set creation
                # However, since we are building an inverse of a one-to-many, 
                # we need to append to existing sets if keys collide.
                
                if f not in inv_map.data:
                    inv_map.data[f] = {particle}
                else:
                    inv_map.data[f].add(particle)

        return inv_map

    def get_cluster(self) -> dict[tuple, list]:
        '''Get Particle <-> Field clusters of mixings. e.g., particles P1, P2, P3
        mixing into fields F1, F2, F3, and particles P4, P5 mixing into F4, F5
        will return {(P1,P2,P3): [F1,F2,F3], (P4, P5): [F4, F5])}'''

        clusters = {}
        for f, particles in self.inverse.items():
            # Keep the owning particle to order the components by mass rank, which
            # is the family-index order of the physical FeynRules class.
            if isinstance(f.rep.spin, WeylSpinor):
                pairs = {(p, c.left_hand) for p in particles for c in p.components
                         if c.left_hand.rep == f.rep.unbroken}
            else:
                pairs = {(p, c) for p in particles for c in p.components
                         if c.rep == f.rep.unbroken}
            ps = tuple(c for _, c in sorted(pairs, key=lambda pc: family_key(pc[0])))
                
            if ps in clusters:
                clusters[ps] = sorted(clusters[ps] + [f])
            else: clusters.update({ps: [f]})

        return clusters

    def get_mixing_matrices(self) -> list[MixingMatrix]:
        '''One MixingMatrix per block (get_cluster), the mass eigenstates in mass order,
        sorted by name.'''
        matrices = []
        for p_components, f_components in self.get_cluster().items():
            particles = sorted({p for f in f_components for p in self.inverse[f]},
                               key=family_key)
            matrices.append(MixingMatrix(p_components, f_components, particles))

        return sorted(matrices, key=lambda m: m.standard_symbol)

