# -*- coding: utf-8 -*-
#
#  Copyright 2017-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
#  This file is part of chython.
#
#  chython is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
from functools import cached_property
from typing import Dict, Iterator, Tuple, Optional, List, Union
from collections import defaultdict

from .bonds import DynamicBond, Bond
from . import cgr_query as query

from ..algorithms.fingerprints import FingerprintsCGR
from ..algorithms.isomorphism import Isomorphism
from ..algorithms.morgan import Morgan
from ..algorithms.rings import Rings
from ..algorithms.smiles import CGRSmiles
from ..algorithms.calculate2d import Calculate2DCGR
from ..algorithms.depict import DepictCGR
from ..algorithms.x3dom import X3domCGR

from ..periodictable import DynamicElement, Element, DynamicQueryElement
from ..exceptions import MappingError, AtomNotFound, BondNotFound


class CGRContainer(CGRSmiles, DepictCGR, Calculate2DCGR, X3domCGR, Morgan, Rings, Isomorphism, FingerprintsCGR):
    __slots__ = ('_atoms', '_bonds', '_conformers', '_charges', '_radicals', '_p_charges', '_p_radicals', '_hybridizations', '_p_hybridizations', '_plane', '__dict__')
    _atoms: Dict[int, DynamicElement]
    _bonds: Dict[int, Dict[int, DynamicBond]]

    def __init__(self):
        self._atoms: Dict[int, DynamicElement] = {}
        self._bonds: Dict[int, Dict[int, DynamicBond]] = {}
        self._conformers: List[Dict[int, Tuple[float, float, float]]] = []
        self._charges: Dict[int, int] = {}
        self._radicals: Dict[int, bool] = {}
        self._p_charges: Dict[int, int] = {}
        self._p_radicals: Dict[int, bool] = {}
        self._hybridizations: Dict[int, int] = {}
        self._p_hybridizations: Dict[int, int] = {}
        self._plane: Dict[int, Tuple[float, float]] = {}
        self.flush_cache()

    # Methods from Graph
    def flush_cache(self):
        self.__dict__.clear()

    def atom(self, n: int) -> DynamicElement:
        return self._atoms[n]

    def has_atom(self, n: int) -> bool:
        return n in self._atoms

    @property
    def atoms_count(self) -> int:
        return len(self._atoms)

    @property
    def atoms_numbers(self) -> Iterator[int]:
        return iter(self._atoms)

    def bond(self, n: int, m: int) -> DynamicBond:
        try:
            return self._bonds[n][m]
        except KeyError as e:
            raise BondNotFound from e

    def has_bond(self, n: int, m: int) -> bool:
        try:
            # check if atom exists and has a bonds dict
            return m in self._bonds[n]
        except KeyError:
            return False # If atom n doesn't exist, it has no bonds.

    def _validate_charge(self, charge: int) -> int:
        # Placeholder: Add actual validation logic if needed
        # For now, just ensure it's an int.
        if not isinstance(charge, int):
            raise TypeError(f"Charge must be an integer, got {type(charge)}")
        # Potentially add range checks, e.g., if charge must be within -4 to 4
        return charge

    def _validate_radical(self, is_radical: bool) -> bool:
        # Placeholder: Add actual validation logic if needed
        if not isinstance(is_radical, bool):
            raise TypeError(f"is_radical must be a boolean, got {type(is_radical)}")
        return is_radical

    def atoms(self) -> Iterator[Tuple[int, DynamicElement]]:
        return iter(self._atoms.items())

    def bonds(self) -> Iterator[Tuple[int, int, DynamicBond]]:
        """
        Iterate other all bonds
        """
        seen = set()
        for n, m_bond in self._bonds.items():
            seen.add(n)
            for m, bond in m_bond.items():
                if m not in seen:
                    yield n, m, bond

    def neighbors(self, n: int) -> Tuple[int, int]:
        """number of neighbors atoms excluding any-bonded"""
        s = p = 0
        for b in self._bonds[n].values():
            if b.order is not None:
                if b.order != 8:
                    s += 1
                if b.p_order is not None and b.p_order != 8:
                    p += 1
            elif b.p_order != 8:
                p += 1
        return s, p

    def add_atom(self, atom: Union[DynamicElement, Element, int, str], n: Optional[int] = None, **kwargs):
        p_charge = self._validate_charge(kwargs.pop('p_charge', 0))
        p_is_radical = self._validate_radical(kwargs.pop('p_is_radical', False))
        xy_coord = kwargs.pop('xy', None)

        if not isinstance(atom, DynamicElement):
            if isinstance(atom, Element):
                atom_cls = DynamicElement.from_atomic_number(atom.atomic_number)
                atom = atom_cls(atom.isotope)
            elif isinstance(atom, str):
                atom_cls = DynamicElement.from_symbol(atom)
                atom = atom_cls()
            elif isinstance(atom, int):
                atom_cls = DynamicElement.from_atomic_number(atom)
                atom = atom_cls()
            else:
                raise TypeError('DynamicElement object expected')

        if n is None:
            n = max(self._atoms, default=0) + 1
        elif not isinstance(n, int):
            raise TypeError('mapping should be integer')
        elif n in self._atoms:
            raise MappingError(f'atom with number {{{n}}} already exists')

        self._atoms[n] = atom
        self._bonds[n] = {}

        if xy_coord is not None:
            self._plane[n] = xy_coord
        self._charges[n] = atom.charge
        self._radicals[n] = atom.is_radical
        self._p_charges[n] = p_charge
        self._p_radicals[n] = p_is_radical
        self._hybridizations[n] = 1
        self._p_hybridizations[n] = 1
        self._conformers.clear()
        self.flush_cache()
        return n

    def add_bond(self, n, m, bond: Union[DynamicBond, Bond, int]):
        if n == m:
            raise MappingError('atom loops impossible')
        if n not in self._bonds or m not in self._bonds:
            raise AtomNotFound('atoms not found')
        if m in self._bonds[n]:
            raise MappingError('atoms already bonded')

        if isinstance(bond, DynamicBond):
            order = bond.order
            p_order = bond.p_order
        elif isinstance(bond, Bond):
            order = p_order = bond.order
            bond = DynamicBond.from_bond(bond)
        else:
            order = p_order = bond
            bond = DynamicBond(order, order)

        self._bonds[n][m] = self._bonds[m][n] = bond
        self._conformers.clear()

        if order != 1 or p_order != 1:
            self._calc_hybridization(n)
            self._calc_hybridization(m)
        self.flush_cache()

    def delete_atom(self, n):
        if n not in self._atoms:
            raise AtomNotFound(f"Atom {n} not found for deletion")

        old_bonds_neighbors = list(self._bonds.get(n, {}).keys())
        for m in old_bonds_neighbors:
            del self._bonds[m][n]

        if n in self._bonds:
            del self._bonds[n]
        del self._atoms[n]

        if n in self._charges: del self._charges[n]
        if n in self._radicals: del self._radicals[n]
        if n in self._p_charges: del self._p_charges[n]
        if n in self._p_radicals: del self._p_radicals[n]
        if n in self._hybridizations: del self._hybridizations[n]
        if n in self._p_hybridizations: del self._p_hybridizations[n]
        if n in self._plane: del self._plane[n]

        for m in old_bonds_neighbors:
            self._calc_hybridization(m)
        self._conformers.clear()
        self.flush_cache()

    def delete_bond(self, n, m):
        if n not in self._bonds or m not in self._bonds.get(n, {}):
            raise BondNotFound(f"Bond between {n} and {m} not found for deletion")

        del self._bonds[n][m]
        del self._bonds[m][n]
        self._conformers.clear()
        self._calc_hybridization(n)
        self._calc_hybridization(m)
        self.flush_cache()

    def remap(self, mapping: Dict[int, int], *, copy=False):
        target = self.copy() if copy else self

        if len(mapping) != len(set(mapping.values())) or \
                not (target._atoms.keys() - mapping.keys()).isdisjoint(mapping.values()):
            raise ValueError('mapping overlap')

        mg = mapping.get
        target._atoms = {mg(n, n): atom for n, atom in target.atoms()}
        target._bonds = {mg(n, n): {mg(m, m): bond for m, bond in m_bond.items()} for n, m_bond in target._bonds.items()}
        target._charges = {mg(n, n): c for n, c in target._charges.items()}
        target._radicals = {mg(n, n): r for n, r in target._radicals.items()}
        target._p_charges = {mg(n, n): c for n, c in target._p_charges.items()}
        target._p_radicals = {mg(n, n): r for n, r in target._p_radicals.items()}
        target._hybridizations = {mg(n, n): h for n, h in target._hybridizations.items()}
        target._p_hybridizations = {mg(n, n): h for n, h in target._p_hybridizations.items()}
        target._plane = {mg(n, n): p for n, p in target._plane.items()}
        target._conformers = [{mg(k, k): v for k, v in conf.items()} for conf in target._conformers]
        target.flush_cache()

        if copy:
            return target

    def copy(self, **kwargs) -> 'CGRContainer':
        copy = object.__new__(self.__class__)
        copy._atoms = {n: atom.copy(full=True) for n, atom in self._atoms.items()}
        copy._bonds = cb = {}
        for n, m_bond in self._bonds.items():
            cb[n] = cbn = {}
            for m, bond in m_bond.items():
                if m in cb:
                    cbn[m] = cb[m][n]
                else:
                    cbn[m] = bond.copy(full=True)

        copy._hybridizations = self._hybridizations.copy()
        copy._conformers = [c.copy() for c in self._conformers]
        copy._p_hybridizations = self._p_hybridizations.copy()
        copy._p_radicals = self._p_radicals.copy()
        copy._p_charges = self._p_charges.copy()
        copy._charges = self._charges.copy()
        copy._radicals = self._radicals.copy()
        copy._plane = self._plane.copy()
        copy.__dict__ = {}
        return copy

    def substructure(self, atoms, *, as_query: bool = False, **kwargs):
        """
        create substructure containing atoms from atoms list

        :param atoms: list of atoms numbers of substructure
        :param as_query: return Query object based on graph substructure
        """
        graph_type = query.QueryCGRContainer if as_query else self.__class__
        atom_type = DynamicQueryElement if as_query else DynamicElement

        atoms_to_include = {n for n in atoms if n in self._atoms}

        sub = graph_type()
        if not hasattr(sub, '_plane'): # QueryCGR may not have plane
            sub._plane = {}

        for n in atoms_to_include:
            new_atom = atom_type.from_atom(self._atoms[n])
            sub.add_atom(new_atom, n, xy=self._plane.get(n))

        for n in atoms_to_include:
            for m, bond in self._bonds[n].items():
                if m in atoms_to_include and n < m:
                    sub.add_bond(n, m, bond.copy())

        if not hasattr(sub, '_charges'):
            sub._charges = {}
            sub._radicals = {}
        sub._charges.update({n: self._charges[n] for n in atoms_to_include if n in self._charges})
        sub._radicals.update({n: self._radicals[n] for n in atoms_to_include if n in self._radicals})
        sub._p_charges = {n: self._p_charges[n] for n in atoms_to_include if n in self._p_charges}
        sub._p_radicals = {n: self._p_radicals[n] for n in atoms_to_include if n in self._p_radicals}

        if as_query:
            sub._hybridizations = {n: (self._hybridizations[n],) for n in atoms_to_include if n in self._hybridizations}
            sub._p_hybridizations = {n: (self._p_hybridizations[n],) for n in atoms_to_include if n in self._p_hybridizations}

            sub._neighbors = cn = {}
            sub._p_neighbors = cpn = {}
            for n in atoms_to_include:
                sn, pn = self.neighbors(n)
                cn[n] = (sn,)
                cpn[n] = (pn,)
        else:
            sub._conformers = [{n: c[n] for n in atoms_to_include if n in c} for c in self._conformers]
            # Hybridizations are already calculated by add_bond, but we can re-ensure.
            for n in sub._atoms:
                sub._calc_hybridization(n)

        return sub

    def augmented_substructure(self, atoms, deep: int = 1):
        atoms_set = set(atoms)
        bonds = self._bonds

        for _ in range(deep):
            n = {y for x in atoms_set for y in bonds.get(x, {})} | atoms_set
            if n == atoms_set:
                break
            atoms_set = n
        return self.substructure(atoms_set)

    def union(self, other, *, remap: bool = False, copy: bool = True):
        from . import molecule
        if isinstance(other, (CGRContainer, molecule.MoleculeContainer)):
            if self._atoms.keys() & other._atoms.keys():
                if not remap:
                    raise MappingError('mapping of graphs is not disjoint')
                other = other.copy()
                other.remap({n: i for i, n in enumerate(other, start=max(self._atoms, default=0) + 1)})
            else:
                other = other.copy()

            u = self.copy() if copy else self
            u._atoms.update(other._atoms)
            u._bonds.update(other._bonds)

            if isinstance(other, CGRContainer):
                u._charges.update(other._charges)
                u._radicals.update(other._radicals)
                u._p_charges.update(other._p_charges)
                u._p_radicals.update(other._p_radicals)
                u._hybridizations.update(other._hybridizations)
                u._p_hybridizations.update(other._p_hybridizations)
                u._plane.update(other._plane)
            elif isinstance(other, molecule.MoleculeContainer):
                oc = {n: a.charge for n, a in other._atoms.items()}
                or_ = {n: a.is_radical for n, a in other._atoms.items()}
                op = {n: getattr(a, 'xy', (0.0, 0.0)) for n, a in other._atoms.items()}
                u._charges.update(oc)
                u._radicals.update(or_)
                u._p_charges.update(oc) # Product state takes reactant state
                u._p_radicals.update(or_)
                u._hybridizations.update(other._hybridizations if hasattr(other, '_hybridizations') else {})
                u._p_hybridizations.update(other._hybridizations if hasattr(other, '_hybridizations') else {})
                u._plane.update(op)

            u._conformers.clear()
            if not copy:
                u.flush_cache()
            return u
        else:
            raise TypeError('CGRContainer or MoleculeContainer expected')

    def compose(self, other: Union['molecule.MoleculeContainer', 'CGRContainer']) -> 'CGRContainer':
        """
        compose 2 graphs to CGR

        :param other: Molecule or CGR Container
        :return: CGRContainer
        """
        from . import molecule
        sa = self._atoms
        spc = self._p_charges
        spr = self._p_radicals
        sp = self._plane
        sb = self._bonds

        bonds = []
        adj: Dict[int, Dict[int, List[Optional[int]]]] = defaultdict(lambda: defaultdict(lambda: [None, None]))
        h = self.__class__()
        atoms = h._atoms

        if isinstance(other, molecule.MoleculeContainer):
            oa = other._atoms
            oc = {n: a.charge for n, a in oa.items()}
            or_ = {n: a.is_radical for n, a in oa.items()}
            op = {n: getattr(a, 'xy', (0.0, 0.0)) for n, a in oa.items()}
            ob = other._bonds
            common = sa.keys() & other

            for n in sa.keys() - common:
                h.add_atom(sa[n].copy(), n, charge=sa[n].charge, is_radical=sa[n].is_radical, xy=sp[n], p_charge=spc.get(n, 0), p_is_radical=spr.get(n, False))
                for m, bond in sb[n].items():
                    if m not in atoms:
                        if m in common:
                            order = bond.order
                            if order:
                                bond = object.__new__(DynamicBond)
                                bond._DynamicBond__order, bond._DynamicBond__p_order = order, None
                                bonds.append((n, m, bond))
                        else:
                            bonds.append((n, m, bond))
            for n in other._atoms.keys() - common:
                h.add_atom(oa[n], n, charge=oc[n], is_radical=or_[n], xy=op[n], p_charge=oc[n], p_is_radical=or_[n])
                for m, bond in ob[n].items():
                    if m not in atoms:
                        if m in common:
                            order = bond.order
                            bond = object.__new__(DynamicBond)
                            bond._DynamicBond__order, bond._DynamicBond__p_order = None, order
                        bonds.append((n, m, bond))
            for n in common:
                an = adj[n]
                for m, bond in ob[n].items():
                    if m in common:
                        an[m][1] = bond.order
                for m, bond in sb[n].items():
                    if m in an or m in common and bond.order:
                        an[m][0] = bond.order
            for n in common:
                san = sa[n]
                if san.atomic_number != oa[n].atomic_number or san.isotope != oa[n].isotope:
                    raise MappingError(f'atoms with number {{{n}}} not equal')
                h.add_atom(san.copy(), n, charge=san.charge, is_radical=san.is_radical, xy=sp[n], p_charge=oc[n], p_is_radical=or_[n])
                for m, (o1, o2) in adj[n].items():
                    if m not in atoms:
                        bond = object.__new__(DynamicBond)
                        bond._DynamicBond__order, bond._DynamicBond__p_order = o1, o2
                        bonds.append((n, m, bond))
        elif isinstance(other, CGRContainer):
            oa = other._atoms
            oc = other._charges
            or_ = other._radicals
            opc = other._p_charges
            opr = other._p_radicals
            op = other._plane
            ob = other._bonds
            common = sa.keys() & other

            for n in sa.keys() - common:
                h.add_atom(sa[n].copy(), n, charge=sa[n].charge, is_radical=sa[n].is_radical, xy=sp[n], p_charge=spc.get(n, 0), p_is_radical=spr.get(n, False))
                for m, bond in sb[n].items():
                    if m not in atoms:
                        if m in common:
                            order = bond.order
                            if order:
                                bond = object.__new__(DynamicBond)
                                bond._DynamicBond__order, bond._DynamicBond__p_order = order, None
                                bonds.append((n, m, bond))
                        else:
                            bonds.append((n, m, bond))
            for n in other._atoms.keys() - common:
                h.add_atom(oa[n].copy(), n, charge=oc[n], is_radical=or_[n], xy=op[n], p_charge=opc[n],
                           p_is_radical=opr[n])
                for m, bond in ob[n].items():
                    if m not in atoms:
                        if m in common:
                            order = bond.p_order
                            if order:
                                bond = object.__new__(DynamicBond)
                                bond._DynamicBond__order, bond._DynamicBond__p_order = None, order
                                bonds.append((n, m, bond))
                        else:
                            bonds.append((n, m, bond))
            for n in common:
                an = adj[n]
                for m, bond in sb[n].items():
                    if m in common and bond.order:
                        an[m][0] = bond.order
                for m, bond in ob[n].items():
                    if m in an or m in common and bond.p_order:
                        an[m][1] = bond.p_order
            for n in common:
                san = sa[n]
                if san.atomic_number != oa[n].atomic_number or san.isotope != oa[n].isotope:
                    raise MappingError(f'atoms with number {{{n}}} not equal')
                h.add_atom(san.copy(), n, charge=san.charge, is_radical=san.is_radical, xy=sp[n], p_charge=opc[n], p_is_radical=opr[n])
                for m, (o1, o2) in adj[n].items():
                    if m not in atoms:
                        bond = object.__new__(DynamicBond)
                        bond._DynamicBond__order, bond._DynamicBond__p_order = o1, o2
                        bonds.append((n, m, bond))
        else:
            raise TypeError('MoleculeContainer or CGRContainer expected')

        for n, m, bond in bonds:
            h.add_bond(n, m, bond)
        return h

    def __xor__(self, other):
        """
        G ^ H is CGR generation
        """
        return self.compose(other)

    def get_mapping(self, other: 'CGRContainer', **kwargs):
        if isinstance(other, CGRContainer):
            return self._get_mapping(other, **kwargs)
        raise TypeError('CGRContainer expected')

    def get_mcs_mapping(self, other: 'CGRContainer', **kwargs):
        if isinstance(other, CGRContainer):
            return super().get_mcs_mapping(other, **kwargs)
        raise TypeError('CGRContainer expected')

    def decompose(self) -> Tuple['molecule.MoleculeContainer', 'molecule.MoleculeContainer']:
        """
        decompose CGR to pair of Molecules, which represents reactants and products state of reaction

        :return: tuple of two molecules
        """
        from . import molecule
        p_charges = self._p_charges
        p_radicals = self._p_radicals
        plane = self._plane

        reactants = molecule.MoleculeContainer()
        products = molecule.MoleculeContainer()

        for n, atom_obj in self._atoms.items():
            atom_template = Element.from_atomic_number(atom_obj.atomic_number)(atom_obj.isotope)
            # Safely get coordinates, defaulting to (0.0, 0.0) if not found
            xy_coords = plane.get(n, (0.0, 0.0)) 
            reactants.add_atom(atom_template, n, charge=atom_obj.charge, is_radical=atom_obj.is_radical, xy=xy_coords)
            products.add_atom(atom_template.copy(), n, charge=p_charges[n], is_radical=p_radicals[n], xy=xy_coords)

        for n, m, bond in self.bonds():
            if bond.order:
                reactants.add_bond(n, m, bond.order)
            if bond.p_order:
                products.add_bond(n, m, bond.p_order)
        return reactants, products

    def __invert__(self):
        """
        decompose CGR
        """
        return self.decompose()

    def _calc_hybridization(self, n: int):
        hybridization = p_hybridization = 1
        for bond in self._bonds[n].values():
            order = bond.order
            p_order = bond.p_order
            if order and hybridization != 4:
                if order == 4:
                    hybridization = 4
                elif order == 3:
                    if hybridization != 3:
                        hybridization = 3
                elif order == 2:
                    if hybridization == 2:
                        hybridization = 3
                    elif hybridization == 1:
                        hybridization = 2
            if p_order and p_hybridization != 4:
                if p_order == 4:
                    p_hybridization = 4
                elif p_order == 3:
                    if p_hybridization != 3:
                        p_hybridization = 3
                elif p_order == 2:
                    if p_hybridization == 2:
                        p_hybridization = 3
                    elif p_hybridization == 1:
                        p_hybridization = 2
        self._hybridizations[n] = hybridization
        self._p_hybridizations[n] = p_hybridization

    def __getstate__(self):
        return {slot: getattr(self, slot) for slot in self.__slots__ if hasattr(self, slot)}

    def __setstate__(self, state):
        for slot, value in state.items():
            setattr(self, slot, value)

    def __iter__(self):
        return iter(self._atoms)

    def __len__(self):
        return len(self._atoms)

    def __copy__(self):
        return self.copy()

    def __or__(self, other):
        return self.union(other, remap=True)

    def __ior__(self, other):
        return self.union(other, remap=True, copy=False)
    
    @cached_property
    def center_atoms(self) -> Tuple[int, ...]:
        """ Get list of atoms of reaction center (atoms with dynamic: bonds, charges, radicals).
        """
        center = {n for n, a in self._atoms.items() if a.is_dynamic}
        center.update(n for n, m_bond in self._bonds.items() if any(bond.is_dynamic for bond in m_bond.values()))
        return tuple(center)


__all__ = ['CGRContainer']
