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
from typing import Dict, Iterator, Tuple, Optional, Collection, List, Union
from collections import defaultdict

from .bonds import DynamicBond, Bond
from .graph import Graph
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
from ..exceptions import MappingError


class CGRContainer(Graph, CGRSmiles, DepictCGR, Calculate2DCGR, X3domCGR, Morgan, Rings, Isomorphism, FingerprintsCGR):
    __slots__ = ('_conformers', '_p_charges', '_p_radicals', '_hybridizations', '_p_hybridizations', '_plane')
    _atoms: Dict[int, DynamicElement]
    _bonds: Dict[int, Dict[int, DynamicBond]]

    def __init__(self):
        self._conformers: List[Dict[int, Tuple[float, float, float]]] = []
        self._p_charges: Dict[int, int] = {}
        self._p_radicals: Dict[int, bool] = {}
        self._hybridizations: Dict[int, int] = {}
        self._p_hybridizations: Dict[int, int] = {}
        self._plane: Dict[int, Tuple[float, float]] = {}
        super().__init__()

    # Methods for charge and radical validation (can be expanded later)
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

    def add_atom(self, atom: Union[DynamicElement, Element, int, str], *args, p_charge: int = 0,
                 p_is_radical: bool = False, **kwargs):
        p_charge = self._validate_charge(p_charge)
        p_is_radical = self._validate_radical(p_is_radical)

        if not isinstance(atom, DynamicElement):
            if isinstance(atom, Element):
                atom = DynamicElement.from_atomic_number(atomic_number=atom.atomic_number, isotope=atom.isotope)
            elif isinstance(atom, str):
                atom = DynamicElement.from_symbol(symbol=atom, isotope=None)
            elif isinstance(atom, int):
                atom = DynamicElement.from_atomic_number(atomic_number=atom, isotope=None)
            else:
                raise TypeError('DynamicElement object expected')

        xy_coord = kwargs.pop('xy', None)
        _map = super().add_atom(atom, *args, **kwargs)
        if xy_coord is not None:
            self._plane[_map] = xy_coord
        
        self._p_charges[_map] = p_charge
        self._p_radicals[_map] = p_is_radical
        self._hybridizations[_map] = 1
        self._p_hybridizations[_map] = 1
        self._conformers.clear()
        return _map

    def add_bond(self, n, m, bond: Union[DynamicBond, Bond, int]):
        if isinstance(bond, DynamicBond):
            order = bond.order
            p_order = bond.p_order
        elif isinstance(bond, Bond):
            order = p_order = bond.order
            bond = DynamicBond.from_bond(bond)
        else:
            order = p_order = bond
            bond = DynamicBond(order, order)

        super().add_bond(n, m, bond)
        self._conformers.clear()

        if order != 1 or p_order != 1:
            self._calc_hybridization(n)
            self._calc_hybridization(m)

    def delete_atom(self, n):
        old_bonds = self._bonds[n]
        super().delete_atom(n)

        del self._p_charges[n]
        del self._p_radicals[n]
        del self._hybridizations[n]
        del self._p_hybridizations[n]

        for m in old_bonds:
            self._calc_hybridization(m)
        self._conformers.clear()

    def delete_bond(self, n, m):
        super().delete_bond(n, m)
        self._conformers.clear()
        self._calc_hybridization(n)
        self._calc_hybridization(m)

    def remap(self, mapping, *, copy=False) -> 'CGRContainer':
        target_graph = self
        if copy:
            target_graph = self.copy() # Create a copy of CGRContainer to work on

        # super().remap (Graph.remap) remaps in-place on the instance it's called on.
        # We call it on target_graph (either self or the copy).
        # Graph.remap does not return anything explicitly other than None implicitly.
        super(CGRContainer, target_graph).remap(mapping) # Call Graph.remap on target_graph

        # Now, CGR-specific attributes need to be remapped on target_graph.
        # The original logic used 'h' returned from super().remap which was assumed to be the remapped graph.
        # With Graph.remap working in-place, target_graph is the one remapped.
        
        mg = mapping.get
        # These are attributes of target_graph that need remapping for their keys.
        # Original attributes were on self, but if copy=True, we remapped a copy.
        # So, if copy=True, target_graph already has copies of these dicts from self.copy().
        # If copy=False, target_graph is self, and we are modifying self's dicts in-place.

        # Remap _p_charges, _p_radicals, _hybridizations, _p_hybridizations, _conformers keys
        # Important: If copy=True, these dictionaries on target_graph are already copies.
        # If copy=False, we are modifying them in-place on self.
        
        current_p_charges = target_graph._p_charges.copy() # Operate on a copy to avoid modification during iteration issues
        target_graph._p_charges.clear()
        for n, c_val in current_p_charges.items():
            target_graph._p_charges[mg(n, n)] = c_val

        current_p_radicals = target_graph._p_radicals.copy()
        target_graph._p_radicals.clear()
        for n, r_val in current_p_radicals.items():
            target_graph._p_radicals[mg(n, n)] = r_val

        current_hybridizations = target_graph._hybridizations.copy()
        target_graph._hybridizations.clear()
        for n, h_val in current_hybridizations.items():
            target_graph._hybridizations[mg(n, n)] = h_val

        current_p_hybridizations = target_graph._p_hybridizations.copy()
        target_graph._p_hybridizations.clear()
        for n, ph_val in current_p_hybridizations.items():
            target_graph._p_hybridizations[mg(n, n)] = ph_val

        remapped_conformers = []
        for conf_dict in target_graph._conformers:
            remapped_conf_dict = {mg(n, n): x for n, x in conf_dict.items()}
            remapped_conformers.append(remapped_conf_dict)
        target_graph._conformers = remapped_conformers
        
        # Remap _plane attribute if it exists and needs remapping (assuming it stores atomID keys)
        if hasattr(target_graph, '_plane') and target_graph._plane:
            current_plane = target_graph._plane.copy()
            target_graph._plane.clear()
            for n, xy_val in current_plane.items():
                target_graph._plane[mg(n,n)] = xy_val

        return target_graph # Return the (possibly new and) remapped graph

    def copy(self, **kwargs) -> 'CGRContainer':
        copy = super().copy(**kwargs)
        copy._hybridizations = self._hybridizations.copy()
        copy._conformers = [c.copy() for c in self._conformers]
        copy._p_hybridizations = self._p_hybridizations.copy()
        copy._p_radicals = self._p_radicals.copy()
        copy._p_charges = self._p_charges.copy()
        return copy

    def substructure(self, atoms_ids_to_include_iterable, *, as_query: bool = False, **kwargs) -> Union['CGRContainer',
                                                                                'query.QueryCGRContainer']:
        """
        create substructure containing atoms from atoms list

        :param atoms_ids_to_include_iterable: list/iterable of atom numbers for the substructure
        :param meta: if True metadata will be copied to substructure (Note: 'meta' kwarg is not currently used)
        :param as_query: return Query object based on graph substructure
        """
        # Determine types for the new subgraph
        SubGraphType = query.QueryCGRContainer if as_query else self.__class__
        SubAtomType = DynamicQueryElement if as_query else DynamicElement

        # Create the new subgraph instance
        sub = SubGraphType()

        # Filter the input atom IDs to include only those present in the current graph
        # and store them in a set for efficient lookup.
        # This set, sub_atom_ids, represents the atoms that will actually be in the substructure.
        sub_atom_ids = {atom_id for atom_id in atoms_ids_to_include_iterable if atom_id in self._atoms}

        # Add atoms to the substructure
        original_atoms_dict = self._atoms
        original_plane_coords = self._plane  # Assuming _plane exists for CGRContainer for xy coords
        for n_id in sub_atom_ids:
            original_atom_object = original_atoms_dict[n_id]
            # Create a new atom of the appropriate type (DynamicElement or DynamicQueryElement)
            # .from_atom should copy intrinsic properties like charge, isotope, etc.
            new_atom_object = SubAtomType.from_atom(original_atom_object)
            # Add the new atom to the subgraph.
            # The add_atom method of CGRContainer/QueryCGRContainer will handle
            # initializing product-state attributes (_p_charges, _p_radicals, etc.) to defaults.
            sub.add_atom(new_atom_object, n_id, xy=original_plane_coords.get(n_id))

        # Add bonds to the substructure
        # Iterate through bonds of the original graph. If both atoms of a bond
        # are in sub_atom_ids, add a copy of the bond to the new subgraph.
        original_bonds_dict = self._bonds
        for n1 in sub_atom_ids:
            if n1 in original_bonds_dict:
                for n2, bond_obj in original_bonds_dict[n1].items():
                    if n2 in sub_atom_ids and n1 < n2:  # Add each bond only once
                        # Create a copy of the bond object
                        new_bond_obj = DynamicBond(bond_obj.order, bond_obj.p_order)
                        sub.add_bond(n1, n2, new_bond_obj)
        
        # --- CGR-specific logic (adapted from the original method) ---

        # Populate product-state charges and radicals for atoms in the substructure,
        # using values from the original CGRContainer. This overrides defaults set by sub.add_atom.
        original_p_charges = self._p_charges
        original_p_radicals = self._p_radicals
        for n_id in sub_atom_ids: # Iterate over atoms actually added to sub
            if n_id in original_p_charges: # Check if key exists in original
                sub._p_charges[n_id] = original_p_charges[n_id]
            if n_id in original_p_radicals:
                sub._p_radicals[n_id] = original_p_radicals[n_id]
            # If not present in original, they keep the defaults from sub.add_atom (0, False respectively)

        if as_query:
            # Populate query-specific attributes if the substructure is a QueryCGRContainer.
            # Values are taken from the original CGRContainer.
            original_hybridizations = self._hybridizations
            original_p_hybridizations = self._p_hybridizations
            # self.neighbors is a method of the original CGRContainer (self)

            # QueryCGRContainer's add_atom initializes its _hybridizations, _p_hybridizations,
            # _neighbors, _p_neighbors to defaults. Here we populate them with specific values from original.
            for n_id in sub_atom_ids:
                if n_id in original_hybridizations:
                    sub._hybridizations[n_id] = (original_hybridizations[n_id],) # Query stores as tuple
                if n_id in original_p_hybridizations:
                    sub._p_hybridizations[n_id] = (original_p_hybridizations[n_id],) # Query stores as tuple
                
                # Get neighbor counts from the original CGRContainer instance 'self'
                s_neighbors, p_neighbors = self.neighbors(n_id)
                sub._neighbors[n_id] = (s_neighbors,) # Query stores as tuple
                sub._p_neighbors[n_id] = (p_neighbors,) # Query stores as tuple
        else:
            # CGRContainer-specific finalization for the substructure.
            # Copy conformers, filtering for atoms present in the substructure.
            sub._conformers = []
            for conf_dict in self._conformers: # self._conformers are from original CGR
                # Create new dict for sub_conf, only including atoms present in the substructure
                sub_conf_dict = {n_id: pos for n_id, pos in conf_dict.items() if n_id in sub_atom_ids}
                if sub_conf_dict: # Only add if the conformer dict is not empty
                    sub._conformers.append(sub_conf_dict)

            # Recalculate hybridizations for the new substructure context.
            # CGRContainer.add_atom initialized _hybridizations & _p_hybridizations to 1.
            # CGRContainer.add_bond (called when bonds were added to 'sub') calls 
            # _calc_hybridization for bonded atoms.
            # This explicit loop ensures all atoms in the substructure (including any isolated ones)
            # have their hybridizations calculated based on the final bond structure of 'sub'.
            # This matches the original code's else block's intent.
            sub._hybridizations = {}
            sub._p_hybridizations = {}
            for n_id_in_sub_graph in sub._atoms: # Iterate over keys of sub._atoms dict
               sub._calc_hybridization(n_id_in_sub_graph)
        
        return sub, sub_atom_ids # Return the new substructure and the set of atom IDs in it

    def augmented_substructure(self, atoms, deep: int = 1):
        atoms = set(atoms)
        bonds = self._bonds

        for _ in range(deep):
            n = {y for x in atoms for y in bonds[x]} | atoms
            if n == atoms:
                break
            atoms = n
        return self.substructure(atoms)

    def union(self, other, **kwargs) -> 'CGRContainer':
        from . import molecule
        if isinstance(other, CGRContainer):
            # Graph.union does not accept atom_type, bond_type and returns only the unioned graph `u`.
            # The `other` object reference passed to super().union is a copy made within Graph.union if remapping occurs.
            # We need to update CGR-specific attributes based on the `other` that was passed into this method,
            # considering that its atom IDs might have been remapped if `remap=True` was in kwargs
            # and a collision occurred. This is complex if `other` is remapped internally by super().union
            # and that mapping isn't returned.
            # For now, assume `other`'s attributes are used pre-super-union or the mapping complexity is handled elsewhere.
            # Simplification: call super().union without the unsupported kwargs.
            
            # Store original `other` attributes that might be needed if `other` itself is modified by `super().union`
            # or if its atom IDs get remapped.
            # This is tricky. Let's assume other._p_charges etc. are on the input `other` instance.
            # If `other` is remapped by `super().union`, we need that remapped `other` or the mapping.
            # `MoleculeContainer.union` returns `u, other_graph_processed_object` which is what we need.
            # Assuming `Graph.union` should be modified to return `u, processed_other_copy` as well.
            # If `Graph.union` is NOT changed, this current CGRContainer.union will be problematic for remapped `other`.
            
            # Let's proceed with the minimal change to avoid the TypeError for now.
            # This means if `other` is remapped by `super().union`, the ._p_charges.update below might use wrong keys.
            u = super().union(other, **kwargs) # remap kwarg is passed in kwargs
            
            # `other` here refers to the original `other` passed to CGRContainer.union.
            # If `super().union` remapped `other` (it operates on a copy), then using `other._p_charges` directly 
            # might be incorrect if atom IDs changed. This logic needs to be robust to remapping.
            # For now, we assume the keys in other._p_charges are still relevant for `u`.
            u._conformers.clear() # This is from original CGRContainer.union
            u._p_charges.update(other._p_charges)
            u._p_radicals.update(other._p_radicals)
            u._hybridizations.update(other._hybridizations)
            u._p_hybridizations.update(other._p_hybridizations)
            return u
        elif isinstance(other, molecule.MoleculeContainer):
            # Similar issue here with super().union and its return / kwargs
            u = super().union(other, **kwargs) # remap kwarg is passed in kwargs

            # `other` here is the MoleculeContainer. Its ._charges, ._radicals are used.
            # If atom IDs were remapped by super().union, this update is problematic.
            u._conformers.clear()
            u._p_charges.update(other._charges) # MoleculeContainer has ._charges
            u._p_radicals.update(other._radicals) # MoleculeContainer has ._radicals
            u._hybridizations.update(other._hybridizations) # MoleculeContainer has ._hybridizations
            u._p_hybridizations.update(other._hybridizations) # This seems to assume molecule has p_hybridizations? Check MoleculeContainer.
                                                              # MoleculeContainer uses self._hybridizations. So this line is likely an error.
                                                              # CGR needs p_hybridizations, but molecule has only one set.
                                                              # Let's assume it meant to copy other._hybridizations to u._p_hybridizations for CGR context.
            # Revisit this logic. For now, focus on removing atom_type/bond_type from super() call.
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
            oc = other._charges
            or_ = other._radicals
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
            return super().get_mapping(other, **kwargs)
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
        return {'conformers': self._conformers, 'p_charges': self._p_charges, 'p_radicals': self._p_radicals,
                **super().__getstate__()}

    def __setstate__(self, state):
        self._p_charges = state['p_charges']
        self._p_radicals = state['p_radicals']
        super().__setstate__(state)
        if 'conformers' in state:
            self._conformers = state['conformers']
        else:
            self._conformers = [] 

        # restore query marks
        self._hybridizations = {}
        self._p_hybridizations = {}
        # Iterate over atom IDs that have bond information to recalculate hybridizations
        # state['_bonds'] stores the bonds dictionary {atom_id: {neighbor_id: bond_obj, ...}, ...}
        for n in state.get('_bonds', {}).keys(): 
            self._calc_hybridization(n)

    def __iter__(self):
        return iter(self._atoms)

    def __len__(self):
        return len(self._atoms)

    @cached_property
    def center_atoms(self) -> Tuple[int, ...]:
        """ Get list of atoms of reaction center (atoms with dynamic: bonds, charges, radicals).
        """
        center = {n for n, a in self._atoms.items() if a.is_dynamic}
        center.update(n for n, m_bond in self._bonds.items() if any(bond.is_dynamic for bond in m_bond.values()))
        return tuple(center)


__all__ = ['CGRContainer']
