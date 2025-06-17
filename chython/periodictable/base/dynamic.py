# -*- coding: utf-8 -*-
#
#  Copyright 2020-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import Type, Optional, Union # Added Union
from .element import Element
from .core import Core
from ...exceptions import IsNotConnectedAtom


class Dynamic(Core):
    __slots__ = ()

    @Core.charge.setter
    def charge(self, charge):
        try:
            g = self._graph()
            g._charges[self._map] = g._validate_charge(charge)
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom

    @Core.is_radical.setter
    def is_radical(self, is_radical):
        try:
            g = self._graph()
            g._radicals[self._map] = g._validate_radical(is_radical)
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def p_charge(self) -> int:
        try:
            return self._graph()._p_charges[self._map]
        except AttributeError:
            raise IsNotConnectedAtom

    @p_charge.setter
    def p_charge(self, charge):
        try:
            g = self._graph()
            g._p_charges[self._map] = g._validate_charge(charge)
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def p_is_radical(self) -> bool:
        try:
            return self._graph()._p_radicals[self._map]
        except AttributeError:
            raise IsNotConnectedAtom

    @p_is_radical.setter
    def p_is_radical(self, is_radical):
        try:
            g = self._graph()
            g._p_radicals[self._map] = g._validate_radical(is_radical)
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def p_hybridization(self):
        """
        Product state hybridization of atom
        """
        try:
            return self._graph()._p_hybridizations[self._map]
        except AttributeError:
            raise IsNotConnectedAtom


class DynamicElement: # Removed ABC
    __slots__ = ('_atomic_number', '_charge', '_is_radical', '_p_charge', '_p_is_radical', '_isotope')

    def __init__(self, atomic_number: int, isotope: Optional[int]):
        self._atomic_number = atomic_number
        self._isotope = isotope
        self._charge = 0 # Default reactant charge
        self._p_charge = 0 # Default product charge
        self._is_radical = False # Default reactant radical state
        self._p_is_radical = False # Default product radical state

    @property
    def isotope(self):
        return self._isotope

    @property
    def atomic_symbol(self) -> str:
        # To get the symbol, we refer to the base Element class
        try:
            # Instantiate the base element class to access its instance properties like atomic_symbol
            base_element_instance = Element.from_atomic_number(self._atomic_number)()
            return base_element_instance.atomic_symbol
        except ValueError:
            return "X" # Fallback for unknown atomic number

    @property
    def atomic_number(self) -> int:
        """
        Element number
        """
        return self._atomic_number

    @classmethod
    def from_symbol(cls, symbol: str, isotope: Optional[int] = None) -> 'DynamicElement':
        """
        Create DynamicElement instance from atomic symbol.
        """
        base_element_class = Element.from_symbol(symbol)
        # Element.from_symbol returns a class, so access its atomic_number property
        # This requires atomic_number to be a classmethod/property on Element subclasses or accessible via an instance.
        # Let's assume Element.from_symbol(symbol)().atomic_number works, or atomic_number is a class var.
        # A safer way: get atomic_number from an instance of the base element class.
        atomic_num = base_element_class().atomic_number # Instantiate to get atomic_number if it's instance property
        return cls(atomic_number=atomic_num, isotope=isotope)

    @classmethod
    def from_atomic_number(cls, atomic_number: int, isotope: Optional[int] = None) -> 'DynamicElement':
        """
        Create DynamicElement instance from atomic number.
        """
        # Basic validation (can be expanded, e.g. checking against known elements)
        if not isinstance(atomic_number, int) or atomic_number < 0: # 0 might be for AnyElement, check conventions
            raise ValueError(f'Invalid atomic number: {atomic_number}')
        return cls(atomic_number=atomic_number, isotope=isotope)

    @classmethod
    def from_atom(cls, atom: Union[Element, 'DynamicElement']) -> 'DynamicElement': # Added DynamicElement to Union
        """
        Create DynamicElement from Element or another DynamicElement.
        Copies reactant state from source, product state defaults.
        """
        if isinstance(atom, DynamicElement):
            new_de = cls(atomic_number=atom.atomic_number, isotope=atom.isotope)
            new_de._charge = atom.charge
            new_de._is_radical = atom.is_radical
            # _p_charge and _p_is_radical will use defaults from __init__ (0, False)
            return new_de
        elif isinstance(atom, Element):
            new_de = cls(atomic_number=atom.atomic_number, isotope=atom.isotope)
            # For Element, reactant and product states mirror the Element's state initially
            new_de._charge = atom.charge
            new_de._is_radical = atom.is_radical
            new_de._p_charge = atom.charge 
            new_de._p_is_radical = atom.is_radical
            return new_de
        else:
            raise TypeError(f'Expected Element or DynamicElement, got {type(atom)}')

    @classmethod
    def from_atoms(cls, atom1: Element, atom2: Element) -> 'DynamicElement': # Kept original Element type hints
        """
        Create DynamicElement from a pair of Element objects, representing reactant and product states.
        """
        if not isinstance(atom1, Element) or not isinstance(atom2, Element):
            raise TypeError('Element instances expected for atom1 and atom2')
        if atom1.atomic_number != atom2.atomic_number:
            raise ValueError('Elements must be of the same atomic type for from_atoms')
        # Isotope consistency check might be desired too, depending on use case.
        # if atom1.isotope != atom2.isotope:
        #     raise ValueError('Elements must be of the same isotope for from_atoms')
        
        new_de = cls(atomic_number=atom1.atomic_number, isotope=atom1.isotope) # or atom2.isotope, should be consistent
        new_de._charge = atom1.charge
        new_de._is_radical = atom1.is_radical
        new_de._p_charge = atom2.charge
        new_de._p_is_radical = atom2.is_radical
        return new_de

    @property
    def charge(self) -> int:
        return self._charge

    @charge.setter
    def charge(self, value: int):
        self._charge = value

    @property
    def is_radical(self) -> bool:
        return self._is_radical

    @is_radical.setter
    def is_radical(self, value: bool):
        self._is_radical = value

    @property
    def p_charge(self) -> int:
        return self._p_charge

    @property
    def p_is_radical(self) -> bool:
        return self._p_is_radical

    def __eq__(self, other):
        if not isinstance(other, DynamicElement):
            return False
        return (self.atomic_number == other.atomic_number and
                self.isotope == other.isotope and
                self.charge == other.charge and self.is_radical == other.is_radical and
                self.p_charge == other.p_charge and self.p_is_radical == other.p_is_radical)

    def __hash__(self):
        return hash((self.atomic_number, self.isotope or 0, # Use 0 for None isotope in hash for consistency
                     self.charge, self.p_charge,
                     self.is_radical, self.p_is_radical))

    @property
    def is_dynamic(self) -> bool:
        """
        Atom has dynamic features (reactant state differs from product state).
        """
        return self.charge != self.p_charge or self.is_radical != self.p_is_radical

    def copy(self, *args, **kwargs): # Added *args, **kwargs to match potential Graph.copy call
        """Create a new copy of this DynamicElement."""
        # Uses the new __init__ that takes atomic_number
        new_copy = self.__class__(atomic_number=self.atomic_number, isotope=self.isotope) 
        new_copy._charge = self.charge
        new_copy._is_radical = self.is_radical
        new_copy._p_charge = self.p_charge
        new_copy._p_is_radical = self.p_is_radical
        return new_copy

    def __copy__(self):
        return self.copy()


__all__ = ['DynamicElement']
