# -*- coding: utf-8 -*-
#
#  Copyright 2018-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import Dict, Tuple, Union
from .bonds import Bond, QueryBond
from .graph import Graph
from ..algorithms.calculate2d import Calculate2DQuery
from ..algorithms.depict import DepictQuery
from ..algorithms.isomorphism import QueryIsomorphism
from ..algorithms.smarts import Smarts
from ..periodictable import Element, QueryElement
from ..periodictable.base import Query


class QueryContainer(Graph[Query, QueryBond], QueryIsomorphism, Smarts, DepictQuery, Calculate2DQuery):
    __slots__ = ('_smarts', '_plane')

    def __init__(self, smarts: str):
        self._plane: Dict[int, Tuple[float, float]] = {}
        super().__init__()
        self._smarts = smarts

    def add_atom(self, atom: Union[Query, Element, int, str], *args, xy=None, **kwargs):
        if not isinstance(atom, Query):
            # set only basic labels: charge, radical, isotope. use Query object directly for the full control.
            if isinstance(atom, Element):
                atom = QueryElement.from_atom(atom)
            elif isinstance(atom, str):
                atom = QueryElement.from_symbol(atom)()
            elif isinstance(atom, int):
                atom = QueryElement.from_atomic_number(atom)()
            else:
                raise TypeError('QueryElement object expected')
        n = super().add_atom(atom, *args, **kwargs)
        if xy is not None:
            self._plane[n] = xy
        return n

    def add_bond(self, n, m, bond: Union[QueryBond, Bond, int, Tuple[int, ...]]):
        if isinstance(bond, Bond):
            bond = QueryBond.from_bond(bond)
        elif not isinstance(bond, QueryBond):
            bond = QueryBond(bond)
        super().add_bond(n, m, bond)

    def delete_atom(self, n: int):
        super().delete_atom(n)
        self._plane.pop(n, None)

    def copy(self, **kwargs) -> 'QueryContainer':
        copy = super().copy(**kwargs)
        copy._plane = self._plane.copy()
        copy._smarts = self._smarts
        return copy

    def union(self, other: 'QueryContainer', *, remap: bool = False, copy: bool = True) -> 'QueryContainer':
        if not isinstance(other, QueryContainer):
            raise TypeError('QueryContainer expected')
        u = super().union(other, remap=remap, copy=copy)
        if copy:
            u._plane = {**self._plane, **other._plane}
        else:
            self._plane.update(other._plane)
        return u


__all__ = ['QueryContainer']
