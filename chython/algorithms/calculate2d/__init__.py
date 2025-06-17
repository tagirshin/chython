# -*- coding: utf-8 -*-
#
#  Copyright 2019-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from .molecule import *
from .reaction import *
from math import sqrt
from importlib.resources import files
from random import random
from ...exceptions import ImplementationError
from ...periodictable import Element


class Calculate2D:
    __slots__ = ()

    def clean2d(self):
        """
        Calculate 2d layout of graph. https://pubs.acs.org/doi/10.1021/acs.jcim.7b00425 JS implementation used.
        """
        if ctx is None:
            raise ImportError('py-mini-racer required for clean2d')
        plane = {}
        for _ in range(5):
            smiles, order = self._clean2d_prepare()
            try:
                xy = ctx.call('$.clean2d', smiles)
            except JSEvalException:
                continue
            break
        else:
            raise ImplementationError

        shift_x, shift_y = xy[0]
        for n, (x, y) in zip(order, xy):
            plane[n] = (x - shift_x, shift_y - y)

        bonds = []
        for n, m, _ in self.bonds():
            xn, yn = plane[n]
            xm, ym = plane[m]
            bonds.append(sqrt((xm - xn) ** 2 + (ym - yn) ** 2))
        if bonds:
            bond_reduce = sum(bonds) / len(bonds) / .825
        else:
            bond_reduce = 1.

        self_plane = self._plane
        for n, (x, y) in plane.items():
            self_plane[n] = (x / bond_reduce, y / bond_reduce)

        self.__dict__.pop('__cached_method__repr_svg_', None)

    def _fix_plane_mean(self, shift_x, shift_y=0, component=None):
        plane = self._plane
        if component is None:
            component = plane

        left_atom = min(component, key=lambda x: plane[x][0])
        right_atom = max(component, key=lambda x: plane[x][0])

        min_x = plane[left_atom][0] - shift_x
        if len(self._atoms[left_atom].atomic_symbol) == 2:
            min_x -= .2

        max_x = plane[right_atom][0] - min_x
        min_y = min(plane[x][1] for x in component)
        max_y = max(plane[x][1] for x in component)
        mean_y = (max_y + min_y) / 2 - shift_y
        for n in component:
            x, y = plane[n]
            plane[n] = (x - min_x, y - mean_y)

        if isinstance(self, molecule.MoleculeContainer):
            if -.18 <= plane[right_atom][1] <= .18:
                factor = self._hydrogens[right_atom]
                if factor == 1:
                    max_x += .15
                elif factor:
                    max_x += .25
        return max_x

    def _fix_plane_min(self, shift_x, shift_y=0, component=None):
        plane = self._plane
        if component is None:
            component = plane

        right_atom = max(component, key=lambda x: plane[x][0])
        min_x = min(plane[x][0] for x in component) - shift_x
        max_x = plane[right_atom][0] - min_x
        min_y = min(plane[x][1] for x in component) - shift_y

        for n in component:
            x, y = plane[n]
            plane[n] = (x - min_x, y - min_y)

        if isinstance(self, molecule.MoleculeContainer):
            if shift_y - .18 <= plane[right_atom][1] <= shift_y + .18:
                factor = self._hydrogens[right_atom]
                if factor == 1:
                    max_x += .15
                elif factor:
                    max_x += .25
        return max_x


class Calculate2DQuery(Calculate2D):
    __slots__ = ()

    def _clean2d_prepare(self):
        mol = molecule.MoleculeContainer()
        for n, atom in self._atoms.items():
            atom = Element.from_atomic_number(atom.atomic_number or 6)()
            mol.add_atom(atom, n)
        for n, m, bond in self.bonds():
            mol.add_bond(n, m, bond.order[0])
        mol._hydrogens = {n: 0 for n in mol._hydrogens}
        smiles, order = mol._smiles(lambda x: random(), _return_order=True)
        return ''.join(smiles), order


class Calculate2DCGR(Calculate2D):
    __slots__ = ()

    def _clean2d_prepare(self):
        mol = molecule.MoleculeContainer()
        for n, atom in self._atoms.items():
            atom = Element.from_atomic_number(atom.atomic_number or 6)()
            mol.add_atom(atom, n)
        for n, m, bond in self.bonds():
            mol.add_bond(n, m, bond.order or 1)
        mol._hydrogens = {n: 0 for n in mol._hydrogens}
        smiles, order = mol._smiles(lambda x: random(), _return_order=True)
        return ''.join(smiles), order


try:
    from py_mini_racer import MiniRacer, JSEvalException

    ctx = MiniRacer()
    ctx.eval('const self = this')
    ctx.eval(files(__package__).joinpath('clean2d.js').read_text())
except RuntimeError:
    ctx = None


__all__ = ['Calculate2DMolecule', 'Calculate2DReaction', 'Calculate2DCGR', 'Calculate2DQuery']
