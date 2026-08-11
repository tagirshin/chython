# -*- coding: utf-8 -*-
#
#  Copyright 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
"""Molecule whose atoms may carry a Synt-On reaction-centre label."""

from math import cos, sin
from re import compile
from typing import Optional, Union

from ..algorithms.depict import _optimal_label_direction, _render_config
from ..periodictable import Element
from ..periodictable.base.element import Element as _Element
from ..periodictable.base.synthon import ROLE_COLOR, Synthon
from .molecule import MoleculeContainer


_label_insert = compile(r'(:[0-9]{1,4})?\]$')


def _glyph_svg(atom, size):
    """Glyph as SVG text content: superscript = bond count, dot = radical mechanism."""
    g = atom.glyph
    if g.endswith('2'):
        return f'{g[:-1]}<tspan dy="-{0.5 * size:.2f}" font-size="{0.7 * size:.2f}">2</tspan>'
    if g.endswith('.'):
        return f'{g[:-1]}&#183;'
    return g


class SynthonContainer(MoleculeContainer):
    """MoleculeContainer whose atoms are Synthon{El} and may carry a reaction-centre label."""

    __slots__ = ()

    def add_atom(self, atom: Union[Synthon, Element, int, str], *args, synthon_label: Optional[str] = None, **kwargs):
        if isinstance(atom, str):
            atom = Synthon.from_symbol(atom)(synthon_label=synthon_label)
        elif isinstance(atom, int):
            atom = Synthon.from_atomic_number(atom)(synthon_label=synthon_label)
        elif isinstance(atom, _Element) and not isinstance(atom, Synthon):
            # add_atom re-assigns charge/is_radical from kwargs, so they travel as kwargs
            kwargs.setdefault('charge', atom.charge)
            kwargs.setdefault('is_radical', atom.is_radical)
            atom = Synthon.from_symbol(atom.atomic_symbol)(
                atom.isotope, implicit_hydrogens=atom.implicit_hydrogens, synthon_label=synthon_label
            )
        elif synthon_label is not None:
            atom.label = synthon_label
        return super().add_atom(atom, *args, **kwargs)

    def _format_atom(self, n, adjacency, **kwargs):
        s = super()._format_atom(n, adjacency, **kwargs)
        # getattr: explicify_hydrogens(), BaseReactor._patcher and union() with a plain molecule
        # all put plain Elements in here.
        label = getattr(self._atoms[n], '_label', None)
        if label is None or not kwargs.get('synthon_label', True):
            return s
        if not s.startswith('['):  # force brackets so the label has somewhere to live
            s = super()._format_atom(n, adjacency, **{**kwargs, 'hydrogens': True})
        return _label_insert.sub(lambda m: f'_{label}{m.group(1) or ""}]', s)

    @property
    def synthon_labels(self) -> dict[int, str]:
        """{atom number: token str} for every labelled atom."""
        return {n: a._label for n, a in self.atoms() if getattr(a, '_label', None) is not None}

    def unlabelled(self) -> 'SynthonContainer':
        """Copy with every label cleared, so matching it is label-blind again.

        The escape hatch for algorithms the reference defines as mark-blind (PAS, and the
        Ketones re-classification hook).
        """
        copy = self.copy()
        for _, a in copy.atoms():
            if isinstance(a, Synthon):
                a._label = None
        copy.flush_cache()
        return copy

    # ponytail: DepictMolecule.__render_atoms is name-mangled, so the override must spell the
    # mangled name. Upgrade path: rename it `_render_atoms` upstream and drop the prefix here.
    def _DepictMolecule__render_atoms(self, uid):
        svg, define, mask = super()._DepictMolecule__render_atoms(uid)
        labelled = [(n, a) for n, a in self.atoms() if getattr(a, '_label', None) is not None]
        if not labelled:
            return svg, define, mask
        size = _render_config['other_size']
        offset = 0.9 * _render_config['font_size']
        mono = _render_config['monochrome']
        define.append(
            f'      <g id="{uid}-synthon" font-size="{size:.2f}" font-family="{_render_config["other_font_style"]}">'
        )
        for n, a in labelled:
            x, y = a.x, -a.y
            # the helper Depict already uses for its neighbour/hybridization labels
            angle = _optimal_label_direction(x, y, [(self._atoms[m].x, -self._atoms[m].y) for m in self._bonds[n]])
            dx, dy = cos(angle) * offset, sin(angle) * offset
            fill = 'black' if mono else ROLE_COLOR[a.role]
            define.append(
                f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx:.2f}" dy="{dy:.2f}" '
                f'text-anchor="middle" fill="{fill}">{_glyph_svg(a, size)}</text>'
            )
        define.append('      </g>')
        svg.append(f'    <use xlink:href="#{uid}-synthon"/>')
        if mask:
            mask.insert(-1, f'          <use xlink:href="#{uid}-synthon" stroke-width="{size * 0.1:.2f}"/>')
        return svg, define, mask


def restore_synthons(structure: SynthonContainer, synthon_labels: Optional[dict[int, str]] = None) -> SynthonContainer:
    """Re-wrap plain Elements left by BaseReactor._patcher as Synthon{El}, then stamp tokens.

    reactor/base.py rebuilds every atom the rule names as a plain Element, copying exactly three
    attributes, so the re-wrap is needed even for a completely unlabelled rule.
    """
    atoms = structure._atoms
    for n, a in atoms.items():
        if not isinstance(a, Synthon):
            new = Synthon.from_symbol(a.atomic_symbol)()
            new.__setstate__(a.__getstate__())  # chython's own slot-transfer protocol
            new._label = None
            atoms[n] = new
    if synthon_labels:
        for n, token in synthon_labels.items():
            atoms[n].label = token
    structure.flush_cache()
    return structure


__all__ = ['SynthonContainer', 'restore_synthons']
