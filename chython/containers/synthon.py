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

from math import cos, hypot, sin
from re import compile
from typing import Optional, Union

from ..algorithms.depict import _optimal_label_direction, _render_config
from ..periodictable import Element
from ..periodictable.base.element import Element as _Element
from ..periodictable.base.synthon import ROLE_COLOR, Synthon
from .molecule import MoleculeContainer


_label_insert = compile(r'(:[0-9]{1,4})?\]$')
# the base renderer's own symbol lines, the only <text id=...> it puts in `define`
_symbol_text = compile(r'<text id="[-0-9a-f]+-([0-9]+)"')
_dx_attr = compile(r'dx="(-?[0-9.]+)"')
_tags = compile(r'<[^>]*>')


def _glyph_svg(atom, size):
    """Glyph as SVG text content: superscript = bond count, dot = radical mechanism."""
    g = atom.glyph
    if g.endswith('2'):
        # typographic superscript: ~1/3 em up, not half, or the 2 reads as a separate glyph
        return f'{g[:-1]}<tspan dy="-{0.35 * size:.2f}" font-size="{0.75 * size:.2f}">2</tspan>'
    if g.endswith('.'):
        return f'{g[:-1]}&#183;'
    return g


def _bond_dot(mol, n, reach, x, y, font_size, fill, lead):
    """Dot marking the atom, always left of the glyph run, or None when nowhere is left.

    The symbol, its hydrogen count and a trailing tag all grow rightwards, so a rightward bond is
    buried under text for most of its length; the left edge is the one side reliably clear. Rides
    a leftward bond when one has room, otherwise sits straight out from the glyph edge - except
    under a leading tag, which holds that edge itself and would push the fallback a whole tag
    width into open space. There the tag is the only marker the atom gets.
    """
    r = .16 * font_size
    need = reach[0] + r  # reach[0] already carries the tag width when the tag leads
    best = None
    for m in mol._bonds[n]:
        px, py = mol._atoms[m].x, -mol._atoms[m].y
        length = hypot(px - x, py - y)
        if px <= x and need < .45 * length and (best is None or px - x < best[0]):
            best = (px - x, (px - x) * need / length, (py - y) * need / length)
    if best is None and lead:
        return None
    dx, dy = best[1:] if best else (-need, 0.)
    return f'    <circle cx="{x + dx:.2f}" cy="{y + dy:.2f}" r="{r:.2f}" fill="{fill}"/>'


def _tag(atom, size, font, fill, up):
    """The label as a raised tspan, to be spliced into the base renderer's symbol <text>."""
    return (f'<tspan dy="-{up:.2f}" font-size="{size:.2f}" font-family="{font}" fill="{fill}">'
            f'{_glyph_svg(atom, size)}</tspan>')


def _crowded_right(mol, n, drawn, x, y):
    """Is another atom's rendered symbol sitting where a trailing tag would land?"""
    return any(px > x and abs(py - y) < .6 and hypot(px - x, py - y) < 1.6
               for m in drawn if m != n for px, py in [(mol._atoms[m].x, -mol._atoms[m].y)])


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
        font_size = _render_config['font_size']
        font = _render_config['other_font_style']
        mono = _render_config['monochrome']
        # which define line carries each atom's element symbol; absent = the base drew a bare vertex
        drawn = {int(m[1]): i for i, line in enumerate(define) if (m := _symbol_text.search(line))}
        labels, dots = [], []
        for n, a in labelled:
            x, y = a.x, -a.y
            fill = 'black' if mono else ROLE_COLOR[a.role]
            i = drawn.get(n)
            if i is not None:
                # ride the symbol, the way the hydrogen count already does: one text run, so the
                # label can neither drift off its atom nor land on the glyphs of that atom
                head, body = define[i].split('>', 1)
                body = body[:body.rindex('</text>')]
                up = .4 * font_size
                # ponytail: 0.66 em per character is the base font's average advance, enough to
                # keep the dot off the run. Upgrade path: real text metrics, as for the tag.
                reach = [.4 * font_size, .66 * font_size * len(_tags.sub('', body)) - .4 * font_size]
                wide = .6 * size * len(a.glyph)
                lead = _crowded_right(self, n, drawn, x, y)
                if lead:
                    # lead with the tag and pull the run left by the tag's own width - exact,
                    # because the label font is monospace even though the symbol font is not
                    head = _dx_attr.sub(lambda m: f'dx="{float(m[1]) - wide:.2f}"', head, 1)
                    body = f'{_tag(a, size, font, fill, up)}<tspan dy="{up:.2f}">{body}</tspan>'
                    reach[0] += wide
                else:  # trailing tag: span_dy first undoes the hydrogen-count subscript
                    body += _tag(a, size, font, fill, up + (_render_config['span_dy'] if '<tspan' in body else 0.))
                    reach[1] += wide
                define[i] = f'{head}>{body}</text>'
                # same anchor as a bare vertex, pushed clear of the text-occupied centre
                if (dot := _bond_dot(self, n, reach, x, y, font_size, fill, lead)) is not None:
                    dots.append(dot)
                continue
            # bare vertex: nothing on screen names this atom, so a dot in the label colour does
            dots.append(f'    <circle cx="{x:.2f}" cy="{y:.2f}" r="{.16 * font_size:.2f}" fill="{fill}"/>')
            taken = [(self._atoms[m].x, -self._atoms[m].y) for m in self._bonds[n]]
            if _render_config['mapping']:  # aam number, anchored end so its box hangs further left
                taken.append((x - _render_config['dx_m'] - .3 * _render_config['mapping_size'] * len(str(n)),
                              y + _render_config['dy_m']))
            # the helper Depict already uses for its neighbour/hybridization labels
            angle = _optimal_label_direction(x, y, taken)
            offset = hypot(_render_config['dx_nh'], _render_config['dy_nh'])
            dx, dy = cos(angle) * offset, sin(angle) * offset
            if dy > 0:  # a baseline sits at the glyph's foot, so downwards has to clear a cap height
                dy += .72 * size
            labels.append(
                f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx:.2f}" dy="{dy:.2f}" '
                f'text-anchor="{"end" if dx < 0 else "start"}" fill="{fill}">{_glyph_svg(a, size)}</text>'
            )
        svg.extend(dots)  # only bare vertices also need a text label beside the dot
        if labels:
            define.append(f'      <g id="{uid}-synthon" font-size="{size:.2f}" font-family="{font}">')
            define.extend(labels)
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
