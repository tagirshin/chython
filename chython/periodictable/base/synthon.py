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
"""Synthon atom family: an Element carrying a Synt-On reaction-centre label."""

from typing import Optional
from .element import Element


# token (sigil stripped) -> (role, attachment points, mechanism, glyph).
# Eight labels, not Synt-On's nine: code 11 ("electrophilic nitrogen") collapses into 'elec',
# because marksCombinations has no N:10 key, so on nitrogen "electrophile" has one meaning already.
LABEL_TABLE = {
    'elec': ('electrophile', 1, 'polar', 'E'),
    'nuc': ('nucleophile', 1, 'polar', 'Nu'),
    'elec2': ('electrophile', 2, 'polar', 'E2'),
    'nuc2': ('nucleophile', 2, 'polar', 'Nu2'),
    'neut2': ('neutral', 2, 'metathesis', 'X2'),
    'elec*': ('electrophile', 1, 'radical', 'E.'),
    'nuc*': ('nucleophile', 1, 'radical', 'Nu.'),
    'elecB': ('electrophile', 1, 'polar', 'E(B)'),
}
SYNTHON_LABELS = frozenset(LABEL_TABLE)
BIVALENT_LABELS = frozenset(k for k, v in LABEL_TABLE.items() if v[1] == 2)
ROLE_COLOR = {'electrophile': '#7D2E8D', 'nucleophile': '#1B7837', 'neutral': '#4A4A4A'}
# __hash__ MUST NOT hash the token string: str hashing is PYTHONHASHSEED-randomised, so the
# canonical SMILES of a labelled molecule would differ between processes. This index is the
# ONE integer at runtime; it never surfaces above __hash__.
LABEL_INDEX = {t: i for i, t in enumerate(LABEL_TABLE)}


class Synthon(Element):
    """Element carrying one Synt-On reaction-centre label.

    Derives Element so the per-element subclasses have a single slot layout - a bare object
    mixin raises "multiple bases have instance lay-out conflict".
    """

    __slots__ = ('_label',)

    def __init__(self, *args, synthon_label: Optional[str] = None, **kwargs):
        super().__init__(*args, **kwargs)
        self.label = synthon_label

    @classmethod
    def from_symbol(cls, symbol: str) -> type['Synthon']:
        """Synthon{El} class by element symbol. Mirrors Element.from_symbol."""
        from .. import synthon_elements

        try:
            return synthon_elements[symbol]
        except KeyError:
            raise ValueError(f'Element with symbol "{symbol}" not found') from None

    @classmethod
    def from_atomic_number(cls, number: int) -> type['Synthon']:
        """Synthon{El} class by atomic number. Mirrors Element.from_atomic_number."""
        return cls.from_symbol(Element.from_atomic_number(number).__name__)

    @property
    def atomic_symbol(self) -> str:
        return self.__class__.__name__[7:]  # mirrors query.py, len('Synthon') == 7

    @property
    def label(self) -> Optional[str]:
        """The synthon token without its sigil: 'elec', 'nuc2', 'elecB'. Never an int."""
        return self._label

    @label.setter
    def label(self, value):
        if value is not None:
            if not isinstance(value, str):
                raise TypeError('synthon label should be the token str, e.g. "elec"')
            if value not in LABEL_TABLE:
                raise ValueError(f'unknown synthon label {value!r}')
        self._label = value

    @property
    def attachment_points(self) -> int:
        """Bond count the label implies. Hard-coded in Synt-On SyntOn.py:539-556."""
        return LABEL_TABLE[self._label][1] if self._label is not None else 0

    @property
    def role(self) -> Optional[str]:
        return LABEL_TABLE[self._label][0] if self._label is not None else None

    @property
    def mechanism(self) -> Optional[str]:
        return LABEL_TABLE[self._label][2] if self._label is not None else None

    @property
    def glyph(self) -> Optional[str]:
        return LABEL_TABLE[self._label][3] if self._label is not None else None

    def copy(self, **kwargs) -> 'Synthon':
        copy = super().copy(**kwargs)
        copy._label = self._label
        return copy

    def __hash__(self):
        # Element.__hash__ seeds Morgan (algorithms/morgan.py), which fixes the canonical
        # traversal order. label=None must hash EXACTLY like a plain Element, or an unlabelled
        # synthon canonicalises differently from the same molecule as a MoleculeContainer.
        # Hash the LABEL_INDEX, never the token string - str hashing is PYTHONHASHSEED-randomised.
        if self._label is None:
            return super().__hash__()
        return hash((super().__hash__(), LABEL_INDEX[self._label]))

    def __eq__(self, other):
        # the matcher evaluates `s_atom == o_atom` (algorithms/isomorphism.py).
        # symmetric wildcard: refuse only when BOTH sides carry a KNOWN and DIFFERENT token.
        # ponytail: this makes atom __eq__ non-transitive; no consumer needs transitivity
        # (is_equal has zero call sites). Upgrade path: move the constraint to the query layer.
        if (
            self._label is not None
            and isinstance(other, Synthon)
            and other._label is not None
            and self._label != other._label
        ):
            return False
        return super().__eq__(other)

    def __repr__(self):
        if self._label is None:
            return super().__repr__()
        return f'{self.__class__.__name__}(synthon_label={self._label!r})'


__all__ = ['Synthon', 'LABEL_TABLE', 'LABEL_INDEX', 'SYNTHON_LABELS', 'BIVALENT_LABELS', 'ROLE_COLOR']
