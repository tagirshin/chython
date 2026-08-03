# -*- coding: utf-8 -*-
#
#  Copyright 2025 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2025 Tagir Akhmetshin <tagirshin@gmail.com>
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
from chython.files.daylight.tokenize import smiles_tokenize, smarts_tokenize
from chython.containers import QueryBond
from chython.exceptions import IncorrectSmarts, IncorrectSmiles
from pytest import mark, raises


def test_smiles_tokenize():
    assert smiles_tokenize('C') == [(0, {'element': 'C'})]
    assert smiles_tokenize('CC') == [(0, {'element': 'C'}), (0, {'element': 'C'})]
    assert smiles_tokenize('C=O') == [(0, {'element': 'C'}), (1, 2), (0, {'element': 'O'})]
    assert smiles_tokenize('C(O)N') == [(0, {'element': 'C'}), (2, None), (0, {'element': 'O'}),
                                        (3, None), (0, {'element': 'N'})]
    assert smiles_tokenize('C2CC2') == [(0, {'element': 'C'}), (6, 2), (0, {'element': 'C'}),
                                        (0, {'element': 'C'}), (6, 2)]


def test_smiles_tokenize_atom():
    assert smiles_tokenize('[NH4+]') == [(0, {'element': 'N', 'isotope': None, 'parsed_mapping': None, 'charge': 1,
                                              'implicit_hydrogens': 4, 'stereo': None})]
    assert smiles_tokenize('[14N]') == [(0, {'element': 'N', 'isotope': 14, 'parsed_mapping': None, 'charge': 0,
                                              'implicit_hydrogens': 0, 'stereo': None})]
    assert smiles_tokenize('[N@H]') == [(0, {'element': 'N', 'isotope': None, 'parsed_mapping': None, 'charge': 0,
                                             'implicit_hydrogens': 1, 'stereo': True})]
    assert smiles_tokenize('[N@@H--]') == [(0, {'element': 'N', 'isotope': None, 'parsed_mapping': None, 'charge': -2,
                                                'implicit_hydrogens': 1, 'stereo': False})]
    assert smiles_tokenize('[N@+3]') == [(0, {'element': 'N', 'isotope': None, 'parsed_mapping': None, 'charge': 3,
                                              'implicit_hydrogens': 0, 'stereo': True})]
    assert smiles_tokenize('[CH2:2]') == [(0, {'element': 'C', 'isotope': None, 'parsed_mapping': 2, 'charge': 0,
                                               'implicit_hydrogens': 2, 'stereo': None})]
    with raises(IncorrectSmiles):
        smiles_tokenize('[@N]')


def test_smarts_tokenize_atom():
    # Per standard SMARTS semantics, uppercase aromaticity-capable symbols
    # ([C], [N], …) carry an implicit aliphatic constraint, encoded as
    # hybridization ∈ {1, 2, 3} so the matcher rejects aromatic atoms.
    # The marker ``_default_aliphatic`` is consumed by the parser to relax
    # the constraint when the atom is connected by an aromatic-only bond
    # (legacy ``[C]:[C]`` chython idiom).
    _aliph = {'hybridization': [1, 2, 3], '_default_aliphatic': True}
    assert smarts_tokenize('[C]') == [(0, {'element': 'C', **_aliph})]
    assert smarts_tokenize('[C,N]') == [(0, {'element': ['C', 'N'], **_aliph})]
    assert smarts_tokenize('[C+]') == [(0, {'charge': 1, 'element': 'C', **_aliph})]
    assert smarts_tokenize('[#1]') == [(0, {'element': 1})]
    assert smarts_tokenize('[C;h1;@]') == [(0, {'element': 'C', 'implicit_hydrogens': [1], 'stereo': True, **_aliph})]
    assert smarts_tokenize('[O;z1,z2;y1]') == [(0, {'element': 'O', 'heteroatoms': [1], 'hybridization': [1, 2]})]
    assert smarts_tokenize('[Se;a;D1,D2;r4,r7:3]') == [(8, {'parsed_mapping': 3, 'element': 'Se', 'hybridization': 4, 'neighbors': [1, 2], 'ring_sizes': [4, 7]})]
    assert smarts_tokenize('[Cl;M]') == [(0, {'element': 'Cl', 'masked': True})]
    # A is any aliphatic atom per Daylight; * is any atom
    assert smarts_tokenize('[A:1]') == [(0, {'parsed_mapping': 1, 'element': 'A',
                                             'hybridization': [1, 2, 3]})]
    assert smarts_tokenize('[*:1]') == [(0, {'parsed_mapping': 1, 'element': 'A'})]
    assert smarts_tokenize('[M]') == [(0, {'element': 'M'})]


def test_smarts_tokenize_bonds():
    _aliph_C = {'element': 'C', 'hybridization': [1, 2, 3], '_default_aliphatic': True}
    assert smarts_tokenize('[C][C]') == [(0, _aliph_C), (0, _aliph_C)]
    assert smarts_tokenize('[C]-[C]') == [(0, _aliph_C), (1, 1), (0, _aliph_C)]
    assert smarts_tokenize('[C]~[C]') == [(0, _aliph_C), (1, 8), (0, _aliph_C)]
    assert smarts_tokenize('[C]!:[C]') == [(0, _aliph_C), (10, [1, 2, 3]), (0, _aliph_C)]
    assert smarts_tokenize('[C]-,=[C]') == [(0, _aliph_C), (10, [1, 2]), (0, _aliph_C)]
    assert smarts_tokenize('[C]-;@[C]') == [(0, _aliph_C), (12, QueryBond(1, True)), (0, _aliph_C)]
    assert smarts_tokenize('[C]!-;!@[C]') == [(0, _aliph_C), (12, QueryBond((2, 3, 4), False)),
                                              (0, _aliph_C)]
    assert smarts_tokenize('[C]-,=;!@[C]') == [(0, _aliph_C), (12, QueryBond((1, 2), False)),
                                               (0, _aliph_C)]
    # ring marks without an explicit order, and orders glued to the mark
    assert smarts_tokenize('[C]@[C]') == [(0, _aliph_C), (12, QueryBond(8, True)), (0, _aliph_C)]
    assert smarts_tokenize('[C]!@[C]') == [(0, _aliph_C), (12, QueryBond(8, False)), (0, _aliph_C)]
    assert smarts_tokenize('[C]-!@[C]') == [(0, _aliph_C), (12, QueryBond(1, False)), (0, _aliph_C)]
    assert smarts_tokenize('[C]-&!@[C]') == [(0, _aliph_C), (12, QueryBond(1, False)), (0, _aliph_C)]
    assert smarts_tokenize('[C]-@[C]') == [(0, _aliph_C), (12, QueryBond(1, True)), (0, _aliph_C)]
    # more than two alternatives
    assert smarts_tokenize('[C]-,=,#[C]') == [(0, _aliph_C), (10, [1, 2, 3]), (0, _aliph_C)]
    assert smarts_tokenize('[C]-,=,#;!@[C]') == [(0, _aliph_C), (12, QueryBond((1, 2, 3), False)), (0, _aliph_C)]


_aliph = {'hybridization': [1, 2, 3], '_default_aliphatic': True}


@mark.parametrize('pattern,expected', [
    # juxtaposition and <&> are the same AND as <;>
    ('[CX3]', {'element': 'C', 'total_connectivity': [3], **_aliph}),
    ('[C&X3]', {'element': 'C', 'total_connectivity': [3], **_aliph}),
    ('[X3C]', {'element': 'C', 'total_connectivity': [3], **_aliph}),
    ('[OX2H1]', {'element': 'O', 'total_connectivity': [2], 'implicit_hydrogens': [1], **_aliph}),
    ('[nH]', {'element': 'N', 'implicit_hydrogens': [1], 'hybridization': 4}),
    ('[C@H]', {'element': 'C', 'implicit_hydrogens': [1], 'stereo': True, **_aliph}),
    ('[!#6!#1]', {'element': 'A', 'excluded_elements': (6, 1)}),
    # brackets without an element
    ('[r5]', {'element': 'A', 'ring_sizes': [5]}),
    ('[D2]', {'element': 'A', 'neighbors': [2]}),
    ('[+]', {'element': 'A', 'charge': 1}),
    ('[N+0]', {'element': 'N', 'charge': 0, **_aliph}),
    ('[X3,X4]', {'element': 'A', 'total_connectivity': [3, 4]}),
    # bare primitives take their Daylight default count, bare <r> means any ring
    ('[r]', {'element': 'A', 'rings_count': list(range(1, 15))}),
    # H is the hydrogen atom only at the head of the bracket
    ('[H]', {'element': 'H'}),
    ('[H;D2]', {'element': 'H', 'neighbors': [2]}),
    ('[CH]', {'element': 'C', 'implicit_hydrogens': [1], **_aliph}),
    # two-letter elements win over D/H/R/X primitives
    ('[Rh]', {'element': 'Rh'}),
])
def test_smarts_bracket_primitives(pattern, expected):
    assert smarts_tokenize(pattern) == [(8 if expected.get('hybridization') == 4 else 0, expected)]


@mark.parametrize('pattern,expected', [
    # a negated count primitive is the complement over its domain
    ('[!H0]', {'element': 'A', 'implicit_hydrogens': [1, 2, 3, 4]}),
    ('[C!D1]', {'element': 'C', 'neighbors': [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], **_aliph}),
    ('[!z1]', {'element': 'A', 'hybridization': [2, 3, 4]}),
    # primitives of the same kind are ANDed, not overwritten
    ('[C;!H0;!H1]', {'element': 'C', 'implicit_hydrogens': [2, 3, 4], **_aliph}),
    # ring membership: r0 and !R mean no ring, !r0 and R mean any ring
    ('[r0]', {'element': 'A', 'ring_sizes': 0}),
    ('[!r]', {'element': 'A', 'ring_sizes': 0}),
    ('[!r0]', {'element': 'A', 'rings_count': list(range(1, 15))}),
    ('[!R0]', {'element': 'A', 'rings_count': list(range(1, 15))}),
    # valence
    ('[Nv3]', {'element': 'N', 'valence': [3], **_aliph}),
    ('[!v4]', {'element': 'A', 'valence': [0, 1, 2, 3, 5, 6, 7, 8]}),
])
def test_smarts_negated_primitives(pattern, expected):
    assert smarts_tokenize(pattern) == [(0, expected)]


def test_smarts_or_across_primitives():
    # needs an expression tree chython's atom attributes cannot express
    with raises(IncorrectSmarts):
        smarts_tokenize('[C,X3]')


def test_smarts_contradictory_primitives():
    with raises(IncorrectSmarts):
        smarts_tokenize('[C;H1;H2]')
