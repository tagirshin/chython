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
import pytest
from chython import smiles
from chython.exceptions import IncorrectSmiles
from chython.files.daylight.parser import parser
from chython.files.daylight.tokenize import smiles_tokenize


def test_daylight_smiles_basic():
    # Test basic SMILES
    result = parser(list(smiles_tokenize('CC')), True)
    assert len(result['atoms']) == 2
    assert len(result['bonds']) == 1
    
    result = parser(list(smiles_tokenize('O')), True)
    assert len(result['atoms']) == 1
    assert not result['bonds']


def test_daylight_smiles_empty():
    # Test empty SMILES
    tokens = list(smiles_tokenize(''))
    assert len(tokens) == 0  # empty string should produce empty token list
    
    # Empty token list should raise error when parsing
    with pytest.raises(IndexError):
        parser(tokens, True)


def test_daylight_smiles_invalid():
    # Test invalid SMILES
    with pytest.raises(IncorrectSmiles):
        parser(list(smiles_tokenize('C1CC')), True)  # unclosed cycle


def test_daylight_smiles_complex():
    # Test complex SMILES
    result = parser(list(smiles_tokenize('C1=CC=CC=C1')), True)
    assert len(result['atoms']) == 6
    assert len(result['bonds']) == 6
    
    result = parser(list(smiles_tokenize('C(=O)O')), True)
    assert len(result['atoms']) == 3
    assert len(result['bonds']) == 2


def test_daylight_smiles_charged():
    # Test charged species
    result = parser(list(smiles_tokenize('[NH4+]')), True)
    assert len(result['atoms']) == 1

    result = parser(list(smiles_tokenize('[OH-]')), True)
    assert len(result['atoms']) == 1


def test_daylight_reaction_normalizes_chromium_trioxide_reagent():
    reaction = smiles('C>[O-2].[O-2].[O-2].[Cr+6]>C')

    assert [str(x) for x in reaction.reagents] == ['O=[Cr](=O)=O']
    assert reaction.meta['chython_parsing_log'] == [
        'normalized oxide fragments: [O-2].[O-2].[O-2].[Cr+6] -> O=[Cr](=O)=O'
    ]


def test_daylight_reaction_normalizes_mapped_chromium_trioxide_bundle():
    reaction = smiles('[O-2:1].[Cr+6:4].[O-2:2].[O-2:3]>>[O:1]=[Cr:4](=[O:2])=[O:3]')

    assert format(reaction.reactants[0], 'm') == '[O:1]=[Cr:4](=[O:2])=[O:3]'


@pytest.mark.parametrize(('source', 'reagent'), (
    ('[O-2].[O-2].[O-2].[U+6]', 'O=[U](=O)=O'),
    ('[O-2].[O-2].[O-2].[O-2].[O-2].[V+5].[V+5]', 'O=[V](=O)O[V](=O)=O'),
    ('[O-2].[O-2].[O-2].[O-2].[O-2].[Nb+5].[Nb+5]', 'O=[Nb](=O)O[Nb](=O)=O'),
    ('[O-2].[O-2].[O-2].[O-2].[O-2].[Ta+5].[Ta+5]', 'O=[Ta](=O)O[Ta](=O)=O'),
))
def test_daylight_reaction_normalizes_high_valent_oxide_bundles(source, reagent):
    reaction = smiles(f'C>{source}>C')

    assert [str(x) for x in reaction.reagents] == [reagent]


def test_daylight_reaction_drops_standalone_phosphorus_pentahydride():
    reaction = smiles('C.[PH5]>O>C')

    assert [str(x) for x in reaction.reactants] == ['C']
    assert reaction.meta['chython_parsing_log'] == ['ignored unsupported standalone fragment: [PH5]']


def test_daylight_reaction_normalizes_phosphorus_oxychloride_fragments():
    reaction = smiles('C.O(Cl)Cl.[P+5]>O>C')

    assert [str(x) for x in reaction.reactants] == ['C', 'ClP(Cl)(Cl)=O']
    assert reaction.meta['chython_parsing_log'] == [
        'normalized phosphorus(V) oxyhalide fragments: O(Cl)Cl.[P+5] -> O=P(Cl)(Cl)Cl'
    ]


def test_daylight_reaction_normalizes_mapped_phosphorus_oxychloride_fragments():
    reaction = smiles('C.O(Cl)[Cl:2].[P+5]>O>C')

    assert format(reaction.reactants[1], 'm') == '[Cl:2][P:5]([Cl:6])([Cl:7])=[O:4]'
    assert reaction.meta['chython_parsing_log'] == [
        'normalized phosphorus(V) oxyhalide fragments: O(Cl)[Cl:2].[P+5] -> O=P(Cl)([Cl:2])Cl'
    ]


def test_daylight_reaction_normalizes_phosphorus_oxybromide_fragments():
    reaction = smiles('C.O(Br)[Br:2].[P+5]>O>C')

    assert [str(x) for x in reaction.reactants] == ['C', 'BrP(Br)(Br)=O']
    assert '[Br:2]' in format(reaction.reactants[1], 'm')


def test_daylight_smiles_aromatic():
    # Test aromatic SMILES
    result = parser(list(smiles_tokenize('c1ccccc1')), True)
    assert len(result['atoms']) == 6
    assert len(result['bonds']) == 6
    
    result = parser(list(smiles_tokenize('n1cccc1')), True)
    assert len(result['atoms']) == 5
    assert len(result['bonds']) == 5


def test_daylight_smiles_isotopes():
    # Test isotope labels
    result = parser(list(smiles_tokenize('[2H]O')), True)
    assert len(result['atoms']) == 2
    assert result['atoms'][0]['isotope'] == 2
    
    result = parser(list(smiles_tokenize('[13C]')), True)
    assert len(result['atoms']) == 1
    assert result['atoms'][0]['isotope'] == 13


def test_daylight_smiles_stereo():
    # Test stereochemistry
    result = parser(list(smiles_tokenize('C/C=C/C')), True)
    assert result['stereo_bonds']

    result = parser(list(smiles_tokenize('C/C=C\\C')), True)
    assert result['stereo_bonds']
