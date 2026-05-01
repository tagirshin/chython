# -*- coding: utf-8 -*-
#
#  Copyright 2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from chython import smiles, smarts, Reactor
from pytest import mark


data = [
    (('[B;D3;x2;z1:4]([O:5])([O:6])-[C;@@;h1:3]1([O;M][C;M]1)', '[Cl,Br,I;D1:1]-[C;a:2]'), ('[A;@:3]-[A:2]',),
     ('CC1O[C@@H]1B(O)O', 'Brc1ccccc1'), ('CC1O[C@H]1c1ccccc1',)),  # inverse stereo check
    (('[B;D3;x2;z1:4]([O:5])([O:6])-[C;@@;h1:3]1([O;M][C;M]1)', '[Cl,Br,I;D1:1]-[C;a:2]'), ('[A;@@:3]-[A:2]',),
     ('CC1O[C@@H]1B(O)O', 'Brc1ccccc1'), ('CC1O[C@@H]1c1ccccc1',)),  # keep stereo on RC
    (('[B;D3;x2;z1:4]([O:5])([O:6])-[C;@@;h1:3]1([O;M][C;M]1)', '[Cl,Br,I;D1:1]-[C;a:2]'), ('[A:3]-[A:2]',),
     ('CC1O[C@@H]1B(O)O', 'Brc1ccccc1'), ('CC1OC1c1ccccc1',)),  # drop stereo on RC
]


@mark.parametrize('patterns, products, source, result', data)
def test_transformer(patterns, products, source, result):
    for q, m in zip(patterns, source):
        assert smarts(q) <= smiles(m)

    reactor = Reactor([smarts(x) for x in patterns], [smarts(x) for x in products])
    out = {format(smiles(x), 'h') for x in result}
    assert {format(x, 'h') for x in next(reactor(*(smiles(x) for x in source))).products} == out


def test_reactor_does_not_match_aliphatic_carbon_to_aromatic_carbon():
    rule_rxn = smarts('[C;D3:1]-[O;D1:2]>>[C;D3:1]-[O;D2:2]-[C;D1:3]')
    reactor = Reactor(patterns=tuple(rule_rxn.reactants), products=tuple(rule_rxn.products), delete_atoms=False)
    reaction = smiles(
        '[cH:18]1[cH:17][c:16]([cH:21][cH:20][cH:19]1)[CH:15]2[c:7]3[cH:8][c:9]([c:11]([c:13]([c:6]3[C:4]([CH:3]2[CH2:2][CH3:1])=[O:5])[Cl:14])[Cl:12])[O:10][CH3:22]'
        '>>'
        '[cH:18]1[cH:17][c:16]([cH:21][cH:20][cH:19]1)[CH:15]2[c:7]3[cH:8][c:9]([c:11]([Cl:12])[c:13]([Cl:14])[c:6]3[C:4](=[O:5])[CH:3]2[CH2:2][CH3:1])[OH:10]'
    )

    assert list(rule_rxn.reactants[0].get_mapping(reaction.products[0])) == []
    assert list(reactor(reaction.products[0])) == []
