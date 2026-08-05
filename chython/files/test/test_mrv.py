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
from io import BytesIO, StringIO

from pytest import mark

from chython import smiles
from chython.files import MRVRead, MRVWrite


MOLECULES = [
    'CCO',
    'c1ccccc1',
    'CC(=O)[O-].[Na+]',
    'CC(=O)Nc1ccc(O)cc1',
    # no tetrahedral centres here: MRV writes wedges against all-zero coordinates,
    # so stereo does not survive the round-trip for coordinate-free molecules
]

REACTIONS = [
    'CC(=O)O.CCO>>CC(=O)OCC.O',
    'Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1',
]


def _roundtrip(obj):
    out = StringIO()
    with MRVWrite(out) as w:
        w.write(obj)
    with MRVRead(BytesIO(out.getvalue().encode())) as r:
        return next(iter(r))


@mark.parametrize('smi', MOLECULES)
def test_molecule_roundtrip(smi):
    m = smiles(smi)
    m.canonicalize()
    back = _roundtrip(m)
    back.canonicalize()
    assert back == m


@mark.parametrize('smi', REACTIONS)
def test_reaction_roundtrip(smi):
    r = smiles(smi)
    r.canonicalize()
    back = _roundtrip(r)
    back.canonicalize()
    assert str(back) == str(r)


def test_metadata_roundtrip():
    m = smiles('CCO')
    m.canonicalize()
    m.meta['source'] = 'chython'
    m.meta['n'] = '42'
    back = _roundtrip(m)
    assert back.meta == {'source': 'chython', 'n': '42'}


def test_multiple_records():
    out = StringIO()
    mols = []
    with MRVWrite(out) as w:
        for smi in MOLECULES:
            m = smiles(smi)
            m.canonicalize()
            mols.append(m)
            w.write(m)
    with MRVRead(BytesIO(out.getvalue().encode())) as r:
        back = list(r)
    assert len(back) == len(mols)
    for b, m in zip(back, mols):
        b.canonicalize()
        assert b == m
