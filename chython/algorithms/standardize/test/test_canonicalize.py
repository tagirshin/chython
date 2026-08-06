# -*- coding: utf-8 -*-
from pytest import mark

from chython import smiles


def test_canonicalize_nitro():
    mol = smiles('c1ccccc1N(=O)=O')
    mol.canonicalize()
    s = str(mol)
    assert '[N+]' in s and '[O-]' in s


def test_canonicalize_benzene():
    mol = smiles('C1=CC=CC=C1')
    mol.canonicalize()
    assert str(mol) == 'c1ccccc1'


def test_kekule_thiele():
    mol = smiles('c1ccccc1')
    mol.kekule()
    s = str(mol)
    assert 'c' not in s
    assert '=' in s
    mol.thiele()
    assert str(mol) == 'c1ccccc1'



def test_neutralize():
    mol = smiles('[NH3+]CC(=O)[O-]')
    mol.neutralize()
    s = str(mol)
    assert '+' not in s or '-' not in s


def test_check_valence_clean():
    mol = smiles('CCO')
    errors = mol.check_valence()
    assert errors == []


def test_canonicalize_idempotent():
    mol = smiles('OCC')
    mol.canonicalize()
    s1 = str(mol)
    mol.canonicalize()
    s2 = str(mol)
    assert s1 == s2


# canonicalize must be a fixed point. standardize_tautomers invalidates the cached
# morgan ordering but used to leave int_adjacency (which it is derived from) stale,
# so amidines and guanidines canonicalized to a different string on the second call.
@mark.parametrize('smi', ['CC(=N)Nc1ccccc1', 'NC(=N)Nc1ccccc1', 'CC(=N)NC', 'NC(=N)N',
                          'CC(=N)Nc1ccc(C)cc1', 'CCSC(N)=N', 'N#CNC(N)=N'])
def test_canonicalize_is_fixed_point(smi):
    m = smiles(smi)
    m.canonicalize()
    once = str(m)
    m.canonicalize()
    assert str(m) == once, 'canonicalize not idempotent'
    reparsed = smiles(once)
    reparsed.canonicalize()
    assert str(reparsed) == once, 'canonical form not stable across a parse round-trip'
