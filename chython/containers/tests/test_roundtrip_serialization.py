import pickle

from chython import smiles
from chython.containers import MoleculeContainer, ReactionContainer
from chython.containers.cgr_query import QueryCGRContainer
from chython.containers.cgr import CGRContainer


def assert_roundtrip_pickle(obj):
    data = pickle.dumps(obj)
    restored = pickle.loads(data)
    return restored


def test_molecule_smiles_pickle_roundtrip():
    mol: MoleculeContainer = smiles('CCO')  # simple ethanol
    assert len(mol) == 3

    mol_smiles = str(mol)
    mol_from_smiles = smiles(mol_smiles)
    assert isinstance(mol_from_smiles, MoleculeContainer)
    assert len(mol_from_smiles) == len(mol)
    assert str(mol_from_smiles) == mol_smiles

    mol_restored = assert_roundtrip_pickle(mol)
    assert isinstance(mol_restored, MoleculeContainer)
    assert len(mol_restored) == len(mol)
    assert str(mol_restored) == str(mol)


def test_reaction_smiles_pickle_roundtrip():
    rxn: ReactionContainer = smiles('[CH3:1][CH2:2][Br:3]>>[CH3:1][CH2:2][OH:3]')
    assert rxn

    rxn_smiles = format(rxn, 'm')
    rxn_from_smiles = smiles(rxn_smiles)
    assert isinstance(rxn_from_smiles, ReactionContainer)
    assert format(rxn_from_smiles, 'm') == rxn_smiles

    rxn_restored = assert_roundtrip_pickle(rxn)
    assert isinstance(rxn_restored, ReactionContainer)
    assert format(rxn_restored, 'm') == format(rxn, 'm')


def test_cgr_and_query_cgr_roundtrip_and_rule():
    rxn: ReactionContainer = smiles('[CH3:1][CH2:2][OH:3]>>[CH2:1]=[CH2:2].[OH2:3]')
    cgr: CGRContainer = ~rxn
    assert len(cgr.center_atoms) >= 1

    cgr_smiles = str(cgr)
    cgr_restored = assert_roundtrip_pickle(cgr)
    assert isinstance(cgr_restored, CGRContainer)
    assert str(cgr_restored) == cgr_smiles

    # Query CGR from substructure
    qcgr: QueryCGRContainer = cgr.substructure(cgr.center_atoms, as_query=True)
    qcgr_smiles = str(qcgr)
    qcgr_restored = assert_roundtrip_pickle(qcgr)
    assert isinstance(qcgr_restored, QueryCGRContainer)
    assert str(qcgr_restored) == qcgr_smiles

    # Reaction rule: center atoms + first environment
    rule = cgr.augmented_substructure(cgr.center_atoms, deep=1)
    assert isinstance(rule, CGRContainer)
    assert len(rule) >= len(cgr.center_atoms)
    rule_smiles = str(rule)
    rule_restored = assert_roundtrip_pickle(rule)
    assert isinstance(rule_restored, CGRContainer)
    assert str(rule_restored) == rule_smiles
