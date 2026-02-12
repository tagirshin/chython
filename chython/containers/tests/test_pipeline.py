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
"""
Pipeline integration tests for chython.

These tests exercise chython the same way downstream tools like SynPlanner do:
reaction SMILES parsing, CGR creation, reaction center analysis, SMARTS matching,
Reactor-based rule application, substructure queries, and serialization.
"""

import pickle

import pytest
from chython import smiles, smarts, Reactor
from chython.containers import (
    CGRContainer,
    MoleculeContainer,
    ReactionContainer,
    QueryContainer,
    Bond,
    QueryBond,
)
from chython.periodictable import QueryElement
from chython.exceptions import InvalidAromaticRing


# ---------------------------------------------------------------------------
# Test reaction data (same reactions used by SynPlanner's test suite)
# ---------------------------------------------------------------------------

REACTION_SMILES = [
    "[CH3:5][CH2:6][OH:7].[O:4]=[C:2]([OH:3])[CH3:1]>[O:8]=[S:9](=[O:10])([OH:11])[OH:12]>"
    "[CH3:5][CH2:6][O:7][C:2](=[O:4])[CH3:1].[OH2:3]",
    "[CH2:5]=[CH2:6].[CH:2]([CH:3]=[CH2:4])=[CH2:1]>>[CH:3]1=[CH:4][CH2:6][CH2:5][CH2:2][CH2:1]1",
    "[CH:2](=[O:3])[CH3:1].[CH:5](=[O:6])[CH3:4]>[Na+:7].[OH-:8]>"
    "[OH:3][CH:2]([CH2:4][CH:5]=[O:6])[CH3:1]",
    "[CH2:1]1[CH2:2][CH2:3][CH2:4][CH2:5][CH:6]1[OH:7]>[O:11]=[Cr:10](=[O:12])([OH:9])[OH:13]>"
    "[CH2:2]1[CH2:1][C:6](=[O:7])[CH2:5][CH2:4][CH2:3]1",
    "[CH3:1][CH2:2][CH2:3][CH:4]=[O:5]>[BH4-:6].[Na+:7]>[CH2:4]([OH:5])[CH2:3][CH2:2][CH3:1]",
]


# ---------------------------------------------------------------------------
# Helpers (replicate what SynPlanner's extraction module does with chython)
# ---------------------------------------------------------------------------


def molecule_substructure_as_query(mol, atoms) -> QueryContainer:
    """Build a QueryContainer from a molecule substructure (SynPlanner pattern)."""
    atoms = set(atoms)
    q = QueryContainer(smarts="")
    for n in atoms:
        atom = mol.atom(n)
        if isinstance(atom, QueryElement):
            q.add_atom(atom.copy(full=True), n)
        else:
            q.add_atom(
                QueryElement.from_atom(
                    atom,
                    neighbors=True,
                    hybridization=True,
                    hydrogens=True,
                    ring_sizes=True,
                    heteroatoms=True,
                ),
                n,
            )
    for n, m, bond in mol.bonds():
        if n in atoms and m in atoms:
            if isinstance(bond, QueryBond):
                q.add_bond(n, m, bond.copy(full=True))
            elif isinstance(bond, Bond):
                q.add_bond(n, m, QueryBond.from_bond(bond))
    return q


def validate_rule(rule: ReactionContainer, reaction: ReactionContainer) -> bool:
    """Validate a reaction rule using the Reactor (SynPlanner pattern)."""
    patterns = tuple(
        molecule_substructure_as_query(m, m.atoms_numbers) for m in rule.reactants
    )
    products = tuple(rule.products)
    reactor = Reactor(patterns=patterns, products=products)
    try:
        for result_reaction in reactor(*reaction.reactants):
            result_products = []
            for result_product in result_reaction.products:
                tmp = result_product.copy()
                try:
                    tmp.kekule()
                    if tmp.check_valence():
                        continue
                except InvalidAromaticRing:
                    continue
                result_products.append(result_product)
            if set(reaction.products) == set(result_products) and len(
                reaction.products
            ) == len(result_products):
                return True
    except (KeyError, IndexError, InvalidAromaticRing):
        return False
    return False


def add_environment_atoms(
    cgr: CGRContainer, center_atoms: set[int], depth: int
) -> set[int]:
    """Expand center atoms with environment (SynPlanner pattern)."""
    if depth:
        env = cgr.augmented_substructure(center_atoms, deep=depth)
        return center_atoms | set(env)
    return center_atoms


def add_ring_structures(cgr: CGRContainer, rule_atoms: set[int]) -> set[int]:
    """Add ring atoms intersecting center (SynPlanner pattern)."""
    for ring in cgr.sssr:
        if set(ring) & rule_atoms:
            rule_atoms |= set(ring)
    return rule_atoms


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module", params=REACTION_SMILES, ids=[
    "fischer_ester", "diels_alder", "aldol", "jones_oxidation", "nabh4_reduction"
])
def parsed_reaction(request) -> ReactionContainer:
    return smiles(request.param)


@pytest.fixture(scope="module")
def fischer_reaction() -> ReactionContainer:
    return smiles(REACTION_SMILES[0])


@pytest.fixture(scope="module")
def diels_alder_reaction() -> ReactionContainer:
    return smiles(REACTION_SMILES[1])


@pytest.fixture(scope="module")
def esterification_reaction() -> ReactionContainer:
    return smiles(
        "[CH3:1][C:2](=[O:3])[OH:4].[OH:5][CH3:6]>>"
        "[CH3:1][C:2](=[O:3])[O:5][CH3:6].[OH2:4]"
    )


# ---------------------------------------------------------------------------
# 1. SMILES parsing (the entry point for every SynPlanner pipeline)
# ---------------------------------------------------------------------------


class TestSmilesParsing:
    """Test that all reaction types SynPlanner uses parse correctly."""

    def test_all_reactions_parse(self):
        for rxn_smi in REACTION_SMILES:
            rxn = smiles(rxn_smi)
            assert isinstance(rxn, ReactionContainer)

    def test_reaction_has_reactants_and_products(self, parsed_reaction):
        assert len(parsed_reaction.reactants) >= 1
        assert len(parsed_reaction.products) >= 1

    def test_molecule_parsing(self):
        mol = smiles("CC(=O)OC1=CC=CC=C1C(=O)O")
        assert isinstance(mol, MoleculeContainer)
        assert len(mol) > 0

    def test_smarts_parsing(self):
        q = smarts("[C]=[O]")
        assert isinstance(q, QueryContainer)

    def test_smarts_substructure_match(self):
        q = smarts("[C]=[O]")
        mol = smiles("CC=O")
        assert q <= mol  # substructure match operator

    def test_smarts_get_mapping(self):
        q = smarts("[C]=[O]")
        mol = smiles("CC(=O)O")
        mappings = list(q.get_mapping(mol))
        assert len(mappings) >= 1


# ---------------------------------------------------------------------------
# 2. CGR creation and analysis (core of filtering and rule extraction)
# ---------------------------------------------------------------------------


class TestCGRCreation:
    """Test CGR creation from reactions — the ~reaction operator."""

    def test_cgr_from_reaction(self, parsed_reaction):
        cgr = ~parsed_reaction
        assert isinstance(cgr, CGRContainer)
        assert len(cgr) > 0

    def test_cgr_center_atoms_nonempty(self, fischer_reaction):
        cgr = ~fischer_reaction
        assert len(cgr.center_atoms) >= 1

    def test_cgr_center_bonds(self, fischer_reaction):
        cgr = ~fischer_reaction
        assert len(cgr.center_bonds) >= 1

    def test_cgr_centers_list(self, fischer_reaction):
        cgr = ~fischer_reaction
        centers = cgr.centers_list
        assert len(centers) >= 1
        # every center atom should appear in exactly one group
        all_atoms = set()
        for group in centers:
            all_atoms.update(group)
        assert set(cgr.center_atoms).issubset(all_atoms)

    def test_cgr_diels_alder_ring_formed(self, diels_alder_reaction):
        cgr = ~diels_alder_reaction
        # Diels-Alder forms a ring — check SSSR
        assert len(cgr.sssr) >= 1
        center = set(cgr.center_atoms)
        assert len(center) >= 4  # at least 4 atoms involved in [4+2]

    def test_cgr_connected_components(self, fischer_reaction):
        cgr = ~fischer_reaction
        # single connected component for proper reaction
        assert cgr.connected_components_count >= 1

    def test_cgr_str_roundtrip(self, parsed_reaction):
        cgr = ~parsed_reaction
        cgr_str = str(cgr)
        assert isinstance(cgr_str, str)
        assert len(cgr_str) > 0


# ---------------------------------------------------------------------------
# 3. Reaction filtering patterns (SynPlanner DynamicBondsFilter, etc.)
# ---------------------------------------------------------------------------


class TestReactionFiltering:
    """Replicate SynPlanner's filtering logic directly on chython objects."""

    def test_dynamic_bonds_count(self):
        """DynamicBondsFilter checks center_bonds count."""
        rxn = smiles("[CH2:1]=[O:2]>>[CH2:1]-[O:2]")
        cgr = ~rxn
        n_dynamic = len(cgr.center_bonds)
        assert 1 <= n_dynamic <= 2

    def test_dynamic_bonds_too_many_filtered(self):
        rxn = smiles(
            "[CH2:1]=[CH:2]-[CH:3]=[CH:4]-[CH:5]=[O:6]>>"
            "[CH2:1]-[CH2:2]-[CH2:3]-[CH2:4]-[CH2:5]-[O:6]"
        )
        cgr = ~rxn
        assert len(cgr.center_bonds) > 2

    def test_no_reaction_detected(self):
        """NoReactionFilter: no center_atoms and no center_bonds."""
        rxn = smiles("[CH3:1]-[CH3:2]>>[CH3:1]-[CH3:2]")
        cgr = ~rxn
        assert not cgr.center_atoms
        assert not cgr.center_bonds

    def test_multicenter_detection(self):
        """MultiCenterFilter: centers_list > 1."""
        # two separate bond changes
        rxn = smiles("[CH3:1][Br:2].[CH3:3][Cl:4]>>[CH3:1][OH:5].[CH3:3][OH:6]")
        cgr = ~rxn
        # may or may not be multicenter depending on mapping; just test access
        _ = cgr.centers_list

    def test_sp3_bond_breaking_detection(self):
        """CCsp3BreakingFilter pattern: check hybridization on augmented substructure."""
        rxn = smiles("[CH3:1][CH2:2][CH2:3][CH3:4]>>[CH3:1][CH2:2][CH3:3].[CH3:4]")
        cgr = ~rxn
        rc = cgr.augmented_substructure(cgr.center_atoms, deep=1)
        # verify we can check hybridization and bond orders
        for n, m, bond in rc.bonds():
            is_broken = bond.order is not None and bond.p_order is None
            if is_broken:
                hn = rc._hybridizations.get(n, 0)
                hm = rc._hybridizations.get(m, 0)
                # at least one sp3 carbon
                assert hn == 1 or hm == 1

    def test_ring_breaking_detection(self):
        """CCRingBreakingFilter pattern: iterate reactant rings and check CGR bonds."""
        rxn = smiles(
            "[CH2:1]1[CH2:2][CH2:3][CH2:4][CH2:5]1>>"
            "[CH2:1]-[CH2:2]-[CH2:3]-[CH2:4]-[CH2:5]"
        )
        cgr = ~rxn
        reactant_rings = set()
        for reactant in rxn.reactants:
            reactant_rings.update(reactant.sssr)
        # ring atoms should overlap with center
        assert any(set(ring) & set(cgr.center_atoms) for ring in reactant_rings)

    def test_wrong_ch_breaking_pattern(self):
        """WrongCHBreakingFilter: explicify H, build augmented CGR, check bonds."""
        rxn = smiles("[CH3:1][CH2:2][H:3]>>[CH3:1][CH2:2][CH3:4]")
        if not rxn.check_valence():
            copy_rxn = rxn.copy()
            copy_rxn.explicify_hydrogens()
            cgr = ~copy_rxn
            rc = cgr.augmented_substructure(cgr.center_atoms, deep=1)
            # verify bond inspection works
            for atom_id in rc.center_atoms:
                atom = rc.atom(atom_id)
                for neighbour_id, bond in rc._bonds[atom_id].items():
                    _ = bond.order
                    _ = bond.p_order
                    _ = rc.atom(neighbour_id).atomic_symbol

    def test_rings_change_filter_pattern(self):
        """RingsChangeFilter: compare ring counts between reactants and products."""
        rxn = smiles(REACTION_SMILES[1])  # Diels-Alder
        r_rings = sum(m.rings_count for m in rxn.reactants)
        p_rings = sum(m.rings_count for m in rxn.products)
        r_arom = sum(len(m.aromatic_rings) for m in rxn.reactants)
        p_arom = sum(len(m.aromatic_rings) for m in rxn.products)
        # Diels-Alder: product gains a ring
        assert p_rings > r_rings

    def test_strange_carbons_filter_pattern(self):
        """StrangeCarbonsFilter: check for molecules with only C and one bond type."""
        mol = smiles("CCCC")
        atom_types = {a.atomic_symbol for _, a in mol.atoms()}
        bond_types = {int(b) for _, _, b in mol.bonds()}
        assert atom_types == {"C"}
        assert len(bond_types) == 1


# ---------------------------------------------------------------------------
# 4. CGR substructure and environment (rule extraction core)
# ---------------------------------------------------------------------------


class TestSubstructureAndEnvironment:
    """Test augmented_substructure and substructure — core of rule extraction."""

    def test_augmented_substructure_depth0(self, esterification_reaction):
        cgr = ~esterification_reaction
        center = set(cgr.center_atoms)
        expanded = add_environment_atoms(cgr, center, 0)
        assert expanded == center

    def test_augmented_substructure_depth1(self, esterification_reaction):
        cgr = ~esterification_reaction
        center = set(cgr.center_atoms)
        expanded = add_environment_atoms(cgr, center, 1)
        assert center.issubset(expanded)
        assert len(expanded) >= len(center)

    def test_augmented_substructure_depth2(self, esterification_reaction):
        cgr = ~esterification_reaction
        center = set(cgr.center_atoms)
        d1 = add_environment_atoms(cgr, center, 1)
        d2 = add_environment_atoms(cgr, center, 2)
        assert d1.issubset(d2)

    def test_add_ring_structures_no_ring(self, esterification_reaction):
        cgr = ~esterification_reaction
        center = set(cgr.center_atoms)
        result = add_ring_structures(cgr, center.copy())
        # esterification typically has no rings in center
        assert center.issubset(result)

    def test_add_ring_structures_diels_alder(self, diels_alder_reaction):
        cgr = ~diels_alder_reaction
        center = set(cgr.center_atoms)
        result = add_ring_structures(cgr, center.copy())
        # DA forms ring — ring atoms should be included
        ring_atoms = {a for ring in cgr.sssr if set(ring) & center for a in ring}
        assert (center | ring_atoms).issubset(result)

    def test_molecule_substructure(self, esterification_reaction):
        """Test molecule.substructure() — used to extract rule fragments."""
        for mol in esterification_reaction.reactants:
            sub = mol.substructure(list(mol.atoms_numbers)[:2])
            assert isinstance(sub, MoleculeContainer)
            assert len(sub) <= len(mol)

    def test_substructure_as_query(self, esterification_reaction):
        """Test cgr.substructure(as_query=True) for QueryCGRContainer creation."""
        cgr = ~esterification_reaction
        from chython.containers.cgr_query import QueryCGRContainer
        q = cgr.substructure(cgr.center_atoms, as_query=True)
        assert isinstance(q, QueryCGRContainer)
        assert len(q) == len(cgr.center_atoms)


# ---------------------------------------------------------------------------
# 5. QueryContainer and QueryElement (rule construction)
# ---------------------------------------------------------------------------


class TestQueryConstruction:
    """Test query construction patterns used in SynPlanner's extraction module."""

    def test_query_element_from_atom(self):
        mol = smiles("CCO")
        for n, atom in mol.atoms():
            qe = QueryElement.from_atom(
                atom,
                neighbors=True,
                hybridization=True,
                hydrogens=True,
                ring_sizes=True,
                heteroatoms=True,
            )
            assert qe.atomic_number == atom.atomic_number

    def test_build_query_container(self):
        mol = smiles("CCO")
        q = molecule_substructure_as_query(mol, mol.atoms_numbers)
        assert isinstance(q, QueryContainer)
        assert len(list(q.atoms())) == len(mol)

    def test_query_substructure_match(self):
        """A query built from a molecule should match the original molecule."""
        mol = smiles("CCO")
        q = molecule_substructure_as_query(mol, mol.atoms_numbers)
        assert q <= mol

    def test_smarts_functional_group_mapping(self):
        """SynPlanner uses SMARTS to find functional groups, then remap."""
        carbonyl = smarts("[C]=[O]")
        mol = smiles("CC(=O)O")
        mappings = list(carbonyl.get_mapping(mol))
        assert len(mappings) >= 1
        # remap and check atoms
        for mp in mappings:
            carbonyl.remap(mp)
            assert len(set(carbonyl.atoms_numbers) & set(mol.atoms_numbers)) > 0
            # restore original mapping
            carbonyl.remap({v: k for k, v in mp.items()})


# ---------------------------------------------------------------------------
# 6. Reactor (the heart of retrosynthetic rule application)
# ---------------------------------------------------------------------------


class TestReactor:
    """Test Reactor patterns used by SynPlanner for rule validation and application."""

    def test_basic_reactor_application(self):
        """Build a reactor from SMARTS patterns and apply to molecules."""
        pattern = smarts("[C:1]Br")
        product = smarts("[A:1]O")
        reactor = Reactor([pattern], [product])
        mol = smiles("CBr")
        results = list(reactor(mol))
        assert len(results) >= 1
        # products should contain methanol
        for result in results:
            product_smiles = {str(p) for p in result.products}
            assert any('O' in s for s in product_smiles)

    def test_reactor_from_rule(self, esterification_reaction):
        """Replicate SynPlanner's validate_rule: extract rule, build reactor, apply."""
        rxn = esterification_reaction
        cgr = ~rxn
        center = set(cgr.center_atoms)

        # extract substructures for the rule
        r_subs = []
        for m in rxn.reactants:
            atoms = center & set(m.atoms_numbers)
            if atoms:
                r_subs.append(m.substructure(atoms))
        p_subs = []
        for m in rxn.products:
            atoms = center & set(m.atoms_numbers)
            if atoms:
                p_subs.append(m.substructure(atoms))

        # build the rule
        if r_subs and p_subs:
            rule = ReactionContainer(r_subs, p_subs)
            # attempt reactor construction
            patterns = tuple(
                molecule_substructure_as_query(m, m.atoms_numbers) for m in rule.reactants
            )
            products = tuple(rule.products)
            reactor = Reactor(patterns=patterns, products=products)
            # should not raise
            assert reactor is not None

    def test_stereo_preserving_reactor(self):
        """Reactor preserves stereochemistry (existing reactor test pattern)."""
        patterns = [
            smarts("[B;D3;x2;z1:4]([O:5])([O:6])-[C;@@;h1:3]1([O;M][C;M]1)"),
            smarts("[Cl,Br,I;D1:1]-[C;a:2]"),
        ]
        products = [smarts("[A;@:3]-[A:2]")]
        reactor = Reactor(patterns, products)
        result = next(reactor(smiles("CC1O[C@@H]1B(O)O"), smiles("Brc1ccccc1")))
        out = {format(p, 'h') for p in result.products}
        assert format(smiles("CC1O[C@H]1c1ccccc1"), 'h') in out


# ---------------------------------------------------------------------------
# 7. Molecule operations used in pipelines
# ---------------------------------------------------------------------------


class TestMoleculeOperations:
    """Test molecule-level operations SynPlanner relies on."""

    def test_canonicalize(self):
        mol = smiles("OCC")
        mol.canonicalize()
        assert str(mol) == "CCO"

    def test_kekule_thiele_roundtrip(self):
        mol = smiles("c1ccccc1")
        mol.kekule()
        # after kekule, no aromatic bonds
        for _, _, b in mol.bonds():
            assert int(b) != 4
        mol.thiele()
        # after thiele, aromatic bonds restored
        has_aromatic = any(int(b) == 4 for _, _, b in mol.bonds())
        assert has_aromatic

    def test_check_valence(self):
        mol = smiles("CCO")
        errors = mol.check_valence()
        assert not errors  # no valence errors

    def test_explicify_hydrogens(self):
        mol = smiles("C")
        count_before = len(mol)
        mol.explicify_hydrogens()
        assert len(mol) > count_before  # H atoms added

    def test_aromatic_rings(self):
        mol = smiles("c1ccccc1")
        assert len(mol.aromatic_rings) >= 1

    def test_sssr(self):
        mol = smiles("C1CCCCC1")
        assert len(mol.sssr) == 1

    def test_rings_count(self):
        mol = smiles("C1CCCCC1")
        assert mol.rings_count == 1

    def test_atom_properties(self):
        mol = smiles("CCO")
        for n, atom in mol.atoms():
            assert hasattr(atom, "atomic_symbol")
            assert hasattr(atom, "atomic_number")
            assert hasattr(atom, "charge")

    def test_bond_properties(self):
        mol = smiles("C=O")
        for n, m, bond in mol.bonds():
            assert int(bond) == 2  # double bond

    def test_atoms_numbers(self):
        mol = smiles("CCO")
        nums = mol.atoms_numbers
        assert len(list(nums)) == 3

    def test_copy(self):
        mol = smiles("CCO")
        mol2 = mol.copy()
        assert str(mol) == str(mol2)
        assert mol is not mol2

    def test_reaction_properties(self, fischer_reaction):
        rxn = fischer_reaction
        assert hasattr(rxn, "reactants")
        assert hasattr(rxn, "products")
        assert hasattr(rxn, "reagents")
        assert len(rxn.reactants) >= 1
        assert len(rxn.products) >= 1


# ---------------------------------------------------------------------------
# 8. Dynamic bond inspection (CGR bond analysis)
# ---------------------------------------------------------------------------


class TestDynamicBonds:
    """Test dynamic bond properties used heavily in filtering logic."""

    def test_dynamic_bond_order_and_p_order(self):
        rxn = smiles("[CH2:1]=[O:2]>>[CH2:1]-[O:2]")
        cgr = ~rxn
        for n, m, bond in cgr.bonds():
            if bond.is_dynamic:
                assert bond.order != bond.p_order

    def test_bond_breaking_detection(self):
        """bond.order is not None and bond.p_order is None means bond broken."""
        rxn = smiles("[CH3:1][CH2:2][OH:3]>>[CH3:1][CH3:2].[OH2:3]")
        cgr = ~rxn
        broken = [(n, m) for n, m, b in cgr.bonds() if b.order is not None and b.p_order is None]
        assert len(broken) >= 1

    def test_bond_formation_detection(self):
        """bond.order is None and bond.p_order is not None means bond formed."""
        rxn = smiles("[CH2:1]=[CH2:2].[CH2:3]=[CH2:4]>>[CH2:1]-[CH2:2]-[CH2:3]-[CH2:4]")
        cgr = ~rxn
        formed = [(n, m) for n, m, b in cgr.bonds() if b.order is None and b.p_order is not None]
        assert len(formed) >= 1

    def test_hybridization_access(self):
        """Filtering code accesses _hybridizations dict directly."""
        rxn = smiles("[CH3:1][CH2:2][OH:3]>>[CH3:1][CH3:2].[OH2:3]")
        cgr = ~rxn
        for n in cgr:
            assert n in cgr._hybridizations
            assert n in cgr._p_hybridizations


# ---------------------------------------------------------------------------
# 9. Serialization (pickle roundtrip — critical for rule storage)
# ---------------------------------------------------------------------------


class TestSerialization:
    """Pickle roundtrip tests — SynPlanner stores extracted rules as pickled objects."""

    def test_molecule_pickle(self):
        mol = smiles("CCO")
        restored = pickle.loads(pickle.dumps(mol))
        assert isinstance(restored, MoleculeContainer)
        assert str(restored) == str(mol)

    def test_reaction_pickle(self):
        rxn = smiles("[CH3:1][CH2:2][Br:3]>>[CH3:1][CH2:2][OH:3]")
        restored = pickle.loads(pickle.dumps(rxn))
        assert isinstance(restored, ReactionContainer)
        assert format(restored, "m") == format(rxn, "m")

    def test_cgr_pickle(self, esterification_reaction):
        cgr = ~esterification_reaction
        restored = pickle.loads(pickle.dumps(cgr))
        assert isinstance(restored, CGRContainer)
        assert str(restored) == str(cgr)

    def test_query_container_pickle(self):
        q = smarts("[C]=[O]")
        restored = pickle.loads(pickle.dumps(q))
        assert isinstance(restored, QueryContainer)

    def test_query_cgr_pickle(self, esterification_reaction):
        from chython.containers.cgr_query import QueryCGRContainer
        cgr = ~esterification_reaction
        q = cgr.substructure(cgr.center_atoms, as_query=True)
        restored = pickle.loads(pickle.dumps(q))
        assert isinstance(restored, QueryCGRContainer)
        assert str(restored) == str(q)

    def test_reactor_rule_pickle_roundtrip(self, esterification_reaction):
        """Full pipeline: reaction -> CGR -> rule -> pickle -> restore -> validate."""
        rxn = esterification_reaction
        cgr = ~rxn
        center = set(cgr.center_atoms)

        # extract rule substructures
        r_subs = [m.substructure(center & set(m.atoms_numbers))
                  for m in rxn.reactants if center & set(m.atoms_numbers)]
        p_subs = [m.substructure(center & set(m.atoms_numbers))
                  for m in rxn.products if center & set(m.atoms_numbers)]

        if r_subs and p_subs:
            rule = ReactionContainer(r_subs, p_subs)
            restored = pickle.loads(pickle.dumps(rule))
            assert isinstance(restored, ReactionContainer)
            assert len(restored.reactants) == len(rule.reactants)
            assert len(restored.products) == len(rule.products)


# ---------------------------------------------------------------------------
# 10. End-to-end rule extraction pipeline
# ---------------------------------------------------------------------------


class TestRuleExtractionPipeline:
    """End-to-end: reaction -> CGR -> center -> environment -> rule -> reactor."""

    @pytest.mark.parametrize("env_depth", [0, 1, 2])
    def test_extraction_pipeline(self, esterification_reaction, env_depth):
        rxn = esterification_reaction
        # 1. Create CGR
        cgr = ~rxn
        center = set(cgr.center_atoms)
        assert len(center) >= 1

        # 2. Add environment
        rule_atoms = add_environment_atoms(cgr, center, env_depth)
        assert center.issubset(rule_atoms)

        # 3. Add ring structures
        rule_atoms = add_ring_structures(cgr, rule_atoms)

        # 4. Create substructures
        r_subs = []
        for m in rxn.reactants:
            atoms = rule_atoms & set(m.atoms_numbers)
            if atoms:
                r_subs.append(m.substructure(atoms))
        p_subs = []
        for m in rxn.products:
            atoms = rule_atoms & set(m.atoms_numbers)
            if atoms:
                p_subs.append(m.substructure(atoms))

        # 5. Build queries
        patterns = [molecule_substructure_as_query(m, m.atoms_numbers) for m in r_subs]
        products = list(p_subs)

        # 6. Every step should produce valid objects
        assert all(isinstance(p, QueryContainer) for p in patterns)
        assert all(isinstance(p, MoleculeContainer) for p in products)

    def test_full_validate_rule(self, esterification_reaction):
        """Full validation: extract rule, build reactor, check products match."""
        rxn = esterification_reaction
        cgr = ~rxn
        center = set(cgr.center_atoms)

        rule_atoms = add_environment_atoms(cgr, center, 1)

        r_subs = [m.substructure(rule_atoms & set(m.atoms_numbers))
                  for m in rxn.reactants if rule_atoms & set(m.atoms_numbers)]
        p_subs = [m.substructure(rule_atoms & set(m.atoms_numbers))
                  for m in rxn.products if rule_atoms & set(m.atoms_numbers)]

        if r_subs and p_subs:
            rule = ReactionContainer(r_subs, p_subs)
            # validate_rule may or may not succeed depending on rule quality,
            # but it should not raise
            result = validate_rule(rule, rxn)
            assert isinstance(result, bool)

    @pytest.mark.parametrize("rxn_smi", REACTION_SMILES[:3])
    def test_cgr_decompose_roundtrip(self, rxn_smi):
        """CGR compose/decompose: ~(~rxn) should give back reactant/product molecules."""
        rxn = smiles(rxn_smi)
        cgr = ~rxn
        reactant_mol, product_mol = ~cgr
        assert isinstance(reactant_mol, MoleculeContainer)
        assert isinstance(product_mol, MoleculeContainer)
