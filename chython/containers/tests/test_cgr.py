import unittest
import pickle

# Assuming chython is in PYTHONPATH or installed
# If running this test file directly, ensure paths are set up correctly.
from chython import smiles
from chython.containers.cgr import CGRContainer
from chython.containers.molecule import MoleculeContainer
from chython.containers.cgr_query import QueryCGRContainer
from chython.periodictable.base.dynamic import DynamicElement
from chython.periodictable.base.element import Element
from chython.periodictable import C
from chython.exceptions import AtomNotFound
# Explicitly import some elements to ensure periodictable __init__ runs fully
from chython.periodictable import C, O, N 
from chython.containers.bonds import DynamicBond

# Helper function to create a CGRContainer for testing
def create_simple_cgr():
    """
    Creates a simple CGR:
    Reactant: C1=C2-C3
    Product:  C1-C2=C3
    Atom maps: 1, 2, 3
    """
    cgr = CGRContainer()
    c1_atom = DynamicElement.from_atomic_number(atomic_number=6)(0) # Carbon
    c2_atom = DynamicElement.from_atomic_number(atomic_number=6)(0)
    c3_atom = DynamicElement.from_atomic_number(atomic_number=6)(0)

    cgr.add_atom(c1_atom, 1, p_charge=0, p_is_radical=False)
    cgr.add_atom(c2_atom, 2, p_charge=0, p_is_radical=False)
    cgr.add_atom(c3_atom, 3, p_charge=0, p_is_radical=False)

    cgr.add_bond(1, 2, DynamicBond(2, 1))  # C1=C2 -> C1-C2
    cgr.add_bond(2, 3, DynamicBond(1, 2))  # C2-C3 -> C2=C3
    return cgr

class TestCGRContainer(unittest.TestCase):

    def test_initial_creation_from_reaction_and_str(self):
        reaction_smiles = "[CH2:5]=[CH2:6].[CH:2]([CH:3]=[CH2:4])=[CH2:1]>>[CH:3]1=[CH:4][CH2:6][CH2:5][CH2:2][CH2:1]1"
        try:
            reaction = smiles(reaction_smiles)
        except Exception as e:
            self.fail(f"SMILES parsing failed for reaction: {e}")
        
        cgr = ~reaction
        
        self.assertIsInstance(cgr, CGRContainer)
        
        cgr_smiles_str = str(cgr)
        self.assertIsInstance(cgr_smiles_str, str)
        self.assertTrue(len(cgr_smiles_str) > 0)

        # Smiles layout order can differ but should remain consistent and non-empty.
        self.assertIn(cgr_smiles_str, (
            'C[.>-]1[=>-]2[->.]C(=C[.>-]C[=>-]C[.>-]1)[.>-]C[=>-]2',
            'C[.>-]1[=>-]C[.>-]C=C[->.]2[.>-]C[=>-]C[.>-]1[->.]2',
        ))

        self.assertEqual(len(cgr), 6)
        atom_ids_in_cgr = set(iter(cgr))
        self.assertEqual(atom_ids_in_cgr, {1, 2, 3, 4, 5, 6})

    def test_substructure(self):
        reaction_smiles = "[CH2:5]=[CH2:6].[CH:2]([CH:3]=[CH2:4])=[CH2:1]>>[CH:3]1=[CH:4][CH2:6][CH2:5][CH2:2][CH2:1]1"
        reaction = smiles(reaction_smiles)
        cgr = ~reaction
        
        sub = cgr.substructure([1, 2, 3, 7]) # Atom 7 is not in CGR
        
        self.assertIsInstance(sub, CGRContainer)
        self.assertEqual(len(sub), 3)
        
        self.assertTrue(all(sub.has_atom(i) for i in [1, 2, 3]))
        self.assertFalse(sub.has_atom(4))
        self.assertFalse(sub.has_atom(7))
        
        for atom_id in [1, 2, 3]:
            self.assertEqual(sub._atoms[atom_id].charge, cgr._atoms[atom_id].charge)
            # Ensure p_charges and p_radicals are correctly copied or defaulted
            self.assertEqual(sub._p_charges.get(atom_id, -99), cgr._p_charges.get(atom_id, -99))
            self.assertEqual(sub._p_radicals.get(atom_id, None), cgr._p_radicals.get(atom_id, None))


        # Check bonds within the substructure
        expected_bonds_in_sub = {(1,2), (1,3), (2,3)}
        actual_bonds_in_sub = set()
        for n1, n2, _ in sub.bonds():
            actual_bonds_in_sub.add(tuple(sorted((n1,n2))))
        self.assertEqual(actual_bonds_in_sub, expected_bonds_in_sub)

        for n1_orig, n2_orig, bond_orig in cgr.bonds():
            if sub.has_atom(n1_orig) and sub.has_atom(n2_orig):
                self.assertTrue(sub.has_bond(n1_orig, n2_orig))
                self.assertEqual(sub.bond(n1_orig,n2_orig).order, bond_orig.order)
                self.assertEqual(sub.bond(n1_orig,n2_orig).p_order, bond_orig.p_order)
            else: # If one atom is outside sub, bond should not be in sub
                # This bond connects an atom in the substructure to one outside.
                # It should NOT exist in the substructure itself.
                if sub.has_atom(n1_orig) != sub.has_atom(n2_orig):
                    try:
                        self.assertFalse(sub.has_bond(n1_orig, n2_orig))
                    except AtomNotFound:
                        pass


    def test_substructure_as_query(self):
        reaction_smiles = "[CH2:5]=[CH2:6].[CH:2]([CH:3]=[CH2:4])=[CH2:1]>>[CH:3]1=[CH:4][CH2:6][CH2:5][CH2:2][CH2:1]1"
        reaction = smiles(reaction_smiles)
        cgr = ~reaction
        
        sub_q = cgr.substructure([1, 2], as_query=True)
        self.assertIsInstance(sub_q, QueryCGRContainer)
        self.assertEqual(len(sub_q), 2)
        
        if 1 in sub_q._hybridizations: self.assertIsInstance(sub_q._hybridizations[1], tuple)
        if 1 in sub_q._neighbors: self.assertIsInstance(sub_q._neighbors[1], tuple)
        if 1 in cgr._hybridizations and 1 in sub_q._hybridizations:
             self.assertEqual(sub_q._hybridizations[1][0], cgr._hybridizations[1])


    def test_augmented_substructure(self):
        reaction_smiles = "[CH2:5]=[CH2:6].[CH:2]([CH:3]=[CH2:4])=[CH2:1]>>[CH:3]1=[CH:4][CH2:6][CH2:5][CH2:2][CH2:1]1"
        reaction = smiles(reaction_smiles)
        cgr = ~reaction

        aug_sub = cgr.augmented_substructure([1], deep=1)
        self.assertEqual(set(aug_sub), {1, 2, 3})

        aug_sub_d2 = cgr.augmented_substructure([1], deep=2)
        self.assertEqual(set(aug_sub_d2), {1, 2, 3, 4, 5})

    def test_add_delete_atom_bond_basic(self):
        cgr = CGRContainer()
        c_atom = DynamicElement.from_atomic_number(atomic_number=6)(0)
        
        m1 = cgr.add_atom(c_atom.copy(), 1, p_charge=1, p_is_radical=True)
        self.assertEqual(m1, 1)
        self.assertTrue(cgr.has_atom(1))
        self.assertEqual(cgr._atoms[1].charge, 0)
        self.assertEqual(cgr._p_charges[1], 1)
        self.assertTrue(cgr._p_radicals[1])

        m2 = cgr.add_atom(c_atom.copy(), 2)
        cgr.add_bond(m1, m2, DynamicBond(1, 2)) # R: 1-2, P: 1=2
        self.assertTrue(cgr.has_bond(1, 2))
        self.assertEqual(cgr.bond(1,2).order, 1)
        self.assertEqual(cgr.bond(1,2).p_order, 2)
        # Check hybridization update
        self.assertEqual(cgr._hybridizations[1], 1) # Based on bond order 1
        self.assertEqual(cgr._p_hybridizations[1], 2) # Based on p_order 2

        cgr.delete_bond(1, 2)
        self.assertFalse(cgr.has_bond(1,2))
        self.assertEqual(cgr._hybridizations[1], 1) # Recalculated to default
        self.assertEqual(cgr._p_hybridizations[1], 1)

        cgr.delete_atom(1)
        self.assertFalse(cgr.has_atom(1))

    def test_atoms_bonds_iterators_properties(self):
        cgr = create_simple_cgr()
        
        atoms_list = list(cgr.atoms())
        self.assertEqual(len(atoms_list), 3)
        self.assertEqual(set(a[0] for a in atoms_list), {1,2,3})

        bonds_list = list(cgr.bonds())
        self.assertEqual(len(bonds_list), 2)
        
        self.assertEqual(cgr.atoms_count, 3)
        self.assertEqual(set(cgr.atoms_numbers), {1,2,3})

    def test_center_atoms(self):
        cgr = CGRContainer()
        c1_atom = DynamicElement.from_atomic_number(atomic_number=6)(0)
        c2_atom_for_test = DynamicElement.from_atomic_number(atomic_number=6)(0)
        c3_atom = DynamicElement.from_atomic_number(atomic_number=6)(0)

        cgr.add_atom(c1_atom, 1)
        cgr.add_atom(c2_atom_for_test, 2, p_charge=1) # This will make atom 2 dynamic if its charge is 0
        cgr.add_atom(c3_atom, 3)

        cgr.add_bond(1, 2, DynamicBond(1, 2)) # Bond 1-2 is dynamic
        cgr.add_bond(2, 3, DynamicBond(1, 1)) # Bond 2-3 is not

        self.assertEqual(set(cgr.center_atoms), {1, 2}) # Atom 2 (dynamic) and Atom 1 (part of dynamic bond)

    def test_neighbors_method(self):
        cgr = create_simple_cgr()
        self.assertEqual(cgr.neighbors(1), (1, 1)) 
        self.assertEqual(cgr.neighbors(2), (2, 2))
        self.assertEqual(cgr.neighbors(3), (1, 1))

    def test_remap(self):
        cgr_orig = create_simple_cgr()
        cgr_orig._conformers.append({1:(0.,0.,0.), 2:(1.,0.,0.)})
        mapping = {1: 10, 2: 20, 3: 30}
        
        cgr_copy = cgr_orig.copy()
        cgr_copy.remap(mapping) 
        
        self.assertEqual(len(cgr_copy), 3)
        self.assertTrue(all(cgr_copy.has_atom(i) for i in [10, 20, 30]))
        self.assertEqual(cgr_copy.bond(10,20).order, 2)
        self.assertEqual(len(cgr_copy._conformers), 1)
        self.assertTrue(10 in cgr_copy._conformers[0])

        cgr_new = cgr_orig.remap(mapping, copy=True)
        self.assertIsNot(cgr_new, cgr_orig)
        self.assertTrue(cgr_orig.has_atom(1)) # Original unchanged
        self.assertTrue(cgr_new.has_atom(10))

    def test_copy(self):
        cgr_orig = create_simple_cgr()
        cgr_orig._p_charges[1] = 5
        cgr_orig._conformers.append({1:(1.0,2.0,3.0)})

        cgr_copy = cgr_orig.copy()

        self.assertIsNot(cgr_orig, cgr_copy)
        self.assertNotEqual(id(cgr_orig._atoms), id(cgr_copy._atoms))
        self.assertNotEqual(id(cgr_orig._p_charges), id(cgr_copy._p_charges))
        self.assertNotEqual(id(cgr_orig._conformers), id(cgr_copy._conformers))
        if cgr_orig._conformers:
            self.assertNotEqual(id(cgr_orig._conformers[0]), id(cgr_copy._conformers[0]))

        self.assertEqual(cgr_copy._p_charges[1], 5)
        
        cgr_copy.add_atom(DynamicElement.from_atomic_number(atomic_number=7)(), 4)
        self.assertNotEqual(len(cgr_orig), len(cgr_copy))

    def test_decompose_invert_operator(self):
        cgr = create_simple_cgr()
        reactants, products = ~cgr

        self.assertIsInstance(reactants, MoleculeContainer)
        self.assertEqual(len(reactants), 3)
        self.assertEqual(reactants.bond(1,2).order, 2)
        self.assertEqual(products.bond(2,3).order, 2)

    def test_compose_xor_operator_mols(self):
        mol_r = MoleculeContainer(); mol_r.add_atom(Element.from_symbol('C')(),1); mol_r.add_atom(Element.from_symbol('C')(),2); mol_r.add_bond(1,2,1)
        mol_p = MoleculeContainer(); mol_p.add_atom(Element.from_symbol('C')(),1); mol_p.add_atom(Element.from_symbol('C')(),2); mol_p.add_bond(1,2,2)
        cgr = mol_r ^ mol_p
        self.assertIsInstance(cgr, CGRContainer)
        self.assertEqual(cgr.bond(1,2).order, 1)
        self.assertEqual(cgr.bond(1,2).p_order, 2)

    def test_union_or_operator_cgrs(self):
        cgr1 = CGRContainer()
        cgr1.add_atom(DynamicElement.from_atomic_number(atomic_number=6)(0), 1)
        cgr2 = CGRContainer()
        cgr2.add_atom(DynamicElement.from_atomic_number(atomic_number=7)(0), 10) # Use distinct atom ID
        
        united_cgr = cgr1 | cgr2 # Default remap=True will handle overlapping if keys were same
        self.assertEqual(len(united_cgr), 2)
        self.assertTrue(united_cgr.has_atom(1) and united_cgr.has_atom(10))

    def test_pickle_unpickle(self):
        cgr_original = create_simple_cgr()
        cgr_original._p_charges[1] = 1
        cgr_original._conformers.append({1:(0.,0.,0.), 2:(1.,0.,0.), 3:(2.,0.,0.)})

        pickled_cgr = pickle.dumps(cgr_original)
        unpickled_cgr = pickle.loads(pickled_cgr)

        self.assertIsInstance(unpickled_cgr, CGRContainer)
        self.assertEqual(len(cgr_original), len(unpickled_cgr))
        self.assertEqual(str(cgr_original), str(unpickled_cgr))
        self.assertEqual(cgr_original._p_charges, unpickled_cgr._p_charges)
        self.assertEqual(cgr_original._hybridizations, unpickled_cgr._hybridizations)
        self.assertEqual(cgr_original._conformers, unpickled_cgr._conformers)


class TestCGRBehavior(unittest.TestCase):
    def test_substructure_returns_cgr(self):
        rxn = smiles('[CH3:1][CH2:2][OH:3]>>[CH3:1][CH3:2].[OH2:3]')
        cgr = ~rxn

        sub = cgr.substructure([1])

        self.assertIsInstance(sub, CGRContainer)
        self.assertNotIsInstance(sub, tuple)
        self.assertEqual(set(sub), {1})

    def test_augmented_substructure_returns_cgr(self):
        rxn = smiles('[CH3:1][CH2:2][OH:3]>>[CH3:1][CH3:2].[OH2:3]')
        cgr = ~rxn

        aug = cgr.augmented_substructure([1], deep=1)

        self.assertIsInstance(aug, CGRContainer)
        self.assertEqual(set(aug), {1, 2})

    def test_dynamic_element_factories_keep_dynamic_subclasses(self):
        dyn_cls = DynamicElement.from_symbol('C')
        self.assertTrue(isinstance(dyn_cls, type))
        self.assertEqual(dyn_cls.__name__, 'DynamicC')

        dyn_cls_num = DynamicElement.from_atomic_number(6)
        self.assertIs(dyn_cls_num, dyn_cls)

        dyn_atom = dyn_cls()
        self.assertEqual(dyn_atom.atomic_symbol, 'C')
        self.assertEqual(dyn_atom.atomic_number, 6)

        elem = C()
        dyn_from_elem = DynamicElement.from_atom(elem)
        self.assertIsInstance(dyn_from_elem, DynamicElement)
        self.assertEqual(dyn_from_elem.atomic_symbol, elem.atomic_symbol)
        self.assertEqual(dyn_from_elem.charge, elem.charge)
        self.assertEqual(dyn_from_elem.p_charge, elem.charge)

if __name__ == '__main__':
    unittest.main() 
