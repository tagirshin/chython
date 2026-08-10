"""SynthonContainer, the `_token` dialect in SMILES and SMARTS, and the reactor post-hook."""

from subprocess import run
from sys import executable

import pytest

from chython import MoleculeContainer, QueryContainer, SynthonContainer, smarts, smiles, synthon_smiles
from chython.containers.synthon import restore_synthons
from chython.exceptions import IncorrectSmarts, IncorrectSmiles
from chython.files.daylight.tokenize import atom_re
from chython.periodictable import LABEL_TABLE, Element, Synthon
from chython.reactor import Reactor, Transformer


ROUND_TRIP = [
    '[NH2_nuc]C',  # bare
    'C[CH3_elec]',  # hcount
    '[13CH2_elec2]C',  # isotope + bivalent
    'CC(=O)[N_nuc](C)Cc1ccccc1',  # the worked example
    'c1ccccc1[c_elecB]1ccccc1',  # aromatic + elecB
    '[N_nuc*]C',  # radical variant token
    '[O-_nuc]C',  # charge before the label
    '[Se_neut2](C)C',  # two-letter symbol + neut2
    'O=[CH_elec]c1ccc([NH2_nuc])cc1',  # two labels, one molecule
    '[13CH2_elec2:7]C',  # isotope + label + atom map
]

# every token shape that parses today, to prove the new optional group re-brackets none of them
SHAPES = [
    'C',
    'CH3',
    'CH',
    'Cl',
    'Se',
    'c',
    'n',
    'nH',
    '13C',
    '13CH2',
    'C+',
    'C-',
    'N+1',
    'O-',
    'C@',
    'C@@',
    'C@H',
    'C:1',
    'CH3:12',
    '13CH2:7',
    'C+:3',
    'C@H:9',
    '13CH+:4',
    'Cl-:2',
    'te',
    'as',
    'C++',
    'CH4',
]

CORPUS = [
    'CC(=O)N(C)Cc1ccccc1',
    'c1ccc2c(c1)[nH]c1ccccc12',
    'OC(=O)c1ccc(N)cc1',
    'C1CC2CCC1CC2',
    'CC(C)(C)OC(=O)N1CCC(N)CC1',
    'Clc1ccc(B(O)O)cc1',
    'O=S(=O)(N)c1ccccc1F',
    'C/C=C/C(=O)OCC',
    'N[C@@H](C)C(=O)O',
    'c1ccsc1',
]

AMIDE_CUT = (
    '[C;$(C([#7])(=[#8])[#6]):1]-!@[#7;A;+0;$([#7;D2]([#6])[#6]),$([#7;D3]([#6])([#6])[#6]):2]>>[#6_elec:1].[#7_nuc:2]'
)


def query_labels(q):
    return {n: a._label for n, a in q.atoms() if getattr(a, '_label', None) is not None}


def apply_rule(rule_smarts, structure):
    """What SynthonTransformer does, inlined: stamp the replacement's tokens on the products."""
    lhs, rhs = rule_smarts.split('>>')
    pattern, replacement = smarts(lhs), smarts(rhs)
    labels = query_labels(replacement)
    tr = Transformer(pattern, replacement)
    out = []
    for mapping in tr._pattern.get_mapping(structure, automorphism_filter=tr._automorphism_filter):
        product = tr._patcher(structure, mapping)
        out.append(restore_synthons(product, {mapping[n]: t for n, t in labels.items() if n in mapping}))
    return out


@pytest.mark.parametrize('text', ROUND_TRIP)
def test_smiles_round_trip(text):
    m = synthon_smiles(text)
    assert m.synthon_labels
    assert str(m) == str(synthon_smiles(str(m)))


@pytest.mark.parametrize('label', sorted(LABEL_TABLE))
def test_every_token_round_trips(label):
    m = synthon_smiles(f'C[CH2_{label}]')
    assert list(m.synthon_labels.values()) == [label]


def test_label_and_map_are_independent_fields():
    m = synthon_smiles('[13CH2_elec2:7]C')
    assert '_elec2:7' in format(m, 'm')


def test_eight_labels_exactly():
    assert set(LABEL_TABLE) == {'elec', 'nuc', 'elec2', 'nuc2', 'neut2', 'elec*', 'nuc*', 'elecB'}


@pytest.mark.parametrize('token', SHAPES)
def test_atom_re_did_not_change_an_existing_token(token):
    g = atom_re.fullmatch(token).groups()
    assert g[5] is None  # the label group never claims an existing shape
    assert str(smiles(f'[{token}]')) == str(synthon_smiles(f'[{token}]'))


def test_labelled_smiles_in_a_plain_container_fails_loudly():
    with pytest.raises(TypeError):
        smiles('[CH3_elec]C')


@pytest.mark.parametrize(
    'text,want',
    [
        ('[#6_elec:1]', {1: 'elec'}),
        ('[#7_nuc:2]', {2: 'nuc'}),
        ('[C_elec2:1]', {1: 'elec2'}),
        ('[c_elecB:3]', {3: 'elecB'}),
        ('[#6;$(C(=O)):1]', {}),
        ('[#6_elec;$(C(=O)):1]', {1: 'elec'}),
        ('[#6;$(C(=O))_elec:1]', {1: 'elec'}),  # order-free, as Daylight primitives are
        ('[N_nuc*:4]', {4: 'nuc*'}),
        ('[#6_elec:1].[#7_nuc:2]', {1: 'elec', 2: 'nuc'}),
    ],
)
def test_smarts_carries_the_token(text, want):
    assert query_labels(smarts(text)) == want


@pytest.mark.parametrize('text', ['[C_Q:1]', '[C_:1]', '[C_nuc3:1]'])
def test_bad_token_raises(text):
    with pytest.raises(IncorrectSmarts):
        smarts(text)


def test_bogus_token_rejected_in_smiles_too():
    with pytest.raises(IncorrectSmiles):
        synthon_smiles('[C_Zz]')


def test_label_does_not_change_what_the_query_matches():
    target = smiles('CC(=O)NC')
    target.canonicalize()
    plain = len(list(smarts('[#6;$(C(=O)[#7]):1]').get_mapping(target)))
    labelled = len(list(smarts('[#6_elec;$(C(=O)[#7]):1]').get_mapping(target)))
    assert plain == labelled > 0


def test_reactant_side_label_is_inert_which_is_why_loaders_must_raise():
    target = smiles('CCO')
    target.canonicalize()
    assert len(list(smarts('[C_elec:1]').get_mapping(target))) == len(list(smarts('[C:1]').get_mapping(target))) == 2


def test_smarts_writer_drops_the_token_so_never_round_trip_a_rule():
    q = smarts('[C_elec2;H2:1]-[#7_nuc:2]')
    assert q.to_smarts() == '[C_elec2;H2:1]-[#7_nuc:2]'  # the cached source string
    assert '_elec2' not in q.to_smarts(mapping_renumber={1: 1, 2: 2})
    assert '_elec' not in str(Reactor.from_smarts('[C:1]-[#7:2]>>[#6_elec:1].[#7_nuc:2]'))


def test_reactor_products_carry_the_labels_through_the_union():
    rx = Reactor.from_smarts(AMIDE_CUT)
    # ExtendedQuery.copy enumerates its slots by hand: without _label there, this is empty
    assert query_labels(rx._replacement) == {1: 'elec', 2: 'nuc'}
    assert all(not query_labels(p) for p in rx._patterns)


def test_patcher_rebuilds_rc_atoms_as_plain_elements():
    target = synthon_smiles('CC(=O)N(C)Cc1ccccc1')
    target.canonicalize()
    products = next(iter(Reactor.from_smarts(AMIDE_CUT)(target))).products
    assert all(isinstance(p, SynthonContainer) for p in products)
    # this is why the post-hook survives the inline product token
    assert any(not isinstance(a, Synthon) for p in products for _, a in p.atoms())
    assert all(not p.synthon_labels for p in products)


def test_restore_synthons_rewraps_and_stamps():
    target = synthon_smiles('CC(=O)N(C)Cc1ccccc1')
    target.canonicalize()
    got = apply_rule(AMIDE_CUT, target)
    assert len(got) == 1
    assert all(isinstance(a, Synthon) for _, a in got[0].atoms())
    assert sorted(got[0].synthon_labels.values()) == ['elec', 'nuc']
    assert len(got[0].split()) == 2


def test_hierarchical_cut_keeps_the_first_label():
    pre = synthon_smiles('C(CO[CH3_nuc])C(=O)Cl')
    pre.canonicalize()
    out = apply_rule('[C;$(C(=O)[Cl]):1]-[Cl:2]>>[#6_elec:1]', pre)
    assert out and set(out[0].synthon_labels.values()) == {'elec', 'nuc'}


def test_a_product_only_atom_gets_its_token():
    # _patcher extends the mapping dict in place, so a map number absent from the pattern
    # still resolves to a product atom
    acyl = synthon_smiles('CCCC(=O)Cl')
    acyl.canonicalize()
    out = apply_rule('[C;$(C(=O)[Cl]):1]-[Cl:2]>>[#6_elec:1]-[#8_nuc:3]', acyl)
    assert out and sorted(out[0].synthon_labels.values()) == ['elec', 'nuc']


@pytest.mark.parametrize('text', CORPUS)
def test_unlabelled_synthon_canonicalises_like_the_molecule(text):
    a = smiles(text)
    a.canonicalize()
    b = synthon_smiles(text)
    b.canonicalize()
    assert isinstance(b, SynthonContainer)
    assert str(a) == str(b)
    assert hash(a) == hash(b)
    for n in a._atoms:
        assert hash(a._atoms[n]) == hash(b._atoms[n])
        assert a._atoms[n] == b._atoms[n]


def test_a_labelled_atom_hashes_differently():
    labelled = synthon_smiles('CC(=O)[N_nuc](C)C')
    plain = smiles('CC(=O)N(C)C')
    assert hash(labelled._atoms[4]) != hash(plain._atoms[4])
    assert str(labelled.unlabelled()) == str(synthon_smiles('CC(=O)N(C)C'))


def test_label_aware_hash_fixes_the_canonical_form():
    a = synthon_smiles('[CH3_elec]c1ccc([CH3_nuc])cc1')
    b = synthon_smiles('[CH3_nuc]c1ccc([CH3_elec])cc1')
    a.canonicalize()
    b.canonicalize()
    assert str(a) == str(b)


def test_canonical_form_is_reproducible_across_processes():
    # hashing the token STRING instead of LABEL_INDEX gives a different answer per process,
    # and nothing inside one process notices.
    code = (
        "from chython import synthon_smiles as s;m = s('Cc1cccc(C)c1[NH_nuc]C(=O)CN(CC)CC');m.canonicalize();print(m)"
    )
    seen = {run([executable, '-c', code], capture_output=True, text=True, check=True).stdout for _ in range(5)}
    assert len(seen) == 1, seen


def test_symmetric_wildcard_eq():
    e = Synthon.from_symbol('C')(synthon_label='elec')
    n = Synthon.from_symbol('C')(synthon_label='nuc')
    bare = Synthon.from_symbol('C')()
    plain = Element.from_symbol('C')()
    assert e != n and n != e
    assert e == bare and bare == e
    assert e == plain and plain == e


def test_a_synthon_container_may_hold_plain_elements():
    s = synthon_smiles('[CH3_elec]CCO')
    s.canonicalize()
    u = s.union(smiles('CC'), remap=True)
    assert isinstance(u, SynthonContainer)
    assert any(not isinstance(a, Synthon) for _, a in u.atoms())
    assert u.synthon_labels  # getattr, not attribute access
    str(u)
    s.explicify_hydrogens()
    str(s)


def test_copy_and_pickle_keep_the_label():
    from pickle import dumps, loads

    m = synthon_smiles('[CH3_elec]CCO')
    assert m.copy().synthon_labels == m.synthon_labels
    assert loads(dumps(m)).synthon_labels == m.synthon_labels
    assert m.substructure(list(m)[:2]).synthon_labels


def test_depict_is_inert_without_labels_and_renders_with_them():
    plain = smiles('O=Cc1ccc(N)cc1')
    plain.clean2d()
    labelled = synthon_smiles('O=[CH_elec]c1ccc([NH2_nuc])cc1')
    labelled.clean2d()
    assert '-synthon' not in synthon_smiles('O=Cc1ccc(N)cc1').depict()
    svg = labelled.depict()
    assert '-synthon' in svg and '#D1495B' in svg and '#00798C' in svg


# --- the !rN blacklist -------------------------------------------------------------------


def test_excluded_ring_sizes_parses():
    q = smarts('[c;!r3;!r4;!r5;!r6;!r7;!r8;!r9;!r10;!r11:1]')
    assert q.atom(1).excluded_ring_sizes == (3, 4, 5, 6, 7, 8, 9, 10, 11)


def test_excluded_ring_sizes_actually_bites():
    # a ring-fusion atom sits in both a small ring and the macrocycle, so r12 accepts it
    # where !r6 rejects it - the whole reason a positive enumeration is not the same predicate
    m = smiles('O=C1CCCCCCCc2ccccc2CN1')
    m.canonicalize()
    assert len(list(smarts('[c;r12:1]').get_mapping(m, automorphism_filter=False))) == 2
    guarded = smarts('[c;!r3;!r4;!r5;!r6;!r7;!r8;!r9;!r10;!r11:1]')
    assert len(list(guarded.get_mapping(m, automorphism_filter=False))) == 0


def test_excluded_ring_sizes_leaves_the_cython_fast_path():
    q = smarts('[c;!r6:1]')
    assert q._has_extended_query  # or the cython matcher would ignore it silently


def test_excluded_ring_sizes_survives_copy():
    q = smarts('[c;!r6:1]')
    assert q.atom(1).copy().excluded_ring_sizes == (6,)
    assert (q | smarts('[C:2]')).atom(1).excluded_ring_sizes == (6,)


def test_excluded_ring_sizes_accepts_a_bigger_ring():
    macro = smiles('O=C1CCCCCCCCCCCCN1')
    macro.canonicalize()
    guarded = smarts('[C;!r3;!r4;!r5;!r6;!r7;!r8;!r9;!r10;!r11:1]')
    assert list(guarded.get_mapping(macro, automorphism_filter=False))


def test_plain_containers_are_untouched():
    q = smarts('[C;$(C(=O)([#6])):1][Cl,F,Br:2]')
    assert type(q) is QueryContainer
    target = smiles('CCC(=O)Cl')
    target.canonicalize()
    assert type(target) is MoleculeContainer
    assert (
        len(list(q.get_mapping(target, _cython=True, automorphism_filter=False)))
        == len(list(q.get_mapping(target, _cython=False, automorphism_filter=False)))
        == 1
    )
