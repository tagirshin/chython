# -*- coding: utf-8 -*-
#
#  Copyright 2022-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from re import compile, match, search
from itertools import permutations
from ...containers.bonds import QueryBond
from ...exceptions import IncorrectSmiles, IncorrectSmarts
from ...periodictable import elements
from ...periodictable.base.synthon import SYNTHON_LABELS


# -,= OR bonds supported
# !: NOT bonds supported
# ~ ANY bonds supported
# @ bond not supported
# @;!: and any other complicated combinations not supported


# tokens structure:
# (type: int, value)
# types:
# 0: atom
# 1: bond
# 2: open chain (
# 3: close chain )
# 4: dot bond .
# 5: in bracket raw data []
# 6: closure number
# 7: raw closure number
# 8: aromatic atom
# 9: up down bond

# 10: query OR bond
# 11: query NOT bond
# 12: in ring bond
# 13: CGR bond
# 14: CGR atom


iso_re = compile(r'^[0-9]+')
chg_re = compile(r'[+-][0-4+-]?')
mpp_re = compile(r':[1-9][0-9]*$')
str_re = compile(r'@[@?]?')
not_charge_re = compile(r'![+-]')

replace_dict = {'-': 1, '=': 2, '#': 3, ':': 4, '~': 8}
not_dict = {'-': [2, 3, 4], '=': [1, 3, 4], '#': [1, 2, 4], ':': [1, 2, 3]}
# the synthon label is a new bracket field between charge and the atom map. `_` cannot occur in a
# bracket today, so the optional group changes nothing that parses now.
synthon_token_re = compile(r'_elec(?:2|\*|B)?|_nuc(?:2|\*)?|_neut2')
atom_re = compile(r'([1-9][0-9]{0,2})?([A-IK-PR-Zacnopsbt][a-ik-pr-vy]?)(@@|@)?(H[1-4]?)?([+-][1-4+-]?)?'
                  r'(_elec(?:2|\*|B)?|_nuc(?:2|\*)?|_neut2)?(:[0-9]{1,4})?')
dyn_atom_re = compile(r'([1-9][0-9]{0,2})?([A-IK-PR-Zacnopsb][a-ik-pr-vy]?)([+-0][1-4+-]?(>[+-0][1-4+-]?)?)?([*^](>[*^])?)?')
charge_dict = {'+': 1, '+0': 0, '+1': 1, '++': 2, '+2': 2, '+3': 3, '+++': 3, '+4': 4, '++++': 4,
               '-': -1, '-0': 0, '-1': -1, '--': -2, '-2': -2, '-3': -3, '---': -3, '-4': -4, '----': -4}
dynamic_bonds = {'.>-': (None, 1), '.>=': (None, 2), '.>#': (None, 3), '.>:': (None, 4), '.>~': (None, 8),
                 '->.': (1, None), '->=': (1, 2), '->#': (1, 3), '->:': (1, 4), '->~': (1, 8),
                 '=>.': (2, None), '=>-': (2, 1), '=>#': (2, 3), '=>:': (2, 4), '=>~': (2, 8),
                 '#>.': (3, None), '#>-': (3, 1), '#>=': (3, 2), '#>:': (3, 4), '#>~': (3, 8),
                 ':>.': (4, None), ':>-': (4, 1), ':>=': (4, 2), ':>#': (4, 3), ':>~': (4, 8),
                 '~>.': (8, None), '~>-': (8, 1), '~>=': (8, 2), '~>#': (8, 3), '~>:': (8, 4)}
dyn_charge_dict = {'-4': -4, '-3': -3, '-2': -2, '--': -2, '-': -1, '0': 0, '+': 1, '++': 2, '+2': 2, '+3': 3, '+4': 4}
tmp = {f'{i}>{j}': (x, ('charge', y - x)) for (i, x), (j, y) in permutations(dyn_charge_dict.items(), 2) if x != y}
dyn_charge_dict = {k: (v,) for k, v in dyn_charge_dict.items()}
dyn_charge_dict.update(tmp)
dyn_radical_dict = {'*': (True,), '*>^': (True, ('radical', None)), '^>*': (False, ('radical', None))}

_aromatic_upper = {'c': 'C', 'n': 'N', 'o': 'O', 'p': 'P', 's': 'S', 'b': 'B', 'as': 'As', 'se': 'Se', 'te': 'Te'}


def _is_element_token(token):
    """Check if a token looks like an element symbol (possibly negated or atomic number)."""
    if token in ('*', 'a', 'A'):
        return True
    raw = token[1:] if token.startswith('!') else token
    if raw.startswith('#'):
        try:
            int(raw[1:])
            return True
        except ValueError:
            return False
    if len(raw) == 1 and raw[0].isupper():
        return True
    if len(raw) == 2 and raw[0].isupper() and raw[1].islower():
        return True
    return False


def _tokenize(smiles):
    token_type = token = None
    bracket_depth = 0
    tokens = []
    for s in smiles:
        if token_type == 12:  # -;, =;, #; or :; found.
            if s == '!':  # not in ring
                if token:  # !! case
                    raise IncorrectSmarts('Invalid ring bond token')
                token = True
            elif s == '@':
                tokens.append((12, QueryBond(tokens.pop(-1)[1], not token)))
                token_type = token = None  # finalize [not]ring bond
            else:
                raise IncorrectSmarts('Invalid ring bond token')
        # [atom block parser]
        elif s == '[':  # open complex token
            if token_type == 5:  # nested [ inside bracket (e.g. $([...]))
                bracket_depth += 1
                token.append(s)
            elif token_type in (10, 11):
                raise IncorrectSmarts('Query bond invalid')
            elif token_type == 7:  # empty closure
                raise IncorrectSmiles('invalid closure')
            else:
                if token:
                    tokens.append((token_type, token))
                token = []
                token_type = 5
                bracket_depth = 0
        elif s == ']':  # close complex token
            if token_type != 5:
                raise IncorrectSmiles(']..]')
            elif bracket_depth > 0:
                bracket_depth -= 1
                token.append(s)
            elif not token:
                raise IncorrectSmiles('empty [] brackets')
            else:
                if '>' in token:
                    s = ''.join(token)
                    if len(s) == 3:  # bond only possible
                        try:
                            tokens.append((13, dynamic_bonds[s]))
                        except KeyError:
                            raise IncorrectSmiles(f'invalid dynamic bond token {s}')
                    else:  # dynamic atom token
                        tokens.append(_dyn_atom_parse(s))
                else:
                    tokens.append((5, ''.join(token)))
                token = None
                token_type = 0  # mark as atom
        elif token_type == 5:  # grow token with brackets. skip validation
            token.append(s)
        # closure parser
        elif s.isnumeric():  # closures
            if token_type in (10, 11):
                raise IncorrectSmarts('Query bond invalid')
            elif token_type == 2:
                raise IncorrectSmiles('(1 case invalid')
            elif token_type == 7:  # % already found. collect number
                if not token and s == '0':
                    raise IncorrectSmiles('number starts with 0')
                token.append(s)
                if len(token) == 2:
                    tokens.append((6, int(''.join(token))))
                    token = None
                    token_type = 6  # mark finished
            else:
                if s == '0':
                    raise IncorrectSmiles('number starts with 0')
                elif token:
                    tokens.append((token_type, token))
                    token = None
                token_type = 6
                tokens.append((6, int(s)))
        elif token_type == 7:
            raise IncorrectSmiles('expected closure number')
        elif s == '%':
            if token_type in (10, 11):
                raise IncorrectSmarts('Query bond invalid')
            elif token_type == 2:
                raise IncorrectSmiles('(%10 case invalid')
            elif token:
                tokens.append((token_type, token))
            token_type = 7
            token = []
        # bonds parser
        elif s in '=#:-~':  # bonds found
            if token_type == 10:
                token.append(replace_dict[s])
                tokens.append((10, token))
                token_type = token = None  # finalize token
            elif token_type == 11:
                token_type = None
                tokens.append((10, not_dict[s]))
            else:
                if token:
                    tokens.append((token_type, token))
                    token = None
                token_type = 1
                tokens.append((1, replace_dict[s]))
        elif s == '@':  # ring bond
            if token_type == 11:
                tokens.append((12, QueryBond(8, False)))  # bare !@ - any non-ring bond
            elif token_type == 1 or token_type is None and tokens and tokens[-1][0] == 10:
                tokens.append((12, QueryBond(tokens.pop(-1)[1], True)))  # order glued to the mark: -@
            else:
                if token:
                    tokens.append((token_type, token))
                tokens.append((12, QueryBond(8, True)))  # bare @ - a ring bond of any order
            token = token_type = None
        elif token_type in (10, 11):  # expected bond symbol
            raise IncorrectSmarts('query bond invalid')
        elif s in r'\/':
            if token:
                tokens.append((token_type, token))
                token = None
            token_type = 9
            tokens.append((9, s == '/'))  # Up is true
        elif s == '.':
            if token:
                tokens.append((token_type, token))
                token = None
            token_type = 4
            tokens.append((4, None))
        elif s in ';&':  # ;@ or ;!@ - in ring or not in ring bond
            if token_type is not None and token_type != 1:  # start of smiles or bond, list of bonds and not bond.
                raise IncorrectSmarts('Ring bond token invalid')
            token_type = 12
        elif s == ',':  # query bond separator
            if token_type == 1:
                token = [tokens.pop(-1)[1]]
            elif token_type is None and tokens and tokens[-1][0] == 10:
                token = tokens.pop(-1)[1]  # third and further alternative: -,=,#
            else:
                raise IncorrectSmarts('Query bond invalid')
            token_type = 10
        elif s == '!':  # query not bond
            if token_type == 1 or token_type is None and tokens and tokens[-1][0] == 10:
                token_type = 12  # ring mark glued to a bond: -!@
                token = True
            elif token_type not in (0, 2, 3, 6, 8):  # closures, brackets or atoms expected
                raise IncorrectSmarts('Query bond invalid')
            else:
                if token:
                    tokens.append((token_type, token))
                    token = None
                token_type = 11
        # brackets
        elif s == '(':
            if token_type == 2:  # barely opened
                raise IncorrectSmiles('((')
            elif token:
                tokens.append((token_type, token))
                token = None
            token_type = 2
            tokens.append((2, None))
        elif s == ')':
            if token_type == 2:  # barely opened
                raise IncorrectSmiles('()')
            elif token:
                tokens.append((token_type, token))
                token = None
            token_type = 3
            tokens.append((3, None))
        # simple atoms
        elif s in 'NOPSFIA':  # organic atoms and the A any-aliphatic wildcard
            if token:
                tokens.append((token_type, token))
                token = None
            token_type = 0
            tokens.append((0, s))
        elif s in 'cnopsba':  # aromatic ring atom (a = any aromatic)
            if token:
                tokens.append((token_type, token))
                token = None
            token_type = 8
            tokens.append((8, s.upper()))
        elif s in 'CB':  # flag possible Cl or Br
            if token:
                tokens.append((token_type, token))
            token_type = 0
            token = s
        elif s == '*':  # any atom wildcard
            if token:
                tokens.append((token_type, token))
                token = None
            token_type = 0
            tokens.append((0, 'A'))
        elif token_type == 0:
            if s == 'l':
                if token == 'C':
                    tokens.append((0, 'Cl'))
                    token = None
                else:
                    raise IncorrectSmiles('invalid element Bl')
            elif s == 'r':
                if token == 'B':
                    tokens.append((0, 'Br'))
                    token = None
                else:
                    raise IncorrectSmiles('invalid smiles for Cr')
            else:
                raise IncorrectSmiles('invalid smiles')
        else:
            raise IncorrectSmiles('invalid smiles')

    if token_type == 5:
        raise IncorrectSmiles('atom description has not finished')
    elif token_type == 7:
        if token:
            tokens.append((6, int(token[0])))
        else:
            raise IncorrectSmiles('invalid %closure')
    elif token:
        tokens.append((token_type, token))  # C or B
    return tokens


def _atom_parse(token):
    # [isotope]Element[element][@[@]][H[n]][+-charge][_synthon][:mapping]
    _match = atom_re.fullmatch(token)
    if _match is None:
        raise IncorrectSmiles(f'atom token invalid {token}')
    isotope, element, stereo, hydrogen, charge, synthon_label, mapping = _match.groups()

    if isotope:
        isotope = int(isotope)

    if stereo:
        stereo = stereo == '@'

    if hydrogen:
        if len(hydrogen) > 1:
            hydrogen = int(hydrogen[1:])
        else:
            hydrogen = 1
    else:
        hydrogen = 0

    if charge:
        try:
            charge = charge_dict[charge]
        except KeyError:
            raise IncorrectSmiles('charge token invalid')
    else:
        charge = 0

    if mapping:
        try:
            mapping = int(mapping[1:])
        except ValueError:
            raise IncorrectSmiles('invalid mapping token')

    if element in ('c', 'n', 'o', 'p', 's', 'as', 'se', 'b', 'te'):
        _type = 8
        element = element.capitalize()
    else:
        _type = 0
    out = {'element': element, 'isotope': isotope, 'parsed_mapping': mapping, 'charge': charge,
           'implicit_hydrogens': hydrogen, 'stereo': stereo}
    if synthon_label:  # ONLY when present: Element.__init__ has no such keyword
        out['synthon_label'] = synthon_label[1:]  # the token WITHOUT its sigil, a str never an int
    return _type, out


def _split_unnested(token, separator):
    """Split token on separator occurrences outside parentheses and brackets."""
    parts = []
    current = []
    depth = 0
    for char in token:
        if char in '([':
            depth += 1
        elif char in ')]':
            depth -= 1
        if char == separator and depth == 0:
            parts.append(''.join(current))
            current = []
        else:
            current.append(char)
    parts.append(''.join(current))
    return parts


def _pull_recursive(token, group):
    """Pull $() and !$() out of one alternative, tagging them with its OR group."""
    recursive = []
    result = []
    i = 0
    n = len(token)
    while i < n:
        if i < n - 1 and token[i] == '$' and token[i + 1] == '(':
            # positive recursive $()
            i += 2  # skip $(
            inner, i = _balanced_extract(token, i, n)
            recursive.append((True, inner, group))
        elif i < n - 2 and token[i] == '!' and token[i + 1] == '$' and token[i + 2] == '(':
            # negated recursive !$()
            i += 3  # skip !$(
            inner, i = _balanced_extract(token, i, n)
            recursive.append((False, inner, group))
        else:
            result.append(token[i])
            i += 1
    return ''.join(result), recursive


def _extract_recursive(token):
    """Extract $() and !$() recursive SMARTS from bracket content.

    ``;`` is the lowest-precedence AND in SMARTS and ``,`` is OR, so
    ``[O;$(a),$(b)]`` means "O and (a or b)". Each recursive constraint is
    tagged with the index of the ``;``-term it came from: constraints sharing a
    tag are alternatives to be ORed, and different tags are ANDed. Without the
    tag every constraint was ANDed, so an alternation like the one above could
    never match.

    Returns (cleaned_token, [(positive, inner_smarts, group), ...]).
    If no recursive found, returns (token, None).
    """
    recursive = []
    kept = []
    for group, term in enumerate(_split_unnested(token, ';')):
        alternatives = []
        for alternative in _split_unnested(term, ','):
            leftover, found = _pull_recursive(alternative, group)
            recursive.extend(found)
            if leftover:
                alternatives.append(leftover)
        if alternatives:
            kept.append(','.join(alternatives))
    if not recursive:
        return token, None
    # clean trailing/leading separators from result
    cleaned = ';'.join(kept)
    # remove trailing semicolons and commas left from extraction
    while cleaned and cleaned[-1] in ';,':
        cleaned = cleaned[:-1]
    # remove leading semicolons and commas
    while cleaned and cleaned[0] in ';,':
        cleaned = cleaned[1:]
    # remove consecutive semicolons
    while ';;' in cleaned:
        cleaned = cleaned.replace(';;', ';')
    # remove consecutive commas
    while ',,' in cleaned:
        cleaned = cleaned.replace(',,', ',')
    # remove semicolon-comma or comma-semicolon combos
    cleaned = cleaned.replace(';,', ';').replace(',;', ';')
    # final trim
    while cleaned and cleaned[-1] in ';,':
        cleaned = cleaned[:-1]
    while cleaned and cleaned[0] in ';,':
        cleaned = cleaned[1:]
    return cleaned, recursive


def _balanced_extract(token, start, n):
    """Extract content from balanced parentheses starting at position start.

    Returns (inner_content, position_after_closing_paren).
    """
    depth = 1
    i = start
    bracket_depth = 0
    while i < n and depth > 0:
        c = token[i]
        if c == '[':
            bracket_depth += 1
        elif c == ']':
            bracket_depth -= 1
        elif bracket_depth == 0:
            if c == '(':
                depth += 1
            elif c == ')':
                depth -= 1
                if depth == 0:
                    return token[start:i], i + 1
        i += 1
    raise IncorrectSmarts('Unbalanced parentheses in recursive SMARTS')


_elements = frozenset(elements)
_two_letter_elements = frozenset(e for e in _elements if len(e) == 2) | {'as', 'se', 'te'}
_extra_symbols = frozenset(('*', 'A', 'a', 'M', 'c', 'n', 'o', 'p', 's', 'b'))
# Daylight defaults for a primitive written without a count. Bare <r> means
# "in a ring of any size", which chython spells <R>.
_bare_primitives = {'D': 'D1', 'X': 'X1', 'H': 'H1', 'h': 'h1', 'r': 'R'}
# domain of each counting primitive: a negated one is the complement over it.
# implicit hydrogens stop at 4 - no atom carries more and the packed matcher has no room.
_count_domains = {'D': ('neighbors', frozenset(range(15))), 'X': ('total_connectivity', frozenset(range(15))),
                  'H': ('implicit_hydrogens', frozenset(range(5))), 'h': ('implicit_hydrogens', frozenset(range(5))),
                  'x': ('ring_connectivity', frozenset(range(15))), 'y': ('heteroatoms', frozenset(range(15))),
                  'R': ('rings_count', frozenset(range(15))),
                  'z': ('hybridization', frozenset(range(1, 5))), 'v': ('valence', frozenset(range(9)))}


def _expand_bare(primitive):
    if primitive.startswith('!'):
        return '!' + _bare_primitives.get(primitive[1:], primitive[1:])
    return _bare_primitives.get(primitive, primitive)


def _lex_primitives(text):
    """Split one AND-group of a bracket body into primitive strings.

    Juxtaposition and <&> are the same high-precedence AND in Daylight, so
    [CX3], [C&X3] and [X3] all have to lex to the same primitives.
    """
    out = []
    i = 0
    n = len(text)
    while i < n:
        if text[i] == '&':
            i += 1
            continue
        start = i
        if text[i] == '!':
            i += 1
            if i == n:
                raise IncorrectSmarts(f'dangling ! in [{text}]')
        c = text[i]
        if c.isdigit():  # isotope
            while i < n and text[i].isdigit():
                i += 1
        elif c in '+-':  # charge
            i += 1
            while i < n and (text[i] == c or text[i].isdigit()):
                i += 1
        elif c == '@':  # chirality
            i += 1
            if i < n and text[i] in '@?':
                i += 1
        elif c in '#:':  # atomic number, atom mapping
            i += 1
            while i < n and text[i].isdigit():
                i += 1
        elif text[i:i + 2] in _two_letter_elements:  # before D/H/R/X: [Dy] [He] [Rh] [Xe]
            i += 2
        elif c in 'DXRHrhxyzva':
            i += 1
            while i < n and text[i].isdigit():
                i += 1
        elif c == '_':  # synthon label. consumed whole by regex, or its tail re-lexes as a primitive
            m = synthon_token_re.match(text, i)
            if m is None:
                raise IncorrectSmarts(f'invalid synthon token in [{text}]')
            i = m.end()
        elif c in _elements or c in _extra_symbols:
            i += 1
        else:
            raise IncorrectSmarts(f'invalid primitive in [{text}]')
        out.append(text[start:i])
    return out


def _is_element_term(term):
    """Whether every alternative of a term is a positive element symbol."""
    for p in term.split(','):
        if p in _extra_symbols or p in _elements or p in _two_letter_elements:
            continue
        elif p.startswith('#') and p[1:].isdigit():
            continue
        elif p.startswith('a') and p[1:].isdigit():  # aromatic with a neighbors count
            continue
        return False
    return True


def _split_primitives(token):
    """Split a bracket body into chython's ';'-separated primitive terms.

    Daylight puts no constraints on primitive order and does not require an
    element at all, while the semantics below expect the element term first, so
    it is hoisted here and defaults to any atom. <,> alternatives of different
    primitive types have no representation and are rejected.
    """
    terms = []
    for chunk in _split_unnested(token, ';'):
        groups = [g for g in (_lex_primitives(x) for x in _split_unnested(chunk, ',')) if g]
        if not groups:
            continue
        elif len(groups) == 1:
            terms.extend(groups[0])
        elif all(len(g) == 1 for g in groups):
            terms.append(','.join(g[0] for g in groups))
        else:
            raise IncorrectSmarts(f'Unsupported OR statement: {chunk}')
    terms = [t for t in terms if t[0] not in '@:']  # stereo and mapping are parsed separately
    element = None
    for i, t in enumerate(terms):
        # a bare H opening the bracket is the hydrogen atom, further in it is a count: [CH], [nH]
        if _is_element_term(t) and (not i or t != 'H'):
            element = terms.pop(i)
            break
    terms = [','.join(_expand_bare(x) for x in t.split(',')) for t in terms]
    # a bracket with no element is any atom, so '*' — not 'A', which is any *aliphatic* atom
    return [element or '*'] + terms


def _query_parse(token):
    out = {}
    # Extract recursive SMARTS $() and !$() BEFORE any regex processing
    token, recursive = _extract_recursive(token)
    if recursive:
        out['recursive_smarts'] = recursive
    if isotope := match(iso_re, token):
        token = token[isotope.end():]  # remove isotope substring
        out['isotope'] = int(isotope.group())
    # Handle negated charges (!+ / !-) before regular charge
    charge_nots = set()
    while nc := search(not_charge_re, token):
        symbol = token[nc.start() + 1]
        token = token[:nc.start()] + token[nc.end():]
        charge_nots.add(symbol)
    if charge_nots:
        if len(charge_nots) == 2:  # both !+ and !- → must be zero
            out['charge'] = 0
        elif '+' in charge_nots:
            out['charge_not'] = 'positive'
        else:
            out['charge_not'] = 'negative'
    elif charge := search(chg_re, token):
        token = token[:charge.start()] + token[charge.end():]  # remove charge substring
        out['charge'] = charge_dict[charge.group()]

    if mapping := search(mpp_re, token):
        token = token[:mapping.start()]
        out['parsed_mapping'] = int(mapping.group()[1:])

    if stereo := search(str_re, token):  # drop stereo mark. unsupported
        token = token[:stereo.start()] + token[stereo.end():]
        out['stereo'] = stereo.group() == '@'

    primitives = _split_primitives(token)
    element = [int(x[1:]) if x.startswith('#') else x for x in primitives[0].split(',')]
    if len(element) == 1:
        element = element[0]

    # Handle lowercase aromatic element symbols (c, n, o, p, s, b, se, te)
    # Also handle * (any atom) and a (any aromatic atom)
    # Handle negated elements: !#6, !C, etc.
    aromatic_from_symbol = False
    aliphatic_from_symbol = False
    if isinstance(element, str):
        if element.startswith('!'):
            # Negated element: !#6 = not carbon, !N = not nitrogen
            negated = element[1:]
            if negated.startswith('#'):
                out['excluded_elements'] = (int(negated[1:]),)
            else:
                upper = _aromatic_upper.get(negated)
                if upper:
                    negated = upper
                out['excluded_elements'] = (negated,)
            element = 'A'
        elif element == '*':
            element = 'A'
        elif element == 'A':  # Daylight: A is any aliphatic atom, * is any atom
            aliphatic_from_symbol = True
        elif element == 'a':
            element = 'A'
            aromatic_from_symbol = True
        elif len(element) >= 2 and element[0] == 'a' and element[1:].isdigit():
            # a3 = any aromatic atom with N neighbors
            out['neighbors'] = [int(element[1:])]
            element = 'A'
            aromatic_from_symbol = True
        else:
            upper = _aromatic_upper.get(element)
            if upper:
                element = upper
                aromatic_from_symbol = True
    elif isinstance(element, list):
        new_elements = []
        has_any = False
        for e in element:
            if isinstance(e, str):
                if e in ('a', 'A', '*'):
                    has_any = True
                    if e == 'a':
                        aromatic_from_symbol = True
                    continue
                upper = _aromatic_upper.get(e)
                if upper:
                    new_elements.append(upper)
                    aromatic_from_symbol = True
                else:
                    new_elements.append(e)
            else:
                new_elements.append(e)
        if has_any and not new_elements:
            element = 'A'
        elif has_any:
            element = 'A'  # any + specifics = any
        elif new_elements:
            element = new_elements
        else:
            element = 'A'

    out['element'] = element

    for p in primitives[1:]:  # parse hydrogens (h), neighbors (D), rings_sizes (r or !R), hybridization == 4 (a)
        if not p:
            continue
        elif p == 'a':  # aromatic atom
            out['hybridization'] = 4
        elif p == 'A':  # aliphatic mark
            out['hybridization'] = [1, 2, 3]
        elif p == '!a':  # not aromatic = aliphatic
            out['hybridization'] = [1, 2, 3]
        elif p == '!R':
            out['ring_sizes'] = 0
        elif p == 'R':  # bare R: in any ring (ring count >= 1)
            out['rings_count'] = list(range(1, 15))
        elif p == 'M':
            out['masked'] = True
        else:
            p = [x for x in p.split(',') if x]
            if not p:
                continue
            negated = p[0][0] == '!'
            t = p[0][1] if negated else p[0][0]
            if t in _count_domains or t == 'r':
                # y and z are chython marks for heteroatom count and hybridization; x is Daylight's ring bonds
                if len({x.lstrip('!')[0] for x in p}) > 1 or any((x[0] == '!') != negated for x in p):
                    raise IncorrectSmarts('Unsupported OR statement')
                try:
                    values = [int(x[2:] if negated else x[1:]) for x in p]
                except ValueError:
                    raise IncorrectSmarts(f'Invalid primitive value')

                if t == 'r':
                    if not negated:  # r0 is the not-in-ring mark, which ring_sizes spells as a bare zero
                        out['ring_sizes'] = 0 if values == [0] else values
                        continue
                    elif values != [0]:  # !r3;!r4;... - a blacklist beside the ring_sizes whitelist
                        excluded = set(out.get('excluded_ring_sizes', ()))
                        excluded.update(values)
                        out['excluded_ring_sizes'] = sorted(excluded)
                        continue
                    key, allowed = 'rings_count', set(range(1, 15))  # !r0 - member of any ring
                else:
                    key, domain = _count_domains[t]
                    allowed = set()
                    for v in values:
                        allowed |= domain - {v} if negated else {v}
                if key in out:  # primitives of the same kind are ANDed: [C!H0!H1]
                    known = out[key]
                    allowed &= {known} if isinstance(known, int) else set(known)
                    if not allowed:
                        raise IncorrectSmarts(f'contradictory {key} constraints')
                out[key] = sorted(allowed)
            elif all(_is_element_token(x) for x in p):
                # Element constraint in primitive position (e.g., [*;O,S,P,N] or [*;!#6])
                extra_elements = []
                extra_excluded = []
                for item in p:
                    negated = item.startswith('!')
                    raw = item[1:] if negated else item
                    if raw.startswith('#'):
                        val = int(raw[1:])
                    elif raw in ('*', 'A'):
                        continue  # any atom — no additional constraint
                    elif raw == 'a':
                        if not negated:
                            out['hybridization'] = 4
                        continue
                    else:
                        upper = _aromatic_upper.get(raw)
                        val = upper if upper else raw
                        if upper and not negated:
                            aromatic_from_symbol = True
                    if negated:
                        extra_excluded.append(val)
                    else:
                        extra_elements.append(val)
                if extra_elements:
                    current = out.get('element', 'A')
                    if current == 'A':
                        out['element'] = extra_elements if len(extra_elements) > 1 else extra_elements[0]
                    elif isinstance(current, list):
                        intersected = [e for e in current if e in extra_elements]
                        if intersected:
                            out['element'] = intersected if len(intersected) > 1 else intersected[0]
                if extra_excluded:
                    existing = list(out.get('excluded_elements', ()))
                    existing.extend(extra_excluded)
                    out['excluded_elements'] = tuple(existing)
            elif len(p) == 1 and p[0][:1] == '_' and p[0][1:] in SYNTHON_LABELS:  # synthon label
                out['synthon_label'] = p[0][1:]  # token minus its sigil; a membership test
            else:
                raise IncorrectSmarts(f'Unsupported SMARTS primitive: {p[0]}')

    # Infer aromatic hybridization from lowercase symbol if not explicitly set
    if aromatic_from_symbol and 'hybridization' not in out:
        out['hybridization'] = 4
    elif aliphatic_from_symbol and 'hybridization' not in out:
        out['hybridization'] = [1, 2, 3]
    elif not aromatic_from_symbol and 'hybridization' not in out:
        # Standard SMARTS: uppercase aromaticity-capable symbol = aliphatic.
        # Marker is relaxed by the parser if the atom has an aromatic bond.
        elem = out.get('element')
        aromatic_capable = set(_aromatic_upper.values())
        if isinstance(elem, str):
            if elem in aromatic_capable:
                out['hybridization'] = [1, 2, 3]
                out['_default_aliphatic'] = True
        elif isinstance(elem, list):
            if elem and all(isinstance(e, str) and e in aromatic_capable for e in elem):
                out['hybridization'] = [1, 2, 3]
                out['_default_aliphatic'] = True

    # Return type 8 for aromatic atoms so the parser assigns aromatic bonds
    _type = 8 if out.get('hybridization') == 4 else 0
    return _type, out


def _dyn_atom_parse(token):
    # [isotope]Element[element][+-charge[>+-charge]][*^[>*^]]
    _match = dyn_atom_re.fullmatch(token)
    if _match is None:
        raise IncorrectSmiles(f'atom token invalid {token}')
    isotope, element, charge, _, is_radical, _ = _match.groups()

    if isotope:
        isotope = int(isotope)

    if charge:
        try:
            charge, *cgr = dyn_charge_dict[charge]
        except KeyError:
            raise IncorrectSmiles('charge token invalid')
    else:
        charge = 0
        cgr = []

    if is_radical:
        try:
            is_radical, *dyn = dyn_radical_dict[is_radical]
        except KeyError:
            raise IncorrectSmiles('invalid dynamic radical token')
        else:
            cgr.extend(dyn)
    else:
        is_radical = False

    if element in ('c', 'n', 'o', 'p', 's', 'as', 'se', 'b'):
        _type = 12
        element = element.capitalize()
    else:
        _type = 11
    return _type, {'element': element, 'isotope': isotope, 'is_radical': is_radical,
                   'parsed_mapping': 0, 'cgr': cgr, 'charge': charge}


def smiles_tokenize(smi):
    tokens = _tokenize(smi)
    out = []
    for token_type, token in tokens:
        if token_type in (0, 8):  # simple atom
            out.append((token_type, {'element': token}))
        elif token_type == 5:
            out.append(_atom_parse(token))
        elif token_type == 10:
            raise IncorrectSmiles('SMARTS detected')
        else:
            out.append((token_type, token))
    return out


def smarts_tokenize(smi):
    tokens = _tokenize(smi)
    out = []
    for token_type, token in tokens:
        if token_type == 0:  # simple non-aromatic atom
            out.append((0, {'element': token}))
        elif token_type == 8:  # simple aromatic atom (already uppercased by _tokenize)
            out.append((8, {'element': token, 'hybridization': 4}))
        elif token_type == 5:
            out.append(_query_parse(token))
        else:
            out.append((token_type, token))
    return out


__all__ = ['smiles_tokenize', 'smarts_tokenize']
