# Chython Code Review Report

**Date:** 2026-02-06
**Reviewer:** Claude (Opus 4.5)
**Branch:** `claude/review-code-style-0qpNj`
**Focus:** Code style alignment, system design consistency, interface coherence

---

## Table of Contents

1. [CGR Container Issues](#1-cgr-container-issues)
2. [Depict Module & Visualization Issues](#2-depict-module--visualization-issues)
3. [Periodictable Module Issues](#3-periodictable-module-issues)
4. [SMARTS Module Impact Analysis](#4-smarts-module-impact-analysis)
5. [Graph/Container Base Issues](#5-graphcontainer-base-issues)
6. [Summary Tables](#6-summary-tables)

---

## 1. CGR Container Issues

### Issue 1: Docstring Style Inconsistencies

**Location:** `cgr.py:621-622` (`center_bonds`), `cgr.py:612-615` (`center_atoms`)

**Current code:**
```python
@cached_property
def center_bonds(self) -> Tuple[Tuple[int, int], ...]:
    """Get list of bonds of reaction center (bonds with dynamic orders or touching center atoms)."""
```

```python
@cached_property
def center_atoms(self) -> Tuple[int, ...]:
    """ Get list of atoms of reaction center (atoms with dynamic: bonds, charges, radicals).
    """
```

**Problem:** The existing codebase uses `:return:` docstring format consistently. Compare with `molecule.py:244-250`:
```python
@cached_property
def aromatic_rings(self) -> Tuple[Tuple[int, ...], ...]:
    """
    Aromatic rings atoms numbers
    """
```

**Recommendation:** Use consistent docstring format with `:return:` when describing return values:
```python
@cached_property
def center_bonds(self) -> Tuple[Tuple[int, int], ...]:
    """
    Bonds of reaction center (dynamic bonds or bonds touching center atoms).

    :return: tuple of atom pairs representing center bonds
    """
```

---

### Issue 2: `aromatic_rings` Property Duplication/Inconsistency

**Location:** `cgr.py:691-699` vs `molecule.py:244-250`

**In `cgr.py`:**
```python
@cached_property
def aromatic_rings(self) -> Tuple[Tuple[int, ...], ...]:
    """Existing or formed aromatic rings atoms numbers."""
    adj = self._bonds
    return tuple(
        ring for ring in self.sssr
        if (adj[ring[0]][ring[-1]].order == 4 and all(adj[n][m].order == 4 for n, m in zip(ring, ring[1:])))
        or (adj[ring[0]][ring[-1]].p_order == 4 and all(adj[n][m].p_order == 4 for n, m in zip(ring, ring[1:])))
    )
```

**In `molecule.py`:**
```python
@cached_property
def aromatic_rings(self) -> Tuple[Tuple[int, ...], ...]:
    """
    Aromatic rings atoms numbers
    """
    bonds = self._bonds
    return tuple(ring for ring in self.sssr if bonds[ring[0]][ring[-1]] == 4
                 and all(bonds[n][m] == 4 for n, m in zip(ring, ring[1:])))
```

**Problems:**
1. Variable naming inconsistency: `adj` vs `bonds` for `self._bonds`
2. The CGR version correctly checks both reactant (`order`) and product (`p_order`) aromatic states, but the logic differs from how other CGR-specific properties work
3. Docstring style differs (single-line vs multi-line)

**Recommendation:** Use `bonds` consistently as the variable name since that's the established pattern.

---

### Issue 3: `centers_list` Method - Variable Shadowing in Loop

**Location:** `cgr.py:655-668`

**Code:**
```python
for r in self.aromatic_rings:
    if not center.isdisjoint(r):
        n = r[0]
        m = r[-1]
        if n in adj and m in adj[n]:
            for n, m in zip(r, r[1:]):  # <-- n and m are shadowed here!
```

**Problem:** The variables `n` and `m` are defined at line 655-656 and then immediately overwritten in the loop at line 658. This is confusing and could lead to bugs.

**Recommendation:** Rename loop variables to avoid shadowing.

---

### Issue 4: `centers_list` Logic Duplication with `center_atoms`

**Location:** `cgr.py:640-648`

The existing `center_atoms` (line 612-618) uses a clean, readable pattern:
```python
@cached_property
def center_atoms(self) -> Tuple[int, ...]:
    center = {n for n, a in self._atoms.items() if a.is_dynamic}
    center.update(n for n, m_bond in self._bonds.items() if any(bond.is_dynamic for bond in m_bond.values()))
    return tuple(center)
```

But `centers_list` reimplements similar logic differently (lines 640-648):
```python
for n, c in charges.items():
    if c != p_charges.get(n, c) or radicals.get(n, False) != p_radicals.get(n, False):
        center.add(n)

for n, m_bond in self._bonds.items():
    for m, bond in m_bond.items():
        if bond.is_dynamic:
            adj[n].add(m)
center.update(adj)
```

**Recommendation:** Reuse `center_atoms` and check `atom.is_dynamic` to avoid logic duplication.

---

### Issue 5: `center_bonds` Definition Semantic Change

**Location:** `cgr.py:620-628`

**Current implementation:**
```python
@cached_property
def center_bonds(self) -> Tuple[Tuple[int, int], ...]:
    """Get list of bonds of reaction center (bonds with dynamic orders or touching center atoms)."""
    center_atoms = set(self.center_atoms)
    bonds = []
    for n, m, bond in self.bonds():
        if bond.is_dynamic or n in center_atoms or m in center_atoms:
            bonds.append((n, m))
    return tuple(bonds)
```

**Issue:** The definition includes **all bonds touching center atoms**, not just dynamic bonds. This is a significant semantic expansion from what the original `1a8809a` commit intended:
```python
# Original (commit 1a8809a):
return tuple((n, m) for n, m, bond in self.bonds() if bond.is_dynamic)
```

This change may be intentional, but it's a behavioral change that:
1. Makes `center_bonds` return different data than expected based on the property name
2. Creates potential confusion - "center bonds" intuitively means bonds that are changing, not all bonds adjacent to changing atoms

**Recommendation:** Either:
- Rename to something like `center_region_bonds` or add a separate property
- Or revert to original semantics and add a new property for the expanded behavior

---

### Issue 6: `ReactionContainer.from_cgr` Method Pattern

**Location:** `reaction.py:155-175`

```python
@classmethod
def from_cgr(cls, cgr: 'CGRContainer') -> 'ReactionContainer':
    """
    Decompose CGR into reaction
    ...
    """
    if not isinstance(cgr, CGRContainer):
        raise TypeError('CGR expected')
    r, p = cgr.decompose()
    reaction = object.__new__(cls)
    reaction._reactants = tuple(r.split())
    reaction._products = tuple(p.split())
    ...
```

**Issues:**
1. Good: Uses `object.__new__` pattern consistent with other `copy()` methods
2. Good: Type checking is present
3. **Missing:** No conformers/plane data transfer consideration
4. **Style:** The method bypasses `__init__` validation - should document this

**Recommendation:** Consider whether metadata from CGR should transfer to the reaction.

---

### Issue 7: Import Organization in `cgr.py`

**Location:** `cgr.py:19-36`

**Current:**
```python
from functools import cached_property
from typing import Dict, Iterator, Tuple, Optional, List, Union
from collections import defaultdict

from .bonds import DynamicBond, Bond
from . import cgr_query as query
...
```

**Issue:** The codebase pattern (seen in `molecule.py:19-45`) is:
1. Standard library imports
2. Third-party imports
3. Local imports

The `defaultdict` import was added but is placed after the other standard library imports. While minor, this breaks the alphabetical ordering within the stdlib group.

---

### Issue 8: Test Coverage Gaps

**Location:** `test_cgr.py:277-284`

```python
def test_centers_list_and_center_bonds(self):
    cgr = create_simple_cgr()
    # whole structure is one center because both bonds are dynamic
    centers = cgr.centers_list
    self.assertEqual(len(centers), 1)
    self.assertEqual(set(centers[0]), {1, 2, 3})
```

**Missing tests:**
1. No test for multiple disconnected centers
2. No test for aromatic ring propagation in `centers_list`
3. No test for edge cases (empty CGR, single atom, etc.)

---

### Issue 9: Hydrogens Fix Pattern

**Location:** `calculate2d/__init__.py`

**Changed from:**
```python
mol._hydrogens = {n: 0 for n in mol._hydrogens}
```

**To (`06a23be`):**
```python
for n in mol._atoms:
    mol._atoms[n]._implicit_hydrogens = 0
```

**Assessment:** This is correct - the code was accessing a non-existent `_hydrogens` dict. The fix properly uses the atom's `_implicit_hydrogens` attribute. However, this pattern of directly setting private attributes could be encapsulated better.

---

## 2. Depict Module & Visualization Issues

### Issue 10: Critical - Dual Configuration Systems

**Location:** `depict.py:54-60` vs `depict.py:595-630`

The most significant issue is that there are **two separate configuration systems** that don't communicate:

**Module-level config (`depict.py:54-60`):**
```python
_render_config = {'carbon': False, 'dashes': (.2, .1), 'span_dy': .15, 'mapping': True, ...}
```
Used by: `DepictMolecule`, `DepictReaction`, `depict_settings()` function

**Class-level config (`depict.py:595-630`):**
```python
class Depict:
    _render_config = {
        'carbon': False,
        'atoms_colors': cpk,
        ...
        'cgr_aromatic_space': .18,  # CGR-specific
        'broken_color': 'red',
        'formed_color': 'green',
    }
```
Used by: `Depict`, `DepictCGR`, `DepictQuery`, `DepictQueryCGR`

**Problems:**
1. **Two different `depict_settings`**: The module-level function (`depict.py:166`) and class method (`Depict.depict_settings` at line 633) have different signatures and behaviors
2. **Configuration keys differ**: Class-level has CGR-specific keys (`broken_color`, `formed_color`, `cgr_aromatic_space`) that module-level lacks
3. **User confusion**: Calling `depict_settings()` affects Molecules/Reactions but not CGR/Query containers

---

### Issue 11: `depict()` Method Signature Inconsistency

**`DepictMolecule.depict()` (`depict.py:240-294`):**
```python
def depict(self: Union['MoleculeContainer', 'DepictMolecule'], *, width=None, height=None, clean2d: bool = True,
           format: Literal['svg', 'png', 'svgz'] = 'svg', png_width=1000, png_heigh=1000, png_scale=1.,
           _embedding=False) -> Union[str, bytes]:
```

**`Depict.depict()` (for CGR/Query, `depict.py:710-735`):**
```python
def depict(self, *, embedding=False):
    # No width, height, clean2d, format, png_* parameters!
```

**Problems:**
1. CGR containers **cannot** specify output dimensions or format
2. CGR containers **cannot** export PNG
3. Parameter name inconsistency: `_embedding` vs `embedding`
4. No `clean2d` option for CGR (though `Calculate2DCGR` exists)

---

### Issue 12: Typo in Parameter Name

**Location:** `depict.py:241`

```python
png_width=1000, png_heigh=1000  # typo: "heigh" should be "height"
```

This typo is in the public API and should be fixed with backwards compatibility.

---

### Issue 13: Config Key Naming Inconsistency (`cgr_aromatic_space`)

**In class `_render_config` (`depict.py:622`):**
```python
'cgr_aromatic_space': .18,  # New
```

But in the fix (`depict.py:1065`):
```python
aromatic_space = config['aromatic_space']  # Changed FROM 'cgr_aromatic_space'
```

This key `cgr_aromatic_space` is defined but never used after the fix. It's dead configuration.

---

### Issue 14: `x3dom.py` Uses `_render_config` Without Defining It

**Location:** `x3dom.py:102-108`

```python
def __render_atoms(self, xyz):
    config = self._render_config  # Where does this come from?

    font = config['font_size']
    carbon = config['carbon']
    radius = config['atom_radius']
    colors = config['atoms_colors']
    mapping_color = config['mapping_color']
```

The `X3dom` class expects `self._render_config` but doesn't inherit from `Depict`. It relies on the class it's mixed into having this attribute.

**Class hierarchy shows:**
- `MoleculeContainer` → uses `DepictMolecule` which uses **module-level** `_render_config`
- `CGRContainer` → uses `DepictCGR(Depict)` which has **class-level** `_render_config`

So `X3domMolecule` works with molecules because it accesses module-level config via `self`, but `X3domCGR` must rely on `Depict._render_config`.

---

### Issue 15: Missing `_render_config` Attribute in `DepictMolecule`

**Location:** `depict.py:237-239`

```python
class DepictMolecule:
    __slots__ = ()
    # No _render_config attribute!
```

The class directly accesses the **module-level** `_render_config` variable:
```python
def __render_bonds(self: Union['MoleculeContainer', 'DepictMolecule']):
    ...
    double_space = _render_config['double_space']  # Module variable, not self._render_config
```

But `X3dom.__render_atoms()` uses `self._render_config`. This works only because Python's name resolution finds the module-level variable, but it's confusing and inconsistent.

---

### Issue 16: `rotate_vector` Export Pattern

**Location:** `depict.py:95-105`

```python
def _rotate_vector(x1, y1, x2, y2):
    """
    rotate x,y vector over x2-x1, y2-y1 angle
    """
    ...

rotate_vector = _rotate_vector  # Public alias to private function
```

**Problems:**
1. Not in `__all__` (line 2174), so intent is unclear
2. Used internally by `Depict` class methods via `rotate_vector` (not `_rotate_vector`)
3. If it's meant to be public, it should be properly documented and exported

---

### Issue 17: Inconsistent Mask Data Structure

**`DepictMolecule.__render_atoms()` returns (`depict.py:361-492`):**
```python
return svg, defines, masks  # masks is a list
```

**`Depict._render_atoms()` returns (`depict.py:1094-1197`):**
```python
return svg, mask  # mask is a defaultdict(list) with keys like 'center', 'symbols', 'aam', 'other'
```

These are incompatible structures, which is why `Depict` has its own `_graph_svg()` method that handles its mask format.

---

### Issue 18: `_render_3d_bonds` Method Missing in Some Classes

**Location:** `x3dom.py:88`

```python
bonds = self._render_3d_bonds(xyz)
```

This method must be implemented by the inheriting class, but there's no abstract declaration or protocol. The `X3dom` base class just assumes it exists.

---

### Issue 19: Hardcoded Values in Depict Settings Logic

**Location:** `Depict.depict_settings()` (`depict.py:660-698`)

```python
if span_dy is not None: config['span_dy'] = span_dy
else: config['span_dy'] = 0.3 * effective_font_size  # Hardcoded factor

if span_size is not None: config['span_size'] = span_size
else: config['span_size'] = 0.7 * effective_font_size  # Hardcoded factor
```

This recalculates values based on font_size even when the user didn't change them, which could cause unexpected side effects.

---

### Issue 20: File Header Inconsistency

**Location:** `x3dom.py:1-18`

```python
#  Copyright 2020 Ramil Nugmanov <nougmanoff@protonmail.com>
#  This file is part of CGRtools.  # <-- Wrong project name!
#
#  CGRtools is free software;  # <-- Should be "chython"
```

The file header still references "CGRtools" instead of "chython".

---

## 3. Periodictable Module Issues

### Issue 21: Critical - Missing File Header in `dynamic_query.py`

**Location:** `dynamic_query.py:1-6`

```python
from typing import Tuple, Dict, Type, Union, Optional
from .dynamic import Dynamic, DynamicElement
from .element import Element
...
```

**Problem:** This file is missing the standard copyright/license header that all other files have. Compare with `query.py:1-18` which has the full LGPL header.

---

### Issue 22: Critical - `DynamicQueryElement.__init__` Calls Wrong Parent

**Location:** `dynamic_query.py:85-88`

```python
def __init__(self, atomic_number: int, isotope: Optional[int], **kwargs):
    super().__init__()  # Calls DynamicQuery.__init__() which doesn't exist!
    self._atomic_number = atomic_number
    self._isotope = isotope
```

**Problem:** `DynamicQuery` inherits from `Dynamic` which inherits from `Core`. Looking at `Core`:
- `Core` has no `__init__` method
- `Dynamic` has no `__init__` method

So `super().__init__()` resolves to `object.__init__()`, which works but is confusing. The `**kwargs` are also never used.

**Recommendation:** Either remove `super().__init__()` or properly initialize parent state.

---

### Issue 23: `DynamicAnyElement.__init__` Ignores `isotope` Parameter

**Location:** `dynamic_query.py:192-193`

```python
def __init__(self, isotope: Optional[int] = None, **kwargs):
    super().__init__()  # isotope parameter is completely ignored!
```

**Problem:** The `isotope` parameter is declared but never stored or used. This is misleading API design.

---

### Issue 24: `copy()` Method Inconsistencies Across Atom Types

| Class | `copy()` Signature | Behavior |
|-------|-------------------|----------|
| `Element` | `copy(full=False, hydrogens=False, stereo=False)` | Copies with options |
| `Query` | `copy(full=False)` | Copies masked if `full=True` |
| `QueryElement` | `copy(full=False)` | Copies isotope |
| `DynamicElement` | `copy(*_, **__)` | **Ignores all parameters!** |
| `DynamicQueryElement` | `copy(*args, **kwargs_copy_call)` | **Ignores all parameters!** |

**`dynamic.py:205-212`:**
```python
def copy(self, *_, **__):  # All args/kwargs silently ignored
    copy = object.__new__(self.__class__)
    copy._isotope = self.isotope
    copy._charge = self.charge
    copy._is_radical = self.is_radical
    copy._p_is_radical = self.p_is_radical
    copy._p_charge = self.p_charge
    return copy
```

**Problems:**
1. Passing `full=True` to `DynamicElement.copy()` does nothing - silent failure
2. Users might expect consistent behavior across atom types
3. `DynamicQueryElement.copy()` creates unused `new_kwargs` dict

---

### Issue 25: `from_atom()` Factory Method Inconsistencies

| Class | `from_atom()` Behavior |
|-------|----------------------|
| `Element` | No `from_atom()` method |
| `QueryElement.from_atom()` | Returns `QueryElement`, takes optional constraint flags |
| `DynamicElement.from_atom()` | Returns `DynamicElement` with same charge/radical |
| `DynamicQueryElement.from_atom()` | Has dead code branches, may return copy of input |

**`dynamic_query.py:155-170`:**
```python
@classmethod
def from_atom(cls, atom: Union[Element, DynamicElement, QueryElement, 'DynamicQueryElement', AnyElement],
              isotope: Optional[int] = None, **kwargs) -> Union['DynamicQueryElement', 'DynamicAnyElement']:
    source_atomic_number = atom.atomic_number
    source_isotope = atom.isotope if isotope is None else isotope

    dqe = cls(atomic_number=source_atomic_number, isotope=source_isotope, **kwargs)

    if isinstance(atom, (QueryElement, DynamicQueryElement, AnyElement)):
        pass  # Dead code - does nothing!
    elif isinstance(atom, (Element, DynamicElement)):
        pass  # Dead code - does nothing!

    if isinstance(atom, (DynamicQueryElement, DynamicAnyElement)) and isotope is None and not kwargs:
        return atom.copy()  # Returns copy instead of new instance!

    return dqe
```

**Problems:**
1. Lines 162-165 are dead code (empty `pass` statements)
2. Inconsistent return: sometimes creates new instance, sometimes returns copy
3. `**kwargs` passed but never used in `__init__`

---

### Issue 26: Property Type Inconsistency - `neighbors`

| Class | `neighbors` Type | Meaning |
|-------|-----------------|---------|
| `Element.neighbors` | `int` | Actual neighbor count |
| `Query.neighbors` | `Tuple[int, ...]` | Allowed values for matching |
| `DynamicQuery.neighbors` | `Tuple[int, ...]` | Allowed values for reactant state |

**Problem:** Same property name with different types makes generic code difficult:

```python
def process_atom(atom):
    n = atom.neighbors  # Is this int or Tuple[int, ...]?
```

**Recommendation:** Consider renaming query properties to `allowed_neighbors` or similar.

---

### Issue 27: `DynamicQueryElement` Missing `is_dynamic` Property

**Location:** `dynamic.py:198-203` (in `DynamicElement`)

```python
@property
def is_dynamic(self) -> bool:
    """
    Atom has dynamic features
    """
    return self.charge != self.p_charge or self.is_radical != self.p_is_radical
```

**Problem:** `DynamicQueryElement` doesn't override this property, but it doesn't have `_charge`, `_is_radical`, `_p_charge`, `_p_is_radical` slots. It relies on container-stored values via `DynamicQuery` parent.

When `DynamicQueryElement` is not attached to a container:
```python
dqe = DynamicQueryElement(atomic_number=6, isotope=None)
dqe.is_dynamic  # Will raise IsNotConnectedAtom!
```

This asymmetry with `DynamicElement` (which works standalone) is confusing.

---

### Issue 28: Container Attachment Architecture Mismatch

**`DynamicElement`** stores state internally:
```python
__slots__ = ('_charge', '_is_radical', '_p_charge', '_p_is_radical', '_isotope')
```

**`DynamicQueryElement`** relies on container for most state:
```python
__slots__ = ('_atomic_number', '_isotope')  # Only these stored locally
# charge, is_radical, p_charge, p_is_radical come from container via DynamicQuery
```

**Problem:** This creates two different paradigms:
1. `DynamicElement` can exist standalone with full state
2. `DynamicQueryElement` is mostly useless without container attachment

**Evidence in `cgr_query.py:70-71`:**
```python
self._charges[_map] = getattr(atom, '_charge', 0)  # Falls back to 0 because DynamicQueryElement has no _charge
self._radicals[_map] = getattr(atom, '_is_radical', False)
```

---

### Issue 29: Validation Logic Duplication

**In `query.py:27-43`:**
```python
def _validate(value, prop):
    if value is None:
        return ()
    elif isinstance(value, int):
        if value < 0 or value > 14:
            raise ValueError(f'{prop} should be in range [0, 14]')
        return (value,)
    ...
```

**In `cgr_query.py:227-242`:**
```python
@staticmethod
def _validate_neighbors(neighbors):
    if neighbors is None:
        neighbors = ()
    elif isinstance(neighbors, int):
        if neighbors < 0 or neighbors > 14:
            raise ValueError('neighbors should be in range [0, 14]')
        neighbors = (neighbors,)
    ...
```

**Problem:** Same validation logic duplicated. Should use shared utility.

---

### Issue 30: `__eq__` and `__hash__` Inconsistencies

**`DynamicElement.__eq__` (`dynamic.py:186-192`):**
```python
def __eq__(self, other):
    return isinstance(other, DynamicElement) and self.atomic_number == other.atomic_number and \
        self.isotope == other.isotope and self.charge == other.charge and self.is_radical == other.is_radical and \
        self.p_charge == other.p_charge and self.p_is_radical == other.p_is_radical
```

**`DynamicQueryElement.__eq__` (`dynamic_query.py:179-183`):**
```python
def __eq__(self, other):
    if not isinstance(other, DynamicQueryElement):
        return False
    return (self.atomic_number == other.atomic_number and
            self.isotope == other.isotope)  # Missing charge/radical comparison!
```

**Problem:** `DynamicQueryElement` equality ignores charge, radical, and product states. Two atoms with same element but different charges are considered equal.

**`DynamicQueryElement.__hash__` (`dynamic_query.py:185-186`):**
```python
def __hash__(self):
    return hash((self.atomic_number, self.isotope or 0))  # Also ignores charge/radical
```

This could cause hash collisions and incorrect dict/set behavior.

---

### Issue 31: Missing `__repr__` in `DynamicQueryElement`

**`DynamicQueryElement`** has no `__repr__` method, so it falls back to `object.__repr__`:
```python
>>> dqe = DynamicQueryElement(6, None)
>>> repr(dqe)
'<chython.periodictable.base.dynamic_query.DynamicQueryElement object at 0x...>'
```

Compare with `Element.__repr__`:
```python
>>> repr(C())
'C()'
```

---

### Issue 32: `DynamicAnyElement` Missing Several Properties

**Location:** `dynamic_query.py:189-219`

Missing compared to `DynamicQueryElement`:
- `mdl_isotope` property
- `isotope` property (takes parameter but doesn't store it)
- `copy()` method
- `from_symbol()`, `from_atomic_number()`, `from_atom()` factory methods

**Problem:** These classes have different interfaces despite serving similar purposes.

---

### Issue 33: Inconsistent Exception Handling for Container Attachment

**`Dynamic.charge` property (`dynamic.py:28-35`):**
```python
@Core.charge.setter
def charge(self, charge):
    try:
        g = self._graph()
        g._charges[self._map] = g._validate_charge(charge)
        g.flush_cache()
    except AttributeError:
        raise IsNotConnectedAtom
```

**`DynamicQueryElement.charge`** - inherits from `DynamicQuery` which inherits from `Dynamic`, but `DynamicQueryElement` stores `_charge` nowhere, so accessing `charge` on unattached atom raises `IsNotConnectedAtom`.

Meanwhile, `DynamicElement.charge` (`dynamic.py:171-172`):
```python
@property
def charge(self) -> int:
    return self._charge  # Works even when unattached
```

---

### Issue 34: Pickling Issues with `DynamicQueryElement`

**Location:** `cgr_query.py:282-312`

```python
def __getstate__(self):
    atomic_numbers = {n: getattr(a, '_atomic_number', getattr(a, 'atomic_number', 0)) for n, a in self._atoms.items()}
    isotopes = {n: getattr(a, '_isotope', getattr(a, 'isotope', None)) for n, a in self._atoms.items()}
    ...

def __setstate__(self, state):
    ...
    for n, atom in self._atoms.items():
        try:
            atom._atomic_number = atomic_numbers.get(n, getattr(atom, '_atomic_number', 0))
            atom._isotope = isotopes.get(n, getattr(atom, '_isotope', None))
            atom._attach_to_graph(self, n)
        except Exception:
            pass  # Silent failure!
```

**Problems:**
1. Bare `except Exception: pass` hides errors
2. Complex workaround to handle atoms that may or may not have `_atomic_number`
3. State reconstruction is fragile

---

### Issue 35: Missing `_attach_to_graph` in Some Atom Types

**Location:** `cgr_query.py:67`

```python
atom._attach_to_graph(self, _map)
```

But `DynamicQueryElement` inherits `_attach_to_graph` from `Core` via `Dynamic` → `DynamicQuery`. However, `DynamicAnyElement` also inherits this.

**`DynamicElement`** does NOT inherit from `Core`, so it has no `_attach_to_graph` method:
```python
class DynamicElement(ABC):  # No Core inheritance!
    __slots__ = ('_charge', '_is_radical', '_p_charge', '_p_is_radical', '_isotope')
```

This means `DynamicElement` cannot be attached to a graph in the same way.

---

## 4. SMARTS Module Impact Analysis

### Issue 36: Critical - No CGR SMARTS Parser

The current `smarts()` function in `files/daylight/smarts.py` only produces `QueryContainer` objects:

**`smarts.py:79`:**
```python
g = QueryContainer(data)
```

**Problem:** There is no `cgr_smarts()` function to parse CGR SMARTS patterns into `QueryCGRContainer`. The tokenizer supports dynamic bond tokens (type 13) and dynamic atom tokens (types 11, 12), but there's no parser that uses them to create query CGR containers.

**Evidence in `tokenize.py:49-50`:**
```python
# 13: CGR bond
# 14: CGR atom
```

And `tokenize.py:64-69` defines `dynamic_bonds` dict, showing CGR bond syntax is tokenizable.

**Impact:** Users cannot create reaction pattern queries using SMARTS notation. They must manually construct `QueryCGRContainer` objects.

**Recommendation:** Add `cgr_smarts()` function:
```python
def cgr_smarts(data: str) -> QueryCGRContainer:
    """Parse CGR SMARTS string into QueryCGRContainer."""
    ...
```

---

### Issue 37: `QueryCGRSmiles` Generates Non-Parseable Output

**Location:** `smiles.py:696-718`

```python
def _format_atom(self, n, adjacency, **kwargs):
    ...
    if kwargs.get('hybridization', True) and hybridization:
        smi.append(';')
        smi.append(''.join(hybridization_str[x] for x in hybridization))
        smi.append('>')  # Dynamic marker
        smi.append(''.join(hybridization_str[x] for x in p_hybridization))
        ...
```

**Problem:** The output format uses `>` separators for dynamic properties, but there's no corresponding parser in `smarts_tokenize()` or `smiles_tokenize()` that can parse this back.

**Round-trip failure:**
```python
cgr = some_cgr.substructure([1, 2], as_query=True)
smiles_str = str(cgr)  # Generates QueryCGR SMILES
# No way to parse smiles_str back into QueryCGRContainer!
```

---

### Issue 38: `_query_parse()` Doesn't Handle Dynamic Properties

**Location:** `tokenize.py:310-370`

```python
def _query_parse(token):
    out = {}
    # ... parses isotope, charge, mapping, stereo, element

    for p in primitives[1:]:
        if p == 'a':  # aromatic atom
            out['hybridization'] = 4
        elif p == '!R':
            out['ring_sizes'] = 0
        elif p == 'M':
            out['masked'] = True
        # ... handles D, h, r, x, z
```

**Problem:** No handling for:
- `p_charge` (product charge)
- `p_hybridization` (product hybridization)
- `p_neighbors` (product neighbors)
- `p_is_radical` (product radical)

These are required by `QueryCGRContainer.add_atom()`.

---

### Issue 39: SMARTS Returns `QueryContainer`, Not Compatible with CGR

**Current flow:**
```
SMARTS string → smarts() → QueryContainer
```

**Expected for reactions:**
```
CGR SMARTS string → cgr_smarts() → QueryCGRContainer
```

**Interface mismatch:** `QueryContainer` and `QueryCGRContainer` have different interfaces:

| Feature | `QueryContainer` | `QueryCGRContainer` |
|---------|-----------------|---------------------|
| Atom type | `QueryElement` | `DynamicQueryElement` |
| Bond type | `QueryBond` | `DynamicBond` |
| Has `p_charge` | No | Yes |
| Has `p_hybridization` | No | Yes |
| Has `neighbors` type | `Tuple[int, ...]` | `Tuple[int, ...]` (same) |

---

### Issue 40: Missing `DynamicQueryBond` Class

**Location:** `cgr_query.py:80-85`

```python
def add_bond(self, n, m, bond: Union[DynamicBond, Bond, int]):
    if isinstance(bond, Bond):
        bond = DynamicBond.from_bond(bond)
    elif not isinstance(bond, DynamicBond):
        bond = DynamicBond(bond, bond)
    super().add_bond(n, m, bond)
```

**Problem:** `QueryCGRContainer` uses `DynamicBond`, but there's no `DynamicQueryBond` that could support bond order lists (like `QueryBond` does for queries).

**`QueryBond` in `containers/bonds.py`:**
```python
class QueryBond:
    def __init__(self, order, ..., not_order=None):
        self._order = (order,) if isinstance(order, int) else tuple(order)
```

But `DynamicBond` only supports single order values:
```python
class DynamicBond:
    def __init__(self, order, p_order):
        self._order = order  # Single int, not tuple
        self._p_order = p_order
```

**Impact:** Cannot express "bond is single OR double in reactant, becomes triple in product" patterns.

---

### Issue 41: Tokenizer Inconsistency - CGR Atom Types

**Location:** `tokenize.py:402-408`

```python
if element in ('c', 'n', 'o', 'p', 's', 'as', 'se', 'b'):
    _type = 12  # aromatic CGR atom
    element = element.capitalize()
else:
    _type = 11  # aliphatic CGR atom
return _type, {'element': element, 'isotope': isotope, 'is_radical': is_radical,
               'parsed_mapping': 0, 'cgr': cgr, 'charge': charge}
```

**Problem:** Types 11 and 12 are defined but processed identically in `parser.py`:

**`parser.py:162-163`:**
```python
if token_type == 11:
    cgr.extend((atom_num, *x) for x in token.pop('cgr'))
```

Type 12 (aromatic CGR atom) has no special handling - it falls through to regular atom processing, losing the CGR information.

---

### Issue 42: `smarts()` Ignores CGR Tokens

**Location:** `smarts.py:79-103`

```python
g = QueryContainer(data)

mapping = {}
for i, a in enumerate(parsed['atoms']):
    mapping[i] = n = a.pop('parsed_mapping', 0) or next(...)
    e = a.pop('element')
    if isinstance(e, int):
        e = QueryElement.from_atomic_number(e)
    elif isinstance(e, str):
        e = QueryElement.from_symbol(e)
    else:
        e = partial(ListElement, e)
    g.add_atom(e(**a), n)
```

**Problem:** If `parsed` contains a `'cgr'` key (set by `parser.py:181-182`), it's completely ignored:

**`parser.py:179-182`:**
```python
mol = {'atoms': atoms, 'bonds': bonds, ...}
if cgr:
    mol['cgr'] = cgr  # This is set but never used by smarts()
return mol
```

---

### Issue 43: `QueryElement.from_atom()` vs `DynamicQueryElement.from_atom()` Signature Mismatch

**`query.py:410-435`:**
```python
@classmethod
def from_atom(cls, atom: 'Element', neighbors=False, hybridization=False, heteroatoms=False,
              hydrogens=False, ring_sizes=False, stereo=False) -> 'QueryElement':
```

**`dynamic_query.py:155-156`:**
```python
@classmethod
def from_atom(cls, atom: Union[Element, DynamicElement, QueryElement, 'DynamicQueryElement', AnyElement],
              isotope: Optional[int] = None, **kwargs) -> Union['DynamicQueryElement', 'DynamicAnyElement']:
```

**Problems:**
1. Different parameter names: `neighbors=False` vs `**kwargs`
2. Different behavior: `QueryElement.from_atom()` uses flags to copy properties, `DynamicQueryElement.from_atom()` doesn't copy any query properties
3. `DynamicQueryElement.from_atom()` ignores `neighbors`, `hybridization`, etc. even if passed in `**kwargs`

---

### Issue 44: `substructure(as_query=True)` Creates Incomplete Query

**Location:** `cgr.py:277-326`

```python
def substructure(self, atoms, *, as_query: bool = False, **kwargs):
    graph_type = query.QueryCGRContainer if as_query else self.__class__
    atom_type = DynamicQueryElement if as_query else DynamicElement

    for n in atoms_to_include:
        new_atom = atom_type.from_atom(self._atoms[n])  # Uses from_atom
        sub.add_atom(new_atom, n, xy=self._plane.get(n))
```

**Problem:** When creating a query from CGR substructure:
1. `DynamicQueryElement.from_atom()` only copies `atomic_number` and `isotope`
2. Charge, radical state, hybridization, neighbors are NOT transferred as query constraints
3. The resulting query is too permissive - it matches any atom of that element type

**Expected behavior:** Should create queries with actual constraints:
```python
# If CGR has carbon with charge=1, hybridization=2
# Query should match only atoms with charge=1, hybridization in (2,)
```

---

### Issue 45: Test Coverage Gaps for CGR SMARTS

**Location:** `test_daylight_smarts.py`

- All tests use molecule SMARTS, none test CGR/reaction SMARTS
- No tests for dynamic bond parsing (`[->]`, `[=>]`, etc.)
- No tests for dynamic charge parsing (`[C+>0]`, etc.)

**`test_cgr.py:326-336`:**
```python
def test_query_cgr_smiles_non_empty_with_dynamic_tokens(self):
    cgr = create_simple_cgr()
    q = cgr.substructure([1, 2], as_query=True)
    smiles_str = str(q)
    self.assertIsInstance(smiles_str, str)
    self.assertTrue(len(smiles_str) > 0)
    # Query CGR SMILES should contain dynamic product markers '>'
    self.assertIn('>', smiles_str)
```

This only tests output generation, not round-trip parsing.

---

## 5. Graph/Container Base Issues

### Issue 46: `Graph.copy()` Calls `atom.copy(full=True)` - But `DynamicElement.copy()` Ignores It

**Location:** `graph.py:157-158`

```python
def copy(self):
    copy._atoms = {n: atom.copy(full=True) for n, atom in self.atoms()}
```

**Problem:** As noted earlier, `DynamicElement.copy()` and `DynamicQueryElement.copy()` ignore all parameters:

```python
# DynamicElement.copy()
def copy(self, *_, **__):  # Ignores full=True!
```

**Impact:** When copying a `CGRContainer`, the `full=True` parameter is silently ignored, potentially losing stereo or other state.

---

### Issue 47: `Graph.has_bond()` Inconsistent Exception Handling

**Location:** `graph.py:66-70`

```python
def has_bond(self, n: int, m: int) -> bool:
    try:
        return m in self._bonds[n]
    except KeyError:
        raise AtomNotFound  # Raises exception instead of returning False
```

**Problem:** `has_bond()` raises `AtomNotFound` when atom `n` doesn't exist, but the method name implies it should return `bool`. Compare with `has_atom()`:

```python
def has_atom(self, n: int) -> bool:
    return n in self._atoms  # Returns False, doesn't raise
```

**Inconsistency:** One returns `False`, the other raises.

---

### Issue 48: `CGRContainer.has_bond()` Different Implementation

**Location:** `cgr.py:81-86`

```python
def has_bond(self, n: int, m: int) -> bool:
    try:
        return m in self._bonds[n]
    except KeyError:
        return False  # Returns False instead of raising
```

**Problem:** `CGRContainer.has_bond()` returns `False` on missing atom, but `Graph.has_bond()` raises `AtomNotFound`. This is inconsistent behavior within the same inheritance hierarchy.

---

### Issue 49: Missing `DynamicBond` in `__all__` Export

**Location:** `containers/__init__.py:48-49`

```python
__all__ = [x for x in locals() if x.endswith('Container')]
__all__.extend(['Bond', 'QueryBond', 'unpack', 'unpach', 'from_rdkit', 'from_rdkit_molecule'])
```

**Problem:** `DynamicBond` is not exported in `__all__`, but it's needed for CGR operations:

```python
from chython.containers import DynamicBond  # Works due to * import
# But not explicitly listed in __all__
```

---

### Issue 50: `DynamicBond.copy()` Ignores Parameters

**Location:** `bonds.py:137-141`

```python
def copy(self, *args, **kwargs) -> 'DynamicBond':
    copy = object.__new__(self.__class__)
    copy._order = self.order
    copy._p_order = self.p_order
    return copy
```

**Problem:** Same pattern as atoms - accepts but ignores parameters. `Bond.copy()` has meaningful `full` and `stereo` parameters that control behavior.

---

### Issue 51: `QueryContainer` Has Hardcoded SMARTS String

**Location:** `query.py:27-38`

```python
class QueryContainer(Graph[Query, QueryBond], QueryIsomorphism):
    __slots__ = ('_smarts',)

    def __init__(self, smarts: str):
        super().__init__()
        self._smarts = smarts

    def __str__(self):
        return self._smarts  # Returns original string, not regenerated
```

**Problem:** If you modify a `QueryContainer` after creation, `str(query)` returns the original SMARTS, not the current state:

```python
q = smarts('[C]')
q.add_atom(QueryElement.from_symbol('N')(), 2)
str(q)  # Still returns '[C]', not '[C].[N]'
```

---

### Issue 52: `QueryCGRContainer` Has No `_smarts` Slot

**Location:** `cgr_query.py:12-13`

```python
class QueryCGRContainer(Graph, QueryCGRSmiles, DepictQueryCGR, Calculate2DCGR):
    __slots__ = ('_p_charges', '_p_radicals', ...)  # No _smarts!
```

**Problem:** Unlike `QueryContainer`, `QueryCGRContainer` doesn't store the original SMARTS string. This means:
1. No `__repr__` that shows the source
2. Can't round-trip serialize via SMARTS

---

### Issue 53: `DynamicAnyElement.isotope` Property Missing

**Location:** `dynamic_query.py:189-219`

```python
class DynamicAnyElement(DynamicQuery):
    __slots__ = ()

    def __init__(self, isotope: Optional[int] = None, **kwargs):
        super().__init__()  # isotope ignored!

    # No isotope property defined!
```

**Problem:** `DynamicAnyElement` accepts `isotope` in `__init__` but:
1. Doesn't store it
2. Has no `isotope` property
3. This will cause `AttributeError` when accessing `.isotope`

---

### Issue 54: `Reactor` Only Accepts `QueryContainer`, Not `QueryCGRContainer`

**Location:** `reactor.py:64-66`

```python
if not all(isinstance(x, QueryContainer) for x in patterns):
    raise TypeError('invalid params')
```

**Problem:** Cannot use `QueryCGRContainer` as reaction patterns. The Reactor was designed before CGR queries existed.

**Impact:** Users cannot use CGR-based queries (with dynamic properties) in the Reactor.

---

### Issue 55: Isomorphism Module Doesn't Import `DynamicQueryElement`

**Location:** `isomorphism.py:27`

```python
from ..periodictable import Element, Query, AnyElement, AnyMetal, ListElement, QueryElement, ExtendedQuery
```

**Problem:** `DynamicQueryElement` and `DynamicAnyElement` not imported. The isomorphism algorithms may not properly handle CGR query atoms.

---

### Issue 56: `_delete_bond_internal` Not In Original API

**Location:** `graph.py:145-151`

```python
def _delete_bond_internal(self, n: int, m: int):
    """Internal bond deletion without flushing cache."""
    if n not in self._bonds or m not in self._bonds.get(n, {}):
        raise BondNotFound(f"Bond between {n} and {m} not found for deletion")
    del self._bonds[n][m]
    del self._bonds[m][n]
```

**Observation:** This is a new helper method that's good practice, but `CGRContainer` doesn't use it - it has its own `delete_bond` implementation that duplicates this logic.

---

### Issue 57: `DynamicQuery` Missing `isotope` Property

**Location:** `dynamic_query.py:8-79`

```python
class DynamicQuery(Dynamic):
    __slots__ = ()
    # No isotope property!
```

**Problem:** `DynamicQuery` is base for `DynamicQueryElement` and `DynamicAnyElement`, but doesn't define `isotope`. Only `DynamicQueryElement` has it, creating inconsistent interfaces.

---

## 6. Summary Tables

### Issues by Module

| Module | Issue Count |
|--------|-------------|
| CGR Container | 9 |
| Depict/Visualization | 11 |
| Periodictable/Atoms | 15 |
| SMARTS Module | 10 |
| Graph/Container Base | 12 |
| **Total** | **57** |

### Issues by Severity

| Severity | Count | Examples |
|----------|-------|----------|
| **Critical** | 6 | No CGR SMARTS parser, dual config systems, `copy()` ignoring params |
| **High** | 12 | `__eq__`/`__hash__` incomplete, semantic changes, missing exports |
| **Medium** | 25 | Variable shadowing, validation duplication, missing properties |
| **Low** | 14 | Docstring style, file headers, naming inconsistencies |

### Most Critical Issues

1. **No `cgr_smarts()` parser** - Blocks CGR query workflows entirely
2. **Dual configuration systems in depict** - Causes user confusion, CGR ignores settings
3. **`copy()` methods ignoring parameters** - Silent data loss risk
4. **`DynamicQueryElement.__eq__` incomplete** - Dict/set bugs possible
5. **`Reactor` doesn't accept `QueryCGRContainer`** - Blocks integration
6. **Container attachment paradigm mismatch** - Inconsistent behavior between atom types

### Recommended Fix Priority

1. **High Priority:**
   - Unify depict configuration systems
   - Implement `cgr_smarts()` parser
   - Fix `copy()` methods to respect parameters
   - Complete `DynamicQueryElement.__eq__` and `__hash__`

2. **Medium Priority:**
   - Add `DynamicQueryBond` class
   - Fix `from_atom()` signature mismatches
   - Update `Reactor` to accept `QueryCGRContainer`
   - Add missing `__repr__` methods

3. **Low Priority:**
   - Fix docstring formats
   - Update file headers
   - Add missing `__all__` exports
   - Clean up dead configuration keys

---

## Appendix: Files Reviewed

- `chython/containers/cgr.py`
- `chython/containers/cgr_query.py`
- `chython/containers/query.py`
- `chython/containers/graph.py`
- `chython/containers/bonds.py`
- `chython/containers/reaction.py`
- `chython/containers/__init__.py`
- `chython/algorithms/depict.py`
- `chython/algorithms/smiles.py`
- `chython/algorithms/isomorphism.py`
- `chython/algorithms/calculate2d/__init__.py`
- `chython/algorithms/x3dom.py`
- `chython/periodictable/__init__.py`
- `chython/periodictable/base/dynamic.py`
- `chython/periodictable/base/dynamic_query.py`
- `chython/periodictable/base/query.py`
- `chython/files/daylight/smarts.py`
- `chython/files/daylight/tokenize.py`
- `chython/files/daylight/parser.py`
- `chython/reactor/reactor.py`
- `chython/containers/tests/test_cgr.py`
- `chython/files/daylight/test/test_daylight_smarts.py`
