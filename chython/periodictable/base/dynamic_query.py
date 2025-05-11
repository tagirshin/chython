from typing import Tuple, Dict, Type, Union, Optional
from .dynamic import Dynamic, DynamicElement
from .element import Element
from .query import QueryElement, AnyElement
from ...exceptions import IsNotConnectedAtom


class DynamicQuery(Dynamic):
    __slots__ = ()

    @property
    def neighbors(self) -> Tuple[int, ...]:
        """
        Number of neighbors of atom in reactant state.
        """
        try:
            return self._graph()._neighbors[self._map]
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def p_neighbors(self) -> Tuple[int, ...]:
        """
        Number of neighbors of atom in product state.
        """
        try:
            return self._graph()._p_neighbors[self._map]
        except AttributeError:
            raise IsNotConnectedAtom

    @neighbors.setter
    def neighbors(self, neighbors):
        try:
            g = self._graph()
            neighbors = g._validate_neighbors(neighbors)
            neighbors, p_neighbors = g._validate_neighbors_pairing(neighbors, g._p_neighbors[self._map])
            g._neighbors[self._map] = neighbors
            g._p_neighbors[self._map] = p_neighbors
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom

    @p_neighbors.setter
    def p_neighbors(self, p_neighbors):
        try:
            g = self._graph()
            p_neighbors = g._validate_neighbors(p_neighbors)
            neighbors, p_neighbors = g._validate_neighbors_pairing(g._neighbors[self._map], p_neighbors)
            g._neighbors[self._map] = neighbors
            g._p_neighbors[self._map] = p_neighbors
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom

    @Dynamic.hybridization.setter
    def hybridization(self, hybridization):
        try:
            g = self._graph()
            hybridization = g._validate_hybridization(hybridization)
            hybridization, p_hybridization = g._validate_hybridization_pairing(hybridization,
                                                                               g._p_hybridizations[self._map])
            g._hybridizations[self._map] = hybridization
            g._p_hybridizations[self._map] = p_hybridization
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom

    @Dynamic.p_hybridization.setter
    def p_hybridization(self, p_hybridization):
        try:
            g = self._graph()
            p_hybridization = g._validate_hybridization(p_hybridization)
            hybridization, p_hybridization = g._validate_hybridization_pairing(g._hybridizations[self._map],
                                                                               p_hybridization)
            g._hybridizations[self._map] = hybridization
            g._p_hybridizations[self._map] = p_hybridization
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom


class DynamicQueryElement(DynamicQuery):
    __slots__ = ('_atomic_number', '_isotope')

    def __init__(self, atomic_number: int, isotope: Optional[int], **kwargs):
        super().__init__()
        self._atomic_number = atomic_number
        self._isotope = isotope

    @property
    def atomic_number(self) -> int:
        return self._atomic_number

    @property
    def isotope(self) -> Optional[int]:
        return self._isotope

    @property
    def _base_element_instance(self) -> Element:
        "Helper to get a base Element instance for delegation."
        return Element.from_atomic_number(self._atomic_number)()

    @property
    def atomic_radius(self) -> float:
        if self._atomic_number == 0:
            return 0.5
        return self._base_element_instance.atomic_radius

    @property
    def isotopes_distribution(self) -> Dict[int, float]:
        if self._atomic_number == 0:
            return {}
        return self._base_element_instance.isotopes_distribution

    @property
    def isotopes_masses(self) -> Dict[int, float]:
        if self._atomic_number == 0:
            return {}
        return self._base_element_instance.isotopes_masses

    @property
    def mdl_isotope(self) -> int:
        if self._atomic_number == 0:
            return 0
        return self._base_element_instance.mdl_isotope

    @property
    def atomic_symbol(self) -> str:
        try:
            if self._atomic_number == 0:
                return 'A'
            return self._base_element_instance.atomic_symbol
        except ValueError:
            return "Xq"

    @classmethod
    def from_symbol(cls, symbol: str, isotope: Optional[int] = None, **kwargs) -> Union['DynamicQueryElement', 'DynamicAnyElement']:
        if symbol == 'A':
            return cls(atomic_number=0, isotope=isotope, **kwargs)
        
        base_element_class = Element.from_symbol(symbol)
        atomic_num = base_element_class().atomic_number
        return cls(atomic_number=atomic_num, isotope=isotope, **kwargs)

    @classmethod
    def from_atomic_number(cls, atomic_number: int, isotope: Optional[int] = None, **kwargs) -> Union['DynamicQueryElement', 'DynamicAnyElement']:
        if atomic_number == 0:
            return cls(atomic_number=0, isotope=isotope, **kwargs)
        
        if not isinstance(atomic_number, int) or atomic_number < 0:
            raise ValueError(f'Invalid query atomic number: {atomic_number}')
        return cls(atomic_number=atomic_number, isotope=isotope, **kwargs)

    @classmethod
    def from_atom(cls, atom: Union[Element, DynamicElement, QueryElement, 'DynamicQueryElement', AnyElement],
                  isotope: Optional[int] = None, **kwargs) -> Union['DynamicQueryElement', 'DynamicAnyElement']:
        source_atomic_number = atom.atomic_number
        source_isotope = atom.isotope if isotope is None else isotope

        dqe = cls(atomic_number=source_atomic_number, isotope=source_isotope, **kwargs)
        
        if isinstance(atom, (QueryElement, DynamicQueryElement, AnyElement)):
            pass
        elif isinstance(atom, (Element, DynamicElement)):
            pass
        
        if isinstance(atom, (DynamicQueryElement, DynamicAnyElement)) and isotope is None and not kwargs:
            return atom.copy()

        return dqe

    def copy(self, *args, **kwargs_copy_call) -> 'DynamicQueryElement':
        new_kwargs = {}
        
        new_copy = self.__class__(atomic_number=self.atomic_number, isotope=self.isotope, **new_kwargs)
        
        return new_copy

    def __eq__(self, other):
        if not isinstance(other, DynamicQueryElement):
            return False
        return (self.atomic_number == other.atomic_number and
                self.isotope == other.isotope)

    def __hash__(self):
        return hash((self.atomic_number, self.isotope or 0))


class DynamicAnyElement(DynamicQuery):
    __slots__ = ()

    def __init__(self, isotope: Optional[int] = None, **kwargs):
        super().__init__()

    @property
    def atomic_symbol(self) -> str:
        return 'A'

    @property
    def atomic_number(self) -> int:
        return 0

    @property
    def isotopes_distribution(self) -> Dict[int, float]:
        return {}

    @property
    def isotopes_masses(self) -> Dict[int, float]:
        return {}

    @property
    def atomic_radius(self):
        return 0.5

    def __eq__(self, other):
        return isinstance(other, DynamicAnyElement)

    def __hash__(self):
        return hash((self.atomic_number,))


__all__ = ['DynamicQueryElement', 'DynamicAnyElement']