from typing import Dict, Iterable, DefaultDict, TypeVar, NewType
from typing_extensions import TypeAlias

T = TypeVar("T")
ResidueID = NewType("ResidueID", int)
LipidName = NewType("LipidName", str)
LipidId = NewType("LipidId", int)
FrameIndex = NewType("FrameIndex", int)

NestedDictT: TypeAlias = DefaultDict[
    ResidueID, DefaultDict[LipidName, Dict[LipidId, T]]
]

NestedFloatDict: TypeAlias = NestedDictT[float]
NestedIterFloatDict: TypeAlias = NestedDictT[Iterable[float]]
NestedIterIntDict: TypeAlias = NestedDictT[Iterable[FrameIndex]]

# TODO: Fix the type definitions above to this:
# NestedDictT: TypeAlias = DefaultDict[ResidueID, DefaultDict[LipidId, T]]
# NestedIterIntDict: TypeAlias = NestedDictT[Iterable[FrameIndex]]
