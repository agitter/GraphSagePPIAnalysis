"""The small part of CPython 2.7 dictionary behavior needed by this project.

GraphSAGE's PPI row order is reproduced when OhmNet Entrez identifiers are kept
as strings, inserted while reading each edgelist, and iterated in the table
order of a 64-bit CPython 2.7 dictionary without hash randomization.

This module is deliberately *not* a general Python 2 emulator.  It implements
only insertion, resizing, lookup, and occupied-slot iteration for dictionaries
that never delete keys.  The implementation follows the historical perturb
probing rule and resize thresholds and is isolated here so it can be audited
and tested independently from the biological transformations.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass
from typing import Generic, TypeVar

Key = TypeVar("Key")
HashFunction = Callable[[Key, int], int]


class LegacyOrderError(ValueError):
    """Raised for unsupported inputs to the legacy-order model."""


@dataclass(frozen=True)
class TableEntry(Generic[Key]):
    key: Key
    hash_value: int


def _word_mask(word_size_bits: int) -> int:
    if word_size_bits not in (32, 64):
        raise LegacyOrderError("word_size_bits must be 32 or 64")
    return (1 << word_size_bits) - 1


def python2_string_hash(value: str | bytes, word_size_bits: int = 64) -> int:
    """Return CPython 2.7's unrandomized hash for a byte string.

    Entrez IDs and GO IDs in the recovered workflow are ASCII.  ``str`` input
    is therefore encoded strictly as ASCII so accidental Unicode normalization
    cannot silently change the modeled bytes.
    """

    raw = value if isinstance(value, bytes) else value.encode("ascii")
    if not raw:
        return 0
    mask = _word_mask(word_size_bits)
    value_hash = (raw[0] << 7) & mask
    for byte in raw:
        value_hash = ((1_000_003 * value_hash) ^ byte) & mask
    value_hash = (value_hash ^ len(raw)) & mask
    sign_bit = 1 << (word_size_bits - 1)
    if value_hash >= sign_bit:
        value_hash -= 1 << word_size_bits
    return -2 if value_hash == -1 else value_hash


class Python2InsertionDict(Generic[Key]):
    """Insertion-only simulation of a CPython 2.7 dictionary table."""

    _PERTURB_SHIFT = 5
    _MIN_TABLE_SIZE = 8

    def __init__(
        self,
        hash_function: HashFunction[Key],
        *,
        word_size_bits: int = 64,
    ) -> None:
        self.word_size_bits = word_size_bits
        self._mask_for_word = _word_mask(word_size_bits)
        self._hash_function = hash_function
        self._table: list[TableEntry[Key] | None] = [None] * self._MIN_TABLE_SIZE
        self._mask = self._MIN_TABLE_SIZE - 1
        self._used = 0
        self._fill = 0

    def __len__(self) -> int:
        return self._used

    @property
    def table_size(self) -> int:
        return len(self._table)

    def _lookup_slot(self, key: Key, hash_value: int) -> int:
        index = hash_value & self._mask
        perturb = hash_value & self._mask_for_word
        while True:
            entry = self._table[index]
            if entry is None or entry.key == key:
                return index
            index = (index * 5 + 1 + perturb) & self._mask
            perturb >>= self._PERTURB_SHIFT

    def _insert_without_resize(self, key: Key, hash_value: int) -> bool:
        slot = self._lookup_slot(key, hash_value)
        if self._table[slot] is not None:
            return False
        self._table[slot] = TableEntry(key=key, hash_value=hash_value)
        self._used += 1
        self._fill += 1
        return True

    def _resize(self, minimum_used: int) -> None:
        new_size = self._MIN_TABLE_SIZE
        while new_size <= minimum_used:
            new_size <<= 1
        occupied = [entry for entry in self._table if entry is not None]
        self._table = [None] * new_size
        self._mask = new_size - 1
        self._used = 0
        self._fill = 0
        # CPython reinserts entries by scanning the old table, not by original
        # insertion time.  This detail is essential to the final key order.
        for entry in occupied:
            self._insert_without_resize(entry.key, entry.hash_value)

    def insert(self, key: Key) -> bool:
        """Insert *key* and return ``True`` only when it was newly added."""

        hash_value = self._hash_function(key, self.word_size_bits)
        inserted = self._insert_without_resize(key, hash_value)
        if inserted and self._fill * 3 >= self.table_size * 2:
            growth = 2 if self._used > 50_000 else 4
            self._resize(growth * self._used)
        return inserted

    def insert_many(self, keys: Iterable[Key]) -> None:
        for key in keys:
            self.insert(key)

    def occupied_entries(self) -> Iterator[tuple[int, TableEntry[Key]]]:
        """Yield ``(slot, entry)`` in historical dictionary iteration order."""

        for slot, entry in enumerate(self._table):
            if entry is not None:
                yield slot, entry

    def keys(self) -> list[Key]:
        return [entry.key for _, entry in self.occupied_entries()]

    def slots_and_keys(self) -> list[tuple[int, Key]]:
        return [(slot, entry.key) for slot, entry in self.occupied_entries()]


def ordered_string_keys(keys: Iterable[str], *, word_size_bits: int = 64) -> list[str]:
    table: Python2InsertionDict[str] = Python2InsertionDict(
        python2_string_hash, word_size_bits=word_size_bits
    )
    table.insert_many(keys)
    return table.keys()
