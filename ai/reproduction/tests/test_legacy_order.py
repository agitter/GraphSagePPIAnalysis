"""Focused tests for the recovered CPython 2.7 ordering behavior."""

from __future__ import annotations

import hashlib

from graphsage_ppi_repro.legacy_order import (
    Python2InsertionDict,
    ordered_string_keys,
    python2_string_hash,
)


def test_known_python2_string_hashes() -> None:
    assert python2_string_hash("") == 0
    assert python2_string_hash("0", 64) == 6_144_018_481
    assert python2_string_hash("42", 64) == 6_656_039_988_060_076
    assert python2_string_hash("GO:0008150", 64) == -2_827_199_014_123_721_004
    assert python2_string_hash("GO:0008150", 32) == -841_863_468


def test_duplicate_insertions_do_not_change_table() -> None:
    table = Python2InsertionDict(python2_string_hash, word_size_bits=64)
    assert table.insert("101")
    before = table.slots_and_keys()
    assert not table.insert("101")
    assert table.slots_and_keys() == before
    assert len(table) == 1


def test_resize_and_table_iteration_are_not_insertion_order() -> None:
    inserted = [str(index) for index in range(20)]
    observed = ordered_string_keys(inserted)
    assert observed == [
        "11",
        "10",
        "13",
        "12",
        "15",
        "14",
        "17",
        "16",
        "19",
        "18",
        "1",
        "0",
        "3",
        "2",
        "5",
        "4",
        "7",
        "6",
        "9",
        "8",
    ]
    assert observed != inserted


def test_full_class_map_key_order_positive_control() -> None:
    """Match the irregular key order serialized in ppi-class_map.json."""

    ordered = ordered_string_keys(str(index) for index in range(56_944))
    digest = hashlib.sha256(
        "".join(f"{key}\n" for key in ordered).encode("ascii")
    ).hexdigest()
    assert digest == "bdbb25c34f37ff939ef6a3e8ca9f8911479b8bdcbdd7c6395c655091ad42f52f"
    assert ordered[:3] == ["50088", "44884", "11542"]
    assert ordered[-3:] == ["38344", "38347", "38346"]


def test_32_and_64_bit_models_are_distinct() -> None:
    keys = [str(index) for index in range(5_000)]
    assert ordered_string_keys(keys, word_size_bits=32) != ordered_string_keys(
        keys, word_size_bits=64
    )
