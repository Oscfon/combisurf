r"""
Tests for the families of words of :mod:`combisurf.lyndon_word_family`.

Each family is checked against a definition that does not share any code with
it: the words of the family are compared with the words of a larger family
filtered by a predicate from :mod:`combisurf.word`, or with a plain Python
computation.
"""
# ****************************************************************************
#  This file is part of combisurf
#
#       Copyright (C) 2026 Vincent Delecroix
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# ****************************************************************************

from array import array

import pytest

from combisurf.word import word_is_cyclically_reduced
from combisurf.word_family import (
    ByLengthWordFamily,
    WordFamily,
    WordFamilyIterator,
)
from combisurf.lyndon_word_family import (
    CyclicallyReducedLyndonWordFamily,
    LexCyclicallyReducedLyndonWordFamily,
    LexLyndonWordFamily,
    LyndonWordFamily,
    UpToInverseWordFamily,
    cyclically_reduced_lyndon_words,
    lyndon_words,
)


# the two sets of words, each as a family and as a naive definition ##########

def lyndon_naive(n, length):
    r"""
    The Lyndon words of the given length, by brute force.
    """
    from itertools import product
    return [w for w in product(range(n), repeat=length)
            if all(w < w[i:] + w[:i] for i in range(1, length))]


def cyclically_reduced_naive(n, length):
    r"""
    The Lyndon words on ``2n`` letters that are cyclically reduced, taken from
    the family of all of them so that the two algorithms share nothing.
    """
    return [tuple(w) for w in LyndonWordFamily(2 * n, length)
            if word_is_cyclically_reduced(w)]


SETS = (
    ("lyndon", lambda n, length: LyndonWordFamily(n, length),
     lambda n, length: LexLyndonWordFamily(n, 0, length + 1),
     lyndon_naive),
    ("cyclically_reduced", lambda n, length: CyclicallyReducedLyndonWordFamily(n, length),
     lambda n, length: LexCyclicallyReducedLyndonWordFamily(n, 0, length + 1),
     cyclically_reduced_naive),
)


def words(family):
    return [tuple(w) for w in family]


# the words themselves ######################################################

@pytest.mark.parametrize("name,fixed,lex,naive", SETS)
def test_matches_naive_definition(name, fixed, lex, naive):
    for n in range(4):
        for length in range(6):
            assert words(fixed(n, length)) == naive(n, length), (name, n, length)


@pytest.mark.parametrize("name,fixed,lex,naive", SETS)
def test_words_are_sorted(name, fixed, lex, naive):
    for n in range(4):
        for length in range(7):
            got = words(fixed(n, length))
            assert got == sorted(got), (name, n, length)


def test_empty_alphabet_and_single_generator():
    # no letter at all, so no word of positive length
    for length in range(1, 5):
        assert words(LyndonWordFamily(0, length)) == []
        assert words(CyclicallyReducedLyndonWordFamily(0, length)) == []
    # one letter: only the word of length one, a^k being a proper power
    assert words(LyndonWordFamily(1, 1)) == [(0,)]
    assert words(LyndonWordFamily(1, 5)) == []
    # one generator: its two letters are inverse of each other
    assert words(CyclicallyReducedLyndonWordFamily(1, 1)) == [(0,), (1,)]
    assert words(CyclicallyReducedLyndonWordFamily(1, 5)) == []
    # the empty word is the only one of length zero
    for n in range(4):
        assert words(LyndonWordFamily(n, 0)) == [()]
        assert words(CyclicallyReducedLyndonWordFamily(n, 0)) == [()]


@pytest.mark.parametrize("cls", [LyndonWordFamily, CyclicallyReducedLyndonWordFamily])
def test_negative_arguments(cls):
    with pytest.raises(ValueError):
        cls(-1, 3)
    with pytest.raises(ValueError):
        cls(2, -1)


# the four ways of getting the words agree ##################################

@pytest.mark.parametrize("name,fixed,lex,naive", SETS)
def test_list_iter_flat_count_agree(name, fixed, lex, naive):
    for n in range(4):
        for length in range(7):
            want = words(fixed(n, length))
            assert [tuple(w) for w in fixed(n, length).list()] == want
            assert fixed(n, length).count() == len(want)
            flat = fixed(n, length).flat()
            if length:
                assert [tuple(flat[length * i: length * (i + 1)])
                        for i in range(len(flat) // length)] == want
            else:
                assert len(flat) == 0


@pytest.mark.parametrize("name,fixed,lex,naive", SETS)
def test_cardinality_against_enumeration(name, fixed, lex, naive):
    r"""
    ``count`` uses a closed form where one is known; it has to agree with
    running through the family.
    """
    for n in range(5):
        for length in range(8):
            family = fixed(n, length)
            assert family.count() == len(family.list()), (name, n, length)


def test_list_returns_distinct_arrays():
    got = CyclicallyReducedLyndonWordFamily(2, 6).list()
    assert len({id(w) for w in got}) == len(got)
    assert len({tuple(w) for w in got}) == len(got)


def test_flat_refuses_a_family_of_several_lengths():
    with pytest.raises(ValueError):
        LexLyndonWordFamily(2, 1, 4).flat()
    with pytest.raises(ValueError):
        ByLengthWordFamily(LyndonWordFamily(2, 0), 1, 4).flat()


# families and iterators are different things ###############################

def test_family_is_not_an_iterator():
    family = LyndonWordFamily(2, 3)
    with pytest.raises(TypeError):
        next(family)


def test_iterators_are_independent():
    family = LyndonWordFamily(3, 4)
    all_words = words(family)
    first, second = iter(family), iter(family)
    assert isinstance(first, WordFamilyIterator)
    assert first is not second
    next(first)
    assert [tuple(w) for w in first] == all_words[1:]
    assert [tuple(w) for w in second] == all_words
    # and the family is none the wiser
    assert words(family) == all_words
    assert family.count() == len(all_words)


def test_iterating_an_iterator_gives_it_back():
    it = iter(LyndonWordFamily(2, 3))
    assert iter(it) is it
    first = next(it)
    assert [tuple(w) for w in it] == words(LyndonWordFamily(2, 3))[1:]
    assert tuple(first) == words(LyndonWordFamily(2, 3))[0]


def test_iterator_needs_a_family():
    with pytest.raises(ValueError):
        WordFamilyIterator(None)


def test_empty_base_family():
    assert WordFamily().list() == []
    assert WordFamily().count() == 0
    assert list(WordFamily()) == []


# the buffer that iteration hands out #######################################

def test_a_kept_word_is_never_overwritten():
    r"""
    The array of a word that the caller let go of is filled again; a word
    that is still held has to be left alone.
    """
    family = CyclicallyReducedLyndonWordFamily(2, 6)
    want = words(family)

    kept = []
    for w in family:
        kept.append(w)
    assert [tuple(w) for w in kept] == want
    assert len({id(w) for w in kept}) == len(kept)

    # holding the previous word while asking for the next one is enough
    previous = None
    for w in family:
        assert w is not previous
        previous = w

    # a memoryview holds the array too
    views = []
    for w in CyclicallyReducedLyndonWordFamily(2, 5):
        views.append(memoryview(w))
        del w
    assert [v.tolist() for v in views] == [list(w) for w in words(
        CyclicallyReducedLyndonWordFamily(2, 5))]


def test_a_released_word_is_recycled():
    family = CyclicallyReducedLyndonWordFamily(2, 6)
    it = iter(family)
    w = next(it)
    address = id(w)
    del w
    # had the array been freed, this would likely take its place
    ballast = [array('i', [9] * 6) for _ in range(50)]
    assert id(next(it)) == address
    del ballast


# the two orders over a range of lengths ####################################

@pytest.mark.parametrize("name,fixed,lex,naive", SETS)
def test_by_length_is_the_union(name, fixed, lex, naive):
    for n in range(4):
        for size_min in range(4):
            for size_max in range(6):
                want = [w for length in range(size_min, size_max)
                        for w in words(fixed(n, length))]
                family = ByLengthWordFamily(fixed(n, 0), size_min, size_max)
                assert words(family) == want, (name, n, size_min, size_max)
                assert family.count() == len(want)


@pytest.mark.parametrize("name,fixed,lex,naive", SETS)
def test_lex_is_the_union_sorted(name, fixed, lex, naive):
    for n in range(4):
        for size_max in range(7):
            want = sorted(w for length in range(size_max)
                          for w in words(fixed(n, length)))
            family = lex(n, size_max - 1)
            got = words(family)
            assert got == want, (name, n, size_max)
            assert got == sorted(got)
            assert family.count() == len(want)


def test_lex_and_by_length_hold_the_same_words():
    for n in range(4):
        for size_min in range(4):
            for size_max in range(7):
                by_length = lyndon_words(n, size_min, size_max)
                lexicographic = lyndon_words(n, size_min, size_max, order="lex")
                assert sorted(words(by_length)) == sorted(words(lexicographic))
                assert by_length.count() == lexicographic.count()


def test_a_single_length_ignores_the_order():
    for order in ("length", "lex"):
        for length in range(6):
            assert words(lyndon_words(2, length, order=order)) == \
                words(LyndonWordFamily(2, length))


def test_unknown_order():
    with pytest.raises(ValueError):
        lyndon_words(2, 1, 4, order="by-weight")


def test_by_length_over_a_family_that_wraps_another():
    r"""
    :class:`UpToInverseWordFamily` does not own its buffer, so running it at
    several lengths goes through the resizing part of the protocol.
    """
    for n in range(4):
        for size_max in range(6):
            base = UpToInverseWordFamily(CyclicallyReducedLyndonWordFamily(n, 0))
            family = ByLengthWordFamily(base, 1, size_max)
            want = [w for length in range(1, size_max)
                    for w in words(UpToInverseWordFamily(
                        CyclicallyReducedLyndonWordFamily(n, length)))]
            assert words(family) == want, (n, size_max)
            assert family.count() == len(want)


# keeping one word out of each inverse pair #################################

def inverse(w):
    return tuple(x ^ 1 for x in reversed(w))


def lyndon_of_class(w):
    return min(w[i:] + w[:i] for i in range(len(w))) if w else ()


def test_up_to_inverse_matches_its_definition():
    for n in range(4):
        for length in range(8):
            base = words(CyclicallyReducedLyndonWordFamily(n, length))
            want = [w for w in base if lyndon_of_class(inverse(w)) > w]
            family = UpToInverseWordFamily(CyclicallyReducedLyndonWordFamily(n, length))
            assert words(family) == want, (n, length)
            assert family.count() == len(want)


def test_up_to_inverse_halves_the_family():
    r"""
    No element of a free group is conjugate to its inverse, so the words pair
    up and exactly half of them are kept.
    """
    for n in range(1, 4):
        for length in range(1, 9):
            whole = CyclicallyReducedLyndonWordFamily(n, length)
            half = UpToInverseWordFamily(CyclicallyReducedLyndonWordFamily(n, length))
            assert 2 * half.count() == whole.count(), (n, length)


def test_up_to_inverse_over_a_range():
    for n in range(4):
        for size_max in range(7):
            base = [w for length in range(size_max)
                    for w in words(CyclicallyReducedLyndonWordFamily(n, length))]
            want = sorted(w for w in base if lyndon_of_class(inverse(w)) > w)
            for order in ("length", "lex"):
                family = cyclically_reduced_lyndon_words(
                    n, 0, size_max, order=order, up_to_inverse=True)
                got = words(family)
                assert sorted(got) == want, (n, size_max, order)
                assert family.count() == len(want)
                if order == "lex":
                    assert got == sorted(got)


def test_up_to_inverse_refuses_other_families():
    with pytest.raises(ValueError):
        UpToInverseWordFamily(LyndonWordFamily(4, 3))
    with pytest.raises(ValueError):
        UpToInverseWordFamily(None)


def test_against_the_naive_definition_over_a_long_range():
    r"""
    The same check as :func:`test_matches_naive_definition`, pushed as far as
    brute force can go; the naive side runs through the ``(2n)^length`` words,
    so the longest length has to shrink as ``n`` grows.
    """
    for n, max_length in ((1, 16), (2, 12), (3, 10), (4, 8)):
        for length in range(max_length + 1):
            got = iter(cyclically_reduced_lyndon_words(n, length))
            want = iter(cyclically_reduced_naive(n, length))
            while True:
                a = next(got, None)
                b = next(want, None)
                assert (a is None) == (b is None), (n, length)
                if a is None:
                    break
                assert tuple(a) == b, (n, length, tuple(a), b)
