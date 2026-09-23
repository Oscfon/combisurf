r"""
Families of Lyndon words

A Lyndon word is a word that is strictly smaller than each of its proper
rotations. This module holds them, and the ones that are cyclically reduced
when the letters ``i`` and ``i ^ 1`` are read as inverse of each other in a
free group, as families in the sense of :mod:`combisurf.word_family`.

The words of a single length come out in lexicographic order::

    sage: from combisurf.lyndon_word_family import lyndon_words
    sage: [tuple(w) for w in lyndon_words(2, 4)]
    [(0, 0, 0, 1), (0, 0, 1, 1), (0, 1, 1, 1)]

Over a range of lengths ``size_min <= n < size_max`` there are two orders to
choose from, and since a family is a set of words together with an order,
they are two different families. By length, one length after the other, or
lexicographically, in which a word comes before every longer word it
begins::

    sage: [tuple(w) for w in lyndon_words(2, 1, 4)]
    [(0,), (1,), (0, 1), (0, 0, 1), (0, 1, 1)]
    sage: [tuple(w) for w in lyndon_words(2, 1, 4, order="lex")]
    [(0,), (0, 0, 1), (0, 1), (0, 1, 1), (1,)]

The cyclically reduced ones are the Lyndon representatives of the conjugacy
classes of the primitive elements of the free group, and
``up_to_inverse`` keeps one word out of each pair of a class and its
inverse::

    sage: from combisurf.lyndon_word_family import cyclically_reduced_lyndon_words
    sage: [tuple(w) for w in cyclically_reduced_lyndon_words(2, 2)]
    [(0, 2), (0, 3), (1, 2), (1, 3)]
    sage: [tuple(w) for w in cyclically_reduced_lyndon_words(2, 2, up_to_inverse=True)]
    [(0, 2), (0, 3)]

The functions :func:`lyndon_words` and :func:`cyclically_reduced_lyndon_words`
build any of these; the classes behind them are also usable directly.
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

from cpython cimport array
from libc.limits cimport LLONG_MAX, LLONG_MIN

from combisurf.word_family cimport WordFamily, ByLengthWordFamily


cdef inline int _moebius(int m) noexcept nogil:
    r"""
    Return the Moebius function of the positive integer ``m``.
    """
    cdef int p = 2, s = 1

    while p <= m // p:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            s = -s
        p += 1
    if m > 1:
        s = -s
    return s


cdef inline long long _pow_checked(long long base, int e) noexcept nogil:
    r"""
    Return ``base ** e`` for ``base >= 1`` and ``e >= 0``, or ``-1`` when the
    result does not fit in a ``long long``.
    """
    cdef long long r = 1
    cdef int _k

    for _k in range(e):
        if r > LLONG_MAX // base:
            return -1
        r *= base
    return r


cdef inline long long _lyndon_cardinality(int n, int length) noexcept nogil:
    r"""
    Return the number of Lyndon words of length ``length`` on ``n`` letters,
    or ``-1`` when the count does not fit in a ``long long``.
    """
    cdef int d, mu
    cdef long long c, total = 0

    if length == 0:
        # the empty word
        return 1
    if n == 0:
        # an empty alphabet has no word of positive length
        return 0

    for d in range(1, length + 1):
        if length % d:
            continue
        mu = _moebius(length // d)
        if mu == 0:
            continue

        # the number of words of length d
        c = _pow_checked(n, d)
        if c < 0:
            return -1

        if mu > 0:
            if total > LLONG_MAX - c:
                return -1
            total += c
        else:
            if total < LLONG_MIN + c:
                return -1
            total -= c

    return total // length


cdef inline long long _cyclically_reduced_cardinality(int n, int length) noexcept nogil:
    r"""
    Return the number of cyclically reduced Lyndon words of length ``length``
    on ``n`` generators, or ``-1`` when it does not fit in a ``long long``.
    """
    cdef int d, mu
    cdef long long base, c, term, total = 0

    if length == 0:
        # the empty word
        return 1
    if n == 0:
        # an empty alphabet has no word of positive length; the formula below
        # would be correct as well, but the power of a negative base is not
        # what _pow_checked expects
        return 0

    base = 2 * n - 1
    for d in range(1, length + 1):
        if length % d:
            continue
        mu = _moebius(length // d)
        if mu == 0:
            continue

        # the number of cyclically reduced words of length d
        c = _pow_checked(base, d)
        if c < 0:
            return -1
        term = 1 if d % 2 else 1 + 2 * (n - 1)
        if c > LLONG_MAX - term:
            return -1
        term += c

        if mu > 0:
            if total > LLONG_MAX - term:
                return -1
            total += term
        else:
            if total < LLONG_MIN + term:
                return -1
            total -= term

    return total // length


cdef class LyndonWordFamily(WordFamily):
    r"""
    Family of the Lyndon words of given length on a given alphabet.

    A Lyndon word is a word that is strictly smaller than each of its proper
    rotations. The words are produced in lexicographic order, by Duval's
    algorithm.

    INPUT:

    - ``n`` -- non-negative integer; the size of the alphabet, whose letters
      are ``0``, ..., ``n - 1``

    - ``length`` -- non-negative integer; the length of the words

    EXAMPLES::

        sage: from combisurf.lyndon_word_family import LyndonWordFamily

        sage: for w in LyndonWordFamily(2, 3):
        ....:     print(w)
        array('i', [0, 0, 1])
        array('i', [0, 1, 1])

        sage: [LyndonWordFamily(2, length).count() for length in range(12)]
        [1, 2, 1, 2, 3, 6, 9, 18, 30, 56, 99, 186]
        sage: [sum(moebius(length // d) * 2**d for d in divisors(length)) // length for length in range(1,10)]
        [2, 1, 2, 3, 6, 9, 18, 30, 56]

        sage: [LyndonWordFamily(3, length).count() for length in range(10)]
        [1, 3, 3, 8, 18, 48, 116, 312, 810, 2184]
        sage: [sum(moebius(length // d) * 3**d for d in divisors(length)) // length for length in range(1,10)]
        [3, 3, 8, 18, 48, 116, 312, 810, 2184]

    The empty word is the only word of length zero, and there is no Lyndon
    word of length at least two on a single letter::

        sage: LyndonWordFamily(5, 0).list()
        [array('i')]
        sage: LyndonWordFamily(1, 1).list()
        [array('i', [0])]
        sage: LyndonWordFamily(1, 7).list()
        []
        sage: LyndonWordFamily(0, 3).list()
        []

    The family knows its own size, so that :meth:`~WordFamily.count` does
    not run through it::

        sage: LyndonWordFamily(2, 60).count()
        19215358392200893

    The number of Lyndon words of length `\ell` on `n` letters is
    `\frac{1}{\ell} \sum_{d | \ell} \mu(\ell/d) n^d`. The computation is
    made with machine integers, and gives up, leaving
    :meth:`~WordFamily.count` to run through the family, once `n^{\ell}`
    no longer fits in one.

    TESTS:

    The comparisons with the formula above are circular, since that formula
    is what :meth:`~WordFamily.count` evaluates; here it is checked against
    running through the family instead::

        sage: for n in range(5):
        ....:     for length in range(9):
        ....:         it = LyndonWordFamily(n, length)
        ....:         assert it.count() == len(it.list()), (n, length)

    ::

        sage: LyndonWordFamily(-1, 3)
        Traceback (most recent call last):
        ...
        ValueError: n (=-1) must be non-negative
        sage: LyndonWordFamily(2, -1)
        Traceback (most recent call last):
        ...
        ValueError: length (=-1) must be non-negative
    """
    def __init__(self, int n, int length):
        r"""
        TESTS::

            sage: from combisurf.lyndon_word_family import LyndonWordFamily
            sage: LyndonWordFamily(2, 5).count()
            6
        """
        if n < 0:
            raise ValueError(f"n (={n}) must be non-negative")
        if length < 0:
            raise ValueError(f"length (={length}) must be non-negative")
        self.n = n
        self._alloc(length)

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from combisurf.lyndon_word_family import LyndonWordFamily
            sage: LyndonWordFamily(3, 5)
            Lyndon words of length 5 on 3 letters
        """
        return f"Lyndon words of length {self.size} on {self.n} letters"

    cdef WordFamily _new(self):
        return LyndonWordFamily(self.n, self.size)

    cdef long long cardinality(self) noexcept nogil:
        return _lyndon_cardinality(self.n, self.size)

    cdef bint first(self) noexcept nogil:
        cdef int k

        for k in range(self.size):
            self.ww[k] = 0
        self.i = -1

        if self.size == 0:
            # the empty word
            return 1
        if self.n <= 0 or (self.n == 1 and self.size >= 2):
            # an empty alphabet has no letter and a single letter a^k is a
            # proper power as soon as k >= 2
            return 0
        if self.size == 1:
            # so that the increment in next() brings the first letter to 0
            self.ww[0] = -1

        self.i = self.size - 1
        return self.next()

    cdef bint next(self) noexcept nogil:
        cdef int i = self.i, j, length = self.size, n = self.n
        cdef int *w = self.ww
        cdef bint done

        if i == -1:
            return 0

        while True:
            # increment at position i and repeat the prefix w[:i+1] periodically
            w[i] += 1
            for j in range(1, length - i):
                w[j + i] = w[j - 1]

            # w is a Lyndon word exactly when its period is its whole length
            done = (i == length - 1)

            # the next position to increment is the last non-maximal letter
            i = length - 1
            while i >= 0 and w[i] == n - 1:
                i -= 1

            if done:
                self.i = i
                return 1
            if i == -1:
                self.i = -1
                return 0


cdef inline int _cyclically_reduced_bump(int *w, int i, int length, int top) noexcept nogil:
    r"""
    Return the smallest letter allowed at position ``i`` of ``w`` that is
    larger than ``w[i]``, or ``-1`` when there is none.

    A letter is allowed at position ``i`` when it is not the inverse of the
    letter before it, and, at the last position, not the inverse of the first
    letter either.
    """
    cdef int c = w[i] + 1
    cdef int f1 = (w[i - 1] ^ 1) if i > 0 else -1
    cdef int f2 = (w[0] ^ 1) if (i == length - 1 and i > 0) else -1

    while c <= top:
        if c != f1 and c != f2:
            return c
        c += 1
    return -1


cdef class CyclicallyReducedLyndonWordFamily(WordFamily):
    r"""
    Family of the cyclically reduced Lyndon words of given length on a given
    number of generators.

    The alphabet has ``2n`` letters, seen as the generators of a free group
    and their inverses, the letters ``i`` and ``i ^ 1`` being inverse of each
    other. A word is cyclically reduced when no letter is followed by its
    inverse, the last letter being followed by the first one. The words are
    produced in lexicographic order.

    Together these words are the Lyndon representatives of the conjugacy
    classes of the primitive elements of the free group.

    INPUT:

    - ``n`` -- non-negative integer; the number of generators, so that the
      alphabet is ``0``, ..., ``2n - 1``

    - ``length`` -- non-negative integer; the length of the words

    EXAMPLES::

        sage: from combisurf.lyndon_word_family import CyclicallyReducedLyndonWordFamily

        sage: for w in CyclicallyReducedLyndonWordFamily(2, 3):
        ....:     print(w)
        array('i', [0, 0, 2])
        array('i', [0, 0, 3])
        array('i', [0, 2, 2])
        array('i', [0, 3, 3])
        array('i', [1, 1, 2])
        array('i', [1, 1, 3])
        array('i', [1, 2, 2])
        array('i', [1, 3, 3])

        sage: [CyclicallyReducedLyndonWordFamily(2, length).count() for length in range(13)]
        [1, 4, 4, 8, 18, 48, 116, 312, 810, 2184, 5880, 16104, 44220]
        sage: [CyclicallyReducedLyndonWordFamily(3, length).count() for length in range(11)]
        [1, 6, 12, 40, 150, 624, 2580, 11160, 48750, 217000, 976248]

    These are exactly the reduced Lyndon words that are cyclically reduced,
    but they are produced without running through the others::

        sage: from combisurf.lyndon_word_family import LyndonWordFamily
        sage: from combisurf.word import word_is_cyclically_reduced
        sage: for length in range(9):
        ....:     naive = [w for w in LyndonWordFamily(4, length) if word_is_cyclically_reduced(w)]
        ....:     assert naive == CyclicallyReducedLyndonWordFamily(2, length).list()

    The empty word is the only word of length zero, and a single generator
    gives nothing beyond length one::

        sage: CyclicallyReducedLyndonWordFamily(5, 0).list()
        [array('i')]
        sage: CyclicallyReducedLyndonWordFamily(1, 1).list()
        [array('i', [0]), array('i', [1])]
        sage: CyclicallyReducedLyndonWordFamily(1, 7).list()
        []
        sage: CyclicallyReducedLyndonWordFamily(0, 3).list()
        []

    ALGORITHM:

    Duval's algorithm with two conditions on the letters: none is the
    inverse of the one before it, and the closing one is not the inverse of
    the first one either. Being cyclically reduced adds the
    single condition that the last letter is not the inverse of the first
    one, which is checked nowhere else than where the last letter is chosen,
    so no word outside the family is ever built.

    Note that the condition is vacuous for half of the words: a Lyndon word
    of length at least two starts with the letter that is strictly smallest
    in it, so its last letter is strictly larger than its first one, and when
    the first letter is odd its inverse is the letter below it, which
    therefore cannot occur at the end. The constraint only bites when the
    first letter is even, and then it forbids exactly one letter.

    Unlike the merely reduced words, the cyclically reduced ones are stable
    under rotation, which is what the usual Moebius inversion counting Lyndon
    words among the words of a family needs. The family therefore knows its
    own size, so that :meth:`~WordFamily.count` does not run through it::

        sage: CyclicallyReducedLyndonWordFamily(2, 39).count()
        103911670590189280

    The number of cyclically reduced words of length `d` on `2n` letters is
    `(2n-1)^d + 1 + (n-1)(1 + (-1)^d)`, and the number of cyclically reduced
    Lyndon words of length `\ell` is `\sum_{d | \ell} \mu(\ell/d)` times
    that, divided by `\ell`. The computation is made with machine integers,
    and gives up, leaving :meth:`~WordFamily.count` to run through the
    family, once `(2n-1)^{\ell}` no longer fits in one.

    TESTS::

        sage: for n in range(4):
        ....:     for length in range(9):
        ....:         it = CyclicallyReducedLyndonWordFamily(n, length)
        ....:         assert it.count() == len(it.list()), (n, length)

        sage: CyclicallyReducedLyndonWordFamily(-1, 3)
        Traceback (most recent call last):
        ...
        ValueError: n (=-1) must be non-negative
        sage: CyclicallyReducedLyndonWordFamily(2, -1)
        Traceback (most recent call last):
        ...
        ValueError: length (=-1) must be non-negative
    """
    def __init__(self, int n, int length):
        r"""
        TESTS::

            sage: from combisurf.lyndon_word_family import CyclicallyReducedLyndonWordFamily
            sage: CyclicallyReducedLyndonWordFamily(2, 5).count()
            48
        """
        if n < 0:
            raise ValueError(f"n (={n}) must be non-negative")
        if length < 0:
            raise ValueError(f"length (={length}) must be non-negative")
        self.n = n
        self._alloc(length)

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from combisurf.lyndon_word_family import CyclicallyReducedLyndonWordFamily
            sage: CyclicallyReducedLyndonWordFamily(3, 5)
            cyclically reduced Lyndon words of length 5 on 3 generators
        """
        return f"cyclically reduced Lyndon words of length {self.size} on {self.n} generators"

    cdef WordFamily _new(self):
        return CyclicallyReducedLyndonWordFamily(self.n, self.size)

    cdef long long cardinality(self) noexcept nogil:
        return _cyclically_reduced_cardinality(self.n, self.size)

    cdef bint first(self) noexcept nogil:
        cdef int k

        for k in range(self.size):
            self.ww[k] = 0
        self.i = -1

        if self.size == 0:
            # the empty word
            return 1
        if self.n <= 0 or (self.n == 1 and self.size >= 2):
            # without generator there is no letter, and on the two letters 0
            # and 1, which are inverse of each other, a Lyndon word of length
            # at least two uses both of them, hence is not reduced
            return 0
        if self.size == 1:
            # so that the increment in next() brings the first letter to 0
            self.ww[0] = -1

        self.i = self.size - 1
        return self.next()

    cdef int next_at(self, int p) noexcept nogil:
        cdef int i, length = self.size, top = 2 * self.n - 1

        if self.i == -1:
            return 0
        if p <= 0:
            # every word shares the empty prefix
            self.i = -1
            return 0
        if p > length:
            p = length

        # forcing the increment at p - 1 or before is what steps over the
        # words that keep the first p letters; when the increment was going
        # to happen before p - 1 anyway there is nothing to force
        i = p - 1
        if i > self.i:
            i = self.i
        while i >= 0 and _cyclically_reduced_bump(self.ww, i, length, top) == -1:
            i -= 1
        self.i = i
        return self.next()

    cdef bint next(self) noexcept nogil:
        cdef int i = self.i, j, length = self.size, top = 2 * self.n - 1, c
        cdef int *w = self.ww
        cdef bint done

        if i == -1:
            return 0

        while True:
            # increment at position i; the position was chosen so that this
            # letter exists
            w[i] = _cyclically_reduced_bump(w, i, length, top)

            # repeating the prefix would put w[0] right after w[i]; when that
            # is not reduced, lengthen the prefix with the smallest letter
            # that is allowed there instead
            while i <= length - 2 and w[i] ^ 1 == w[0]:
                i += 1
                c = w[0] + 1
                if i == length - 1 and c == (w[0] ^ 1):
                    # that letter closes the word, so it must not be the
                    # inverse of the first one either
                    c += 1
                w[i] = c

            # repeat the prefix w[:i+1] periodically on w[i+1:]
            for j in range(1, length - i):
                w[j + i] = w[j - 1]

            # w is a Lyndon word exactly when its period is its whole length
            done = (i == length - 1)

            # the next position to increment is the last one holding a letter
            # that can be replaced by a larger allowed one
            i = length - 1
            while i >= 0 and _cyclically_reduced_bump(w, i, length, top) == -1:
                i -= 1

            if done:
                self.i = i
                return 1
            if i == -1:
                self.i = -1
                return 0


cdef class LexWordFamily(WordFamily):
    r"""
    Base class for the families of Lyndon words of a range of lengths taken
    in lexicographic order.

    Words of different lengths are compared as words, so a word comes before
    every longer word it begins. That is what sets these families apart from
    :class:`ByLengthWordFamily`, which takes one length after the other.

    A subclass says which words it holds through two C methods: ``_bump(i,
    c)``, returning the smallest letter allowed at position ``i`` that is at
    least ``c``, or ``-1``, and ``_emit()``, telling whether the word now in
    the buffer belongs to the family.

    ALGORITHM:

    The words of a family of this module are the nodes of the tree of Lyndon
    prefixes whose Lyndon prefix is the whole word, and a pre-order walk of
    that tree visits the nodes in lexicographic order. Each node carries the
    length ``p`` of its Lyndon prefix: the child continuing the period keeps
    it, a child raising a letter above the period gets ``p`` equal to its own
    length, and a node is a Lyndon word exactly when ``p`` is its length.

    Note that the walk has to emit as it descends, and not only after a
    backtracking step as the fixed length algorithm does. When the letter
    that would continue the period is not allowed, raising it to the next
    allowed one makes a Lyndon word on the spot, and such a word would
    otherwise be stepped over.
    """
    def _init_range(self, int n, int top, int size_min, int size_max):
        r"""
        Set the parameters shared by the subclasses.

        TESTS::

            sage: from combisurf.lyndon_word_family import LexLyndonWordFamily
            sage: LexLyndonWordFamily(2, 1, 4).count()
            5
        """
        if n < 0:
            raise ValueError(f"n (={n}) must be non-negative")
        if size_min < 0:
            raise ValueError(f"size_min (={size_min}) must be non-negative")
        if size_max < 0:
            raise ValueError(f"size_max (={size_max}) must be non-negative")

        self.n = n
        self.top = top
        self.size_min = size_min
        self.size_max = size_max
        self._alloc(size_max - 1 if size_max > size_min else 0)
        self.varying_size = True
        self.has_empty_word = (size_min == 0 and size_max > size_min)

    cdef int _bump(self, int i, int c) noexcept nogil:
        return c if c <= self.top else -1

    cdef bint _emit(self) noexcept nogil:
        return 1

    cdef bint first(self) noexcept nogil:
        self.size = 0
        self.p = 0
        self.phase = 0
        return self.next()

    cdef bint next(self) noexcept nogil:
        cdef int c, c0
        cdef bint descended
        cdef int *w = self.ww

        if self.size_max <= self.size_min:
            return 0

        if self.phase == 0:
            # the empty word is a prefix of every other, hence the smallest
            self.phase = 1
            if self.size_min == 0:
                self.size = 0
                return 1

        while True:
            descended = False
            if self.size < self.capacity:
                c0 = 0 if self.size == 0 else w[self.size - self.p]
                c = self._bump(self.size, c0)
                if c >= 0:
                    w[self.size] = c
                    if self.size == 0 or c != c0:
                        # the period is broken, the word is its own prefix
                        self.p = self.size + 1
                    self.size += 1
                    descended = True

            if not descended:
                # no child, so move on to the next brother, going up through
                # the letters that cannot be raised any further
                while True:
                    if self.size == 0:
                        return 0
                    c = self._bump(self.size - 1, w[self.size - 1] + 1)
                    if c >= 0:
                        w[self.size - 1] = c
                        self.p = self.size
                        break
                    self.size -= 1

            if self.p == self.size and self.size >= self.size_min and self._emit():
                return 1


cdef class LexLyndonWordFamily(LexWordFamily):
    r"""
    Family of the Lyndon words whose length lies in a range, in
    lexicographic order.

    INPUT:

    - ``n`` -- non-negative integer; the size of the alphabet

    - ``size_min`` -- non-negative integer; the smallest length

    - ``size_max`` -- non-negative integer; one past the largest length

    EXAMPLES::

        sage: from combisurf.lyndon_word_family import LexLyndonWordFamily
        sage: F = LexLyndonWordFamily(2, 1, 4)
        sage: F
        Lyndon words of length 1 <= n < 4 on 2 letters, in lexicographic order
        sage: F.list()
        [array('i', [0]), array('i', [0, 0, 1]), array('i', [0, 1]),
         array('i', [0, 1, 1]), array('i', [1])]

    The same words taken length by length come out differently ordered::

        sage: from combisurf.word_family import ByLengthWordFamily
        sage: from combisurf.lyndon_word_family import LyndonWordFamily
        sage: ByLengthWordFamily(LyndonWordFamily(2, 0), 1, 4).list()
        [array('i', [0]), array('i', [1]), array('i', [0, 1]),
         array('i', [0, 0, 1]), array('i', [0, 1, 1])]
        sage: sorted(tuple(w) for w in F) == sorted(
        ....:     tuple(w) for w in ByLengthWordFamily(LyndonWordFamily(2, 0), 1, 4))
        True

    The empty word is smaller than every other, so it comes first when the
    range holds it::

        sage: LexLyndonWordFamily(2, 0, 3).list()
        [array('i'), array('i', [0]), array('i', [0, 1]), array('i', [1])]

    TESTS::

        sage: LexLyndonWordFamily(2, 3, 3).list()
        []
        sage: LexLyndonWordFamily(2, -1, 3)
        Traceback (most recent call last):
        ...
        ValueError: size_min (=-1) must be non-negative
    """
    def __init__(self, int n, int size_min, int size_max):
        r"""
        TESTS::

            sage: from combisurf.lyndon_word_family import LexLyndonWordFamily
            sage: LexLyndonWordFamily(3, 0, 5).count()
            33
        """
        self._init_range(n, n - 1, size_min, size_max)

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from combisurf.lyndon_word_family import LexLyndonWordFamily
            sage: LexLyndonWordFamily(3, 2, 7)
            Lyndon words of length 2 <= n < 7 on 3 letters, in lexicographic order
        """
        return (f"Lyndon words of length {self.size_min} <= n < {self.size_max} "
                f"on {self.n} letters, in lexicographic order")

    cdef WordFamily _new(self):
        return LexLyndonWordFamily(self.n, self.size_min, self.size_max)

    cdef long long cardinality(self) noexcept nogil:
        cdef long long total = 0, c
        cdef int length

        for length in range(self.size_min, self.size_max):
            c = _lyndon_cardinality(self.n, length)
            if c < 0:
                return -1
            if total > LLONG_MAX - c:
                return -1
            total += c
        return total


cdef class LexCyclicallyReducedLyndonWordFamily(LexWordFamily):
    r"""
    Family of the cyclically reduced Lyndon words whose length lies in a
    range, in lexicographic order.

    INPUT:

    - ``n`` -- non-negative integer; the number of generators, so that the
      alphabet is ``0``, ..., ``2n - 1``

    - ``size_min`` -- non-negative integer; the smallest length

    - ``size_max`` -- non-negative integer; one past the largest length

    EXAMPLES::

        sage: from combisurf.lyndon_word_family import LexCyclicallyReducedLyndonWordFamily
        sage: LexCyclicallyReducedLyndonWordFamily(2, 1, 4).list()
        [array('i', [0]), array('i', [0, 0, 2]), array('i', [0, 0, 3]),
         array('i', [0, 2]), array('i', [0, 2, 2]), array('i', [0, 3]),
         array('i', [0, 3, 3]), array('i', [1]), array('i', [1, 1, 2]),
         array('i', [1, 1, 3]), array('i', [1, 2]), array('i', [1, 2, 2]),
         array('i', [1, 3]), array('i', [1, 3, 3]), array('i', [2]),
         array('i', [3])]

        sage: from combisurf.lyndon_word_family import CyclicallyReducedLyndonWordFamily
        sage: sorted(tuple(w) for w in LexCyclicallyReducedLyndonWordFamily(2, 1, 6)) == sorted(
        ....:     tuple(w) for length in range(1, 6)
        ....:     for w in CyclicallyReducedLyndonWordFamily(2, length))
        True

    .. NOTE::

        Unlike :class:`CyclicallyReducedLyndonWordFamily`, this family cannot
        prune on the closing letter. Over a range of lengths every prefix is
        itself a word of the family to be, and a prefix whose last letter is
        the inverse of its first one still extends to words that are
        cyclically reduced, ``(0, 2, 1)`` to ``(0, 2, 1, 3)`` for instance.
        The condition is therefore tested when a word is emitted rather than
        used to cut the walk short.
    """
    def __init__(self, int n, int size_min, int size_max):
        r"""
        TESTS::

            sage: from combisurf.lyndon_word_family import LexCyclicallyReducedLyndonWordFamily
            sage: LexCyclicallyReducedLyndonWordFamily(2, 0, 5).count()
            35
        """
        self._init_range(n, 2 * n - 1, size_min, size_max)

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from combisurf.lyndon_word_family import LexCyclicallyReducedLyndonWordFamily
            sage: LexCyclicallyReducedLyndonWordFamily(3, 2, 7)
            cyclically reduced Lyndon words of length 2 <= n < 7 on 3 generators, in lexicographic order
        """
        return (f"cyclically reduced Lyndon words of length {self.size_min} <= n < "
                f"{self.size_max} on {self.n} generators, in lexicographic order")

    cdef WordFamily _new(self):
        return LexCyclicallyReducedLyndonWordFamily(self.n, self.size_min, self.size_max)

    cdef long long cardinality(self) noexcept nogil:
        cdef long long total = 0, c
        cdef int length

        for length in range(self.size_min, self.size_max):
            c = _cyclically_reduced_cardinality(self.n, length)
            if c < 0:
                return -1
            if total > LLONG_MAX - c:
                return -1
            total += c
        return total

    cdef int _bump(self, int i, int c) noexcept nogil:
        cdef int f

        if i > 0:
            f = self.ww[i - 1] ^ 1
            while c <= self.top and c == f:
                c += 1
        return c if c <= self.top else -1

    cdef bint _emit(self) noexcept nogil:
        return self.size <= 1 or (self.ww[self.size - 1] ^ 1) != self.ww[0]


cdef bint _is_cyclically_reduced(family):
    r"""
    Return whether every word of ``family`` is cyclically reduced.

    Only the families of this module are, so a type test answers it; running
    it up through :class:`~combisurf.word_family.ByLengthWordFamily` is what
    lets a range of lengths be filtered as well.
    """
    cdef ByLengthWordFamily by_length

    if isinstance(family, (CyclicallyReducedLyndonWordFamily,
                           LexCyclicallyReducedLyndonWordFamily)):
        return True
    if isinstance(family, ByLengthWordFamily):
        by_length = <ByLengthWordFamily>family
        return _is_cyclically_reduced(by_length.base)
    return False


cdef inline int _least_rotation(int *u, int length) noexcept nogil:
    r"""
    Return the starting index of the lexicographically least rotation of the
    word ``u`` of length ``length``, by Duval's algorithm.
    """
    cdef int i = 0, j = 1, k = 0
    cdef int a, b, x, y

    # the bound on k, rather than on i + k and j + k, is what lets the
    # comparison of two rotations run past the end of the word
    while i < length and j < length and k < length:
        # i + k and j + k stay below twice the length, so one subtraction
        # does what a remainder would
        x = i + k
        if x >= length:
            x -= length
        y = j + k
        if y >= length:
            y -= length
        a = u[x]
        b = u[y]
        if a == b:
            k += 1
            continue
        if a > b:
            i = i + k + 1
        else:
            j = j + k + 1
        if i == j:
            j += 1
        k = 0

    return i if i < j else j


cdef class UpToInverseWordFamily(WordFamily):
    r"""
    Family of the cyclically reduced Lyndon words that are smaller than the
    inverse of their own conjugacy class.

    The words are read as free group elements, the letters ``i`` and
    ``i ^ 1`` being inverse of each other. Writing `\bar{w}` for the Lyndon
    word of the conjugacy class of `w^{-1}`, that is the lexicographically
    least rotation of the reversed and inverted word, a word ``w`` is kept
    exactly when ``w`` is lexicographically smaller than `\bar{w}`.

    Since the cyclically reduced Lyndon words are the Lyndon representatives
    of the conjugacy classes of the primitive elements of the free group,
    this runs through those classes up to inversion.

    INPUT:

    - ``base`` -- a family of cyclically reduced Lyndon words

    EXAMPLES::

        sage: from combisurf.lyndon_word_family import UpToInverseWordFamily
        sage: from combisurf.lyndon_word_family import CyclicallyReducedLyndonWordFamily

        sage: UpToInverseWordFamily(CyclicallyReducedLyndonWordFamily(2, 1)).list()
        [array('i', [0]), array('i', [2])]

        sage: for w in UpToInverseWordFamily(CyclicallyReducedLyndonWordFamily(2, 3)):
        ....:     print(w)
        array('i', [0, 0, 2])
        array('i', [0, 0, 3])
        array('i', [0, 2, 2])
        array('i', [0, 3, 3])

    This splits the family in two exactly. No element of a free group is
    conjugate to its own inverse, so `w \mapsto \bar{w}` is an involution
    without fixed point on the cyclically reduced Lyndon words, and taking
    the smaller of each pair keeps one word out of two::

        sage: for n in range(1, 4):
        ....:     for length in range(1, 9):
        ....:         it = CyclicallyReducedLyndonWordFamily(n, length)
        ....:         assert 2 * UpToInverseWordFamily(it).count() == it.count()

    That halving is where :meth:`~WordFamily.count` gets its answer, so it
    inherits the closed form of the underlying family and does not run
    through anything::

        sage: UpToInverseWordFamily(CyclicallyReducedLyndonWordFamily(3, 24)).count()
        1241763427726250

    A family whose words are not all cyclically reduced is refused, since
    the criterion is meaningless on a family that `w \mapsto \bar{w}` does
    not send to itself, or that it fixes a word of::

        sage: from combisurf.lyndon_word_family import LyndonWordFamily
        sage: UpToInverseWordFamily(LyndonWordFamily(4, 4))
        Traceback (most recent call last):
        ...
        ValueError: base must be a family of cyclically reduced Lyndon words

    ALGORITHM:

    Let ``a`` be the first letter of ``w``, which for a Lyndon word is the
    smallest letter occurring in it. The first letter of `\bar{w}` is the
    smallest letter of `w^{-1}`, which is readily read off ``w``, and settles
    most of the comparison without computing `\bar{w}` at all:

    - when ``a`` is odd, its inverse ``a - 1`` occurs in `w^{-1}`, and no
      letter of `w^{-1}` is smaller, so `\bar{w}` starts below ``w`` and the
      word is dropped. Being a condition on the first letter alone, it holds
      for the whole block of consecutive words starting with ``a``, which is
      therefore stepped over at once through ``next_at(1)`` rather than
      generated and rejected word by word;

    - when ``a`` is even and its inverse ``a + 1`` does not occur in ``w``,
      every letter of `w^{-1}` is at least ``a + 1``, so `\bar{w}` starts
      above ``w`` and the word is kept;

    - otherwise both start with ``a``, and `\bar{w}` is computed by Duval's
      least rotation algorithm and compared with ``w``.

    Only the first of these is a condition on a prefix, so it is the only one
    that prunes; it keeps between 7 and 25 percent of the family from being
    generated at all, over the range of parameters that can be run through.

    TESTS:

    The empty word is its own inverse, so it is not kept either::

        sage: UpToInverseWordFamily(CyclicallyReducedLyndonWordFamily(2, 0)).list()
        []
        sage: UpToInverseWordFamily(CyclicallyReducedLyndonWordFamily(2, 0)).count()
        0

    The count is checked against running through the family::

        sage: for n in range(4):
        ....:     for length in range(9):
        ....:         it = UpToInverseWordFamily(CyclicallyReducedLyndonWordFamily(n, length))
        ....:         assert it.count() == len(it.list()), (n, length)
    """
    def __init__(self, WordFamily base):
        r"""
        TESTS::

            sage: from combisurf.lyndon_word_family import UpToInverseWordFamily
            sage: from combisurf.lyndon_word_family import CyclicallyReducedLyndonWordFamily
            sage: UpToInverseWordFamily(CyclicallyReducedLyndonWordFamily(2, 4)).count()
            9
            sage: UpToInverseWordFamily(None)
            Traceback (most recent call last):
            ...
            ValueError: base must be a family of cyclically reduced Lyndon words
        """
        if not _is_cyclically_reduced(base):
            raise ValueError("base must be a family of cyclically reduced Lyndon words")
        self.base = base
        # share the buffer of the underlying family, there is nothing to copy
        self.w = base.w
        self.ww = base.ww
        self.capacity = base.capacity
        self.size = base.size
        self.u = array.array('i', [0] * base.capacity)
        self.uu = self.u.data.as_ints
        self.varying_size = base.varying_size
        self.has_empty_word = False

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from combisurf.lyndon_word_family import UpToInverseWordFamily
            sage: from combisurf.lyndon_word_family import CyclicallyReducedLyndonWordFamily
            sage: UpToInverseWordFamily(CyclicallyReducedLyndonWordFamily(3, 5))
            cyclically reduced Lyndon words of length 5 on 3 generators, up to inverse
        """
        return f"{self.base!r}, up to inverse"

    cdef WordFamily _new(self):
        return UpToInverseWordFamily(self.base._new())

    cdef int _resize(self, int capacity) except -1:
        self.base._resize(capacity)
        # the buffer is the one of the underlying family, keep sharing it
        self.w = self.base.w
        self.ww = self.base.ww
        self.capacity = self.base.capacity
        self.size = self.base.size
        self.out = None
        self.out_size = -1
        self.u = array.array('i', [0] * self.capacity)
        self.uu = self.u.data.as_ints
        return 0

    cdef void _set_size(self, int size) noexcept nogil:
        self.base._set_size(size)
        self.size = size

    cdef bint _keep(self) noexcept nogil:
        cdef int *w = self.ww
        cdef int *u = self.uu
        cdef int length = self.base.size
        cdef int a, k, s, x, y

        if length == 0:
            # the empty word is its own inverse
            return 0

        a = w[0]
        if a & 1:
            # the inverse of a is a - 1, the smallest letter of the inverse
            # word, so the latter starts below w
            return 0

        # does the inverse of a occur in w?
        for k in range(1, length):
            if w[k] == a + 1:
                break
        else:
            # every letter of the inverse word is at least a + 1
            return 1

        # both start with a, so compare in full
        for k in range(length):
            u[k] = w[length - 1 - k] ^ 1
        s = _least_rotation(u, length)
        for k in range(length):
            x = u[s]
            y = w[k]
            if x != y:
                return x > y
            s += 1
            if s == length:
                s = 0
        return 0

    cdef long long cardinality(self) noexcept nogil:
        cdef long long c = self.base.cardinality()

        if c < 0:
            return -1

        if self.base.has_empty_word:
            # the empty word is its own inverse and is dropped
            c -= 1

        # on the rest the map is an involution without fixed point, so it
        # pairs the words two by two and exactly half of them are kept
        return c // 2

    cdef bint _advance(self) noexcept nogil:
        r"""
        Move the underlying family to the next word that may be kept.

        A word whose first letter is odd is dropped, and so is every word
        that follows it with the same first letter, so the whole block is
        stepped over at once.
        """
        cdef int r

        if self.base.size and (self.ww[0] & 1):
            r = self.base.next_at(1)
            if r >= 0:
                return r
        return self.base.next()

    cdef bint first(self) noexcept nogil:
        if not self.base.first():
            return 0
        self.size = self.base.size
        if self._keep():
            return 1
        return self.next()

    cdef bint next(self) noexcept nogil:
        while self._advance():
            self.size = self.base.size
            if self._keep():
                return 1
        return 0


def _range_family(kind, int n, size_min, size_max, order):
    r"""
    Build the family of ``kind`` words of the given length or range.

    TESTS::

        sage: from combisurf.lyndon_word_family import lyndon_words
        sage: lyndon_words(2, 1, 4, order="nonsense")
        Traceback (most recent call last):
        ...
        ValueError: order must be "length" or "lex", not 'nonsense'
    """
    fixed, lex = kind

    if size_max is None:
        # a single length, where the two orders agree
        return fixed(n, size_min)
    if order == "length":
        return ByLengthWordFamily(fixed(n, 0), size_min, size_max)
    if order == "lex":
        return lex(n, size_min, size_max)
    raise ValueError(f'order must be "length" or "lex", not {order!r}')


def lyndon_words(int n, size_min, size_max=None, order="length"):
    r"""
    Return the family of Lyndon words on ``n`` letters, of a given length or
    of a range of lengths.

    INPUT:

    - ``n`` -- non-negative integer; the size of the alphabet

    - ``size_min`` -- non-negative integer; the length of the words, or the
      smallest of them when ``size_max`` is given

    - ``size_max`` -- non-negative integer (default: ``None``); one past the
      largest length, as for ``range``

    - ``order`` -- ``"length"`` (default) or ``"lex"``; whether the words are
      taken one length after the other, or in lexicographic order, in which a
      word comes before every longer word it begins. The two agree on a
      single length, so this is only read when ``size_max`` is given

    EXAMPLES::

        sage: from combisurf.lyndon_word_family import lyndon_words
        sage: for w in lyndon_words(2, 3):
        ....:     print(w)
        array('i', [0, 0, 1])
        array('i', [0, 1, 1])

    Over a range of lengths, the two orders hold the same words::

        sage: [tuple(w) for w in lyndon_words(2, 1, 4)]
        [(0,), (1,), (0, 1), (0, 0, 1), (0, 1, 1)]
        sage: [tuple(w) for w in lyndon_words(2, 1, 4, order="lex")]
        [(0,), (0, 0, 1), (0, 1), (0, 1, 1), (1,)]
        sage: (lyndon_words(2, 1, 4).count(), lyndon_words(2, 1, 4, order="lex").count())
        (5, 5)
    """
    return _range_family((LyndonWordFamily, LexLyndonWordFamily),
                         n, size_min, size_max, order)


def cyclically_reduced_lyndon_words(int n, size_min, size_max=None,
                                    order="length", bint up_to_inverse=False):
    r"""
    Return the family of cyclically reduced Lyndon words on ``n`` generators,
    that is on the ``2n`` letters ``0``, ..., ``2n - 1``, of a given length or
    of a range of lengths.

    INPUT:

    - ``n`` -- non-negative integer; the number of generators

    - ``size_min`` -- non-negative integer; the length of the words, or the
      smallest of them when ``size_max`` is given

    - ``size_max`` -- non-negative integer (default: ``None``); one past the
      largest length

    - ``order`` -- ``"length"`` (default) or ``"lex"``; see
      :func:`lyndon_words`

    - ``up_to_inverse`` -- boolean (default: ``False``); whether to keep only
      one word out of each pair of a word and the Lyndon word of the
      conjugacy class of its inverse

    EXAMPLES::

        sage: from combisurf.lyndon_word_family import cyclically_reduced_lyndon_words
        sage: for w in cyclically_reduced_lyndon_words(2, 2):
        ....:     print(w)
        array('i', [0, 2])
        array('i', [0, 3])
        array('i', [1, 2])
        array('i', [1, 3])

    With ``up_to_inverse`` the words left out are the ones whose inverse is
    kept, here ``(1, 2)`` for ``(0, 3)`` and ``(1, 3)`` for ``(0, 2)``::

        sage: for w in cyclically_reduced_lyndon_words(2, 2, up_to_inverse=True):
        ....:     print(w)
        array('i', [0, 2])
        array('i', [0, 3])

    Over a range of lengths::

        sage: [tuple(w) for w in cyclically_reduced_lyndon_words(2, 1, 4)]
        [(0,), (1,), (2,), (3,), (0, 2), (0, 3), (1, 2), (1, 3),
         (0, 0, 2), (0, 0, 3), (0, 2, 2), (0, 3, 3), (1, 1, 2), (1, 1, 3), (1, 2, 2), (1, 3, 3)]
        sage: [tuple(w) for w in cyclically_reduced_lyndon_words(2, 1, 4, order="lex")]
        [(0,), (0, 0, 2), (0, 0, 3), (0, 2), (0, 2, 2), (0, 3), (0, 3, 3),
         (1,), (1, 1, 2), (1, 1, 3), (1, 2), (1, 2, 2), (1, 3), (1, 3, 3), (2,), (3,)]

    TESTS:

    Keeping one word out of each inverse pair is available in both orders and
    at a single length::

        sage: for kwds in ({}, {"order": "lex"}):
        ....:     F = cyclically_reduced_lyndon_words(2, 1, 6, up_to_inverse=True, **kwds)
        ....:     G = cyclically_reduced_lyndon_words(2, 1, 6, **kwds)
        ....:     assert 2 * F.count() == G.count()
        ....:     assert sorted(tuple(w) for w in F) == sorted(
        ....:         tuple(w) for length in range(1, 6)
        ....:         for w in cyclically_reduced_lyndon_words(2, length, up_to_inverse=True))
    """
    F = _range_family((CyclicallyReducedLyndonWordFamily,
                       LexCyclicallyReducedLyndonWordFamily),
                      n, size_min, size_max, order)
    if up_to_inverse:
        return UpToInverseWordFamily(F)
    return F
