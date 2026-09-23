r"""
Families of words on non-negative integers

A family is a set of words together with an order on them, so the same words
taken in another order make another family. This module holds what every
family has in common; the families themselves are in
:mod:`combisurf.lyndon_word_family`.

The first half of what follows is what there is to know to run through a
family; the second half is what a new family has to provide.

**Using a family**

A family is not an iterator, no more than a ``range`` is. It is iterated any
number of times, each ``iter`` handing out a :class:`WordFamilyIterator` of
its own, and it also answers :meth:`~WordFamily.list`,
:meth:`~WordFamily.flat` and :meth:`~WordFamily.count` without any of them
disturbing an iteration under way::

    sage: from combisurf.lyndon_word_family import lyndon_words
    sage: f = lyndon_words(2, 3)
    sage: f.list()
    [array('i', [0, 0, 1]), array('i', [0, 1, 1])]
    sage: f.count()
    2
    sage: f.flat()
    array('i', [0, 0, 1, 0, 1, 1])

:meth:`~WordFamily.count` uses a closed form when the family knows one,
rather than running through anything, which is why counting words that could
never be enumerated is instant::

    sage: from combisurf.lyndon_word_family import cyclically_reduced_lyndon_words
    sage: cyclically_reduced_lyndon_words(3, 24).count()
    2483526855452500

:meth:`~WordFamily.flat` lays the words end to end in a single array instead
of building one object per word, and is therefore only available when they
all have the same length.

:class:`ByLengthWordFamily` runs any family at one length after the other,
so that the words of a range of lengths come out shortest first::

    sage: from combisurf.word_family import ByLengthWordFamily
    sage: from combisurf.lyndon_word_family import LyndonWordFamily
    sage: [tuple(w) for w in ByLengthWordFamily(LyndonWordFamily(2, 0), 1, 4)]
    [(0,), (1,), (0, 1), (0, 0, 1), (0, 1, 1)]

Iteration hands out one array per word, and refills the array handed out
last instead of allocating a new one, but only while nothing else holds it,
which its reference count tells; see :meth:`~WordFamily.__next__`. Getting
a run through a family down to that speed is a matter of how the loop is
written, not of the family. The steps below all run through the 179756976
words of ``cyclically_reduced_lyndon_words(4, 11)``, timed on one
machine.

Written at the top level of a session, a plain loop takes 7.0 seconds::

    sage: for w in cyclically_reduced_lyndon_words(2, 6):
    ....:     pass

Dropping the word at the end of the body brings it to 6.1 seconds. Without
the ``del`` the loop variable still refers to the previous word at the
moment the next one is asked for, so that array is not free and a new one is
allocated for every word::

    sage: for w in cyclically_reduced_lyndon_words(2, 6):
    ....:     del w

Moving the very same loop inside a function brings it to 2.6 seconds. At the
top level ``w`` is a global, so binding it and deleting it go through a
dictionary, which costs about as much as the allocation the ``del`` saves;
as a local both are array slots::

    sage: def local_iter(n, length):
    ....:     for w in cyclically_reduced_lyndon_words(n, length):
    ....:         del w
    sage: local_iter(2, 6)

Finally, no Python object need be built at all: a Cython caller can
``cimport`` the class and drive ``first`` and ``next`` itself, which is what
``word_family.pxd`` is for. That takes 1.2 seconds, counting the words as
it goes for good measure. In a session, with the source tree on the include
path so that the declarations are found::

    %%cython --include-dirs /path/to/combisurf

    from combisurf.lyndon_word_family cimport CyclicallyReducedLyndonWordFamily

    def cython_iter(int n, int length):
        cdef CyclicallyReducedLyndonWordFamily f
        cdef long total = 0
        f = CyclicallyReducedLyndonWordFamily(n, length)
        if not f.first():
            return total
        total = 1
        while f.next():
            total += 1
        return total

which agrees with :meth:`~WordFamily.count`, ``cython_iter(2, 6)`` giving
116 and ``cython_iter(4, 11)`` giving 179756976. This last step is shown
rather than doctested, since compiling Cython from a doctest needs
``sage.misc.cython``, which is not available everywhere.

**Adding a family**

A family is a subclass of :class:`WordFamily`. Everything above comes for
free from the C protocol below, declared in ``word_family.pxd``; a subclass
is only ever a few methods long.

The words live in a buffer of C integers, ``ww[0]``, ..., ``ww[size - 1]``,
allocated once by ``_alloc(capacity)`` and never reallocated while the
family is run. ``size`` may change from one word to the next, as it does for
the families of a range of lengths, as long as it stays at most
``capacity``.

Two methods are mandatory and move that buffer:

- ``first()`` sets the buffer to the first word and returns 1, or returns 0
  when the family is empty;

- ``next()`` moves it to the word after the current one and returns 1, or
  returns 0 when there is none. It may only be called when the buffer holds
  a word of the family, that is after a ``first()`` or a ``next()`` that
  returned 1; once either returned 0 the content of the buffer is
  unspecified.

Both are ``noexcept nogil``, which the compiler enforces: they work on
``ww`` and on C attributes, never on a Python object. That is what lets a
Cython caller run through a family at the speed measured above, and what
makes the ``.pxd`` worth cimporting.

One more method is mandatory, and is the only one that builds an object:

- ``_new()`` returns a fresh family of the same words. It is what
  :meth:`~WordFamily.__iter__` hands out, so that iterators are independent
  of one another and of the family they come from. Its default returns the
  empty family, so a subclass that forgets it yields nothing at all.

Two are optional refinements. Each has a default returning ``-1``, which
stands for "not supported by this family"; a caller falls back on running
through the words:

- ``cardinality()`` returns the number of words, or ``-1`` when that number
  does not fit in a ``long long``. :meth:`~WordFamily.count` uses it, and so
  does :meth:`~WordFamily.list`, which then fills a list of the right size
  rather than growing one;

- ``next_at(p)`` is the pruning primitive: it moves the buffer to the next
  word whose first ``p`` letters are not those of the current word, stepping
  over the whole block of words sharing that prefix at once. A family whose
  words are chosen by a condition on a prefix can be run through much faster
  with it.

Two more are only needed by a family that is run at several lengths, or that
wraps another one. :class:`ByLengthWordFamily` calls them, so a family that
has to work underneath it implements them:

- ``_resize(capacity)`` gives the family a buffer of that many letters. The
  default calls ``_alloc``; a family that does not own its buffer, because
  it shares the one of the family it wraps, passes the request on;

- ``_set_size(size)`` says at which length the family is about to be run.
  The default sets ``size``; again a wrapping family passes it on, since
  that is where the words come from.

Two flags say what the words are like, and are read by the machinery rather
than by the algorithm:

- ``varying_size``, whether the words have different lengths, which
  :meth:`~WordFamily.flat` refuses;

- ``has_empty_word``, whether the empty word is one of them.
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

import gc

from cpython cimport array
from cpython.ref cimport PyObject
from libc.limits cimport LLONG_MAX
from libc.string cimport memcpy

cdef extern from "Python.h":
    # declared on a borrowed PyObject * rather than on an object, so that
    # Cython does not create the temporary reference that would make every
    # count one too many
    Py_ssize_t _refcount "Py_REFCNT" (PyObject *o) noexcept nogil


cdef class WordFamily:
    r"""
    Base class for families of words.

    A subclass provides the words by overriding a handful of C methods,
    ``first`` and ``next`` at the very least. They are declared in
    ``word_family.pxd`` and spelled out under "Adding a family" in the
    documentation of this module, :mod:`combisurf.word_family`.

    EXAMPLES::

        sage: from combisurf.lyndon_word_family import LyndonWordFamily
        sage: f = LyndonWordFamily(2, 3)
        sage: f.list()
        [array('i', [0, 0, 1]), array('i', [0, 1, 1])]
        sage: f.count()
        2

    The base class is itself the empty family::

        sage: from combisurf.word_family import WordFamily
        sage: WordFamily().list()
        []
        sage: WordFamily().count()
        0
    """
    def __cinit__(self, *args, **kwds):
        r"""
        Set up an empty buffer.

        The arguments are ignored; every class in the hierarchy receives the
        arguments of the constructor, so they have to be accepted here.

        TESTS::

            sage: from combisurf.lyndon_word_family import LyndonWordFamily
            sage: LyndonWordFamily(2, 3).count()
            2
        """
        self.w = array.array('i', [])
        self.ww = self.w.data.as_ints
        self.capacity = 0
        self.size = 0
        self.out = None
        self.out_size = -1
        self.varying_size = False
        self.has_empty_word = False

    cdef int _alloc(self, int capacity) except -1:
        r"""
        Allocate a buffer of ``capacity`` letters and rewind the family.
        """
        if capacity < 0:
            raise ValueError(f"capacity (={capacity}) must be non-negative")
        self.w = array.array('i', [0] * capacity)
        self.ww = self.w.data.as_ints
        self.capacity = capacity
        self.size = capacity
        self.has_empty_word = (capacity == 0)
        self.out = None
        self.out_size = -1
        self.varying_size = False
        self.has_empty_word = False
        return 0

    cdef int _resize(self, int capacity) except -1:
        r"""
        Give this family a buffer holding ``capacity`` letters, so that it
        can be run at any length up to that.

        A family that does not own its buffer, because it wraps another one,
        has to override this and pass the request on.
        """
        return self._alloc(capacity)

    cdef void _set_size(self, int size) noexcept nogil:
        r"""
        Set the length at which this family is about to be run.

        A family that wraps another one has to override this and pass the
        length on, since that is where the words come from.
        """
        self.size = size

    cdef array.array _copy(self):
        r"""
        Return the current word as an array of its own.

        The array handed out by the previous call is filled again instead of
        a new one being allocated, but only when nothing else holds it any
        more, which its reference count tells. A caller that keeps the words
        it is given, or merely holds the previous one while asking for the
        next, therefore always gets distinct arrays.
        """
        cdef array.array out

        if (self.out is not None and self.out_size == self.size
                and _refcount(<PyObject *>self.out) == 1):
            memcpy(self.out.data.as_ints, self.ww, self.size * sizeof(int))
            return self.out

        out = array.clone(self.w, self.size, False)
        memcpy(out.data.as_ints, self.ww, self.size * sizeof(int))
        self.out = out
        self.out_size = self.size
        return out

    cdef WordFamily _new(self):
        r"""
        Return a new family of the same words, rewound.

        Every subclass has to override this, since this is what makes
        :meth:`__iter__` hand out iterators that are independent of one
        another and of the family they come from.
        """
        return WordFamily()

    cdef bint first(self) noexcept nogil:
        return 0

    cdef bint next(self) noexcept nogil:
        return 0

    cdef int next_at(self, int p) noexcept nogil:
        return -1

    cdef long long cardinality(self) noexcept nogil:
        return -1

    def __iter__(self):
        r"""
        Return a new iterator through this family.

        A family is not itself an iterator, no more than a ``range`` is: it
        is iterated any number of times, and each ``iter`` hands out a
        :class:`WordFamilyIterator` of its own that the others know nothing
        about.

        EXAMPLES::

            sage: from combisurf.lyndon_word_family import LyndonWordFamily
            sage: f = LyndonWordFamily(2, 3)
            sage: [tuple(w) for w in f]
            [(0, 0, 1), (0, 1, 1)]
            sage: [tuple(w) for w in f]
            [(0, 0, 1), (0, 1, 1)]

            sage: it1 = iter(f)
            sage: it2 = iter(f)
            sage: next(it1)
            array('i', [0, 0, 1])
            sage: next(it2)
            array('i', [0, 0, 1])

        A family is not an iterator, so it has no ``__next__`` of its own::

            sage: next(f)
            Traceback (most recent call last):
            ...
            TypeError: 'combisurf.lyndon_word_family.LyndonWordFamily' object is not an iterator
        """
        return WordFamilyIterator(self._new())

    def list(self):
        r"""
        Return the list of the words of this family.

        EXAMPLES::

            sage: from combisurf.lyndon_word_family import LyndonWordFamily
            sage: LyndonWordFamily(3, 2).list()
            [array('i', [0, 1]), array('i', [0, 2]), array('i', [1, 2])]

        .. NOTE::

            The garbage collector is turned off while the list is built and
            restored afterwards. Arrays are tracked by the cyclic collector,
            so building a long list of them otherwise makes the collector
            walk a growing population over and over, which turns a linear
            construction into a superlinear one. Reference counting is not
            affected, hence neither is the memory actually used.
        """
        cdef bint was_enabled = gc.isenabled()
        cdef long long c = self.cardinality()
        cdef Py_ssize_t i = 0
        ans = []
        gc.disable()
        try:
            if not self.first():
                return ans
            if c > 0:
                # the size is known, so fill a list of that size rather than
                # let it grow; the bound keeps a cardinality that claims too
                # much from running past the end
                ans = [None] * c
                ans[0] = self._copy()
                i = 1
                while i < c and self.next():
                    ans[i] = self._copy()
                    i += 1
                if i < c:
                    del ans[i:]
            else:
                ans.append(self._copy())
                while self.next():
                    ans.append(self._copy())
        finally:
            if was_enabled:
                gc.enable()
        return ans

    def flat(self):
        r"""
        Return all the words of this family concatenated in a single array.

        This is the cheap way of getting the whole family at once: it builds
        one object instead of one per word. It is only available for families
        whose words all have the same length.

        EXAMPLES::

            sage: from combisurf.lyndon_word_family import LyndonWordFamily
            sage: LyndonWordFamily(3, 2).flat()
            array('i', [0, 1, 0, 2, 1, 2])

        The words are read back by slicing::

            sage: f = LyndonWordFamily(2, 4).flat()
            sage: [tuple(f[4*i: 4*i + 4]) for i in range(len(f) // 4)]
            [(0, 0, 0, 1), (0, 0, 1, 1), (0, 1, 1, 1)]

        Slicing a :class:`memoryview` of the result reads the words without
        copying them, and ``toreadonly`` makes such a view read only::

            sage: m = memoryview(f).toreadonly()
            sage: [m[4*i: 4*i + 4].tolist() for i in range(len(m) // 4)]
            [[0, 0, 0, 1], [0, 0, 1, 1], [0, 1, 1, 1]]
            sage: m[0] = 3
            Traceback (most recent call last):
            ...
            TypeError: cannot modify read-only memory
        """
        cdef long long c
        cdef Py_ssize_t sz, total, k = 0
        cdef array.array out
        cdef int *o

        if self.varying_size:
            raise ValueError("flat() needs all the words to have the same length")

        c = self.count()
        sz = self.size
        total = c * sz
        out = array.clone(self.w, total, False)

        if total == 0:
            return out

        o = out.data.as_ints
        with nogil:
            if self.first():
                while k < total:
                    memcpy(o + k, self.ww, sz * sizeof(int))
                    k += sz
                    if not self.next():
                        break
        return out

    def count(self):
        r"""
        Return the number of words of this family.

        The family is run through unless it provides a ``cardinality``, in
        which case that count is returned directly.

        EXAMPLES::

            sage: from combisurf.lyndon_word_family import LyndonWordFamily
            sage: [LyndonWordFamily(2, length).count() for length in range(12)]
            [1, 2, 1, 2, 3, 6, 9, 18, 30, 56, 99, 186]

        :class:`CyclicallyReducedLyndonWordFamily` is counted without being
        run through::

            sage: from combisurf.lyndon_word_family import CyclicallyReducedLyndonWordFamily
            sage: CyclicallyReducedLyndonWordFamily(3, 24).count()
            2483526855452500
        """
        cdef long long c = self.cardinality()
        if c >= 0:
            return c
        c = 0
        with nogil:
            if self.first():
                c = 1
                while self.next():
                    c += 1
        return c


cdef class WordFamilyIterator:
    r"""
    Iterator through a family of words.

    This is what :meth:`WordFamily.__iter__` hands out, in the way a
    ``range`` hands out a ``range_iterator``; it is not meant to be built
    directly. It owns a family of its own, so that several runs through the
    same family, and the family itself, never interfere.

    EXAMPLES::

        sage: from combisurf.lyndon_word_family import LyndonWordFamily
        sage: it = iter(LyndonWordFamily(2, 3))
        sage: it
        iterator through Lyndon words of length 3 on 2 letters
        sage: next(it)
        array('i', [0, 0, 1])

    Iterating an iterator gives it back, so a partial run is resumed rather
    than started again::

        sage: iter(it) is it
        True
        sage: next(it)
        array('i', [0, 1, 1])
        sage: next(it)
        Traceback (most recent call last):
        ...
        StopIteration

    TESTS::

        sage: from combisurf.word_family import WordFamilyIterator
        sage: WordFamilyIterator(None)
        Traceback (most recent call last):
        ...
        ValueError: family must be a WordFamily
    """
    def __init__(self, WordFamily family):
        r"""
        TESTS::

            sage: from combisurf.lyndon_word_family import LyndonWordFamily
            sage: len(list(iter(LyndonWordFamily(3, 3))))
            8
        """
        if family is None:
            raise ValueError("family must be a WordFamily")
        self.family = family
        self.state = 0

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from combisurf.lyndon_word_family import CyclicallyReducedLyndonWordFamily
            sage: iter(CyclicallyReducedLyndonWordFamily(3, 5))
            iterator through cyclically reduced Lyndon words of length 5 on 3 generators
        """
        return f"iterator through {self.family!r}"

    def __iter__(self):
        r"""
        Return this iterator.

        EXAMPLES::

            sage: from combisurf.lyndon_word_family import LyndonWordFamily
            sage: it = iter(LyndonWordFamily(2, 3))
            sage: iter(it) is it
            True
        """
        return self

    def __next__(self):
        r"""
        Return the next word, as an array.

        EXAMPLES::

            sage: from combisurf.lyndon_word_family import LyndonWordFamily
            sage: it = iter(LyndonWordFamily(2, 3))
            sage: next(it)
            array('i', [0, 0, 1])
            sage: next(it)
            array('i', [0, 1, 1])
            sage: next(it)
            Traceback (most recent call last):
            ...
            StopIteration

        Words that are kept are arrays of their own, so that accumulating
        them works as expected::

            sage: list(LyndonWordFamily(2, 4))
            [array('i', [0, 0, 0, 1]), array('i', [0, 0, 1, 1]), array('i', [0, 1, 1, 1])]

            sage: it = iter(LyndonWordFamily(2, 4))
            sage: u = next(it)
            sage: v = next(it)
            sage: u is v
            False
            sage: u, v
            (array('i', [0, 0, 0, 1]), array('i', [0, 0, 1, 1]))

        A word that the caller has let go of is not allocated again, its
        array is filled anew instead. A loop that ends by dropping the word
        it was given therefore allocates one array in all rather than one per
        word::

            sage: it = iter(LyndonWordFamily(2, 4))
            sage: w = next(it)
            sage: address = id(w)
            sage: del w
            sage: id(next(it)) == address
            True
        """
        if self.state == 0:
            if not self.family.first():
                self.state = 2
                raise StopIteration
            self.state = 1
        elif self.state == 1:
            if not self.family.next():
                self.state = 2
                raise StopIteration
        else:
            raise StopIteration
        return self.family._copy()


cdef class ByLengthWordFamily(WordFamily):
    r"""
    Family of the words of another family whose length lies in a range,
    ordered by length first and by the order of that family within a length.

    INPUT:

    - ``base`` -- a :class:`WordFamily`; it is used at every length in turn,
      so its own length is immaterial and its buffer has to be large enough,
      which the constructor sees to

    - ``size_min`` -- non-negative integer; the smallest length

    - ``size_max`` -- non-negative integer; one past the largest length, as
      for ``range``

    EXAMPLES::

        sage: from combisurf.word_family import ByLengthWordFamily
        sage: from combisurf.lyndon_word_family import LyndonWordFamily
        sage: F = ByLengthWordFamily(LyndonWordFamily(2, 0), 1, 4)
        sage: F
        Lyndon words of length 1 <= n < 4 on 2 letters
        sage: F.list()
        [array('i', [0]), array('i', [1]), array('i', [0, 1]),
         array('i', [0, 0, 1]), array('i', [0, 1, 1])]

    The words of one length come out in the order of the underlying family,
    and the lengths in increasing order::

        sage: [len(w) for w in F]
        [1, 1, 2, 3, 3]

    An empty range gives nothing, and a range reduced to one length gives
    that length alone::

        sage: ByLengthWordFamily(LyndonWordFamily(2, 0), 3, 3).list()
        []
        sage: ByLengthWordFamily(LyndonWordFamily(2, 0), 3, 4).list()
        [array('i', [0, 0, 1]), array('i', [0, 1, 1])]

    Since the words have different lengths, they cannot be laid end to end::

        sage: F.flat()
        Traceback (most recent call last):
        ...
        ValueError: flat() needs all the words to have the same length

    TESTS::

        sage: ByLengthWordFamily(LyndonWordFamily(2, 0), -1, 4)
        Traceback (most recent call last):
        ...
        ValueError: size_min (=-1) must be non-negative
        sage: ByLengthWordFamily(None, 1, 4)
        Traceback (most recent call last):
        ...
        ValueError: base must be a WordFamily
    """
    def __init__(self, WordFamily base, int size_min, int size_max):
        r"""
        TESTS::

            sage: from combisurf.word_family import ByLengthWordFamily
            sage: from combisurf.lyndon_word_family import LyndonWordFamily
            sage: ByLengthWordFamily(LyndonWordFamily(3, 0), 0, 5).count()
            33
        """
        if base is None:
            raise ValueError("base must be a WordFamily")
        if size_min < 0:
            raise ValueError(f"size_min (={size_min}) must be non-negative")
        if size_max < 0:
            raise ValueError(f"size_max (={size_max}) must be non-negative")

        self.base = base
        self.size_min = size_min
        self.size_max = size_max
        self.varying_size = True
        self.has_empty_word = (size_min == 0 and size_max > size_min)

        # the base runs at every length in turn, so give it a buffer that
        # holds the longest of them and move its size rather than build one
        # family per length
        base._resize(size_max - 1 if size_max > size_min else 0)
        self.w = base.w
        self.ww = base.ww
        self.capacity = base.capacity
        self.size = size_min

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from combisurf.word_family import ByLengthWordFamily
            sage: from combisurf.lyndon_word_family import CyclicallyReducedLyndonWordFamily
            sage: ByLengthWordFamily(CyclicallyReducedLyndonWordFamily(2, 0), 2, 7)
            cyclically reduced Lyndon words of length 2 <= n < 7 on 2 generators
        """
        base = repr(self.base)
        old = f"of length {self.base.size}"
        new = f"of length {self.size_min} <= n < {self.size_max}"
        return base.replace(old, new) if old in base else f"{base}, of length {self.size_min} <= n < {self.size_max}"

    cdef WordFamily _new(self):
        return ByLengthWordFamily(self.base._new(), self.size_min, self.size_max)

    cdef long long cardinality(self) noexcept nogil:
        cdef long long total = 0, c
        cdef int length

        for length in range(self.size_min, self.size_max):
            self.base._set_size(length)
            c = self.base.cardinality()
            if c < 0:
                return -1
            if total > LLONG_MAX - c:
                return -1
            total += c
        return total

    cdef bint first(self) noexcept nogil:
        self.size = self.size_min - 1
        return self.next()

    cdef bint next(self) noexcept nogil:
        # still words to come at the current length?
        if self.size >= self.size_min and self.base.next():
            return 1

        # otherwise move on to the next length that has any
        while self.size + 1 < self.size_max:
            self.size += 1
            self.base._set_size(self.size)
            if self.base.first():
                return 1

        self.size = self.size_min
        return 0
