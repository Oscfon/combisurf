from cpython cimport array

from combisurf.word_family cimport WordFamily


cdef class LyndonWordFamily(WordFamily):
    cdef int n              # size of the alphabet
    cdef int i              # position of the letter to be incremented next, -1 when exhausted


cdef class CyclicallyReducedLyndonWordFamily(WordFamily):
    cdef int n              # number of generators, the alphabet has 2n letters
    cdef int i              # position of the letter to be incremented next, -1 when exhausted


cdef class LexWordFamily(WordFamily):
    cdef int n              # letters, or generators for the cyclically reduced words
    cdef int top            # largest letter
    cdef int size_min
    cdef int size_max
    cdef int p              # length of the Lyndon prefix of the current node
    cdef int phase          # 0 before the empty word, 1 running

    cdef int _bump(self, int i, int c) noexcept nogil
    cdef bint _emit(self) noexcept nogil


cdef class LexLyndonWordFamily(LexWordFamily):
    pass


cdef class LexCyclicallyReducedLyndonWordFamily(LexWordFamily):
    pass


cdef class UpToInverseWordFamily(WordFamily):
    cdef WordFamily base    # the family whose words it filters
    cdef array.array u      # scratch holding the inverse of the current word
    cdef int *uu            # C view on it

    cdef bint _keep(self) noexcept nogil
    cdef bint _advance(self) noexcept nogil
