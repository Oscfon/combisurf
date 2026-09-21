from cpython cimport array


cdef class WordFamily:
    cdef array.array w      # the buffer, handed out to Python as copies
    cdef int *ww            # C view on the buffer; the only thing the C protocol touches
    cdef int capacity       # allocated number of letters
    cdef int size           # length of the current word, at most capacity
    cdef array.array out    # the array handed out last, reused when free
    cdef int out_size       # its length
    cdef bint varying_size  # whether the words have different lengths
    cdef bint has_empty_word      # whether the empty word is one of them

    cdef int _alloc(self, int capacity) except -1
    cdef int _resize(self, int capacity) except -1
    cdef void _set_size(self, int size) noexcept nogil
    cdef array.array _copy(self)
    cdef WordFamily _new(self)

    # mandatory C protocol
    cdef bint first(self) noexcept nogil
    cdef bint next(self) noexcept nogil

    # optional C protocol; -1 means "not supported by this family"
    cdef int next_at(self, int p) noexcept nogil
    cdef long long cardinality(self) noexcept nogil


cdef class WordFamilyIterator:
    cdef WordFamily family  # a family of its own, nothing else iterates it
    cdef int state          # 0 fresh, 1 running, 2 exhausted


cdef class ByLengthWordFamily(WordFamily):
    cdef WordFamily base    # the family whose size this one moves
    cdef int size_min
    cdef int size_max
