import pytest

def test_str_to_cycles():
    from combisurf.permutation import str_to_cycles

    assert str_to_cycles("(0,1)") == [[0, 1]]
    assert str_to_cycles("(0,1)(3,2)") == [[0, 1], [3, 2]]
    assert str_to_cycles("()(0,1)()(2,3)") == [[0, 1], [2, 3]]
    assert str_to_cycles("(0,1,2)(~0,~1,~2)") == [[0, 1, 2], [-1, -2, -3]]

    with pytest.raises(TypeError):
        str_to_cycles(2)


def test_order():
    from combisurf.permutation import perm_init, perm_order

    p = perm_init("()")
    assert perm_order(p) == 1

    p = perm_init("(1)")
    assert perm_order(p) == 1

    p = perm_init("(1,2)")
    assert perm_order(p) == 2

    p = perm_init("(1,2)(3,4)(6,7,8)")
    assert perm_order(p) == 6


def test_trim():
    from combisurf.permutation import perm_init, perm_trim

    for data in [[], [0], [1, 0], [2, 1, 0]]:
        for n in range(4):
            p = perm_init(data + [-1] * n)
            perm_trim(p)
            assert p == perm_init(data)


def test_cycles():
    from combisurf.permutation import perm_random, perm_cycles, perm_are_in_same_orbit, perm_dense_cycles, perm_dense_cycle_positions

    for cycle_length in [1, 2, 3, 5, 10, 50]:
        p = perm_random(10)
        cycles = perm_cycles(p)
        dense_cycles = perm_dense_cycles(p)
        dense_positions = perm_dense_cycle_positions(p)
        for i0, c0 in enumerate(cycles):

            for pos0, j0 in enumerate(c0):
                assert dense_positions[j0] == pos0
                assert dense_cycles[j0] == i0

            for i1, c1 in enumerate(cycles):
                for j0 in c0:
                    for j1 in c1:
                        assert perm_are_in_same_orbit(p, j0, j1) == (i0 == i1)


def test_perm_dense_cycles():
    from array import array
    from itertools import permutations
    from combisurf.permutation import perm_cycles, perm_dense_cycles

    # the labels are consecutive and constant along the cycles
    for n in range(1, 8):
        for q in permutations(range(n)):
            p = array('i', q)
            res = perm_dense_cycles(p)
            cycles = perm_cycles(p)
            assert max(res) + 1 == len(cycles), (q, res)
            for k, c in enumerate(cycles):
                assert all(res[j] == k for j in c), (q, res)

    # same with inactive points
    for n in range(1, 7):
        for q in permutations(range(n)):
            for mask in range(1 << n):
                p = array('i', [-1 if (mask >> i) & 1 or (mask >> q[i]) & 1 else q[i]
                                for i in range(n)])
                if any(p[i] != -1 and p[p[i]] == -1 for i in range(n)):
                    continue
                res = perm_dense_cycles(p)
                active = [i for i in range(n) if p[i] != -1]
                assert all(res[i] == -1 for i in range(n) if p[i] == -1), (q, mask, res)
                labels = sorted(set(res[i] for i in active))
                assert labels == list(range(len(labels))), (q, mask, res)
                for i in active:
                    assert res[i] == res[p[i]], (q, mask, res)


def reference_cycle_positions(p, n=None):
    r"""
    Return the position of each point in its cycle, built from ``perm_cycles``.

    This is the construction the docstring of ``perm_dense_cycle_positions``
    gives as its definition.
    """
    from array import array
    from combisurf.permutation import perm_cycles

    if n is None:
        n = len(p)
    ans = array('i', [-1] * n)
    for c in perm_cycles(p, n=n):
        for pos, j in enumerate(c):
            ans[j] = pos
    return ans


def test_perm_dense_cycle_positions():
    from array import array
    from itertools import permutations
    from combisurf.permutation import perm_dense_cycle_positions

    # every permutation, against the construction from perm_cycles. The two
    # agree because both start a cycle at the smallest point it contains: they
    # scan 0, 1, ..., n - 1 and open a cycle at the first point not yet seen.
    for n in range(1, 8):
        for q in permutations(range(n)):
            p = array('i', q)
            assert perm_dense_cycle_positions(p) == reference_cycle_positions(p), q

    # the empty permutation
    assert perm_dense_cycle_positions(array('i', [])) == array('i', [])


def test_perm_dense_cycle_positions_inactive():
    from array import array
    from itertools import permutations
    from combisurf.permutation import perm_dense_cycle_positions

    # same with inactive points, which must get -1 and must not shift the
    # positions of the others
    tested = 0
    for n in range(1, 7):
        for q in permutations(range(n)):
            for mask in range(1 << n):
                p = array('i', [-1 if (mask >> i) & 1 or (mask >> q[i]) & 1 else q[i]
                                for i in range(n)])
                if any(p[i] != -1 and p[p[i]] == -1 for i in range(n)):
                    continue
                res = perm_dense_cycle_positions(p)
                assert res == reference_cycle_positions(p), (q, mask, list(res))

                # the positions run 0, 1, ... along each cycle and wrap to 0
                for i in range(n):
                    if p[i] == -1:
                        assert res[i] == -1, (q, mask, list(res))
                    else:
                        assert res[p[i]] in (res[i] + 1, 0), (q, mask, list(res))
                tested += 1
    assert tested


def test_perm_dense_cycle_positions_n():
    from array import array
    from itertools import permutations
    from combisurf.permutation import perm_dense_cycle_positions, perm_dense_cycles

    # with n given, only the first n points are scanned. Build p as the direct
    # sum of a permutation of [0, k) and one of [k, n), so that [0, k) is
    # stable and n = k is legitimate.
    for k in range(4):
        for m in range(4):
            for a in permutations(range(k)):
                for b in permutations(range(m)):
                    p = array('i', list(a) + [k + x for x in b])
                    head = array('i', list(a))
                    assert perm_dense_cycle_positions(p, k) == reference_cycle_positions(head)
                    assert perm_dense_cycles(p, k) == perm_dense_cycles(head)

    # n = -1 is the sentinel for "use all of p", not an out of range value
    p = array('i', [1, 0, 3, 2])
    assert perm_dense_cycle_positions(p, -1) == perm_dense_cycle_positions(p)
    assert perm_dense_cycles(p, -1) == perm_dense_cycles(p)

    # n outside [0, len(p)]
    for bad in (5, 9, -2, -3):
        with pytest.raises(ValueError, match=r"must be between 0 and len\(p\)"):
            perm_dense_cycle_positions(p, bad)
        with pytest.raises(ValueError, match=r"must be between 0 and len\(p\)"):
            perm_dense_cycles(p, bad)

    # n in range, but a cycle leaves [0, n)
    q = array('i', [3, 1, 2, 0])
    for bad in (1, 2, 3):
        with pytest.raises(ValueError, match=r"does not map \[0, \d+\) to itself"):
            perm_dense_cycle_positions(q, bad)
        with pytest.raises(ValueError, match=r"does not map \[0, \d+\) to itself"):
            perm_dense_cycles(q, bad)
