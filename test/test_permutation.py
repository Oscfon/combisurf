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
