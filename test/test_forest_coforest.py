import pytest


def some_roots(n):
    r"""
    Return a sample of random subsets of {0, 1, ..., n - 1}
    """
    if n <= 0:
        raise ValueError

    yield tuple(range(n))
    if n == 1:
        return

    for i in range(min(n, 3)):
        yield (i,)

    if n == 2:
        return

    yield (0,1)
    yield tuple(i for i in range(n) if i % 7 == 2 or i % 13 == 3)
    yield (0, n - 1)


def test_forest_coforest():
    from sage.all import Graph
    from combisurf import OrientedMap

    sample = [
        OrientedMap(vp="(1)(~1,17,13,6,~15)(~13,15,4,~2,~6)(~17,7,~4,8)(~7,9,10,11,~9)(2,~8,~10,~11)"),
        OrientedMap(vp="(2,17,13,5,6,~5)(~13,4,~2,~6)(~17,7,~4,8)(~7,9,10,11,~9)(~8,~10,~11)")
    ]

    for m in sample:
        nv = m.num_vertices()
        nf = m.num_faces()
        ne = m.num_edges()
        for root_vertices in some_roots(nv):
            for root_faces in some_roots(nf):
                forest, coforest, comp_edges = m.forest_coforest_decomposition(root_vertices, root_faces)

                assert len(forest) == nv, m
                assert sum(h != -1 for h in forest) == nv - len(root_vertices)
                assert all(forest[v] == -1 for v in root_vertices)

                assert len(coforest) == nf, m
                assert sum(h != -1 for h in coforest) == nf - len(root_faces)
                assert all(coforest[f] == -1 for f in root_faces)

                assert len(comp_edges) == ne - nv - nf + len(root_vertices) + len(root_faces), (m, root_vertices, root_faces, comp_edges)

                # check that we indeed get an edge partition
                edges = ([h // 2 for h in forest if h != -1] +
                         [h // 2 for h in coforest if h != -1] +
                         list(comp_edges))
                assert set(edges) == set(m.edge_indices())

                # check the forest
                F = Graph(nv)
                h2v = m.half_edge_to_vertex()
                for h in forest:
                    if h == -1:
                        continue
                    F.add_edge(h2v[h], h2v[h ^ 1])
                assert F.is_forest()
                for cc in F.connected_components(sort=False):
                    assert sum(v in cc for v in root_vertices) == 1

                # check the coforest
                C = Graph(nf)
                h2f = m.half_edge_to_face()
                for h in coforest:
                    if h == -1:
                        continue
                    C.add_edge(h2f[h], h2f[h ^ 1])
                assert C.is_forest()
                for cc in C.connected_components(sort=False):
                    assert sum(f in cc for f in root_faces) == 1


def test_forest_coforest_roots():
    # a tree never leaves its connected component. Left to itself the
    # decomposition roots one tree per component; given roots that miss one it
    # says so rather than leaving a vertex or a face out.
    from combisurf import OrientedMap
    from test_fold import maps_with_a_folded_edge

    m = OrientedMap(fp="(0,1,3)(~0,~1,~3)(2,4,5)(~2,~4,~5)")
    assert not m.is_connected()
    assert len(m.connected_components()) == 2
    assert m.num_vertices() == 2 and m.num_faces() == 4

    forest, coforest, comp = m.forest_coforest_decomposition()
    assert sum(1 for x in forest if x == -1) == 2, list(forest)
    assert sum(1 for x in coforest if x == -1) == 2, list(coforest)
    assert all(x >= -1 for x in forest), list(forest)
    assert all(x >= -1 for x in coforest), list(coforest)

    with pytest.raises(ValueError, match="root_vertices must contain a vertex"):
        m.forest_coforest_decomposition((0,), (0, 2))
    with pytest.raises(ValueError, match="root_faces must contain a face"):
        m.forest_coforest_decomposition((0, 1), (0,))

    # a single tree can not span a map that is not connected
    with pytest.raises(ValueError, match="requires a connected map"):
        m.tree_cotree_decomposition()

    # the quad system, on the other hand, is the disjoint union of the quad
    # systems of the components
    q = m.quad_system()
    assert len(q.connected_components()) == 2
    assert set(q.face_profile()) == {4}

    # an empty root list fails even on a connected map, whichever of the two
    # it is, since the single component then holds no root
    c = OrientedMap(vp="(0,1,2)(~0,~1,~2)")
    assert c.is_connected()
    with pytest.raises(ValueError, match="root_vertices must contain a vertex"):
        c.forest_coforest_decomposition((), (0,))
    with pytest.raises(ValueError, match="root_faces must contain a face"):
        c.forest_coforest_decomposition((0,), ())
    with pytest.raises(ValueError, match="root_vertices must contain a vertex"):
        c.forest_coforest_decomposition((), ())

    # on a connected map the automatic roots are the former defaults
    assert c.forest_coforest_decomposition() == c.forest_coforest_decomposition((0,), (0,))

    # in general one tree and one cotree per component, rooted at the
    # smallest vertex resp. face index that the component contains
    for mm in maps_with_a_folded_edge():
        h2v = mm.half_edge_to_vertex()
        h2f = mm.half_edge_to_face()
        vmins, fmins = set(), set()
        for cc in mm.connected_components():
            hs = [h for e in cc for h in (2 * e, 2 * e + 1)
                  if h < len(h2v) and h2v[h] != -1]
            vmins.add(min(h2v[h] for h in hs))
            fmins.add(min(h2f[h] for h in hs))

        forest, coforest, comp = mm.forest_coforest_decomposition()
        assert all(x >= -1 for x in forest), (mm, list(forest))
        assert all(x >= -1 for x in coforest), (mm, list(coforest))
        assert {v for v, x in enumerate(forest) if x == -1} == vmins, mm
        assert {f for f, x in enumerate(coforest) if x == -1} == fmins, mm


def test_forest_coforest_empty_map():
    # the empty map has one vertex, one face and no edge, so its decomposition
    # is the two roots and nothing else
    from combisurf import OrientedMap

    z = OrientedMap()
    assert z.num_vertices() == 1 and z.num_faces() == 1 and z.num_edges() == 0

    forest, coforest, comp = z.forest_coforest_decomposition()
    assert list(forest) == [-1]
    assert list(coforest) == [-1]
    assert list(comp) == []
    assert z.tree_cotree_decomposition() == z.forest_coforest_decomposition()

    # and the quad system of the empty map is the empty map
    assert z.quad_system() == z
    q, proj = z.quad_system(mapping=True)
    assert q == z and proj == []


def test_forest_coforest_with_folded_edges():
    # a folded edge is a loop, so the search skips it: it is never a forest or
    # a coforest edge and it always ends up complementary. The decomposition
    # used to name its inactive half-edge instead.
    from combisurf import OrientedMap
    from test_fold import maps_with_a_folded_edge

    def check(m):
        vp = m.vertex_permutation(copy=False)
        forest, coforest, comp = m.forest_coforest_decomposition()
        active = set(m.half_edges())
        folded = {e for e in m.edge_indices() if vp[(2 * e) ^ 1] == -1}
        assert folded, m
        for h in list(forest) + list(coforest):
            if h == -1:
                continue
            # the surviving half-edge of a folded edge is active, so being
            # active is not enough: the edge itself must not be folded
            assert h in active, (m, h)
            assert h // 2 not in folded, (m, h)
        assert folded <= set(comp), (m, folded, list(comp))

    f = OrientedMap("(0,2,~2)")
    assert f.has_folded_edge()
    check(f)

    tested = 0
    for m in maps_with_a_folded_edge():
        if not m.is_connected():
            continue
        check(m)
        tested += 1
        # and once more with a second edge folded
        for h in list(m.half_edges()):
            r = m.copy(mutable=True)
            try:
                r.fold_half_edge(h)
            except ValueError:
                continue
            if r.is_connected():
                check(r)
                tested += 1
    assert tested
