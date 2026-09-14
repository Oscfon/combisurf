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
