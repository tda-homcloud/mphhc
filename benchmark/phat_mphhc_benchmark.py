import numpy as np
import time
import msgpack

import homcloud.alpha_shape3_ext as alpha_shape3_ext
from homcloud.alpha_filtration import boundary_of_simplex
import homcloud.phat_ext as phat
import mphhc
from homcloud.cubical_ext import CubicalFiltrationExt


def run(f):
    t = time.perf_counter()
    f()
    return time.perf_counter() - t


def alpha_filtration_from_uniform_random_points(rng, n):
    def sortkey(simplex_alpha_pair):
        simplex, alpha = simplex_alpha_pair
        return (alpha, len(simplex))

    pc = rng.uniform(-1, 1, size=(n, 3))
    alpha_shape = alpha_shape3_ext.compute(pc, False)
    index_to_simplex = []
    simplex_to_index = {}

    for index, (simplex, alpha) in enumerate(sorted(alpha_shape, key=sortkey)):
        index_to_simplex.append(simplex)
        simplex_to_index[simplex] = index

    return index_to_simplex, simplex_to_index


def benchmark_alpha_filtration(repeat, num_points):
    rng = np.random.default_rng(58340232)
    times = [[], [], [], [], [], []]
    for _ in range(repeat):
        index_to_simplex, simplex_to_index = (
            alpha_filtration_from_uniform_random_points(rng, num_points)
        )
        matrices = [
            phat.Matrix(len(index_to_simplex)),
            mphhc.Matrix(3),
            mphhc.Matrix(3),
            mphhc.Matrix(3, True),
            mphhc.FlatMatrix(3),
            mphhc.FlatMatrix(3, True),
        ]
        for index, simplex in enumerate(index_to_simplex):
            dim = len(simplex) - 1
            b = [simplex_to_index[t] for t in boundary_of_simplex(simplex)]
            for m in matrices:
                m.set_dim_col(index, dim, b)

        times[0].append(run(matrices[0].reduce_twist))
        times[1].append(run(matrices[1].reduce_twist))
        times[2].append(run(matrices[2].reduce_standard))
        times[3].append(run(matrices[3].reduce_standard))
        times[4].append(run(matrices[4].reduce))
        times[5].append(run(matrices[5].reduce))
        # print(set(m1.birth_death_pairs()) ^ set(m2.birth_death_pairs()))
        # print(set(m1.birth_death_pairs()) ^ set(m3.birth_death_pairs()))
        # assert set(m1.birth_death_pairs()) == set(m2.birth_death_pairs())

    print(
        "phat:",
        np.mean(times[0]),
        "mphhc-twist:",
        np.mean(times[1]),
        "mphhc-standard:",
        np.mean(times[2]),
        "mphhc-standard-with-basis:",
        np.mean(times[3]),
        "mphhc-flat:",
        np.mean(times[4]),
        "mphhc-flat-with-basis:",
        np.mean(times[5]),
    )


def cubical_filtration(rng, size):
    array = np.random.uniform(0, 1, size=size)
    filtration = CubicalFiltrationExt(array, [False, False, False], True)
    return filtration.build_phat_matrix()


def benchmark_random_3d_cubical(repeat, size):
    rng = np.random.default_rng(58340232)

    ptimes_1 = []
    ptimes_2 = []

    for _ in range(repeat):
        m1 = cubical_filtration(rng, size)
        m2 = mphhc.Matrix(3)
        bm = msgpack.unpackb(m1.boundary_map_byte_sequence())
        for index, (dim, col) in enumerate(bm["map"]):
            m2.set_dim_col(index, dim, col)

        ptimes_1.append(run(m1.reduce_twist))
        ptimes_2.append(run(m2.reduce_twist))
        # print(set(m1.birth_death_pairs()) ^ set(m2.birth_death_pairs()))

    print("phat:", np.mean(ptimes_1), "mphhc:", np.mean(ptimes_2))


benchmark_alpha_filtration(repeat=5, num_points=10000)
# phat: 0.021183716505765914 mphhc-twist: 0.012243387848138809 mphhc-standard: 0.1491922608576715
# benchmark_alpha_filtration(repeat=5, num_points=100000)
# benchmark_random_3d_cubical(repeat=5, size=(50, 50, 50))
# phat: 0.13374500460922717 mphhc: 0.09117477238178254
# benchmark_random_3d_cubical(repeat=5, size=(100, 100, 100))
