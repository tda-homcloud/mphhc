import mphhc
import pytest

def test_Matrix_reduce():
    bm = mphhc.Matrix(2)
    # C++ indices: 0, 1, 2(edge), 3, 4, 5(edge), 6(edge), 7(edge), 8(edge), 9(face), 10(face)
    bm.set_dim_col(0, []) # 0
    bm.set_dim_col(0, []) # 1
    bm.set_dim_col(1, [0, 1]) # 2
    bm.set_dim_col(0, []) # 3
    bm.set_dim_col(0, []) # 4
    bm.set_dim_col(1, [3, 4]) # 5
    bm.set_dim_col(1, [1, 3]) # 6
    bm.set_dim_col(1, [0, 4]) # 7
    bm.set_dim_col(1, [1, 4]) # 8
    bm.set_dim_col(2, [5, 6, 8]) # 9
    bm.set_dim_col(2, [2, 7, 8]) # 10

    bm.reduce_twist()
    assert bm.is_reduced() is True
    
    pairs = bm.birth_death_pairs()
    # birth_death_pair in C++ is (dim, birth, death)
    # Python returns list of (dim, birth, death) tuples
    expected = [
        (0, 1, 2),
        (0, 4, 5),
        (0, 3, 6),
        (1, 8, 9),
        (1, 7, 10),
        (0, 0, None),
    ]
    
    # Sort for comparison
    pairs.sort()
    expected.sort()
    
    assert pairs == expected
