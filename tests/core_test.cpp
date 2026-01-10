#include <gtest/gtest.h>
#define MPHHC_UNITTEST
#include "mphhc/core.hpp"
#include <vector>
#include <random>
#include <iterator>

using namespace mphhc;

TEST(BoundaryMatrixTest, AddDimCol) {
  using namespace mphhc;
  boundary_matrix boundary_matrix(2);

  boundary_matrix.set_dim_col(0, 0, column{});
  boundary_matrix.set_dim_col(1, 0, column{});
  boundary_matrix.set_dim_col(2, 1, column{0, 1});

  ASSERT_EQ(boundary_matrix.num_simplices(), 3);
}

TEST(BitTreeColumnStaticFunctionTest, ComputeHeight) {
  using namespace mphhc;

  ASSERT_EQ(bit_tree_column::compute_height(64), 1);
  ASSERT_EQ(bit_tree_column::compute_height(65), 2);
  ASSERT_EQ(bit_tree_column::compute_height(64 * 64), 2);
  ASSERT_EQ(bit_tree_column::compute_height(64 * 64 + 1), 3);
  ASSERT_EQ(bit_tree_column::compute_height(64 * 64 * 64 * 64 * 64), 5);
  ASSERT_EQ(bit_tree_column::compute_height(64 * 64 * 64 * 64 * 64 + 1), -1);
}

TEST(BitTreeColumnStaticFunctionTest, ComputeDataSize) {
  using namespace mphhc;

  ASSERT_EQ(bit_tree_column::compute_data_size(1, 64), 1);
  ASSERT_EQ(bit_tree_column::compute_data_size(2, 65), 3);
  ASSERT_EQ(bit_tree_column::compute_data_size(2, 128), 3);
  ASSERT_EQ(bit_tree_column::compute_data_size(2, 191), 4);
  ASSERT_EQ(bit_tree_column::compute_data_size(5, 64 * 64 * 64 * 64 * 63),
            1 + 64 + (64 * 64) + (64 * 64 * 64) + (64 * 64 * 64 * 63));
            
}

struct BitTreeColumnTestParams {
  int num_simplices;

  friend void PrintTo(const BitTreeColumnTestParams& p, std::ostream* os) {
    *os << "{ num_simplices: " << p.num_simplices << " }"; 
  }
};

class BitTreeColumnTest : public ::testing::TestWithParam<BitTreeColumnTestParams> {
};

INSTANTIATE_TEST_SUITE_P(
    BitTreeColumnTestParameters,
    BitTreeColumnTest,
    ::testing::Values(
         BitTreeColumnTestParams{48},
         BitTreeColumnTestParams{64 * 48},
         BitTreeColumnTestParams{64 * 64 * 48},
         BitTreeColumnTestParams{64 * 64 * 64 * 48},
         BitTreeColumnTestParams{64 * 64 * 64 * 64 * 48}
    )
);

TEST_P(BitTreeColumnTest, SetAndMax) {
  using namespace mphhc;
  BitTreeColumnTestParams const &p = GetParam();
  bit_tree_column column(p.num_simplices);

  ASSERT_EQ(column.max(), -1);
  
  column.set(29);
  column.set(17);
  ASSERT_EQ(column.max(), 29);

  if (p.num_simplices > 64 * 64) {
    column.set(64 * 64 + 34);
    column.set(64 * 64 + 127);
    column.set(64  + 130);
    ASSERT_EQ(column.max(), 64 * 64 + 127);
  }
}

TEST_P(BitTreeColumnTest, Set) {
  using namespace mphhc;
  using C = mphhc::column;

  BitTreeColumnTestParams const &p = GetParam();
  bit_tree_column column(p.num_simplices);

  ASSERT_EQ(column.debug_export_column(), (C{}));
  column.set(29);
  column.set(63);
  column.set(17);
  ASSERT_EQ(column.debug_export_column(), (C{17, 29, 63}));

  if (p.num_simplices <= 64) return;

  column.set((32<<6) + 12);
  column.set((32<<6) + 191);
  ASSERT_EQ(column.debug_export_column(), (C{17, 29, 63, (32<<6) + 12, (32<<6) + 191}));
}

TEST_P(BitTreeColumnTest, SetXor) {
  using namespace mphhc;
  using C = mphhc::column;
  
  BitTreeColumnTestParams const &p = GetParam();
  bit_tree_column column(p.num_simplices);

  column.set_xor(29);
  ASSERT_EQ(column.debug_export_column(), (C{29}));
  column.set_xor(17);
  ASSERT_EQ(column.debug_export_column(), (C{17, 29}));

  column.set_xor(29);
  ASSERT_EQ(column.debug_export_column(), (C{17}));

  column.set_xor(17);
  ASSERT_EQ(column.debug_export_column(), (C{}));

  if (p.num_simplices > 64 * 64) {
    column.set_xor(17);
    ASSERT_EQ(column.debug_export_column(), (C{17}));

    column.set_xor(64 * 64 + 21);
    ASSERT_EQ(column.debug_export_column(), (C{17, 64 * 64 + 21}));
    
    column.set_xor(64 * 32 + 31);
    ASSERT_EQ(column.debug_export_column(), (C{17, 64 * 32 + 31, 64 * 64 + 21}));
    column.set_xor(64 * 64 + 191);
    ASSERT_EQ(column.debug_export_column(), (C{17, 64 * 32 + 31, 64 * 64 + 21, 64 * 64 + 191}));
    column.set_xor(64 * 32 + 31);
    column.set_xor(64 * 64 + 191);
    ASSERT_EQ(column.debug_export_column(), (C{17, 64 * 64 + 21}));

    column.set_xor(64 * 64 + 20);
    column.set_xor(64 * 64 + 21);
    ASSERT_EQ(column.debug_export_column(), (C{17, 64 * 64 + 20}));

    column.set_xor(64 * 64 + 20);
    column.set_xor(17);
    ASSERT_EQ(column.debug_export_column(), (C{}));
  }
}

using rng = std::minstd_rand0;

column generate_random_column(rng* rng, int num_simplices, int num_samples) {
  std::set<mphhc::index> indices;
  std::uniform_int_distribution<> dist(0, num_simplices - 1);

  for (int i = 0; i < num_samples; ++i) {
    indices.insert(dist(*rng));
  }

  column column(indices.begin(), indices.end());
  return column;
}

TEST_P(BitTreeColumnTest, ImportColumnAndExportColumn) {
  BitTreeColumnTestParams const &p = GetParam();
  std::minstd_rand0 rng(7438911);
  column column = generate_random_column(&rng, p.num_simplices, 400);
  bit_tree_column bt_column(p.num_simplices);

  bt_column.import_column(column);
  ASSERT_EQ(bt_column.debug_export_column(), column);
}

TEST_P(BitTreeColumnTest, ExportAndClearColumn) {
  BitTreeColumnTestParams const &p = GetParam();
  std::minstd_rand0 rng(7438911);
  column random_column = generate_random_column(&rng, p.num_simplices, 400);
  bit_tree_column bt_column(p.num_simplices);
  column exported_column;

  bt_column.import_column(random_column);
  bt_column.export_and_clear_column(&exported_column);
  ASSERT_EQ(exported_column, random_column);
  ASSERT_TRUE(bt_column.none());
}

TEST_P(BitTreeColumnTest, Add) {
  BitTreeColumnTestParams const &p = GetParam();
  std::minstd_rand0 rng(7438911);

  column column_1 = generate_random_column(&rng, p.num_simplices, 400);
  column column_2 = generate_random_column(&rng, p.num_simplices, 400);
  column expected;
  std::set_symmetric_difference(column_1.begin(), column_1.end(), column_2.begin(), column_2.end(),
                                std::back_inserter(expected));
  bit_tree_column bt_column(p.num_simplices);
  bt_column.import_column(column_1);
  bt_column.add(column_2);
  
  ASSERT_EQ(bt_column.debug_export_column(), expected);
}


TEST(BoundaryMatrixTest, ReduceStandard) {
  using C = mphhc::column;
  boundary_matrix bm(2);

  bm.set_dim_col(0, 0, C{});
  bm.set_dim_col(1, 0, C{});
  bm.set_dim_col(2, 1, C{0, 1});
  bm.set_dim_col(3, 0, C{});
  bm.set_dim_col(4, 0, C{});
  bm.set_dim_col(5, 1, C{3, 4});
  bm.set_dim_col(6, 1, C{1, 3});
  bm.set_dim_col(7, 1, C{0, 4});
  bm.set_dim_col(8, 1, C{1, 4});
  bm.set_dim_col(9, 2, C{5, 6, 8});
  bm.set_dim_col(10, 2, C{2, 7, 8});
  bm.reduce_standard();

  std::vector<birth_death_pair> pairs = bm.birth_death_pairs();
  std::vector<birth_death_pair> expected = {
    {0, 1, 2},
    {0, 4, 5},
    {0, 3, 6},
    {1, 8, 9},
    {1, 7, 10},
    {0, 0, -1},
  };

  std::sort(pairs.begin(), pairs.end());
  std::sort(expected.begin(), expected.end());
  ASSERT_EQ(pairs, expected);
}

TEST(BoundaryMatrixTest, ReduceTwist) {
  using C = mphhc::column;
  boundary_matrix bm(2);

  bm.set_dim_col(0, 0, C{});
  bm.set_dim_col(1, 0, C{});
  bm.set_dim_col(2, 1, C{0, 1});
  bm.set_dim_col(3, 0, C{});
  bm.set_dim_col(4, 0, C{});
  bm.set_dim_col(5, 1, C{3, 4});
  bm.set_dim_col(6, 1, C{1, 3});
  bm.set_dim_col(7, 1, C{0, 4});
  bm.set_dim_col(8, 1, C{1, 4});
  bm.set_dim_col(9, 2, C{5, 6, 8});
  bm.set_dim_col(10, 2, C{2, 7, 8});
  bm.reduce_twist();

  std::vector<birth_death_pair> pairs = bm.birth_death_pairs();
  std::vector<birth_death_pair> expected = {
    {0, 1, 2},
    {0, 4, 5},
    {0, 3, 6},
    {1, 8, 9},
    {1, 7, 10},
    {0, 0, -1},
  };

  std::sort(pairs.begin(), pairs.end());
  std::sort(expected.begin(), expected.end());
  ASSERT_EQ(pairs, expected);
}
