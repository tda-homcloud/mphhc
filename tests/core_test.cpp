#include <gtest/gtest.h>
#define MPHHC_UNITTEST
#include "mphhc/core.hpp"
#include <vector>

TEST(CoreTest, VersionString) {
    EXPECT_EQ(mphhc::get_version(), "0.1.0");
}

TEST(BoundaryMatrixTest, AddDimCol) {
  using namespace mphhc;
  boundary_matrix boundary_matrix(2);

  boundary_matrix.add_dim_col(0, column{});
  boundary_matrix.add_dim_col(0, column{});
  boundary_matrix.add_dim_col(1, column{0, 1});

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
  ASSERT_EQ(bit_tree_column::compute_data_size(2, 129), 4);
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
         BitTreeColumnTestParams{64 * 64 * 64 * 64 * 48})
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
    column.set(64 * 64 + 129);
    column.set(64  + 130);
    ASSERT_EQ(column.max(), 64 * 64 + 129);
  }
}

TEST_P(BitTreeColumnTest, SetAndToVector) {
  using namespace mphhc;
  using V = std::vector<mphhc::index>;

  BitTreeColumnTestParams const &p = GetParam();
  bit_tree_column column(p.num_simplices);

  ASSERT_EQ(column.to_vector(), (V{}));
  column.set(29);
  column.set(63);
  column.set(17);
  ASSERT_EQ(column.to_vector(), (V{17, 29, 63}));

  if (p.num_simplices <= 64) return;

  column.set((32<<6) + 12);

  ASSERT_EQ(column.max(), (32<<6) + 12);
  ASSERT_EQ(column.to_vector(), (V{17, 29, 63, (32<<6) + 12}));
}

TEST_P(BitTreeColumnTest, SetXor) {
  using namespace mphhc;
  using mphhc::index;
  using V = std::vector<index>;
  
  BitTreeColumnTestParams const &p = GetParam();
  bit_tree_column column(p.num_simplices);

  column.set_xor(29);
  ASSERT_EQ(column.to_vector(), (V{29}));
  column.set_xor(17);
  ASSERT_EQ(column.to_vector(), (V{17, 29}));

  column.set_xor(29);
  ASSERT_EQ(column.to_vector(), (V{17}));

  column.set_xor(17);
  ASSERT_EQ(column.to_vector(), (V{}));

  if (p.num_simplices > 64 * 64) {
    column.set_xor(17);
    ASSERT_EQ(column.to_vector(), (V{17}));

    column.set_xor(64 * 64 + 21);
    ASSERT_EQ(column.to_vector(), (V{17, 64 * 64 + 21}));
    
    column.set_xor(64 * 32 + 31);
    ASSERT_EQ(column.to_vector(), (V{17, 64 * 32 + 31, 64 * 64 + 21}));
    column.set_xor(64 * 64 + 129);
    column.set_xor(64 * 32 + 31);
    column.set_xor(64 * 64 + 129);
    ASSERT_EQ(column.to_vector(), (V{17, 64 * 64 + 21}));

    column.set_xor(64 * 64 + 20);
    column.set_xor(64 * 64 + 21);
    ASSERT_EQ(column.to_vector(), (V{17, 64 * 64 + 20}));

    column.set_xor(64 * 64 + 20);
    column.set_xor(17);
    ASSERT_EQ(column.to_vector(), (V{}));
  }
}
