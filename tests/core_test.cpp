#include <gtest/gtest.h>
#include "mphhc/core.hpp"

TEST(CoreTest, AddFunction) {
    EXPECT_EQ(mphhc::add(1, 2), 3);
    EXPECT_NE(mphhc::add(1, 2), 4);
}

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

TEST(BitTreeColumnTest, ComputeHeight) {
  using namespace mphhc;

  ASSERT_EQ(bit_tree_column::compute_height(64), 1);
  ASSERT_EQ(bit_tree_column::compute_height(65), 2);
  ASSERT_EQ(bit_tree_column::compute_height(64 * 64), 2);
  ASSERT_EQ(bit_tree_column::compute_height(64 * 64 + 1), 3);
  ASSERT_EQ(bit_tree_column::compute_height(64 * 64 * 64 * 64 * 64), 5);
  ASSERT_EQ(bit_tree_column::compute_height(64 * 64 * 64 * 64 * 64 + 1), -1);
}

TEST(BitTreeColumnTest, ComputeDataSize) {
  using namespace mphhc;

  ASSERT_EQ(bit_tree_column::compute_data_size(1, 64), 1);
  ASSERT_EQ(bit_tree_column::compute_data_size(2, 65), 3);
  ASSERT_EQ(bit_tree_column::compute_data_size(2, 128), 3);
  ASSERT_EQ(bit_tree_column::compute_data_size(2, 129), 4);
  ASSERT_EQ(bit_tree_column::compute_data_size(5, 64 * 64 * 64 * 64 * 63),
            1 + 64 + (64 * 64) + (64 * 64 * 64) + (64 * 64 * 64 * 63));
            
}
