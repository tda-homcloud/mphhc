#include <gtest/gtest.h>
#define MPHHC_UNITTEST
#include <iterator>
#include <random>
#include <vector>

#include "mphhc/core.hpp"

using namespace mphhc;

TEST(BoundaryMatrixTest, AddDimCol) {
  using namespace mphhc;
  BoundaryMatrix boundary_matrix(2);

  boundary_matrix.SetMimCol(0, 0, Column{});
  boundary_matrix.SetMimCol(1, 0, Column{});
  boundary_matrix.SetMimCol(2, 1, Column{0, 1});

  ASSERT_EQ(boundary_matrix.NumSimplices(), 3);
}

TEST(BitTreeColumnStaticFunctionTest, ComputeHeight) {
  using namespace mphhc;

  ASSERT_EQ(BitTreeColumn::ComputeHeight(64), 1);
  ASSERT_EQ(BitTreeColumn::ComputeHeight(65), 2);
  ASSERT_EQ(BitTreeColumn::ComputeHeight(64 * 64), 2);
  ASSERT_EQ(BitTreeColumn::ComputeHeight(64 * 64 + 1), 3);
  ASSERT_EQ(BitTreeColumn::ComputeHeight(64 * 64 * 64 * 64 * 64), 5);
  ASSERT_EQ(BitTreeColumn::ComputeHeight(64 * 64 * 64 * 64 * 64 + 1), -1);
}

TEST(BitTreeColumnStaticFunctionTest, ComputeDataSize) {
  using namespace mphhc;

  ASSERT_EQ(BitTreeColumn::ComputeDataSize(1, 64), 1);
  ASSERT_EQ(BitTreeColumn::ComputeDataSize(2, 65), 3);
  ASSERT_EQ(BitTreeColumn::ComputeDataSize(2, 128), 3);
  ASSERT_EQ(BitTreeColumn::ComputeDataSize(2, 191), 4);
  ASSERT_EQ(BitTreeColumn::ComputeDataSize(5, 64 * 64 * 64 * 64 * 63),
            1 + 64 + (64 * 64) + (64 * 64 * 64) + (64 * 64 * 64 * 63));
}

struct BitTreeColumnTestParams {
  int num_simplices;

  friend void PrintTo(const BitTreeColumnTestParams &p, std::ostream *os) {
    *os << "{ num_simplices: " << p.num_simplices << " }";
  }
};

class BitTreeColumnTest
    : public ::testing::TestWithParam<BitTreeColumnTestParams> {};

INSTANTIATE_TEST_SUITE_P(
    BitTreeColumnTestParameters, BitTreeColumnTest,
    ::testing::Values(BitTreeColumnTestParams{48},
                      BitTreeColumnTestParams{64 * 48},
                      BitTreeColumnTestParams{64 * 64 * 48},
                      BitTreeColumnTestParams{64 * 64 * 64 * 48},
                      BitTreeColumnTestParams{64 * 64 * 64 * 64 * 48}));

TEST_P(BitTreeColumnTest, SetAndMax) {
  using namespace mphhc;
  BitTreeColumnTestParams const &p = GetParam();
  BitTreeColumn column(p.num_simplices);

  ASSERT_EQ(column.Max(), -1);

  column.Set(29);
  column.Set(17);
  ASSERT_EQ(column.Max(), 29);

  if (p.num_simplices > 64 * 64) {
    column.Set(64 * 64 + 34);
    column.Set(64 * 64 + 127);
    column.Set(64 + 130);
    ASSERT_EQ(column.Max(), 64 * 64 + 127);
  }
}

TEST_P(BitTreeColumnTest, Set) {
  using namespace mphhc;
  using C = mphhc::Column;

  BitTreeColumnTestParams const &p = GetParam();
  BitTreeColumn column(p.num_simplices);

  ASSERT_EQ(column.DebugExportColumn(), (C{}));
  column.Set(29);
  column.Set(63);
  column.Set(17);
  ASSERT_EQ(column.DebugExportColumn(), (C{17, 29, 63}));

  if (p.num_simplices <= 64) return;

  column.Set((32 << 6) + 12);
  column.Set((32 << 6) + 191);
  ASSERT_EQ(column.DebugExportColumn(),
            (C{17, 29, 63, (32 << 6) + 12, (32 << 6) + 191}));
}

TEST_P(BitTreeColumnTest, SetXor) {
  using namespace mphhc;
  using C = mphhc::Column;

  BitTreeColumnTestParams const &p = GetParam();
  BitTreeColumn column(p.num_simplices);

  column.SetXor(29);
  ASSERT_EQ(column.DebugExportColumn(), (C{29}));
  column.SetXor(17);
  ASSERT_EQ(column.DebugExportColumn(), (C{17, 29}));

  column.SetXor(29);
  ASSERT_EQ(column.DebugExportColumn(), (C{17}));

  column.SetXor(17);
  ASSERT_EQ(column.DebugExportColumn(), (C{}));

  if (p.num_simplices > 64 * 64) {
    column.SetXor(17);
    ASSERT_EQ(column.DebugExportColumn(), (C{17}));

    column.SetXor(64 * 64 + 21);
    ASSERT_EQ(column.DebugExportColumn(), (C{17, 64 * 64 + 21}));

    column.SetXor(64 * 32 + 31);
    ASSERT_EQ(column.DebugExportColumn(), (C{17, 64 * 32 + 31, 64 * 64 + 21}));
    column.SetXor(64 * 64 + 191);
    ASSERT_EQ(column.DebugExportColumn(),
              (C{17, 64 * 32 + 31, 64 * 64 + 21, 64 * 64 + 191}));
    column.SetXor(64 * 32 + 31);
    column.SetXor(64 * 64 + 191);
    ASSERT_EQ(column.DebugExportColumn(), (C{17, 64 * 64 + 21}));

    column.SetXor(64 * 64 + 20);
    column.SetXor(64 * 64 + 21);
    ASSERT_EQ(column.DebugExportColumn(), (C{17, 64 * 64 + 20}));

    column.SetXor(64 * 64 + 20);
    column.SetXor(17);
    ASSERT_EQ(column.DebugExportColumn(), (C{}));
  }
}

using rng = std::minstd_rand0;

Column generate_random_column(rng *rng, int num_simplices, int num_samples) {
  std::set<mphhc::Index> indices;
  std::uniform_int_distribution<> dist(0, num_simplices - 1);

  for (int i = 0; i < num_samples; ++i) {
    indices.insert(dist(*rng));
  }

  Column column(indices.begin(), indices.end());
  return column;
}

TEST_P(BitTreeColumnTest, ImportColumnAndExportColumn) {
  BitTreeColumnTestParams const &p = GetParam();
  std::minstd_rand0 rng(7438911);
  Column column = generate_random_column(&rng, p.num_simplices, 400);
  BitTreeColumn bt_column(p.num_simplices);

  bt_column.ImportColumn(column);
  ASSERT_EQ(bt_column.DebugExportColumn(), column);
}

TEST_P(BitTreeColumnTest, ExportAndClearColumn) {
  BitTreeColumnTestParams const &p = GetParam();
  std::minstd_rand0 rng(7438911);
  Column random_column = generate_random_column(&rng, p.num_simplices, 400);
  BitTreeColumn bt_column(p.num_simplices);
  Column exported_column;

  bt_column.ImportColumn(random_column);
  bt_column.ExportAndClearColumn(&exported_column);
  ASSERT_EQ(exported_column, random_column);
  ASSERT_TRUE(bt_column.None());
}

TEST_P(BitTreeColumnTest, Add) {
  BitTreeColumnTestParams const &p = GetParam();
  std::minstd_rand0 rng(7438911);

  Column column_1 = generate_random_column(&rng, p.num_simplices, 400);
  Column column_2 = generate_random_column(&rng, p.num_simplices, 400);
  Column expected;
  std::set_symmetric_difference(column_1.begin(), column_1.end(),
                                column_2.begin(), column_2.end(),
                                std::back_inserter(expected));
  BitTreeColumn bt_column(p.num_simplices);
  bt_column.ImportColumn(column_1);
  bt_column.Add(column_2);

  ASSERT_EQ(bt_column.DebugExportColumn(), expected);
}

TEST(BoundaryMatrixTest, ReduceStandard) {
  using C = mphhc::Column;
  BoundaryMatrix bm(2);

  bm.SetMimCol(0, 0, C{});
  bm.SetMimCol(1, 0, C{});
  bm.SetMimCol(2, 1, C{0, 1});
  bm.SetMimCol(3, 0, C{});
  bm.SetMimCol(4, 0, C{});
  bm.SetMimCol(5, 1, C{3, 4});
  bm.SetMimCol(6, 1, C{1, 3});
  bm.SetMimCol(7, 1, C{0, 4});
  bm.SetMimCol(8, 1, C{1, 4});
  bm.SetMimCol(9, 2, C{5, 6, 8});
  bm.SetMimCol(10, 2, C{2, 7, 8});
  bm.ReduceStandard();

  std::vector<BirthDeathPair> pairs = bm.BirthDeathPairs();
  std::vector<BirthDeathPair> expected = {
      {0, 1, 2}, {0, 4, 5}, {0, 3, 6}, {1, 8, 9}, {1, 7, 10}, {0, 0, -1},
  };

  std::sort(pairs.begin(), pairs.end());
  std::sort(expected.begin(), expected.end());
  ASSERT_EQ(pairs, expected);
}

TEST(BoundaryMatrixTest, ReduceTwist) {
  using C = mphhc::Column;
  BoundaryMatrix bm(2);

  bm.SetMimCol(0, 0, C{});
  bm.SetMimCol(1, 0, C{});
  bm.SetMimCol(2, 1, C{0, 1});
  bm.SetMimCol(3, 0, C{});
  bm.SetMimCol(4, 0, C{});
  bm.SetMimCol(5, 1, C{3, 4});
  bm.SetMimCol(6, 1, C{1, 3});
  bm.SetMimCol(7, 1, C{0, 4});
  bm.SetMimCol(8, 1, C{1, 4});
  bm.SetMimCol(9, 2, C{5, 6, 8});
  bm.SetMimCol(10, 2, C{2, 7, 8});
  bm.ReduceTwist();

  std::vector<BirthDeathPair> pairs = bm.BirthDeathPairs();
  std::vector<BirthDeathPair> expected = {
      {0, 1, 2}, {0, 4, 5}, {0, 3, 6}, {1, 8, 9}, {1, 7, 10}, {0, 0, -1},
  };

  std::sort(pairs.begin(), pairs.end());
  std::sort(expected.begin(), expected.end());
  ASSERT_EQ(pairs, expected);
}

TEST(BoundaryMatrixTest, ReduceStandardWithBasis) {
  using C = mphhc::Column;
  BoundaryMatrix bm(2, true);

  bm.SetMimCol(0, 0, C{});
  bm.SetMimCol(1, 0, C{});
  bm.SetMimCol(2, 1, C{0, 1});
  bm.SetMimCol(3, 0, C{});
  bm.SetMimCol(4, 0, C{});
  bm.SetMimCol(5, 1, C{3, 4});
  bm.SetMimCol(6, 1, C{1, 3});
  bm.SetMimCol(7, 1, C{0, 4});
  bm.SetMimCol(8, 1, C{1, 4});
  bm.SetMimCol(9, 2, C{5, 6, 8});
  bm.SetMimCol(10, 2, C{2, 7, 8});
  bm.ReduceStandard();

  auto basis = bm.Basis();

  std::vector<Column> expected = {
      C{0}, C{1},          C{2},       C{3}, C{4},     C{5},
      C{6}, C{2, 5, 6, 7}, C{5, 6, 8}, C{9}, C{9, 10},
  };

  ASSERT_EQ(basis, expected);
}

TEST(FlatBoundaryMatrixTest, AddDimCol) {
  using namespace mphhc;
  FlatBoundaryMatrix boundary_matrix(2);

  boundary_matrix.SetMimCol(0, 0, Column{});
  boundary_matrix.SetMimCol(1, 0, Column{});
  boundary_matrix.SetMimCol(2, 1, Column{0, 1});

  ASSERT_EQ(boundary_matrix.NumSimplices(), 3);
}

TEST(FlatBoundaryMatrixTest, Reduce) {
  using C = mphhc::Column;
  FlatBoundaryMatrix bm(2);

  bm.SetMimCol(0, 0, C{});
  bm.SetMimCol(1, 0, C{});
  bm.SetMimCol(2, 1, C{0, 1});
  bm.SetMimCol(3, 0, C{});
  bm.SetMimCol(4, 0, C{});
  bm.SetMimCol(5, 1, C{3, 4});
  bm.SetMimCol(6, 1, C{1, 3});
  bm.SetMimCol(7, 1, C{0, 4});
  bm.SetMimCol(8, 1, C{1, 4});
  bm.SetMimCol(9, 2, C{5, 6, 8});
  bm.SetMimCol(10, 2, C{2, 7, 8});
  bm.Reduce();

  std::vector<BirthDeathPair> pairs = bm.BirthDeathPairs();
  std::vector<BirthDeathPair> expected = {
      {0, 1, 2}, {0, 4, 5}, {0, 3, 6}, {1, 8, 9}, {1, 7, 10}, {0, 0, -1},
  };

  std::sort(pairs.begin(), pairs.end());
  std::sort(expected.begin(), expected.end());
  ASSERT_EQ(pairs, expected);
}

TEST(FlatBoundaryMatrixTest, Reduce_2) {
  using C = mphhc::Column;
  FlatBoundaryMatrix bm(2);

  bm.SetMimCol(0, 0, C{});
  bm.SetMimCol(1, 0, C{});
  bm.SetMimCol(2, 1, C{5, 8});
  bm.SetMimCol(3, 1, C{6, 8});
  bm.Reduce();

  std::vector<BirthDeathPair> pairs = bm.BirthDeathPairs();
  std::vector<BirthDeathPair> expected = {
      {0, 0, -1}, {0, 1, -1}, {0, 8, 2}, {0, 6, 3},
  };
  std::sort(pairs.begin(), pairs.end());
  std::sort(expected.begin(), expected.end());
  ASSERT_EQ(pairs, expected);

}

TEST(FlatBoundaryMatrixTest, Basis) {
  using C = mphhc::Column;
  FlatBoundaryMatrix bm(2, true);

  bm.SetMimCol(0, 0, C{});
  bm.SetMimCol(1, 0, C{});
  bm.SetMimCol(2, 1, C{0, 1});
  bm.SetMimCol(3, 0, C{});
  bm.SetMimCol(4, 0, C{});
  bm.SetMimCol(5, 1, C{3, 4});
  bm.SetMimCol(6, 1, C{1, 3});
  bm.SetMimCol(7, 1, C{0, 4});
  bm.SetMimCol(8, 1, C{1, 4});
  bm.SetMimCol(9, 2, C{5, 6, 8});
  bm.SetMimCol(10, 2, C{2, 7, 8});
  bm.Reduce();

  auto basis = bm.Basis();

  std::vector<Column> expected = {
      C{0}, C{1},          C{2},       C{3}, C{4},     C{5},
      C{6}, C{2, 5, 6, 7}, C{5, 6, 8}, C{9}, C{9, 10},
  };

  ASSERT_EQ(basis, expected);
}
