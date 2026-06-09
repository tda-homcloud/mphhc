#pragma once
#include <boost/core/bit.hpp>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace mphhc {

std::string GetVersion();

using Index = int32_t;
using Column = std::vector<Index>;
struct BirthDeathPair {
  int dim;
  Index birth, death;

  BirthDeathPair(int dim, Index birth, Index death);
  bool operator==(const BirthDeathPair&) const;
  bool operator<(const BirthDeathPair&) const;
  friend void PrintTo(const BirthDeathPair& pair, std::ostream* os);
};

class BoundaryMatrix {
  struct IndexInfo {
    int dim;
    Index nth;
  };

  std::vector<std::vector<Column>> columns_;
  std::vector<std::vector<Column>> basis_;
  std::vector<std::vector<Index>> local_to_global_index_;
  std::vector<IndexInfo> global_to_local_index_;
  bool reduced_;
  bool save_basis_;

  inline void RecordPivotTwist(int d, int i, std::vector<Index>& pivot_table) {
    Index L = columns_[d][i].back();
    pivot_table[L] = i;
    columns_[d - 1][L].clear();
  }

 public:
  BoundaryMatrix(int maxdim, bool save_basis = false);
  int MaxDim() const;
  Index SetMimCol(Index i, int dim, Column&& col);
  Index SetDimCol(Index i, int dim, const Column& col);
  int NumSimplices() const;
  bool IsReduced() const;
  bool IsSaveBasis() const;

  void ReduceStandard();
  void ReduceTwist();
  std::vector<BirthDeathPair> BirthDeathPairs() const;
  std::vector<Column> Basis() const;
};

class FlatBoundaryMatrix {
  std::vector<Column> columns_;
  std::vector<int8_t> dims_;
  std::vector<Column> basis_;
  bool reduced_;
  bool save_basis_;

 public:
  FlatBoundaryMatrix(int maxdim, bool save_basis = false);
  Index SetMimCol(Index i, int dim, Column&& col);
  Index SetDimCol(Index i, int dim, const Column& col);
  int NumSimplices() const;
  bool IsReduced() const;
  bool IsSaveBasis() const;

  void Reduce();
  std::vector<BirthDeathPair> BirthDeathPairs() const;
  std::vector<Column> Basis() const;
};

class Bitset64 {
 public:
  uint64_t data;

  inline Bitset64() : data(0ull) {}

  inline void Clear() { data = 0ull; }
  inline void Flip(int pos) { data ^= (1ull << pos); }
  inline void Reset(int pos) { data &= ~(1ull << pos); }
  inline void Set(int pos) { data |= (1ull << pos); }

  inline bool Test(int pos) const { return (1ull << pos) & data; }
  inline bool None() const { return data == 0ull; }
  inline bool Any() const { return data != 0ull; }
  inline int Max() const {
    return static_cast<int>(boost::core::bit_width(data)) - 1;
  }
  inline int Min() const { return boost::core::countr_zero(data); }
};

class BitTreeColumn {
  static const uint64_t MASK = (1 << 6) - 1;

  static constexpr int NODE_BLOCK_SIZE_TABLE[6] = {
      -1,
      0,
      1,
      1 + (1 << 6),
      1 + (1 << 6) + (1 << 12),
      1 + (1 << 6) + (1 << 12) + (1 << 18)};

  std::vector<Bitset64> data_;
  int num_index_;
  int height_;
  int node_block_size_;

  void Retrieve(int h, int r, int i, std::vector<Index>* indices) const;

 public:
  static int ComputeHeight(int num_index);
  static int ComputeDataSize(int height, int num_level);

  explicit BitTreeColumn(int num_index);
  void ImportColumn(const Column& vec);

  inline void SetNodes(Index i, int r) {
    for (int h = 1; h <= height_ - 1; ++h) {
      r = (r - 1) >> 6;
      int k = (i >> (6 * h)) & MASK;
      if (data_[r].Test(k))
        return;
      else
        data_[r].Set(k);
    }
  }

  inline void Set(Index i) {
    int r = (i >> 6) + node_block_size_;
    int k = i & MASK;
    data_[r].Set(k);
    SetNodes(i, r);
  }

  inline void ResetNodes(Index i, int r) {
    for (int h = 1; h <= height_ - 1; ++h) {
      if (data_[r].Any()) return;
      r = (r - 1) >> 6;
      int k = (i >> (6 * h)) & MASK;
      data_[r].Reset(k);
    }
  }

  inline void SetXor(Index i) {
    int r = (i >> 6) + node_block_size_;
    int k = i & MASK;
    data_[r].Flip(k);

    if (data_[r].Test(k)) {
      SetNodes(i, r);
    } else {
      ResetNodes(i, r);
    }
  }

  inline Index Max() const {
    using boost::core::bit_width;

    if (data_[0].None()) return -1;

    int r = 0;
    Index i = 0;

    for (int h = height_ - 1; h >= 1; --h) {
      int k = data_[r].Max();
      r = (r << 6) + k + 1;
      i = (i << 6) + k;
    }
    return (i << 6) + data_[r].Max();
  }

  inline bool None() const { return data_[0].None(); }

  inline void Add(const Column& other) {
    for (Index i : other) SetXor(i);
  }

  Column DebugExportColumn();
  void ExportAndClearColumn(Column* col);
};

}  // namespace mphhc
