#include "mphhc/core.hpp"

#include <algorithm>
#include <boost/core/bit.hpp>
#include <cassert>
#include <iostream>
#include <tuple>

namespace mphhc {

std::string GetVersion() { return "0.1.0"; }

BirthDeathPair::BirthDeathPair(int d, Index bi, Index de)
    : dim(d), birth(bi), death(de) {}

bool BirthDeathPair::operator==(const BirthDeathPair& other) const {
  return dim == other.dim && birth == other.birth && death == other.death;
}

bool BirthDeathPair::operator<(const BirthDeathPair& other) const {
  using T = std::tuple<int, Index, Index>;
  return T{dim, birth, death} < T{other.dim, other.birth, other.death};
}

void PrintTo(const BirthDeathPair& pair, std::ostream* os) {
  *os << "{ dim: " << pair.dim << ", birth: " << pair.birth
      << ", death: " << pair.death << " }";
}

BoundaryMatrix::BoundaryMatrix(int maxdim, bool save_basis)
    : columns_(maxdim + 1),
      basis_(maxdim + 1),
      local_to_global_index_(maxdim + 1),
      reduced_(false),
      save_basis_(save_basis) {}

int BoundaryMatrix::MaxDim() const {
  return static_cast<int>(columns_.size()) - 1;
}

Index BoundaryMatrix::SetMimCol(Index i, int dim, Column&& col) {
  Index new_index = global_to_local_index_.size();
  assert(new_index == i);
  IndexInfo new_local_index = {dim, static_cast<Index>(columns_[dim].size())};

  auto& new_column = columns_[dim].emplace_back(col);
  for (std::size_t n = 0; n < new_column.size(); ++n)
    new_column[n] = global_to_local_index_[new_column[n]].nth;

  std::sort(new_column.begin(), new_column.end());

  local_to_global_index_[dim].push_back(new_index);
  global_to_local_index_.push_back(new_local_index);

  return new_index;
}

Index BoundaryMatrix::SetDimCol(Index i, int dim, const Column& col) {
  Column copied_col(col);
  return SetMimCol(i, dim, std::move(copied_col));
}

int BoundaryMatrix::NumSimplices() const {
  return global_to_local_index_.size();
}

bool BoundaryMatrix::IsReduced() const { return reduced_; }

bool BoundaryMatrix::IsSaveBasis() const { return save_basis_; }

class StandardAlgorithm {
  int max_dim_;
  BitTreeColumn bt_column_;
  std::vector<std::vector<Column>>& columns_;
  bool save_basis_;
  BitTreeColumn basis_column_;
  std::vector<std::vector<Column>>& basis_;

 public:
  StandardAlgorithm(int num_simplices,
                    std::vector<std::vector<Column>>& columns, bool save_basis,
                    std::vector<std::vector<Column>>& basis)
      : max_dim_(columns.size() - 1),
        bt_column_(num_simplices),
        columns_(columns),
        save_basis_(save_basis),
        basis_column_(save_basis ? num_simplices : 0),
        basis_(basis) {}

  inline void RecordPivot(int d, int i, std::vector<Index>& pivot_table) {
    Index L = columns_[d][i].back();
    pivot_table[L] = i;
  }

  inline void RecordBasisVector(int d) {
    if (save_basis_) {
      basis_[d].emplace_back();
      basis_column_.ExportAndClearColumn(&basis_[d].back());
    }
  }

  inline void RecordSimpleBasisVector(int d, int i) {
    if (save_basis_) basis_[d].push_back(std::vector<Index>{i});
  }

  inline bool IsReduced(int d, int i, const std::vector<Index>& pivot_table) {
    return pivot_table[columns_[d][i].back()] == -1;
  }

  inline void InitBasisColumn(int i) {
    if (save_basis_) basis_column_.Set(i);
  }

  void Run() {
    for (int d = max_dim_; d >= 1; --d) {
      std::vector<Index> pivot_table(columns_[d - 1].size(), -1);

      for (Index i = 0; i < columns_[d].size(); ++i) {
        if (columns_[d][i].empty()) {
          RecordSimpleBasisVector(d, i);
          continue;
        }
        if (IsReduced(d, i, pivot_table)) {
          RecordPivot(d, i, pivot_table);
          RecordSimpleBasisVector(d, i);
          continue;
        }

        bt_column_.ImportColumn(columns_[d][i]);
        InitBasisColumn(i);

        Index m, mx;
        while ((mx = bt_column_.Max()) != -1 && (m = pivot_table[mx]) != -1) {
          bt_column_.Add(columns_[d][m]);
          if (save_basis_) basis_column_.Add(basis_[d][m]);
        }

        bt_column_.ExportAndClearColumn(&columns_[d][i]);
        RecordBasisVector(d);

        if (!columns_[d][i].empty()) {
          RecordPivot(d, i, pivot_table);
        }
      }
    }
  }
};

void BoundaryMatrix::ReduceStandard() {
  if (reduced_) return;

  StandardAlgorithm algorithm(NumSimplices(), columns_, save_basis_, basis_);
  algorithm.Run();

  reduced_ = true;
}

class TwistAlgorithm {
  int max_dim_;
  BitTreeColumn bt_column_;
  std::vector<std::vector<Column>>& columns_;

 public:
  TwistAlgorithm(int num_simplices, std::vector<std::vector<Column>>& columns)
      : max_dim_(columns.size() - 1),
        bt_column_(num_simplices),
        columns_(columns) {}

  inline void RecordPivot(int d, int i, std::vector<Index>& pivot_table) {
    Index L = columns_[d][i].back();
    pivot_table[L] = i;
    columns_[d - 1][L].clear();
  }

  void Run() {
    for (int d = max_dim_; d >= 1; --d) {
      std::vector<Index> pivot_table(columns_[d - 1].size(), -1);

      for (Index i = 0; i < columns_[d].size(); ++i) {
        if (columns_[d][i].empty()) continue;
        if (pivot_table[columns_[d][i].back()] == -1) {
          RecordPivot(d, i, pivot_table);
          continue;
        }

        bt_column_.ImportColumn(columns_[d][i]);

        Index m, mx;
        while ((mx = bt_column_.Max()) != -1 && (m = pivot_table[mx]) != -1) {
          bt_column_.Add(columns_[d][m]);
        }

        bt_column_.ExportAndClearColumn(&columns_[d][i]);
        if (!columns_[d][i].empty()) {
          RecordPivot(d, i, pivot_table);
        }
      }
    }
  }
};

void BoundaryMatrix::ReduceTwist() {
  if (reduced_) return;

  assert(!save_basis_);

  TwistAlgorithm algorithm(NumSimplices(), columns_);
  algorithm.Run();

  reduced_ = true;
}

std::vector<BirthDeathPair> BoundaryMatrix::BirthDeathPairs() const {
  assert(reduced_);

  std::vector<BirthDeathPair> pairs;
  std::vector<bool> used(NumSimplices(), false);

  for (int d = 1; d <= MaxDim(); ++d) {
    for (Index i = 0; i < columns_[d].size(); ++i) {
      if (!columns_[d][i].empty()) {
        Index L = columns_[d][i].back();
        Index birth = local_to_global_index_[d - 1][L];
        Index death = local_to_global_index_[d][i];
        pairs.emplace_back(d - 1, birth, death);
        used[birth] = used[death] = true;
      }
    }
  }

  for (Index i = 0; i < used.size(); ++i) {
    if (!used[i]) pairs.emplace_back(global_to_local_index_[i].dim, i, -1);
  }
  return pairs;
}

std::vector<Column> BoundaryMatrix::Basis() const {
  assert(reduced_ && save_basis_);

  std::vector<Column> ret;

  for (Index i = 0; i < NumSimplices(); ++i) {
    IndexInfo t = global_to_local_index_[i];
    if (t.dim == 0) {
      ret.push_back(Column{i});
      continue;
    }

    Column basis_col;
    for (Index local_index : basis_[t.dim][t.nth]) {
      basis_col.push_back(local_to_global_index_[t.dim][local_index]);
    }
    ret.push_back(std::move(basis_col));
  }

  return ret;
}

FlatBoundaryMatrix::FlatBoundaryMatrix(int maxdim, bool save_basis)
    : reduced_(false),
      save_basis_(save_basis) {}

Index FlatBoundaryMatrix::SetMimCol(Index i, int dim, Column&& col) {
  Index new_index = columns_.size();
  assert(new_index == i);

  dims_.push_back(dim);
  auto& new_column = columns_.emplace_back(col);
  std::sort(new_column.begin(), new_column.end());
  return new_index;
}

Index FlatBoundaryMatrix::SetDimCol(Index i, int dim, const Column& col) {
  Column copied_col(col);
  return SetMimCol(i, dim, std::move(copied_col));
}

int FlatBoundaryMatrix::NumSimplices() const {
  return columns_.size();
}

bool FlatBoundaryMatrix::IsReduced() const { return reduced_; }

bool FlatBoundaryMatrix::IsSaveBasis() const { return save_basis_; }

class FlatAlgorithm {
  BitTreeColumn bt_column_;
  std::vector<Column>& columns_;
  bool save_basis_;
  BitTreeColumn basis_column_;
  std::vector<Column>& basis_;
  std::vector<Index> pivot_table_;
  
 public:
  FlatAlgorithm(std::vector<Column>& columns,
                bool save_basis,
                std::vector<Column>& basis)
      : bt_column_(columns.size()),
        columns_(columns),
        save_basis_(save_basis),
        basis_column_(save_basis ? columns.size() : 0),
        basis_(basis),
        pivot_table_(columns_.size(), -1) {}

  inline void RecordPivot(int i) {
    Index L = columns_[i].back();
    pivot_table_[L] = i;
  }

  inline void RecordBasisVector() {
    if (save_basis_) {
      basis_.emplace_back();
      basis_column_.ExportAndClearColumn(&basis_.back());
    }
  }
  
  inline bool IsReduced(Index i) const {
    return pivot_table_[columns_[i].back()] == -1;
  }

  inline void RecordSimpleBasisVector(Index i) {
    if (save_basis_) basis_.push_back(std::vector<Index>{i});
  }

  inline void InitBasisColumn(int i) {
    if (save_basis_) basis_column_.Set(i);
  }
  
  void Run () {
    for (Index i = 0; i < columns_.size(); ++i) {
      if (columns_[i].empty()) {
        RecordSimpleBasisVector(i);
        continue;
      }

      if (IsReduced(i)) {
        RecordPivot(i);
        RecordSimpleBasisVector(i);
        continue;
      }

      bt_column_.ImportColumn(columns_[i]);
      InitBasisColumn(i);
      
      Index m, mx;
      while ((mx = bt_column_.Max()) != -1 && (m = pivot_table_[mx]) != -1) {
        bt_column_.Add(columns_[m]);
        if (save_basis_) { basis_column_.Add(basis_[m]); }
      }

      bt_column_.ExportAndClearColumn(&columns_[i]);
      RecordBasisVector();
      
      if (!columns_[i].empty()) {
        RecordPivot(i);
      }
    }
  }
};

void FlatBoundaryMatrix::Reduce() {
  if (reduced_) return;

  FlatAlgorithm algorithm(columns_, save_basis_, basis_);
  algorithm.Run();

  reduced_ = true;
}

std::vector<BirthDeathPair> FlatBoundaryMatrix::BirthDeathPairs() const {
  assert(reduced_);

  std::vector<BirthDeathPair> pairs;
  std::vector<bool> used(NumSimplices(), false);

  for (Index i = 0; i < columns_.size(); ++i) {
    if (!columns_[i].empty()) {
      Index birth = columns_[i].back();
      pairs.emplace_back(dims_[i] - 1, birth, i);
      used[birth] = used[i] = true;
    }
  }

  for (Index i = 0; i < NumSimplices(); ++i) {
    if (!used[i]) {
      pairs.emplace_back(dims_[i], i, -1);
    }
  }

  return pairs;
}

int BitTreeColumn::ComputeHeight(int num_index) {
  int64_t s;
  for (int h = 1, s = 64; h <= 5; ++h) {
    if (num_index <= s) return h;
    s *= 64;
  }
  return -1;
}

std::vector<Column> FlatBoundaryMatrix::Basis() const {
  assert(reduced_ && save_basis_);

  return basis_;
}

int BitTreeColumn::ComputeDataSize(int height, int num_level) {
  assert(height >= 1 && height <= 5);
  return NODE_BLOCK_SIZE_TABLE[height] + (num_level + 63) / 64;
}

BitTreeColumn::BitTreeColumn(int num_index) {
  assert(num_index < (1 << 30));
  num_index_ = num_index;
  height_ = ComputeHeight(num_index);
  data_.resize(ComputeDataSize(height_, num_index));
  node_block_size_ = NODE_BLOCK_SIZE_TABLE[height_];
}

void BitTreeColumn::ImportColumn(const Column& column) {
  for (Index i : column) Set(i);
}

void BitTreeColumn::ExportAndClearColumn(Column* col) {
  col->clear();

  for (Index m = Max(); m != -1; m = Max()) {
    SetXor(m);
    col->push_back(m);
  }
  std::reverse(col->begin(), col->end());
}

Column BitTreeColumn::DebugExportColumn() {
  if (None()) return Column{};

  Column indices;
  Retrieve(height_ - 1, 0, 0, &indices);
  return indices;
}

void BitTreeColumn::Retrieve(int h, int r, Index i, Column* indices) const {
  using boost::core::bit_width;
  using boost::core::countr_zero;

  Bitset64 d = data_[r];
  assert(d.Any());

  int max = d.Max();
  int min = d.Min();

  for (int k = min; k <= max; ++k) {
    if (d.Test(k)) {
      if (h == 0) {
        indices->push_back((i << 6) + k);
      } else {
        Retrieve(h - 1, (r << 6) + k + 1, (i << 6) + k, indices);
      }
    }
  }
}

}  // namespace mphhc
