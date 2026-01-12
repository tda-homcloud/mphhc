#include "mphhc/core.hpp"

#include <algorithm>
#include <boost/core/bit.hpp>
#include <cassert>
#include <iostream>
#include <tuple>

namespace mphhc {

std::string get_version() { return "0.1.0"; }

birth_death_pair::birth_death_pair(int d, index bi, index de)
    : dim(d), birth(bi), death(de) {}

bool birth_death_pair::operator==(const birth_death_pair& other) const {
  return dim == other.dim && birth == other.birth && death == other.death;
}

bool birth_death_pair::operator<(const birth_death_pair& other) const {
  using T = std::tuple<int, index, index>;
  return T{dim, birth, death} < T{other.dim, other.birth, other.death};
}

void PrintTo(const birth_death_pair& pair, std::ostream* os) {
  *os << "{ dim: " << pair.dim << ", birth: " << pair.birth
      << ", death: " << pair.death << " }";
}

boundary_matrix::boundary_matrix(int maxdim, bool save_basis)
    : columns_(maxdim + 1),
      basis_(maxdim + 1),
      local_to_global_index_(maxdim + 1),
      reduced_(false),
      save_basis_(save_basis) {}

int boundary_matrix::max_dim() const {
  return static_cast<int>(columns_.size()) - 1;
}

index boundary_matrix::set_dim_col(index i, int dim, column&& col) {
  index new_index = global_to_local_index_.size();
  assert(new_index == i);
  index_info new_local_index = {dim, static_cast<index>(columns_[dim].size())};

  auto& new_column = columns_[dim].emplace_back(col);
  for (std::size_t n = 0; n < new_column.size(); ++n)
    new_column[n] = global_to_local_index_[new_column[n]].nth;

  std::sort(new_column.begin(), new_column.end());

  local_to_global_index_[dim].push_back(new_index);
  global_to_local_index_.push_back(new_local_index);

  return new_index;
}

index boundary_matrix::set_dim_col(index i, int dim, const column& col) {
  column copied_col(col);
  return set_dim_col(i, dim, std::move(copied_col));
}

int boundary_matrix::num_simplices() const {
  return global_to_local_index_.size();
}

bool boundary_matrix::is_reduced() const { return reduced_; }

bool boundary_matrix::is_save_basis() const { return save_basis_; }

class standard_algorithm {
  int max_dim_;
  bit_tree_column bt_column_;
  std::vector<std::vector<column>>& columns_;
  bool save_basis_;
  bit_tree_column basis_column_;
  std::vector<std::vector<column>>& basis_;

 public:
  standard_algorithm(int num_simplices,
                     std::vector<std::vector<column>>& columns, bool save_basis,
                     std::vector<std::vector<column>>& basis)
      : max_dim_(columns.size() - 1),
        bt_column_(num_simplices),
        columns_(columns),
        save_basis_(save_basis),
        basis_column_(save_basis ? num_simplices : 0),
        basis_(basis) {}

  inline void record_pivot(int d, int i, std::vector<index>& pivot_table) {
    index L = columns_[d][i].back();
    pivot_table[L] = i;
  }

  inline void record_basis_vector(int d) {
    if (save_basis_) {
      basis_[d].emplace_back();
      basis_column_.export_and_clear_column(&basis_[d].back());
    }
  }

  inline void record_simple_basis_vector(int d, int i) {
    if (save_basis_) basis_[d].push_back(std::vector<index>{i});
  }

  inline bool reduced(int d, int i, const std::vector<index>& pivot_table) {
    return pivot_table[columns_[d][i].back()] == -1;
  }

  inline void init_basis_column(int i) {
    if (save_basis_) basis_column_.set(i);
  }

  void run() {
    for (int d = max_dim_; d >= 1; --d) {
      std::vector<index> pivot_table(columns_[d - 1].size(), -1);

      for (index i = 0; i < columns_[d].size(); ++i) {
        if (columns_[d][i].empty()) {
          record_simple_basis_vector(d, i);
          continue;
        }
        if (reduced(d, i, pivot_table)) {
          record_pivot(d, i, pivot_table);
          record_simple_basis_vector(d, i);
          continue;
        }

        bt_column_.import_column(columns_[d][i]);
        init_basis_column(i);

        index m, mx;
        while ((mx = bt_column_.max()) != -1 && (m = pivot_table[mx]) != -1) {
          bt_column_.add(columns_[d][m]);
          if (save_basis_) basis_column_.add(basis_[d][m]);
        }

        bt_column_.export_and_clear_column(&columns_[d][i]);
        record_basis_vector(d);

        if (!columns_[d][i].empty()) {
          record_pivot(d, i, pivot_table);
        }
      }
    }
  }
};

void boundary_matrix::reduce_standard() {
  if (reduced_) return;

  standard_algorithm algorithm(num_simplices(), columns_, save_basis_, basis_);
  algorithm.run();

  reduced_ = true;
}

class twist_algorithm {
  int max_dim_;
  bit_tree_column bt_column_;
  std::vector<std::vector<column>>& columns_;

 public:
  twist_algorithm(int num_simplices, std::vector<std::vector<column>>& columns)
      : max_dim_(columns.size() - 1),
        bt_column_(num_simplices),
        columns_(columns) {}

  inline void record_pivot(int d, int i, std::vector<index>& pivot_table) {
    index L = columns_[d][i].back();
    pivot_table[L] = i;
    columns_[d - 1][L].clear();
  }

  void run() {
    for (int d = max_dim_; d >= 1; --d) {
      std::vector<index> pivot_table(columns_[d - 1].size(), -1);

      for (index i = 0; i < columns_[d].size(); ++i) {
        if (columns_[d][i].empty()) continue;
        if (pivot_table[columns_[d][i].back()] == -1) {
          record_pivot(d, i, pivot_table);
          continue;
        }

        bt_column_.import_column(columns_[d][i]);

        index m, mx;
        while ((mx = bt_column_.max()) != -1 && (m = pivot_table[mx]) != -1) {
          bt_column_.add(columns_[d][m]);
        }

        bt_column_.export_and_clear_column(&columns_[d][i]);
        if (!columns_[d][i].empty()) {
          record_pivot(d, i, pivot_table);
        }
      }
    }
  }
};

void boundary_matrix::reduce_twist() {
  if (reduced_) return;

  assert(!save_basis_);

  twist_algorithm algorithm(num_simplices(), columns_);
  algorithm.run();

  reduced_ = true;
}

std::vector<birth_death_pair> boundary_matrix::birth_death_pairs() const {
  assert(reduced_);

  std::vector<birth_death_pair> pairs;
  std::vector<bool> used(num_simplices(), false);

  for (int d = 1; d <= max_dim(); ++d) {
    for (index i = 0; i < columns_[d].size(); ++i) {
      if (!columns_[d][i].empty()) {
        index L = columns_[d][i].back();
        index birth = local_to_global_index_[d - 1][L];
        index death = local_to_global_index_[d][i];
        pairs.emplace_back(d - 1, birth, death);
        used[birth] = used[death] = true;
      }
    }
  }

  for (index i = 0; i < used.size(); ++i) {
    if (!used[i]) pairs.emplace_back(global_to_local_index_[i].dim, i, -1);
  }
  return pairs;
}

std::vector<column> boundary_matrix::basis() const {
  assert(reduced_ && save_basis_);

  std::vector<column> ret;

  for (index i = 0; i < num_simplices(); ++i) {
    index_info t = global_to_local_index_[i];
    if (t.dim == 0) {
      ret.push_back(column{i});
      continue;
    }

    column basis_col;
    for (index local_index : basis_[t.dim][t.nth]) {
      basis_col.push_back(local_to_global_index_[t.dim][local_index]);
    }
    ret.push_back(std::move(basis_col));
  }

  return ret;
}

int bit_tree_column::compute_height(int num_index) {
  int64_t s;
  for (int h = 1, s = 64; h <= 5; ++h) {
    if (num_index <= s) return h;
    s *= 64;
  }
  return -1;
}

int bit_tree_column::compute_data_size(int height, int num_level) {
  assert(height >= 1 && height <= 5);
  return NODE_BLOCK_SIZE_TABLE[height] + (num_level + 63) / 64;
}

bit_tree_column::bit_tree_column(int num_index) {
  assert(num_index < (1 << 30));
  num_index_ = num_index;
  height_ = compute_height(num_index);
  data_.resize(compute_data_size(height_, num_index));
  node_block_size_ = NODE_BLOCK_SIZE_TABLE[height_];
}

void bit_tree_column::import_column(const column& column) {
  for (index i : column) set(i);
}

void bit_tree_column::export_and_clear_column(column* col) {
  col->clear();

  for (index m = max(); m != -1; m = max()) {
    set_xor(m);
    col->push_back(m);
  }
  std::reverse(col->begin(), col->end());
}

column bit_tree_column::debug_export_column() {
  if (none()) return column{};

  column indices;
  retrieve(height_ - 1, 0, 0, &indices);
  return indices;
}

void bit_tree_column::retrieve(int h, int r, index i, column* indices) const {
  using boost::core::bit_width;
  using boost::core::countr_zero;

  bitset64 d = data_[r];
  assert(d.any());

  int max = d.max();
  int min = d.min();

  for (int k = min; k <= max; ++k) {
    if (d.test(k)) {
      if (h == 0) {
        indices->push_back((i << 6) + k);
      } else {
        retrieve(h - 1, (r << 6) + k + 1, (i << 6) + k, indices);
      }
    }
  }
}

}  // namespace mphhc
