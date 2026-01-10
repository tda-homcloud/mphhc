#include "mphhc/core.hpp"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <boost/core/bit.hpp>
#include <tuple>
#include <unordered_map>
#include <deque>

namespace mphhc {

std::string get_version() { return "0.1.0"; }

birth_death_pair::birth_death_pair(int d, index bi, index de):
    dim(d), birth(bi), death(de)
{
}

bool birth_death_pair::operator==(const birth_death_pair& other) const {
  return dim == other.dim && birth == other.birth && death == other.death;
}

bool birth_death_pair::operator<(const birth_death_pair& other) const {
  using T = std::tuple<int, index, index>;
  return T{ dim, birth, death } < T{ other.dim, other.birth, other.death };
}

void PrintTo(const birth_death_pair& pair, std::ostream* os) {
  *os << "{ dim: " << pair.dim << ", birth: " << pair.birth << ", death: " << pair.death << " }"; 
}


boundary_matrix::boundary_matrix(int maxdim):
    columns_(maxdim + 1),
    local_to_global_index_(maxdim + 1),
    reduced_(false)
{
}

int boundary_matrix::max_dim() const {
  return static_cast<int>(columns_.size()) - 1;
}

index boundary_matrix::set_dim_col(index i, int dim, column&& col) {
  index new_index = global_to_local_index_.size();
  assert(new_index == i);
  index_info new_local_index = { dim, static_cast<index>(columns_[dim].size()) };
    
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

bool boundary_matrix::is_reduced() const {
  return reduced_;
}

void boundary_matrix::reduce_standard() {
  if (reduced_)
    return;

  bit_tree_column bt_column(num_simplices());
  
  for (int d = 1; d <= max_dim(); ++d) {
    std::unordered_map<index, index> pivot_table;
    for (index i = 0; i < columns_[d].size(); ++i) {
      if (columns_[d][i].empty())
        continue;
      if (pivot_table.count(columns_[d][i].back()) == 0) {
        pivot_table.insert(std::make_pair(columns_[d][i].back(), i));
        continue;
      }
      
      bt_column.import_column(columns_[d][i]);

      std::unordered_map<index, index>::const_iterator it;
      while ((it = pivot_table.find(bt_column.max())) != pivot_table.end()) {
        bt_column.add(columns_[d][it->second]);
      }
      bt_column.export_and_clear_column(&columns_[d][i]);
      if (!columns_[d][i].empty())
        pivot_table.insert(std::make_pair(columns_[d][i].back(), i));
    }
  }
  
  reduced_ = true;
}

void boundary_matrix::reduce_twist() {
  if (reduced_)
    return;

  bit_tree_column bt_column(num_simplices());
  
  for (int d = max_dim(); d >= 1; --d) {
    // std::unordered_map<index, index> pivot_table;
    std::vector<index> pivot_table(columns_[d - 1].size(), -1);
    
    for (index i = 0; i < columns_[d].size(); ++i) {
      if (columns_[d][i].empty())
        continue;
      if (pivot_table[columns_[d][i].back()] == -1) {
        index L = columns_[d][i].back();
        pivot_table[L] = i;
        columns_[d - 1][L].clear();
        continue;
      }
      
      bt_column.import_column(columns_[d][i]);

      index m;
      while ((m = pivot_table[bt_column.max()]) != -1) {
        bt_column.add(columns_[d][m]);
      }
      bt_column.export_and_clear_column(&columns_[d][i]);
      if (!columns_[d][i].empty()) {
        index L = columns_[d][i].back();
        pivot_table[L] = i;
        columns_[d - 1][L].clear();
      }
    }
  }
  
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
        index birth = local_to_global_index_[d -1][L];
        index death = local_to_global_index_[d][i];
        pairs.emplace_back(d - 1, birth, death);
        used[birth] = used[death] = true;
      }
    }
  }

  for (index i = 0; i < used.size(); ++i) {
    if (!used[i])
      pairs.emplace_back(global_to_local_index_[i].dim, i, -1);
  }
  return pairs;
}

int bit_tree_column::compute_height(int num_index) {
  int64_t s;
  for (int h = 1, s = 64; h <= 5; ++h) {
    if (num_index <= s)
      return h;
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
  for (index i: column)
    set(i);
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
  if (none())
    return column{};

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

}
