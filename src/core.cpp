#include "mphhc/core.hpp"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <boost/core/bit.hpp>
#include <bitset>

namespace mphhc {

std::string get_version() { return "0.1.0"; }

boundary_matrix::boundary_matrix(int maxdim):
    columns_(maxdim + 1),
    local_to_global_index_(maxdim + 1)
{
}

index boundary_matrix::add_dim_col(int dim, column&& col) {
  index new_index = global_to_local_index_.size();
  index_info new_local_index = { dim, static_cast<index>(columns_[dim].size()) };
    
  auto& new_column = columns_[dim].emplace_back(col);
  std::sort(new_column.begin(), new_column.end());

  local_to_global_index_[dim].push_back(new_index);
  global_to_local_index_.push_back(new_local_index);

  return new_index;
}

index boundary_matrix::add_dim_col(int dim, const column& col) {
  column copied_col(col);
  return add_dim_col(dim, std::move(copied_col));
}

int boundary_matrix::num_simplices() const {
  int sum = 0;
  for (const auto& v: columns_) {
    sum += v.size();
  }
  return sum;
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
  static int node_size[6] = {
    0,
    1,
    1 + (1<<6),
    1 + (1<<6) + (1<<12),
    1 + (1<<6) + (1<<12) + (1<<18),
  };
  return node_size[height - 1] + (num_level + 63) / 64;
}

bit_tree_column::bit_tree_column(int num_index) {
  assert(num_index < (1 << 30));
  num_index_ = num_index;
  height_ = compute_height(num_index);
  data_.resize(compute_data_size(height_, num_index));
  clear();
}

void bit_tree_column::clear() {
  data_[0] = 0;
}

void bit_tree_column::set(index i) {
  int r = 0;
  for (int h = height_ - 1; h >= 1; --h) {
    int k = (i >> 6 * h) & MASK;
    int next_r = (r << 6) + k + 1;
    if ((data_[r] & (1ull<<k)) == 0) {
      data_[next_r] = 0;
      data_[r] |= 1ull<<k;
    }
    r = next_r;
  }

  int k = i & MASK;
  data_[r] |= (1ull<<k);
}

index bit_tree_column::max() const {
  using boost::core::bit_width;
  
  if (data_[0] == 0)
    return -1;
  
  int r = 0;
  index i = 0;
  
  for (int h = height_ - 1; h >= 1; --h) {
    int k = bit_width(data_[r]) - 1;
    r = (r << 6) + k + 1;
    i = (i << 6) + k;
  }
  return (i << 6) + bit_width(data_[r]) - 1;
}

void bit_tree_column::set_xor(index i) {
  int r = 0;

  for (int h = height_ - 1; h >= 1; --h) {
    int k = (i >> 6 * h) & MASK;
    int next_r = (r << 6) + k + 1;
    if ((data_[r] & (1ull<<k)) == 0) {
      data_[next_r] = 0;
      data_[r] |= 1ull<<k;
    }
    r = next_r;
  }

  data_[r] ^= 1ull << (i & MASK);

  for (int h = 1; h < height_; ++h) {
    if (data_[r] != 0)
      break;
    int next_r = (r - 1) >> 6;
    int k = (i >> 6 * h) & MASK;
    data_[next_r] ^= (1ull << k);
    r = next_r;
  }
}

std::vector<index> bit_tree_column::to_vector() const {
  std::vector<index> indices;

  if (data_[0] == 0)
    return indices;
  
  retrieve(height_ - 1, 0, 0, &indices);
  return indices;
}

void bit_tree_column::retrieve(int h, int r, index i, std::vector<index>* indices) const {
  using boost::core::bit_width;
  using boost::core::countr_zero;

  uint64_t d = data_[r];
  assert(d != 0);

  uint64_t high = bit_width(d);
  uint64_t low = countr_zero(d);

  for (uint64_t k = low; k < high; ++k) {
    if (d & (1ull << k)) {
      if (h == 0) {
        indices->push_back((i << 6) + k);
      } else {
        retrieve(h - 1, (r << 6) + k + 1, (i << 6) + k, indices);
      }
    }
  }
}

#ifdef MPHHC_UNITTEST
const std::vector<uint64_t>& bit_tree_column::data() const {
  return data_;
}
#endif

}
