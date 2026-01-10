#pragma once
#include <string>
#include <vector>
#include <cstdint>
#include <boost/core/bit.hpp>
#include <iostream>

namespace mphhc {

std::string get_version();

using index = int32_t;
using column = std::vector<index>;
struct birth_death_pair {
  int dim;
  index birth, death;

  birth_death_pair(int dim, index birth, index death);
  bool operator==(const birth_death_pair&) const;
  bool operator<(const birth_death_pair&) const;
  friend void PrintTo(const birth_death_pair& pair, std::ostream* os); 
};

class boundary_matrix {
  using local_index = int32_t;
  
  struct index_info {
    int dim;
    local_index nth;
  };
  
  std::vector<std::vector<column>> columns_;
  std::vector<std::vector<index>> local_to_global_index_;
  std::vector<index_info> global_to_local_index_;
  bool reduced_;
  
 public:
  boundary_matrix(int maxdim);
  int max_dim() const;
  index set_dim_col(index i, int dim, column&& col);
  index set_dim_col(index i, int dim, const column& col);
  int num_simplices() const;
  bool is_reduced() const;

  void reduce_standard();
  void reduce_twist();
  std::vector<birth_death_pair> birth_death_pairs() const;
};

class bitset64 {
 public:
  uint64_t data;

  inline bitset64(): data(0ull) {}
  
  inline void clear() { data = 0ull; }
  inline void flip(int pos) { data ^= (1ull << pos); }
  inline void set(int pos) { data |= (1ull << pos); }

  inline bool test(int pos) const { return (1ull << pos) & data; }
  inline bool none() const { return data == 0ull; }
  inline bool any() const { return data != 0ull; }
  inline int max() const { return static_cast<int>(boost::core::bit_width(data)) - 1; }
  inline int min() const { return boost::core::countr_zero(data); }
};
  
class bit_tree_column {
  static const uint64_t MASK = (1 << 6) - 1;

  static constexpr int NODE_BLOCK_SIZE_TABLE[6] = {
    -1,
    0,
    1,
    1 + (1<<6),
    1 + (1<<6) + (1<<12),
    1 + (1<<6) + (1<<12) + (1<<18)
  };
  
  std::vector<bitset64> data_;
  int num_index_;
  int height_;
  int node_block_size_;

  void retrieve(int h, int r, int i, std::vector<index>* indices) const;

 public:
  static int compute_height(int num_index);
  static int compute_data_size(int height, int num_level);
  
  explicit bit_tree_column(int num_index);
  void import_column(const column& vec);
  
  inline void set_nodes(index i, int r) {
    for (int h = 1; h <= height_ - 1; ++h) {
      r = (r - 1) >> 6;
      int k = (i >> (6 * h)) & MASK;
      if (data_[r].test(k))
        return;
      else 
        data_[r].set(k);
    }
  }

  inline void set(index i) {
    int r = (i >> 6) + node_block_size_;
    int k = i & MASK;
    data_[r].set(k);
    set_nodes(i, r);
  }

  inline void unset_nodes(index i, int r) {
    for (int h = 1; h <= height_ - 1; ++h) {
      if (data_[r].any())
        return;
      r = (r - 1) >> 6;
      int k = (i >> (6 * h)) & MASK;
      data_[r].flip(k);
    }
  }

  inline void set_xor(index i) {
    int r = (i >> 6) + node_block_size_;
    int k = i & MASK;
    data_[r].flip(k);
    
    if (data_[r].test(k)) {
      set_nodes(i, r);
    } else {
      unset_nodes(i, r);
    }
  }

  inline index max() const {
    using boost::core::bit_width;
  
    if (data_[0].none())
      return -1;
  
    int r = 0;
    index i = 0;
  
    for (int h = height_ - 1; h >= 1; --h) {
      int k = data_[r].max();
      r = (r << 6) + k + 1;
      i = (i << 6) + k;
    }
    return (i << 6) + data_[r].max();
  }

  inline bool none() const { return data_[0].none(); }

  inline void add(const column& other) {
    for (index i: other)
      set_xor(i);
  }

  column debug_export_column();
  void export_and_clear_column(column* col);
};

}
