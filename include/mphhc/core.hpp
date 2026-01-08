#pragma once
#include <string>
#include <vector>
#include <cstdint>
#include <boost/core/bit.hpp>

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
  index add_dim_col(int dim, column&& col);
  index add_dim_col(int dim, const column& col);
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
  
  std::vector<bitset64> data_;
  int num_index_;
  int height_;

  void retrieve(int h, int r, int i, std::vector<index>* indices) const;

 public:
  static int compute_height(int num_index);
  static int compute_data_size(int height, int num_level);
  
  explicit bit_tree_column(int num_index);
  void import_column(const column& vec);
  void clear();
  void set(index i);
  void set_xor(index i);
  index max() const;
  bool none() const;
  void add(const column& other);
  column export_column() const;
};

}
