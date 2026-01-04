#pragma once
#include <string>
#include <vector>
#include <cstdint>

namespace mphhc {

std::string get_version();
int add(int a, int b);

using index = int32_t;
using column = std::vector<index>;


class boundary_matrix {
  using local_index = int32_t;
  
  struct index_info {
    int dim_;
    local_index nth_;
  };
  
  std::vector<std::vector<column>> columns_;
  std::vector<std::vector<index>> local_to_global_index_;
  std::vector<index_info> global_to_local_index_;

 public:
  boundary_matrix(int maxdim);
  index add_dim_col(int dim, column&& col);
  index add_dim_col(int dim, const column& col);
  int num_simplices() const;
  
};


class bit_tree_column {
  std::vector<uint64_t> data_;
  int num_index_;
  int height_;

 public:
  static int compute_height(int num_index);
  static int compute_data_size(int height, int num_level);
  
  bit_tree_column(int num_index);
  void clear();
  void set(index i);
};



}
