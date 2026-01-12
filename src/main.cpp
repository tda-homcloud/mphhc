#include <iostream>

#include "mphhc/core.hpp"

int main(int argc, char** argv) {
  using namespace mphhc;

  std::cout << GetVersion() << std::endl;

  return 0;
}
