//
// Vortexje -- std::vector using eigen's aligned_allocator
//

#ifndef __ALIGNED_VECTOR_HPP__
#define __ALIGNED_VECTOR_HPP__

#include <vector>
#include <Eigen/StdVector>

namespace Vortexje
{
  template <typename T>
  using vector_aligned = std::vector<T, Eigen::aligned_allocator<T> >;
};

#endif // __ALIGNED_VECTOR_HPP__
