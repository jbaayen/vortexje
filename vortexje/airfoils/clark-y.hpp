//
// Vortexje -- Clark-Y airfoil generator.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __CLARK_Y_HPP__
#define __CLARK_Y_HPP__

#include <Eigen/Core>
#include <Eigen/StdVector>

namespace Vortexje
{

namespace Airfoils
{

/**
   Class for generation of Clark-Y airfoils.
   
   @brief ClarkY airfoil generation.
*/
class ClarkY
{
public:
    static std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > generate(double chord, int n_points, int &trailing_edge_point_id);
};

};

};

#endif // __CLARK_Y_HPP__
