//
// Vortexje -- NACA 4-series airfoil generator.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __NACA4_GENERATOR_HPP__
#define __NACA4_GENERATOR_HPP__

#include <Eigen/Core>
#include <Eigen/StdVector>

namespace Vortexje
{

namespace ShapeGenerators
{

namespace Airfoils
{

/**
   Class for generation of NACA 4-series airfoils.
   
   @brief NACA-4 series airfoil generation.
*/
class NACA4Generator
{
public:
    static std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > generate(double max_camber, double max_camber_dist, double max_thickness, bool finite_te_thickness, double chord, int n_points, int &trailing_edge_point_id);
};

};

};

};

#endif // __NACA4_GENERATOR_HPP__
