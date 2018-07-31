//
// Vortexje -- NACA 4-series airfoil generator.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __NACA4_AIRFOIL_GENERATOR_HPP__
#define __NACA4_AIRFOIL_GENERATOR_HPP__

#include <Eigen/Core>

#include <vortexje/vector-aligned.hpp>

namespace Vortexje
{

/**
   Class for generation of NACA 4-digit series airfoils.
   
   @brief NACA 4-digit series airfoil generation.
*/
class NACA4AirfoilGenerator
{
public:
    static vector_aligned<Eigen::Vector3d> generate(double max_camber, double max_camber_dist, double max_thickness, bool finite_te_thickness, double chord, int n_points, int &trailing_edge_point_id);
};

};

#endif // __NACA4_AIRFOIL_GENERATOR_HPP__
