//
// Vortexje -- Ellipse generator.
//
// Copyright (C) 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __ELLIPSE_GENERATOR_HPP__
#define __ELLIPSE_GENERATOR_HPP__

#include <Eigen/Core>

#include <vortexje/vector-aligned.hpp>

namespace Vortexje
{

/**
   Class for generation of ellipses.
   
   @brief Ellipse generation.
*/
class EllipseGenerator
{
public:
    static vector_aligned<Eigen::Vector3d> generate(double a, double b, int n_points);
};

};

#endif // __ELLIPSE_GENERATOR_HPP__
