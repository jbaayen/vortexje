//
// Vortexje -- Lifting surface builder.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __LIFTING_SURFACE_BUILDER_HPP__
#define __LIFTING_SURFACE_BUILDER_HPP__

#include <Eigen/Core>
#include <Eigen/StdVector>

#include <vortexje/lifting-surface.hpp>
#include <vortexje/surface-builder.hpp>

namespace Vortexje
{

/**
   Class for construction of lifting surfaces by joining 2-D airfoils together.
   
   @brief Lifting surface construction helper.
*/
class LiftingSurfaceBuilder : protected SurfaceBuilder
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    LiftingSurfaceBuilder(LiftingSurface &lifting_surface);
    
    /**
       LiftingSurface under construction.
    */
    LiftingSurface &lifting_surface;
    
    std::vector<int> create_nodes(const std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > &airfoil_points, int trailing_edge_point_id);
    
    std::vector<int> create_panels_between(const std::vector<int> &first_airfoil_nodes, const std::vector<int> &second_airfoil_nodes, int trailing_edge_point_id);
    
    std::vector<int> create_panels_inside(const std::vector<int> &airfoil_nodes, int trailing_edge_point_id, int z_sign);
};

};

#endif // __LIFTING_SURFACE_BUILDER_HPP__
