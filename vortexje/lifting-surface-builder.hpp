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

#include <vortexje/lifting-surface.hpp>
#include <vortexje/surface-builder.hpp>

namespace Vortexje
{

/**
   Class for construction of lifting surfaces by joining 2-D airfoils together.
   
   @brief Lifting surface construction helper.
*/
class LiftingSurfaceBuilder : public SurfaceBuilder
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    LiftingSurfaceBuilder(LiftingSurface &lifting_surface);
    
    /**
       LiftingSurface under construction.
    */
    LiftingSurface &lifting_surface;
    
    std::vector<int> create_panels_inside_airfoil(const std::vector<int> &airfoil_nodes, int trailing_edge_point_id, int z_sign);
    
    void finish(const std::vector<std::vector<int> > &node_strips, const std::vector<std::vector<int> > &panel_strips, int trailing_edge_point_id);
};

};

#endif // __LIFTING_SURFACE_BUILDER_HPP__
