//
// Vortexje -- Lifting surface.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __LIFTING_SURFACE_HPP__
#define __LIFTING_SURFACE_HPP__

#include <Eigen/Core>

#include <vortexje/surface.hpp>

namespace Vortexje
{

/**
   Representation of lifting surface.
   
   @brief Lifting surface representation.
*/
class LiftingSurface : public Surface
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    LiftingSurface();
         
    /**
       Nodes on the upper side of the surface.  The first dimension is the chordwise direction;  the second the spanwise direction.
    */
    Eigen::MatrixXi upper_nodes;
    
    /**
       Nodes on the lower side of the surface.  The first dimension is the chordwise direction;  the second the spanwise direction.
    */
    Eigen::MatrixXi lower_nodes;
    
    /**
       Panels on the upper side of the surface.  The first dimension is the chordwise direction;  the second the spanwise direction.
    */
    Eigen::MatrixXi upper_panels;
    
    /**
       Panels on the lower side of the surface.  The first dimension is the chordwise direction;  the second the spanwise direction.
    */
    Eigen::MatrixXi lower_panels;
    
    int n_chordwise_nodes() const;
    int n_chordwise_panels() const;
    
    int n_spanwise_nodes() const;
    int n_spanwise_panels() const;
    
    int trailing_edge_node(int index) const;
    int trailing_edge_upper_panel(int index) const;
    int trailing_edge_lower_panel(int index) const;
    
    Eigen::Vector3d trailing_edge_bisector(int node_index) const;
    
    virtual Eigen::Vector3d wake_emission_velocity(const Eigen::Vector3d &apparent_velocity, int node_index) const;
};

};

#endif // __LIFTING_SURFACE_HPP__
