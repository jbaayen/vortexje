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
#include <Eigen/StdVector>

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

    void sort_strips();
    
    /**
       Unit vector pointing along the LiftingSurface chord.
    */
    Eigen::Vector3d chord_direction;
    
    /**
       Unit vector pointing along the LiftingSurface span.
    */
    Eigen::Vector3d span_direction;
         
    /**
       Chord-wise strips of panels on the upper side of the surface.
    */
    std::vector<std::vector<int> > upper_panel_strips;
    
    /**
       Chord-wise strips of panels on the lower side of the surface.
    */
    std::vector<std::vector<int> > lower_panel_strips;
    
    /**
       List of panel numbers that are located above and adjecent to the trailing edge.
     */
    std::vector<int> upper_trailing_edge_panels;

    /**
       List of panel numbers that are located below and adjecent to the trailing edge   
    */
    std::vector<int> lower_trailing_edge_panels; 
    
    /**
       List of node numbers that form the trailing edge.
    */
    std::vector<int> trailing_edge_nodes;
    
    void transform(const Eigen::Matrix3d &transformation);
    
    bool closest_panel(const Eigen::Vector3d &x, int &panel, double &distance) const;
    
    Eigen::Vector3d close_to_body_point(int node) const;
};

};

#endif // __LIFTING_SURFACE_HPP__
