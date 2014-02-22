//
// Vortexje -- Wake.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __WAKE_HPP__
#define __WAKE_HPP__

#include <Eigen/Geometry>
#include <Eigen/StdVector>

#include <vortexje/lifting-surface.hpp>

namespace Vortexje
{

/**
   Representation of a wake, i.e., a vortex sheet.
   
   @brief Wake representation.
*/
class Wake : public Surface
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    Wake(LiftingSurface &lifting_surface);
    
    /**
       Associated lifting surface.
    */
    LiftingSurface &lifting_surface;
    
    void add_layer();
    
    void translate_trailing_edge(const Eigen::Vector3d &translation);
    void transform_trailing_edge(const Eigen::Matrix3d &transformation);
    void transform_trailing_edge(const Eigen::Transform<double, 3, Eigen::Affine> &transformation);
    
    void update_ramasamy_leishman_vortex_core_radii(int panel, double dt);
    
    /**
       Strengths of the doublet, or vortex ring, panels.
    */
    std::vector<double> doublet_coefficients; 
    
    /**
       Radii of the vortex filaments forming the vortex rings.
    */ 
    std::vector<std::vector<double> > vortex_core_radii;
  
private:  
    /**
       Initial lengths of the vortex filaments forming the vortex rings.
    */
    std::vector<std::vector<double> > base_edge_lengths;
};

};

#endif // __WAKE_HPP__
