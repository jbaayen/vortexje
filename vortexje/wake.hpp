//
// Vortexje -- Wake.
//
// Copyright (C) 2012 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __WAKE_HPP__
#define __WAKE_HPP__

#include <vector>

#include <vortexje/mesh.hpp>
#include <vortexje/wing.hpp>

namespace Vortexje
{

/**
   Representation of a wake, i.e., a vortex sheet.
   
   @brief Wake representation.
*/
class Wake : public Mesh
{
public:
    Wake(Wing &wing);
    
    /**
       Associated Wing.
    */
    Wing &wing;
    
    void add_layer(std::vector<Mesh*> &other_meshes);
    
    void translate_trailing_edge(Eigen::Vector3d translation);
    void transform_trailing_edge(Eigen::Matrix3d transformation);
    
    void update_ramasamy_leishman_vortex_core_radii(int panel, double dt);
    
    /**
       Strengths of the doublet, or vortex ring, panels.
    */
    std::vector<double> doublet_coefficients; 
    
    /**
       Radii of the vortex filaments forming the vortex rings.
    */ 
    std::vector<std::vector<double> > vortex_core_radii;
    
    /**
       Initial lengths of the vortex filaments forming the vortex rings.
    */
    std::vector<std::vector<double> > base_edge_lengths;
};

};

#endif // __WAKE_HPP__
