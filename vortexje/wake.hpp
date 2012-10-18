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

class Wake : public Mesh
{
public:
    // Constructor:
    Wing &wing;
    Wake(Wing &wing);
    
    // Add layer of wake panels:
    void add_layer(std::vector<Mesh*> &meshes_without_wakes);
    
    // Geometry of trailing edge nodes:
    void translate_trailing_edge(Eigen::Vector3d translation);
    void transform_trailing_edge(Eigen::Matrix3d transformation);
    
    // Update Ramasamy-Leishman vortex core radii:
    void update_ramasamy_leishman_vortex_core_radii(int panel, double dt);
    
    // Wake panel states:
    std::vector<double> doublet_coefficients;    
    std::vector<std::vector<double> > vortex_core_radii;
    std::vector<std::vector<double> > base_edge_lengths;
};

};

#endif // __WAKE_HPP__
