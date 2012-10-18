//
// Vortexje -- Wing.
//
// Copyright (C) 2012 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __WING_HPP__
#define __WING_HPP__

#include <vector>

#include <Eigen/Core>

#include <vortexje/mesh.hpp>

namespace Vortexje
{

class Wing : public Mesh
{
public:
    // Constructor:
    Wing();
    Wing(Mesh &mesh, Eigen::Vector3d location, Eigen::Vector3d chord_direction, Eigen::Vector3d top_direction, Eigen::Vector3d span_direction);

    // Wing construction:
    void sort_trailing_edge();
    
    // Wing-fixed frame:
    Eigen::Vector3d location;
    
    Eigen::Vector3d chord_direction;
    Eigen::Vector3d top_direction;
    Eigen::Vector3d span_direction;
    
    // Trailing edge:
    std::vector<int> trailing_edge_nodes;
    std::vector<int> trailing_edge_top_panels;
    std::vector<int> trailing_edge_bottom_panels;
    
    // Re-implemented virtual methods:
    void translate(Eigen::Vector3d translation);
    void translate(Eigen::Vector3d translation, std::vector<Mesh*> &corotating_meshes);
    
    void transform(Eigen::Matrix3d transformation);
    void transform(Eigen::Matrix3d transformation, std::vector<Mesh*> &corotating_meshes);
    
    bool closest_panel(Eigen::Vector3d x, int &panel, double &distance);
    
    Eigen::Vector3d close_to_body_point(int node);
};

};

#endif // __WING_HPP__
