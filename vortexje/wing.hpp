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

/**
   Representation of a wake-emitting, or source-doublet, mesh.
   
   This class would typically be used for representing a blade or a wing, whence the name.
   
   @brief Wing representation.
*/
class Wing : public Mesh
{
public:
    Wing();
    Wing(Mesh &mesh, Eigen::Vector3d &location, Eigen::Vector3d &chord_direction, Eigen::Vector3d &top_direction, Eigen::Vector3d &span_direction);

    void sort_trailing_edge();
    
    /**
       Reference location.
    */
    Eigen::Vector3d location;
    
    /**
       Unit vector pointing along the wing chord.
    */
    Eigen::Vector3d chord_direction;
    
    /**
       Unit vector pointing upwards, out of, and perpendicular to the wing.
    */
    Eigen::Vector3d top_direction;
    
    /**
       Unit vector pointing along the wing span.
    */
    Eigen::Vector3d span_direction;
    
    /**
       List of node numbers that form the trailing edge.
    */
    std::vector<int> trailing_edge_nodes;
    
    /**
       List of panel numbers that are located above and adjecent to the trailing edge.
    */
    std::vector<int> trailing_edge_top_panels;
    
    /**
       List of panel numbers that are located below and adjecent to the trailing edge.
    */
    std::vector<int> trailing_edge_bottom_panels;
    
    void translate(Eigen::Vector3d &translation);
    void translate(Eigen::Vector3d &translation, std::vector<Mesh*> &corotating_meshes);
    
    void transform(Eigen::Matrix3d &transformation);
    void transform(Eigen::Matrix3d &transformation, std::vector<Mesh*> &corotating_meshes);
    
    bool closest_panel(Eigen::Vector3d &x, int &panel, double &distance);
    
    Eigen::Vector3d close_to_body_point(int node);
};

};

#endif // __WING_HPP__
