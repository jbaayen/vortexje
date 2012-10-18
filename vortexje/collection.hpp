//
// Vortexje -- Collection.
//
// Copyright (C) 2012 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __COLLECTION_HPP__
#define __COLLECTION_HPP__

#include <vector>
#include <string>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <vortexje/mesh.hpp>
#include <vortexje/wing.hpp>
#include <vortexje/wake.hpp>

namespace Vortexje
{

class Collection
{
public:
    // Constructor:
    std::string        id;
    Mesh              &nolift_mesh;
    std::vector<Wing*> wings;
    std::vector<Wake*> wakes;
    
    Collection(std::string id,
               Mesh        &nolift_mesh);
         
    // Destructor:
    ~Collection();
         
    // Wing management:
    virtual void add_wing(Wing *wing);
    
    // Kinematics:
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
    
    Eigen::Quaterniond attitude;
    Eigen::Vector3d rotational_velocity;

    void set_position(Eigen::Vector3d position);
    void set_attitude(Eigen::Quaterniond attitude);
    
    void set_velocity(Eigen::Vector3d velocity);
    void set_rotational_velocity(Eigen::Vector3d rotational_velocity);
    
    // Panel velocity:
    Eigen::Vector3d panel_kinematic_velocity(Mesh &mesh, int panel);
    
    // Node velocity:
    Eigen::Vector3d node_kinematic_velocity(Mesh &mesh, int node);
    
protected:
    // All meshes without wakes:
    std::vector<Mesh*> meshes_without_wakes;
};

};

#endif // __COLLECTION_HPP__
