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

/**
   Container for wake-emitting Wing meshes, as well as for an obstacle-only mesh.
   
   This class is designed for the application of kinematic operations to a collection of meshes at once.
   A typical application would be an airplane consisting of a fuselage and wake-emitting wings.
   
   @brief Mesh container.
*/
class Collection
{
public:
    /**
       Collection name.
    */
    std::string        id;
    
    /**
       Drag-only mesh.
    */
    Mesh              &nolift_mesh;
    
    /**
       List of wings.
    */
    std::vector<Wing*> wings;
    
    /**
       List of wakes associated with the wings, in the same order.
    */
    std::vector<Wake*> wakes;
    
    Collection(std::string id,
               Mesh        &nolift_mesh);
         
    ~Collection();
         
    virtual void add_wing(Wing *wing);
    
    /**
       Linear position of the entire collection.
    */
    Eigen::Vector3d position;
    
    /**
       Linear velocity of the entire collection.
    */
    Eigen::Vector3d velocity;
    
    /**
       Attitude (orientation) of the entire collection.
    */
    Eigen::Quaterniond attitude;
    
    /**
       Rotational velocity of the entire collection.
    */
    Eigen::Vector3d rotational_velocity;

    void set_position(Eigen::Vector3d position);
    void set_attitude(Eigen::Quaterniond attitude);
    
    void set_velocity(Eigen::Vector3d velocity);
    void set_rotational_velocity(Eigen::Vector3d rotational_velocity);
    
    Eigen::Vector3d panel_kinematic_velocity(Mesh &mesh, int panel);
    
    Eigen::Vector3d node_kinematic_velocity(Mesh &mesh, int node);
    
protected:
    /**
       List of all non-wake meshes.  For internal use.
    */
    std::vector<Mesh*> meshes_without_wakes;
};

};

#endif // __COLLECTION_HPP__
