//
// Vortexje -- Collection.
//
// Copyright (C) 2012 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <vortexje/collection.hpp>

#include <iostream>

using namespace std;
using namespace Eigen;
using namespace Vortexje;
    
/**
   Constructs a new Collection.
   
   @param[in]   id              Name for this collection.
   @param[in]   nolift_mesh     Obstacle-only mesh.
*/
Collection::Collection(string id, Mesh &nolift_mesh) : id(id), nolift_mesh(nolift_mesh)
{
    // Initialize mesh list:
    meshes_without_wakes.push_back(&nolift_mesh);
    
    // Initialize kinematics:
    position = Vector3d(0, 0, 0);
    velocity = Vector3d(0, 0, 0);
    
    attitude = Quaterniond(1, 0, 0, 0);
    rotational_velocity = Vector3d(0, 0, 0);
}

/**
   Collection destructor.
*/
Collection::~Collection()
{
}

/**
   Adds a Wing to this collection.
   
   @param[in]   wing    Wing.
*/
void
Collection::add_wing(Wing *wing)
{
    wings.push_back(wing);
    meshes_without_wakes.push_back(wing);
        
    Wake *wake = new Wake(*wing);
    wakes.push_back(wake);
}

/**
   Sets the linear position of this collection.
   
   @param[in]   position    Linear position.
*/
void
Collection::set_position(Vector3d position)
{
    // Compute differential:
    Vector3d dposition = position - this->position;
       
    // Apply for all non-wake meshes:
    for (int i = 0; i < (int) meshes_without_wakes.size(); i++) {
        Mesh *mesh = meshes_without_wakes[i];
        
        // Translate:
        mesh->translate(dposition, meshes_without_wakes);
    }
    
    // Apply for trailing edge wake nodes:
    for (int i = 0; i < (int) wakes.size(); i++) {
        Wake *wake = wakes[i];
        
        // Translate:
        wake->translate_trailing_edge(dposition);
    }
    
    // Update state:
    this->position = position;
}

/**
   Sets the attitude (orientation) of this collection.
   
   @param[in]   attitude    Attitude (orientation) of this collection, as normalized quaternion.
*/
void
Collection::set_attitude(Quaterniond attitude)
{   
    // Apply for all non-wake meshes:
    for (int i = 0; i < (int) meshes_without_wakes.size(); i++) {
        Mesh *mesh = meshes_without_wakes[i];
        
        // Translate to origin:
        mesh->translate(-position, meshes_without_wakes);
        
        // Transform to canonical orientation:
        mesh->transform(this->attitude.inverse().toRotationMatrix(), meshes_without_wakes);
        
        // Transform to new orientation:
        mesh->transform(attitude.toRotationMatrix(), meshes_without_wakes);
        
        // Translate back:
        mesh->translate(position, meshes_without_wakes);
    }
    
    // Apply for trailing edge wake nodes:
    for (int i = 0; i < (int) wakes.size(); i++) {
        Wake *wake = wakes[i];
        
        // Translate:
        wake->translate_trailing_edge(-position);
        
        // Transform to canonical orientation:
        wake->transform_trailing_edge(this->attitude.inverse().toRotationMatrix());
        
        // Transform to new orientation:
        wake->transform_trailing_edge(attitude.toRotationMatrix());
        
        // Translate back:
        wake->translate_trailing_edge(position);
    }
    
    // Update state:
    this->attitude = attitude;
}

/**
   Sets the linear velocity of this collection.
   
   @param[in]   velocity    Linear velocity.
*/
void 
Collection::set_velocity(Vector3d velocity)
{
    this->velocity = velocity;
}

/**
   Sets the rotational velocity of this collection.
   
   @param[in]   rotational_velocity     Rotational velocity.
*/
void
Collection::set_rotational_velocity(Vector3d rotational_velocity)
{
    this->rotational_velocity = rotational_velocity;
}

/**
   Computes the kinematic velocity of the given panel.
   
   @param[in]   mesh   Mesh, belonging to this collection. 
   @param[in]   panel  Panel, belonging to this mesh.
   
   @return The kinematic velocity.
*/
Vector3d
Collection::panel_kinematic_velocity(Mesh &mesh, int panel)
{
    Vector3d panel_position = mesh.panel_collocation_point(panel, false);
    Vector3d r = panel_position - position;
    return velocity + rotational_velocity.cross(r);
}

/**
   Computes the kinematic velocity of the given node.
   
   @param[in]   mesh   Mesh, belonging to this collection. 
   @param[in]   node   Node, belonging to this mesh.
   
   @return The kinematic velocity.
*/
Vector3d
Collection::node_kinematic_velocity(Mesh &mesh, int node)
{
    Vector3d r = mesh.nodes[node] - position;
    return velocity + rotational_velocity.cross(r);
}
