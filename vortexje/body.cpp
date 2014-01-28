//
// Vortexje -- Body.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <vortexje/body.hpp>

#include <iostream>

using namespace std;
using namespace Eigen;
using namespace Vortexje;
    
/**
   Constructs a new Body.
   
   @param[in]   id   Name for this body.
*/
Body::Body(string id) : id(id)
{
    // Initialize kinematics:
    position = Vector3d(0, 0, 0);
    velocity = Vector3d(0, 0, 0);
    
    attitude = Quaterniond(1, 0, 0, 0);
    rotational_velocity = Vector3d(0, 0, 0);
}

/**
   Body destructor.
*/
Body::~Body()
{
    for (int i = 0; i < (int) wakes.size(); i++)
        delete wakes[i];
}

/**
   Adds a non-lifting surface to this body.
   
   @param[in]   surface    Non-lifting surface.
*/
void
Body::add_non_lifting_surface(Surface *non_lifting_surface)
{
    // Non-lifting surfaces go in the front.
    non_lifting_surfaces.push_back(non_lifting_surface);
    non_wake_surfaces.push_back(non_lifting_surface);
}

/**
   Adds a lifting surface to this body.
   
   @param[in]   lifting_surface    Lifting surface.
*/
void
Body::add_lifting_surface(LiftingSurface *lifting_surface)
{
    // Lifting surfaces go in the back.
    lifting_surfaces.push_back(lifting_surface);
    non_wake_surfaces.push_back(lifting_surface);
        
    Wake *wake = new Wake(*lifting_surface);
    wakes.push_back(wake);
}

/**
   Sets the linear position of this body.
   
   @param[in]   position    Linear position.
*/
void
Body::set_position(const Vector3d &position)
{
    // Compute differential:
    Vector3d dposition = position - this->position;
       
    // Apply for all non-wake surfacees:
    for (int i = 0; i < (int) non_wake_surfaces.size(); i++) {
        Surface *surface = non_wake_surfaces[i];
        
        // Translate:
        surface->translate(dposition);
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
   Sets the attitude (orientation) of this body.
   
   @param[in]   attitude    Attitude (orientation) of this body, as normalized quaternion.
*/
void
Body::set_attitude(const Quaterniond &attitude)
{   
    Eigen::Vector3d translation;
    Eigen::Matrix3d transformation;
    
    // Apply for all non-wake surfacees:
    for (int i = 0; i < (int) non_wake_surfaces.size(); i++) {
        Surface *surface = non_wake_surfaces[i];
        
        // Translate to origin:
        translation = -position;
        surface->translate(translation);
        
        // Transform to canonical orientation:
        transformation = this->attitude.inverse().toRotationMatrix();
        surface->transform(transformation);
        
        // Transform to new orientation:
        transformation = attitude.toRotationMatrix();
        surface->transform(transformation);
        
        // Translate back:
        surface->translate(position);
    }
    
    // Apply for trailing edge wake nodes:
    for (int i = 0; i < (int) wakes.size(); i++) {
        Wake *wake = wakes[i];
        
        // Translate:
        translation = -position;
        wake->translate_trailing_edge(translation);
        
        // Transform to canonical orientation:
        transformation = this->attitude.inverse().toRotationMatrix();
        wake->transform_trailing_edge(transformation);
        
        // Transform to new orientation:
        transformation = attitude.toRotationMatrix();
        wake->transform_trailing_edge(transformation);
        
        // Translate back:
        wake->translate_trailing_edge(position);
    }
    
    // Update state:
    this->attitude = attitude;
}

/**
   Sets the linear velocity of this body.
   
   @param[in]   velocity    Linear velocity.
*/
void 
Body::set_velocity(const Vector3d &velocity)
{
    this->velocity = velocity;
}

/**
   Sets the rotational velocity of this body.
   
   @param[in]   rotational_velocity     Rotational velocity.
*/
void
Body::set_rotational_velocity(const Vector3d &rotational_velocity)
{
    this->rotational_velocity = rotational_velocity;
}

/**
   Computes the kinematic velocity of the given panel.
   
   @param[in]   surface   Surface, belonging to this body. 
   @param[in]   panel  Panel, belonging to this surface.
   
   @return The kinematic velocity.
*/
Vector3d
Body::panel_kinematic_velocity(const Surface &surface, int panel) const
{
    Vector3d panel_position = surface.panel_collocation_point(panel, false);
    Vector3d r = panel_position - position;
    return velocity + rotational_velocity.cross(r);
}

/**
   Computes the kinematic velocity of the given node.
   
   @param[in]   surface   Surface, belonging to this body. 
   @param[in]   node   Node, belonging to this surface.
   
   @return The kinematic velocity.
*/
Vector3d
Body::node_kinematic_velocity(const Surface &surface, int node) const
{
    Vector3d r = surface.nodes[node] - position;
    return velocity + rotational_velocity.cross(r);
}
