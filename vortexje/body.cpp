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
    vector<Wake*>::iterator wi;
    for (wi = wakes.begin(); wi != wakes.end(); wi++)
        delete (*wi);
}

/**
   Adds a non-lifting surface to this body.
   
   @param[in]   surface   Non-lifting surface.
*/
void
Body::add_non_lifting_surface(Surface *non_lifting_surface)
{
    non_lifting_surfaces.push_back(non_lifting_surface);
    non_wake_surfaces.push_back(non_lifting_surface);
}

/**
   Adds a lifting surface to this body.
   
   @param[in]   lifting_surface   Lifting surface.
*/
void
Body::add_lifting_surface(LiftingSurface *lifting_surface)
{
    lifting_surfaces.push_back(lifting_surface);
    non_wake_surfaces.push_back(lifting_surface);
        
    Wake *wake = new Wake(*lifting_surface);
    wakes.push_back(wake);
}

/**
   Sets the linear position of this body.
   
   @param[in]   position   Linear position.
*/
void
Body::set_position(const Vector3d &position)
{
    // Compute differential translation:
    Vector3d translation = position - this->position;
       
    // Apply for all non-wake surfacees:
    vector<Surface*>::iterator si;
    for (si = non_wake_surfaces.begin(); si != non_wake_surfaces.end(); si++)
        (*si)->translate(translation);

    // Apply to trailing edge wake nodes:
    vector<Wake*>::iterator wi;
    for (wi = wakes.begin(); wi != wakes.end(); wi++)
        (*wi)->translate_trailing_edge(translation);
    
    // Update state:
    this->position = position;
}

/**
   Sets the attitude (orientation) of this body.
   
   @param[in]   attitude   Attitude (orientation) of this body, as normalized quaternion.
*/
void
Body::set_attitude(const Quaterniond &attitude)
{   
    // Compute differential transformation:
    Transform<double, 3, Affine> transformation = Translation<double, 3>(position) * attitude * this->attitude.inverse() * Translation<double, 3>(-position);
    
    // Apply for all non-wake surfacees:
    vector<Surface*>::iterator si;
    for (si = non_wake_surfaces.begin(); si != non_wake_surfaces.end(); si++)
        (*si)->transform(transformation);
    
    // Apply for trailing edge wake nodes:
    vector<Wake*>::iterator wi;
    for (wi = wakes.begin(); wi != wakes.end(); wi++)
        (*wi)->transform_trailing_edge(transformation);
    
    // Update state:
    this->attitude = attitude;
}

/**
   Sets the linear velocity of this body.
   
   @param[in]   velocity   Linear velocity.
*/
void 
Body::set_velocity(const Vector3d &velocity)
{
    this->velocity = velocity;
}

/**
   Sets the rotational velocity of this body.
   
   @param[in]   rotational_velocity   Rotational velocity.
*/
void
Body::set_rotational_velocity(const Vector3d &rotational_velocity)
{
    this->rotational_velocity = rotational_velocity;
}

/**
   Computes the kinematic velocity of the given panel.
   
   @param[in]   surface   Surface, belonging to this body. 
   @param[in]   panel     Panel, belonging to this surface.
   
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
   @param[in]   node      Node, belonging to this surface.
   
   @return The kinematic velocity.
*/
Vector3d
Body::node_kinematic_velocity(const Surface &surface, int node) const
{
    Vector3d r = surface.nodes[node] - position;
    return velocity + rotational_velocity.cross(r);
}
