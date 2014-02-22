//
// Vortexje -- Body.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __BODY_HPP__
#define __BODY_HPP__

#include <string>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>

#include <vortexje/surface.hpp>
#include <vortexje/lifting-surface.hpp>
#include <vortexje/wake.hpp>
#include <vortexje/boundary-layer.hpp>

namespace Vortexje
{

/**
   Container for surfaces.
   
   This class is designed for the application of kinematic operations to a body consisting of multiple surfaces.
   A typical application would be an airplane consisting of a fuselage and lifting surfaces.
   
   @brief Surface container.
*/
class Body
{
public:
    /**
       Body name.
    */
    std::string id;
    
    /**
       Data structure grouping a non-lifting surface with its boundary layer.
    */
    class SurfaceData {
    public:
        SurfaceData(Surface &surface, BoundaryLayer &boundary_layer) :
            surface(surface), boundary_layer(boundary_layer) {}
        
        Surface &surface;
        
        BoundaryLayer &boundary_layer;
    };
    
    /**
       Data structure grouping a lifting surface with its boundary layer, and with its wake.
    */
    class LiftingSurfaceData : public SurfaceData {
    public:
        LiftingSurfaceData(LiftingSurface &lifting_surface, BoundaryLayer &boundary_layer) :
            SurfaceData(lifting_surface, boundary_layer), lifting_surface(lifting_surface), wake(Wake(lifting_surface)) { }
        
        LiftingSurface &lifting_surface;
        
        Wake wake;
    };
    
    /**
       List of non-lifting surfaces.
    */
    std::vector<SurfaceData*> non_lifting_surfaces;
    
    /**
       List of lifting surfaces.
    */
    std::vector<LiftingSurfaceData*> lifting_surfaces;
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    Body(const std::string &id);
    
    ~Body();

    void add_non_lifting_surface(Surface &surface);
    void add_non_lifting_surface(Surface &surface, BoundaryLayer &boundary_layer);
    
    void add_lifting_surface(LiftingSurface &lifting_surface);    
    void add_lifting_surface(LiftingSurface &lifting_surface, BoundaryLayer &boundary_layer);
    
    /**
       Linear position of the entire body.
    */
    Eigen::Vector3d position;
    
    /**
       Linear velocity of the entire body.
    */
    Eigen::Vector3d velocity;
    
    /**
       Attitude (orientation) of the entire body.
    */
    Eigen::Quaterniond attitude;
    
    /**
       Rotational velocity of the entire body.
    */
    Eigen::Vector3d rotational_velocity;

    void set_position(const Eigen::Vector3d &position);
    void set_attitude(const Eigen::Quaterniond &attitude);
    
    void set_velocity(const Eigen::Vector3d &velocity);
    void set_rotational_velocity(const Eigen::Vector3d &rotational_velocity);
    
    Eigen::Vector3d panel_kinematic_velocity(const Surface &surface, int panel) const;
    
    Eigen::Vector3d node_kinematic_velocity(const Surface &surface, int node) const;
};

};

#endif // __BODY_HPP__
