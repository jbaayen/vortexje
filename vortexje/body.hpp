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
       List of non-lifting surfaces.
    */
    std::vector<Surface*> non_lifting_surfaces;
    
    /**
       List of surfaces.
    */
    std::vector<LiftingSurface*> lifting_surfaces;
    
    /**
       List of wakes associated with the lifting surfaces, in the same order.
    */
    std::vector<Wake*> wakes;
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    Body(std::string id);
         
    ~Body();

    virtual void add_non_lifting_surface(Surface *surface);
    
    virtual void add_lifting_surface(LiftingSurface *lifting_surface);
    
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
    
protected:
    /**
       List of all non-wake surfaces.  For internal use.
    */
    std::vector<Surface*> non_wake_surfaces;
};

};

#endif // __BODY_HPP__
