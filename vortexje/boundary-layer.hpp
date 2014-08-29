//
// Vortexje -- Boundary layer base class.
//
// Copyright (C) 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __BOUNDARY_LAYER_HPP__
#define __BOUNDARY_LAYER_HPP__

#include <memory>

#include <Eigen/Core>

#include <vortexje/surface.hpp>

namespace Vortexje
{

/**
   Per-surface boundary layer base class.
   
   @brief Boundary layer base class.
*/
class BoundaryLayer
{
public:
    /**
       Destructor.
    */
    virtual ~BoundaryLayer() {};
    
    /**
       Normally, this function solves the relevant boundary layer equations.  Here, it does nothing.
   
       @param[in]   surface_velocities   (n x 3)-matrix of surface velocities.
       
       @returns true on success.
    */    
    virtual bool recalculate(const Eigen::MatrixXd &surface_velocities) = 0;
    
    /**
       Returns the blowing velocity for the given panel.
   
       @param[in]   surface   Reference surface.
       @param[in]   panel     Reference panel.
   
       @returns Blowing velocity for the given panel.
    */
    virtual double blowing_velocity(const std::shared_ptr<Surface> &surface, int panel) const = 0;
    
    /**
       Returns the friction force acting on the given panel.
   
       @param[in]   surface   Reference surface.
       @param[in]   panel     Reference panel.
   
       @returns Friction force acting on the given panel.
    */
    virtual Eigen::Vector3d friction(const std::shared_ptr<Surface> &surface, int panel) const = 0;
};

};

#endif // __BOUNDARY_LAYER_HPP__
