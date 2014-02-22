//
// Vortexje -- Boundary layer base class.
//
// Copyright (C) 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __BOUNDARY_LAYER_HPP__
#define __BOUNDARY_LAYER_HPP__

#include <Eigen/Core>

namespace Vortexje
{

/**
   Per-surface boundary layer base class.
   
   @brief Boundary layer base class.
*/
class BoundaryLayer
{
public:        
    virtual void recalculate(const Eigen::MatrixXd &surface_velocities) = 0;
    
    virtual double blowing_velocity(int panel) const = 0;
    
    virtual Eigen::Vector3d friction(int panel) const = 0;
};

};

#endif // __BOUNDARY_LAYER_HPP__
