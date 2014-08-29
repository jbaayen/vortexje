//
// Vortexje -- Dummy boundary layer class.
//
// Copyright (C) 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __DUMMY_BOUNDARY_LAYER_HPP__
#define __DUMMY_BOUNDARY_LAYER_HPP__

#include <vortexje/boundary-layer.hpp>

namespace Vortexje
{

/**
   Dummy boundary layer class.
   
   @brief Dummy boundary layer.
*/
class DummyBoundaryLayer : public BoundaryLayer
{
public:        
    bool recalculate(const Eigen::MatrixXd &surface_velocities);
    
    double blowing_velocity(const std::shared_ptr<Surface> &surface, int panel) const;
    
    Eigen::Vector3d friction(const std::shared_ptr<Surface> &surface, int panel) const;
};

};

#endif // __DUMMY_BOUNDARY_LAYER_HPP__
