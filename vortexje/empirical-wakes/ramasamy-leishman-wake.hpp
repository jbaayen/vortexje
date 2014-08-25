//
// Vortexje -- Ramasamy-Leishman wake model.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __RAMASAMY_LEISHMAN_WAKE_HPP__
#define __RAMASAMY_LEISHMAN_WAKE_HPP__

#include <vortexje/wake.hpp>

namespace Vortexje
{

/**
   Representation of a wake, i.e., a vortex sheet.
   
   @brief Wake representation.
   
   @note See M. Ramasamy and J. G. Leishman, Reynolds Number Based Blade Tip Vortex Model, University of Maryland, 2005.
*/
class RamasamyLeishmanWake : public Wake
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    RamasamyLeishmanWake(std::shared_ptr<LiftingSurface> lifting_surface);
    
    void add_layer();
    
    void update_properties(double dt);
    
    Eigen::Vector3d vortex_ring_unit_velocity(const Eigen::Vector3d &x, int this_panel) const;

    /**
       Radii of the vortex filaments forming the vortex rings.
    */ 
    std::vector<std::vector<double> > vortex_core_radii;
    
    /**
       Ramasamy-Leishman wake model parameters.
       
       @brief Wake model parameters.
    */
    class Parameters {
    public:
        /**
           Kinematic viscosity of the fluid.
        */
        static double fluid_kinematic_viscosity;
        
        /**
           Initial vortex filament radius, when using the Ramasamy-Leishman vortex sheet model.
           
           @note See M. Ramasamy and J. G. Leishman, Reynolds Number Based Blade Tip Vortex Model, University of Maryland, 2005.
        */
        static double initial_vortex_core_radius;
        
        /**
           Minimum vortex filament radius, when using the Ramasamy-Leishman vortex sheet model.
           
           @note See M. Ramasamy and J. G. Leishman, Reynolds Number Based Blade Tip Vortex Model, University of Maryland, 2005.
        */
        static double min_vortex_core_radius;
    
        /**
           Lamb's constant from the Ramasamy-Leishman vortex sheet model.
           
           @note See M. Ramasamy and J. G. Leishman, Reynolds Number Based Blade Tip Vortex Model, University of Maryland, 2005.
        */
        static double lambs_constant;
        
        /**
           a' constant from the Ramasamy-Leishman vortex sheet model.
           
           @note See M. Ramasamy and J. G. Leishman, Reynolds Number Based Blade Tip Vortex Model, University of Maryland, 2005.
        */
        static double a_prime;
    };
  
private:  
    /**
       Initial lengths of the vortex filaments forming the vortex rings.
    */
    std::vector<std::vector<double> > base_edge_lengths;
    
    void update_vortex_ring_radii(int panel, double dt);
};

};

#endif // __RAMASAMY_LEISHMAN_WAKE_HPP__
