//
// Vortexje -- Numerical simulation parameters.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __PARAMETERS_HPP__
#define __PARAMETERS_HPP__

namespace Vortexje
{

/**
   Parameter settings.
   
   The static parameters of this class control the behavior of the solver.
   
   @brief Parameter settings.
*/
class Parameters
{
public:
    /**
       Maximum number of iterations of BiCGSTAB linear solver.
    */
    static int    linear_solver_max_iterations;
    
    /**
       Tolerance of BiCGSTAB linear solver.
    */
    static double linear_solver_tolerance;
    
    /**
       Whether or not to apply the unsteady Bernoulli equation.
    */
    static bool   unsteady_bernoulli;
    
    /**
       Whether or not to enable wake relaxation/convection.
    */
    static bool   convect_wake;
    
    /**
       Length of the wake, in case of no wake convection.
    */
    static double static_wake_length;
    
    /**
       Kinematic viscosity of the fluid.
    */
    static double fluid_kinematic_viscosity;
    
    /**
       Whether or not to use the Ramasamy-Leishman vortex sheet model.
       
       @note See M. Ramasamy and J. G. Leishman, Reynolds Number Based Blade Tip Vortex Model, University of Maryland, 2005.
    */
    static bool   use_ramasamy_leishman_vortex_sheet;
    
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
    
    /**
       Quantities below this threshold will be treated as nil.
    */
    static double inversion_tolerance;
    
    /**
       Distance to below-surface collocation points (along normal).
    */
    static double collocation_point_delta;
    
    /**
       Use N. Marcov's formula for computing the surface velocities.
       
       @note See L. Dragoş, Mathematical Methods in Aerodynamics, Springer, 2003.
    */
    static bool marcov_surface_velocity;
};

};

#endif // __PARAMETERS_HPP__
