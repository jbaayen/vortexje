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
       Whether or not to enable the convection of wake nodes.
    */
    static bool   convect_wake;
    
    /**
       Whether to emit new wake panels into the direction of the trailing edge bisector, rather than following the
       apparent velocity.
    */
    static bool   wake_emission_follow_bisector;
    
    /**
       Multiplied with the trailing edge velocity to obtain the distance by which the first wake vortex is placed away
       from the trailing edge.
       
       @note See J. Katz and A, Plotkin, Low-Speed Aerodynamics, 2nd Edition, Cambridge University Press, 2001. 
    */
    static double wake_emission_distance_factor;
    
    /**
       Wake vortex filament core radius (Rankine model).  Within the core, the velocity decreases linearly to zero.
       
       For wake-wake interaction problems, set this to an suitably small number greater than zero.  Care must be taken, however, that
       the radius is less than the distance between two subsequent layers of wake panels.  I.e., it must hold that
       
       \f[
         \mbox{wake\_vortex\_core\_radius} < \mbox{wake\_emission\_distance\_factor} \cdot v_{\infty} \cdot dt 
       \f]
    */
    static double wake_vortex_core_radius;
    
    /**
       Length of the wake, in case of no wake convection.
    */
    static double static_wake_length;
    
    /**
       Quantities below this threshold will be treated as nil.
    */
    static double inversion_tolerance;
    
    /**
       Distance to below-surface collocation points (along normal).
    */
    static double collocation_point_delta;
    
    /**
       Thickness of the velocity interpolation layer around the body and its boundary layer.
       
       Inside the interpolation layer, velocities will be interpolated linearly between surface (or boundary layer) velocities, and the 
       velocities a layer thickness away from the body and its boundary layer.
       
       For wake-body interaction problems, set this to an suitably small number greater than zero.  Care must be taken, however, that
       the outer edge of the first row of wake panels lies outside of the interpolation layer.  I.e., it must hold that
       
       \f[
         \mbox{interpolation\_layer\_thickness} < \mbox{wake\_emission\_distance\_factor} \cdot v_{\infty} \cdot dt 
       \f]
    */
    static double interpolation_layer_thickness;
    
    /**
       Maximum number of potential and boundary layer computation iterations.
    */
    static int    max_boundary_layer_iterations;
    
    /**
       Boundary layer iteration tolerance.
    */
    static double boundary_layer_iteration_tolerance;
};

};

#endif // __PARAMETERS_HPP__
