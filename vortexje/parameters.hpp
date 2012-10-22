//
// Vortexje -- Numerical simulation parameters.
//
// Copyright (C) 2012 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __PARAMETERS_HPP__
#define __PARAMETERS_HPP__

namespace Vortexje
{

class Parameters
{
public:
    // BiCGSTAB linear solver parameters.
    static int    linear_solver_max_iterations;
    
    static double linear_solver_tolerance;
    
    // Unsteady Bernoulli equation.
    static bool   unsteady_bernoulli;
    
    // Wake convection.
    static bool   convect_wake;
    
    static double static_wake_length;
    
    // Separation detection.
    static double min_pressure_coefficient;
    
    // Kinematic viscosity of fluid.
    static double fluid_kinematic_viscosity;
    
    // Whether or not to use the Ramasamy-Leishman vortex sheet model, see
    //    M. Ramasamy and J. G. Leishman, Reynolds Number Based Blade Tip Vortex Model, University of Maryland, 2005.
    static bool   use_ramasamy_leishman_vortex_sheet;
    
    // Ramasamy-Leishman vortex sheet model parameters:
    static double initial_vortex_core_radius;
    
    static double min_vortex_core_radius;

    static double lambs_constant;
    
    static double a_prime;
    
    // Close-to-body interpolation layer parameters, see
    //    K. Dixon, C. S. Ferreira, C. Hofemann, G. van Brussel, G. van Kuik,
    //    A 3D Unsteady Panel Method for Vertical Axis Wind Turbines, DUWIND, 2008.
    static double interpolation_layer_thickness;
    
    static double interpolation_layer_notch_angle;
    
    // Inversion threshold.  Values below this threshold will not be inverted.
    static double inversion_tolerance;
    
    // If the inner product of normals of neighbouring panels is more than sharp_edge_tolerance in absolute value,
    // then the edge shared by the two panels is considered "sharp" and treated accordingly.
    // (e.g., for trailing edge detection.)
    static double sharp_edge_threshold;
    
    // Factor applied to the smallest panel edge, to obtain the normal distance for below-surface collocation points:
    //    collocation_point_delta = collocation_point_delta_factor * min_edge.
    static double collocation_point_delta_factor;
};

};

#endif // __PARAMETERS_HPP__
