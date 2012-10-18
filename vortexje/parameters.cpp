//
// Vortexje -- Default parameters
//
// Copyright (C) 2012 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <cmath>

#include <vortexje/parameters.hpp>

using namespace Vortexje;

// Default values:
int    Parameters::linear_solver_max_iterations       = 20000;

double Parameters::linear_solver_tolerance            = 1e-10;

bool   Parameters::unsteady_bernoulli                 = true;

bool   Parameters::convect_wake                       = true;

double Parameters::static_wake_length                 = 100.0;

double Parameters::min_pressure_coefficient           = -12;

double Parameters::fluid_kinematic_viscosity          = 15.68e-6;

bool   Parameters::use_ramasamy_leishman_vortex_sheet = true;

double Parameters::initial_vortex_core_radius         = 0.07;

double Parameters::min_vortex_core_radius             = 0.07;

double Parameters::lambs_constant                     = 1.25643;

double Parameters::a_prime                            = 6.5e-5;

bool   Parameters::use_markov_surface_velocity        = false;

double Parameters::interpolation_layer_thickness      = 0.075;

double Parameters::interpolation_layer_notch_angle    = 2 * M_PI / 3;

double Parameters::inversion_tolerance                = 1e-6;

double Parameters::sharp_edge_threshold               = 0.1;

double Parameters::collocation_point_delta_factor     = 1e-2;
