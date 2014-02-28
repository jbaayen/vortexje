//
// Vortexje -- Default parameters
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
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

bool   Parameters::convect_wake                       = false;

double Parameters::static_wake_length                 = 100.0;

double Parameters::inversion_tolerance                = 1e-7;

double Parameters::collocation_point_delta            = 1e-7;

bool   Parameters::marcov_surface_velocity            = false;

int    Parameters::max_boundary_layer_iterations      = 100;

double Parameters::boundary_layer_iteration_tolerance = 1e-3;
