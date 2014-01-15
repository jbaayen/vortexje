//
// Vortexje -- Simple sphere example.
//
// Copyright (C) 2013 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <iostream>

#include <vortexje/solver.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

int
main (int argc, char **argv)
{
    // For this rather fine mesh, we need to adjust the distance between collocation points and mesh surface.
    Parameters::collocation_point_delta_factor = 5e-2;
    
    // Load sphere mesh:
    Mesh sphere(string("sphere.msh"));
   
    // Create mesh collection:
    Collection collection(string("sphere"), sphere);
    
    // Set up solver:
    Solver solver("sphere-log");
    solver.add_collection(collection);
    
    Vector3d freestream_velocity(1.5, 0, 0);
    solver.set_freestream_velocity(freestream_velocity);
    
    double fluid_density = 1.2;
    solver.set_fluid_density(fluid_density);
    
    // Run simulation:
    solver.update_coefficients(0);
    solver.log_coefficients(0);
    
    // Done:
    return 0;
}
