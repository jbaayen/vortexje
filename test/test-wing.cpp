//
// Vortexje -- Simple wing test case.
//
// Copyright (C) 2012 Baayen & Heinz GmbH.
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
    // Load meshes:
    Mesh nolift_mesh;
    
    // Create wing:
    Vector3d location(0, 0, 0);
    Vector3d chord_direction(1, 0, 0);
    Vector3d top_direction(0, 1, 0);
    Vector3d span_direction(0, 0, 1);
    
    Mesh wing_mesh(string("naca0012.msh"));
    
    Wing wing(wing_mesh, location, chord_direction, top_direction, span_direction);
    
    // Prescribe angle of attack:
    double alpha = 5.0 / 180.0 * M_PI;
    wing.rotate(span_direction, -alpha);
    
    // Create mesh collection:
    Collection collection(string("test-wing"),
                          nolift_mesh);
    collection.add_wing(&wing);
    
    // Set up solver:
    Solver solver("test-wing-log");
    solver.add_collection(collection);
    solver.set_freestream_velocity(Vector3d(30, 0, 0));
    solver.set_fluid_density(1.2);
    
    // Run simulation:
    double t = 0.0;
    double dt = 0.01;
    int step_number = 0;
    
    solver.initialize_wakes(dt);
    while (t < 60) {
        solver.update_coefficients(dt);
        solver.log_coefficients(step_number);
        solver.update_wakes(dt);
        
        t += dt;
        step_number++;
    }
    
    // Done:
    return 0;
}
