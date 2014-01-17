//
// Vortexje -- Simple sphere example.
//
// Copyright (C) 2013 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <iostream>
#include <fstream>

#include <vortexje/solver.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

#define TEST_TOLERANCE 5e-2

int
main (int argc, char **argv)
{
    // Load sphere mesh:
    Mesh sphere(string("sphere.msh"));
   
    // Create mesh collection:
    Collection collection(string("test-sphere"), sphere);
    
    // Set up solver:
    Solver solver("test-sphere-log");
    solver.add_collection(collection);
    
    Vector3d freestream_velocity(30.0, 0, 0);
    solver.set_freestream_velocity(freestream_velocity);
    
    double fluid_density = 1.2;
    solver.set_fluid_density(fluid_density);
    
    // Run simulation:
    solver.update_coefficients(0);
    
    // Check pressure coefficients:
    for (int i = 0; i < sphere.n_panels(); i++) {
        Vector3d x = sphere.panel_collocation_point(i, false);
        
        // Angle between flow and point on sphere:
        double theta = acos(x.dot(freestream_velocity) / (x.norm() * freestream_velocity.norm()));
        
        // Analytical solution for pressure coefficient.  See J. Katz and A. Plotkin, Low-Speed Aerodynamics, Cambridge University Press, 2001.
        double C_p_ref = 1.0 - 9.0 / 4.0 * pow(sin(theta), 2);
        
        // Computed pressure coefficient:
        double C_p = solver.pressure_coefficient(sphere, i);
        
        // Compare:
        if (fabs(C_p - C_p_ref) > TEST_TOLERANCE) {
            cerr << " *** TEST FAILED *** " << endl;
            cerr << " theta = " << theta << " rad" << endl;
            cerr << " C_p(ref) = " << C_p_ref << endl;
            cerr << " C_p = " << C_p << endl;
            cerr << " ******************* " << endl;
            
            exit(1);
        }
    }
    
    // Done:
    return 0;
}
