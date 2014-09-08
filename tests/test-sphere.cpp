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
#include <vortexje/surface-loaders/gmsh-surface-loader.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

#define TEST_TOLERANCE 5e-2

int
main (int argc, char **argv)
{
    // Load sphere surface:
    GmshSurfaceLoader surface_loader;
    
    shared_ptr<Surface> sphere(new Surface());
    surface_loader.load(sphere, string("sphere.msh"));
   
    // Create surface body:
    shared_ptr<Body> body(new Body(string("test-sphere")));
    body->add_non_lifting_surface(sphere);
    
    // Set up solver:
    Solver solver("test-sphere-log");
    solver.add_body(body);
    
    Vector3d freestream_velocity(30.0, 0, 0);
    solver.set_freestream_velocity(freestream_velocity);
    
    double fluid_density = 1.2;
    solver.set_fluid_density(fluid_density);
    
    // Run simulation:
    solver.solve();
    
    // Check pressure coefficients and surface potential values:
    for (int i = 0; i < sphere->n_panels(); i++) {
        Vector3d x = sphere->panel_collocation_point(i, false);
        
        // Angle between flow and point on sphere:
        double theta = acos(x.dot(freestream_velocity) / (x.norm() * freestream_velocity.norm()));
        
        // Analytical solution for pressure coefficient and surface potential.
        // See J. Katz and A. Plotkin, Low-Speed Aerodynamics, Cambridge University Press, 2001.
        double C_p_ref = 1.0 - 9.0 / 4.0 * pow(sin(theta), 2);
        
        double R = x.norm();
        double phi_ref = freestream_velocity.norm() * cos(theta) * (R + pow(R, 3) / (2 * pow(R, 2)));
        
        // Computed pressure coefficient:
        double C_p = solver.pressure_coefficient(sphere, i);
        
        double phi = solver.surface_velocity_potential(sphere, i);
        
        // Compare:
        if (fabs(C_p - C_p_ref) > TEST_TOLERANCE) {
            cerr << " *** TEST FAILED *** " << endl;
            cerr << " theta = " << theta << " rad" << endl;
            cerr << " C_p(ref) = " << C_p_ref << endl;
            cerr << " C_p = " << C_p << endl;
            cerr << " ******************* " << endl;
            
            exit(1);
        }
        
        if (fabs(phi - phi_ref) / (R * freestream_velocity.norm()) > TEST_TOLERANCE) {
            cerr << " *** TEST FAILED *** " << endl;
            cerr << " theta = " << theta << " rad" << endl;
            cerr << " phi(ref) = " << phi_ref << endl;
            cerr << " phi = " << phi << endl;
            cerr << " ******************* " << endl;
            
            exit(1);
        }
    }
    
    // Done:
    return 0;
}
