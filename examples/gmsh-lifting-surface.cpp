//
// Vortexje -- Gmsh wing section construction example.
//
// Copyright (C) 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <iostream>

#include <vortexje/solver.hpp>
#include <vortexje/parameters.hpp>
#include <vortexje/surface-loaders/gmsh-surface-loader.hpp>
#include <vortexje/surface-writers/vtk-surface-writer.hpp>
#include <vortexje/field-writers/vtk-field-writer.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

// Main:
int
main (int argc, char **argv)
{
    // Enable wake convection:
    Parameters::convect_wake = true;
    
    // Create lifting surface object:
    LiftingSurface wing;
    
    // Load Gmsh mesh file:
    GmshSurfaceLoader loader;
    loader.load(wing, "gmsh-lifting-surface.msh");

    // Set lifting surface metadata:
    // N.B.  The node and panel numbers must match those in the Gmsh file,
    // *except* that Gmsh counts from 1, whereas Vortexje counts from 0. 
    // The Gmsh node and panel numbers are one higher than the Vortexje numbers.
    int n_nodes_per_airfoil = 32;
    int n_airfoils = 21;
    
    wing.upper_nodes.resize(n_nodes_per_airfoil / 2 + 1, n_airfoils);
    wing.lower_nodes.resize(n_nodes_per_airfoil / 2 + 1, n_airfoils);
    for (int j = 0; j < n_airfoils; j++) {
        for (int i = 0; i < n_nodes_per_airfoil / 2 + 1; i++) { 
            wing.upper_nodes(i, j) = j * n_nodes_per_airfoil + i;
            
            if (i == 0)
                wing.lower_nodes(i, j) = j * n_nodes_per_airfoil;
            else
                wing.lower_nodes(i, j) = j * n_nodes_per_airfoil + (n_nodes_per_airfoil - i);
        }
    }
    
    wing.upper_panels.resize(n_nodes_per_airfoil / 2, n_airfoils - 1);
    wing.lower_panels.resize(n_nodes_per_airfoil / 2, n_airfoils - 1);
    for (int j = 0; j < n_airfoils - 1; j++) {
        for (int i = 0; i < n_nodes_per_airfoil / 2; i++) {
            wing.upper_panels(i, j) = j * n_nodes_per_airfoil + i;
            
            wing.lower_panels(i, j) = j * n_nodes_per_airfoil + (n_nodes_per_airfoil - 1 - i);
        }
    }
    
    // Terminate neighbor relationships across trailing edge.
    for (int i = 0; i < wing.n_spanwise_panels(); i++)
        wing.cut_panels(wing.trailing_edge_upper_panel(i), wing.trailing_edge_lower_panel(i));
    
    // Prescribe angle of attack:
    double alpha = 5.0 / 180.0 * M_PI;
    wing.rotate(Vector3d::UnitZ(), -alpha);
    
    // Create surface body:
    Body body(string("section"));
    body.add_lifting_surface(wing);
    
    // Set up solver:
    Solver solver("gmsh-lifting-surface-log");
    solver.add_body(body);
    
    Vector3d freestream_velocity(30, 0, 0);
    solver.set_freestream_velocity(freestream_velocity);
    
    double fluid_density = 1.2;
    solver.set_fluid_density(fluid_density);
    
    // Set up surface writer:
    VTKSurfaceWriter surface_writer;
    
    // Set up field writer:
    VTKFieldWriter field_writer;
    
    // Run simulation:
    double t = 0.0;
    double dt = 0.01;
    int step_number = 0;
    
    solver.initialize_wakes(dt);
    while (t < 60) {
        // Solve:
        solver.solve(dt);
        
        // Log coefficients:
        solver.log(step_number, surface_writer);
        
        // Enable below to log the velocity field:
        // field_writer.write_velocity_field(solver, "velocity-field.vtk", 0.1, 0.1, 0.1, 0.2, 0.2, 0.2);
        
        // Update wake:
        solver.update_wakes(dt);
        
        // Step time:
        t += dt;
        step_number++;
    }
    
    // Done:
    return 0;
}
