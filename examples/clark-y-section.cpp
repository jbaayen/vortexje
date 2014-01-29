//
// Vortexje -- Clark-Y wing section construction example.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <iostream>

#include <vortexje/solver.hpp>
#include <vortexje/parameters.hpp>
#include <vortexje/lifting-surface-builder.hpp>
#include <vortexje/shape-generators/airfoils/clark-y-airfoil-generator.hpp>
#include <vortexje/surface-writers/vtk-surface-writer.hpp>
#include <vortexje/field-writers/vtk-field-writer.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

int
main (int argc, char **argv)
{
    // Enable wake convection:
    Parameters::convect_wake = true;
    
    // Create lifting surface object:
    LiftingSurface wing;

    // Construct wing section:
    LiftingSurfaceBuilder surface_builder(wing);

    const int n_points_per_airfoil = 32;
    const int n_airfoils = 21;
    
    const double chord = 1.0;
    const double span = 5.0;
    
    int trailing_edge_point_id;
    vector<int> prev_airfoil_nodes;
    
    vector<vector<int> > node_strips;
    vector<vector<int> > panel_strips;
    
    for (int i = 0; i < n_airfoils; i++) {
        vector<Vector3d, Eigen::aligned_allocator<Vector3d> > airfoil_points =
            ClarkYAirfoilGenerator::generate(chord, n_points_per_airfoil, trailing_edge_point_id);
        for (int j = 0; j < (int) airfoil_points.size(); j++)
            airfoil_points[j](2) += i * span / (double) (n_airfoils - 1);
             
        vector<int> airfoil_nodes = surface_builder.create_nodes_for_points(airfoil_points);
        node_strips.push_back(airfoil_nodes);
        
        if (i > 0) {
            vector<int> airfoil_panels = surface_builder.create_panels_between_shapes(airfoil_nodes, prev_airfoil_nodes, trailing_edge_point_id);
            panel_strips.push_back(airfoil_panels);
        }
            
        prev_airfoil_nodes = airfoil_nodes;
    }

    surface_builder.finish(node_strips, panel_strips, trailing_edge_point_id);
    
    // Translate into the canonical coordinate system:
    Vector3d translation(-chord / 3.0, 0.0, -span / 2.0);
    wing.translate(translation);
    
    // Prescribe angle of attack:
    double alpha = 5.0 / 180.0 * M_PI;
    wing.rotate(Vector3d::UnitZ(), -alpha);
    
    // Create surface body:
    Body body(string("wing"));
    body.add_lifting_surface(&wing);
    
    // Set up solver:
    Solver solver("clark-y-section-log");
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
        solver.update_coefficients(dt);
        
        // Log coefficients:
        solver.log_coefficients(step_number, surface_writer);
        
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
