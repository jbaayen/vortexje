//
// Vortexje -- Simple wing construction example.
//
// Copyright (C) 2012 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <iostream>

#include <vortexje/solver.hpp>
#include <vortexje/wing-builder.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

int
main (int argc, char **argv)
{
    // Load meshes:
    Mesh nolift_mesh;
    
    // Create wing:
    Wing wing;
    
    wing.location = Vector3d(0, 0, 0);
    
    wing.chord_direction = Vector3d(1, 0, 0);
    wing.top_direction   = Vector3d(0, 1, 0);
    wing.span_direction  = Vector3d(0, 0, 1);
    
    WingBuilder wing_builder(wing);
    
    int trailing_edge_point_id;
    vector<int> prev_airfoil_nodes;
    
    const int n_points_per_airfoil = 32;
    const int n_airfoils = 21;
    
    const double chord = 1.0;
    const double span = 5.0;
    
    for (int i = 0; i < n_airfoils; i++) {
        vector<Vector3d, Eigen::aligned_allocator<Vector3d> > airfoil_points = wing_builder.generate_clarky_airfoil(chord, n_points_per_airfoil, trailing_edge_point_id);
        for (int j = 0; j < (int) airfoil_points.size(); j++)
            airfoil_points[j](2) += i * span / (double) (n_airfoils - 1);
            
        vector<int> airfoil_nodes = wing_builder.add_points(airfoil_points, trailing_edge_point_id);
        
        if (i > 0) {
            int trailing_edge_top_panel_id, trailing_edge_bottom_panel_id;
            
            wing_builder.connect_nodes(airfoil_nodes, prev_airfoil_nodes,
                                       trailing_edge_point_id, trailing_edge_top_panel_id, trailing_edge_bottom_panel_id,
                                       true, WingBuilder::QUADRANGLES);
                                       
            wing.trailing_edge_top_panels.push_back(trailing_edge_top_panel_id);
            wing.trailing_edge_bottom_panels.push_back(trailing_edge_bottom_panel_id);
        }
            
        prev_airfoil_nodes = airfoil_nodes;
    }
    
    wing.sort_trailing_edge();
    wing.compute_panel_neighbors();
    
    Vector3d translation(-chord / 3.0, 0.0, -span / 2.0);
    wing.translate(translation);
    
    // Prescribe angle of attack:
    double alpha = 5.0 / 180.0 * M_PI;
    wing.rotate(wing.span_direction, -alpha);
    
    // Create mesh collection:
    Collection collection(string("wing"),
                          nolift_mesh);
    collection.add_wing(&wing);
    
    // Set up solver:
    Solver solver("wing-builder-log");
    solver.add_collection(collection);
    
    Vector3d freestream_velocity(30, 0, 0);
    solver.set_freestream_velocity(freestream_velocity);
    
    double fluid_density = 1.2;
    solver.set_fluid_density(fluid_density);
    
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
