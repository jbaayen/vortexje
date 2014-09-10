//
// Vortexje -- NACA0025 airfoil section with rounded tips, following the setup from
//             K. Bullivant, Tests of the NACA 0025 and 0035 Airfoils in the Full-Scale Wind Tunnel, NACA Report No. 708, 1941.
//
// Copyright (C) 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <cmath>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <sys/stat.h>

#include <vortexje/solver.hpp>
#include <vortexje/lifting-surface-builder.hpp>
#include <vortexje/shape-generators/airfoils/naca4-airfoil-generator.hpp>
#include <vortexje/surface-writers/vtk-surface-writer.hpp>
#include <vortexje/surface-writers/gmsh-surface-writer.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

static const double pi = 3.141592653589793238462643383279502884;

#define ROUNDED_TIPS

#define DELTA_CONVERGENCE 1e-1
#define DELTA_T           5e-2

void
simulate (double alpha_deg)
{ 
    // Set up parameters for unsteady simulation:
    Parameters::unsteady_bernoulli = true;
    Parameters::convect_wake       = true;
    
    // Create wing:
    shared_ptr<LiftingSurface> wing(new LiftingSurface());
    
    LiftingSurfaceBuilder surface_builder(*wing);

    const double chord = 1.8288; // 6 feet
    const double span = 10.9728; // 36 feet
    
    const int n_points_per_airfoil = 32;
    const int n_airfoils = 21;
#ifdef ROUNDED_TIPS
    const int n_airfoils_per_tip = 3;
#endif
    
    vector<Vector3d, Eigen::aligned_allocator<Vector3d> > airfoil_points;
    vector<int> airfoil_nodes;
    vector<int> airfoil_panels;
    
    int trailing_edge_point_id;
    vector<int> first_airfoil_nodes;
    vector<int> last_airfoil_nodes;
    vector<int> prev_airfoil_nodes; 
    
    vector<vector<int> > node_strips;
    vector<vector<int> > panel_strips;
    
    // Rectangular airfoil section.
    for (int i = 0; i < n_airfoils; i++) {
        airfoil_points = NACA4AirfoilGenerator::generate(0, 0, 0.25, true, chord, n_points_per_airfoil, trailing_edge_point_id);

        for (int j = 0; j < (int) airfoil_points.size(); j++)
            airfoil_points[j](2) += i * span / (double) (n_airfoils - 1);
             
        airfoil_nodes = surface_builder.create_nodes_for_points(airfoil_points);
        node_strips.push_back(airfoil_nodes);
        
        if (i > 0)
            airfoil_panels = surface_builder.create_panels_between_shapes(airfoil_nodes, prev_airfoil_nodes, trailing_edge_point_id);
        
        if (i > 0)
            panel_strips.push_back(airfoil_panels);
            
        prev_airfoil_nodes = airfoil_nodes;
        
        if (i == 0)
            first_airfoil_nodes = airfoil_nodes;
    }
    
    last_airfoil_nodes = prev_airfoil_nodes;

#ifdef ROUNDED_TIPS
    // Negative half of a surface of revolution.
    for (int i = 0; i < n_airfoils_per_tip; i++) {
        airfoil_points = NACA4AirfoilGenerator::generate(0, 0, 0.25, true, chord, n_points_per_airfoil, trailing_edge_point_id);
        for (int j = 0; j < (int) airfoil_points.size(); j++) {
            airfoil_points[j] = AngleAxis<double>((n_airfoils_per_tip - i) * pi / 2.0 / (double) (n_airfoils_per_tip), Vector3d::UnitX()) * airfoil_points[j];
            
            if (airfoil_points[j](2) > 0)
                airfoil_points[j](2) *= -1;
        }
        
        airfoil_nodes = surface_builder.create_nodes_for_points(airfoil_points);
        
        airfoil_nodes[0]                      = first_airfoil_nodes[0];
        airfoil_nodes[trailing_edge_point_id] = first_airfoil_nodes[trailing_edge_point_id];
        
        if (i > 0)
            airfoil_panels = surface_builder.create_panels_between_shapes(airfoil_nodes, prev_airfoil_nodes, trailing_edge_point_id);
        
        prev_airfoil_nodes = airfoil_nodes;
    }
    
    surface_builder.create_panels_between_shapes(prev_airfoil_nodes, first_airfoil_nodes, trailing_edge_point_id);
    
    // Positive half of a surface of revolution.
    prev_airfoil_nodes = last_airfoil_nodes;
    
    for (int i = 0; i < n_airfoils_per_tip; i++) {
        airfoil_points = NACA4AirfoilGenerator::generate(0, 0, 0.25, true, chord, n_points_per_airfoil, trailing_edge_point_id);
        for (int j = 0; j < (int) airfoil_points.size(); j++) {
            airfoil_points[j] = AngleAxis<double>((i + 1) * pi / 2.0 / (double) (n_airfoils_per_tip), Vector3d::UnitX()) * airfoil_points[j];
            
            if (airfoil_points[j](2) < 0)
                airfoil_points[j](2) *= -1;
                
            airfoil_points[j](2) += span;
        }
        
        airfoil_nodes = surface_builder.create_nodes_for_points(airfoil_points);
        
        airfoil_nodes[0]                      = last_airfoil_nodes[0];
        airfoil_nodes[trailing_edge_point_id] = last_airfoil_nodes[trailing_edge_point_id];
        
        airfoil_panels = surface_builder.create_panels_between_shapes(airfoil_nodes, prev_airfoil_nodes, trailing_edge_point_id);
        
        prev_airfoil_nodes = airfoil_nodes;
    }
#endif

    surface_builder.finish(node_strips, panel_strips, trailing_edge_point_id);
    
    // Create body:
    shared_ptr<Body> body(new Body(string("section")));
    body->add_lifting_surface(wing);
    
    // Set up orientation: 
    double alpha = pi * alpha_deg / 180.0;
    
    Quaterniond attitude = AngleAxis<double>(-alpha, Vector3d::UnitZ()) * Quaterniond(1, 0, 0, 0);
    body->set_attitude(attitude);
    
    // Set up solver:
    Solver solver("naca0025-rounded-tips-log");
    solver.add_body(body);
    
    Vector3d freestream_velocity(24.6888, 0, 0); // 81 feet/second.
    solver.set_freestream_velocity(freestream_velocity);
    
    double fluid_density = 1.2;
    solver.set_fluid_density(fluid_density);
    
    // Run simulation:
    double t = 0.0;
    double dt = DELTA_T;
    int step_number = 0;
    
    solver.initialize_wakes(dt);
    
    // Iterate until converence:
    Vector3d F_prev(0, 0, 0);
    while (true) {
        // Solve:
        solver.solve(dt);

        // Compute force:
        Vector3d F = solver.force(body);
        
        // Check convergence:
        double delta = (F - F_prev).norm();
        cout << "Force delta = " << delta << " N" << endl;
        if (delta < DELTA_CONVERGENCE)
            break;
        F_prev = F;
        
        // Update wakes:
        solver.update_wakes(dt);
        
        // Step time:
        t += dt;
        step_number++;
    }
    
    // Log output:
    VTKSurfaceWriter surface_writer;
    solver.log(0, surface_writer);
    
    Vector3d F = solver.force(body);
    double S = chord * span;
    
    cout << "C_L = " << F(1) / (0.5 * fluid_density * freestream_velocity.squaredNorm() * S) << endl;
}

int
main (int argc, char **argv)
{
    simulate(6.9);
    
    return 0;
}
