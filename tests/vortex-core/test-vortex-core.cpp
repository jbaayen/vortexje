//
// Vortexje -- Test Rankine wake vortex model.
//
// Copyright (C) 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <iostream>
#include <fstream>

#include <vortexje/solver.hpp>
#include <vortexje/lifting-surface-builder.hpp>
#include <vortexje/shape-generators/airfoils/naca4-airfoil-generator.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

static const double pi = 3.141592653589793238462643383279502884;

#define TEST_TOLERANCE 1e-3

int
main (int argc, char **argv)
{
    // Set parameters:
    Parameters::wake_vortex_core_radius = 1e1;
    Parameters::convect_wake            = false;
    Parameters::static_wake_length      = 1e5;
    
    // Create wing:
    shared_ptr<LiftingSurface> wing(new LiftingSurface("main"));
    
    LiftingSurfaceBuilder surface_builder(*wing);

    const double chord = 1.0;
    const double span = 1e4;
    
    const int n_points_per_airfoil = 10;
    const int n_airfoils = 2;
    
    int trailing_edge_point_id;
    vector<int> prev_airfoil_nodes;
    
    vector<vector<int> > node_strips;
    vector<vector<int> > panel_strips;
    
    for (int i = 0; i < n_airfoils; i++) {
        vector_aligned<Vector3d> airfoil_points =
            NACA4AirfoilGenerator::generate(0, 0, 0.12, true, chord, n_points_per_airfoil, trailing_edge_point_id);
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

    // Rotate by angle of attack:
    double alpha = pi * 5.0 / 180.0;
    wing->rotate(Vector3d::UnitZ(), -alpha);
    
    // Create wake:
    shared_ptr<Wake> wake(new Wake(wing));
    
    // Create surface body:
    shared_ptr<Body> body(new Body(string("wing-section")));
    body->add_lifting_surface(wing, wake);
    
    // Set up solver:
    Solver solver("test-vortex-core-log");
    solver.add_body(body);
    
    Vector3d freestream_velocity(30.0, 0, 0);
    solver.set_freestream_velocity(freestream_velocity);
    
    double fluid_density = 1.2;
    solver.set_fluid_density(fluid_density);
    
    // Run simulation:
    solver.initialize_wakes();
    solver.solve();
    
    // Check wake-induced velocities:
    Vector3d te_middle = 0.5 * wing->nodes[wing->trailing_edge_node(0)] + 0.5 * wing->nodes[wing->trailing_edge_node(1)];
    double r;
    
    r = 2.0 * Parameters::wake_vortex_core_radius;
    Vector3d far_point(te_middle + (Parameters::static_wake_length + r) * Vector3d(1, 0, 0));
    
    Vector3d reference_far_point_velocity(0, -wake->doublet_coefficients[0] / (2 * pi * r), 0);
    
    Vector3d far_point_velocity = solver.velocity(far_point) - freestream_velocity;
    if ((far_point_velocity - reference_far_point_velocity).norm() > TEST_TOLERANCE) {
        cerr << " *** FAR POINT TEST FAILED *** " << endl;
        cerr << " |V_ref| = " << reference_far_point_velocity.norm() << endl;
        cerr << " |V| = " << far_point_velocity.norm() << endl;
        cerr << " ******************* " << endl;
        
        exit(1);
    }
    
    r = 0.5 * Parameters::wake_vortex_core_radius;
    Vector3d close_point(te_middle + (Parameters::static_wake_length + r) * Vector3d(1, 0, 0));
 
    Vector3d reference_close_point_velocity(0, -wake->doublet_coefficients[0] * r / (2 * pi * pow(Parameters::wake_vortex_core_radius, 2)), 0);
    //Vector3d reference_close_point_velocity(0, -wake->doublet_coefficients[0] / (2 * pi * r), 0);
    
    Vector3d close_point_velocity = solver.velocity(close_point) - freestream_velocity;
    if ((close_point_velocity - reference_close_point_velocity).norm() > TEST_TOLERANCE) {
        cerr << " *** CLOSE POINT TEST FAILED *** " << endl;
        cerr << " |V_ref| = " << reference_close_point_velocity.norm() << endl;
        cerr << " |V| = " << close_point_velocity.norm() << endl;
        cerr << " ******************* " << endl;
        
        exit(1);
    }
    
    r = 0.0;
    Vector3d locus_point(te_middle + (Parameters::static_wake_length + r) * Vector3d(1, 0, 0));
 
    Vector3d reference_locus_point_velocity(0, -wake->doublet_coefficients[0] * r / (2 * pi * pow(Parameters::wake_vortex_core_radius, 2)), 0);
    //Vector3d reference_locus_point_velocity(0, -wake->doublet_coefficients[0] / (2 * pi * r), 0);
    
    Vector3d locus_point_velocity = solver.velocity(locus_point) - freestream_velocity;
    if ((locus_point_velocity - reference_locus_point_velocity).norm() > TEST_TOLERANCE) {
        cerr << " *** LOCUS POINT TEST FAILED *** " << endl;
        cerr << " |V_ref| = " << reference_locus_point_velocity.norm() << endl;
        cerr << " |V| = " << locus_point_velocity.norm() << endl;
        cerr << " ******************* " << endl;
        
        exit(1);
    }
    
    // Done:
    return 0;
}
