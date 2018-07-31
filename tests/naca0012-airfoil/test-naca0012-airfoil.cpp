//
// Vortexje -- Test NACA0012 lift- and drag coefficients against a reference.
//
// Copyright (C) 2013 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <cmath>
#include <iostream>

#include <vortexje/solver.hpp>
#include <vortexje/lifting-surface-builder.hpp>
#include <vortexje/shape-generators/airfoils/naca4-airfoil-generator.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

static const double pi = 3.141592653589793238462643383279502884;

#define TEST_TOLERANCE 2e-2

Vector2d
run_test(double alpha)
{
    // Set up parameters for simplest possible simulation:
    Parameters::unsteady_bernoulli = false;
    Parameters::convect_wake       = false;
    
    // Create wing:
    shared_ptr<LiftingSurface> wing(new LiftingSurface("main"));
    
    LiftingSurfaceBuilder surface_builder(*wing);

    const double chord = 0.75;
    const double span = 4.5;
    
    const int n_points_per_airfoil = 32;
    const int n_airfoils = 21;
    
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
    wing->rotate(Vector3d::UnitZ(), -alpha);
    
    // Create surface body:
    shared_ptr<Body> body(new Body(string("wing-section")));
    body->add_lifting_surface(wing);
    
    // Set up solver:
    Solver solver("test-naca0012-log");
    solver.add_body(body);
    
    Vector3d freestream_velocity(30, 0, 0);
    solver.set_freestream_velocity(freestream_velocity);
    
    double fluid_density = 1.2;
    solver.set_fluid_density(fluid_density);
    
    // Compute:
    solver.initialize_wakes();
    solver.solve();
    
    // Output lift and drag coefficients:
    Vector3d F_a = solver.force(body);
    
    double q = 0.5 * fluid_density * freestream_velocity.squaredNorm();
    
    double S = chord * span;
    
    double C_L = F_a(1) / (q * S);
    double C_D = F_a(0) / (q * S);
    
    return Vector2d(C_L, C_D);
}

int
main (int argc, char **argv)
{
    // Load reference values.
    vector_aligned<Vector3d> reference_results;
    
    ifstream f;
    f.open("naca0012-reference-data.txt");
    
    while (f.good()) {
        string line;
        getline(f, line);
        
        if (line.length() == 0)
            break;
        
        istringstream tokens(line);
        
        Vector3d reference_result;
        tokens >> reference_result(0) >> reference_result(1) >> reference_result(2);
        
        reference_results.push_back(reference_result);
    }
    
    f.close();
    
    // Compute new values, and compare.
    for (unsigned int i = 0; i < reference_results.size(); i++) {        
        Vector3d &reference_result = reference_results[i];
        
        Vector2d res = run_test(reference_result(0) / 180.0 * pi);
        
        if (fabs(res[0] - reference_result(1)) > TEST_TOLERANCE) {
            cerr << " *** TEST FAILED *** " << endl;
            cerr << " alpha = " << reference_result(0) << " deg" << endl;
            cerr << " C_L(ref) = " << reference_result(1) << endl;
            cerr << " C_L = " << res(0) << endl;
            cerr << " ******************* " << endl;
            
            exit(1);
        }
       
        if (fabs(res[1] - reference_result(2)) > TEST_TOLERANCE) {
            cerr << " *** TEST FAILED *** " << endl;
            cerr << " alpha = " << reference_result(0) << " deg" << endl;
            cerr << " C_D(ref) = " << reference_result(2) << endl;
            cerr << " C_D = " << res(1) << endl;
            cerr << " ******************* " << endl;
            
            exit(1);
        }
    }
    
    // Done.
    return 0;
}
