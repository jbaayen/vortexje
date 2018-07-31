//
// Vortexje -- Elliptic wing with NACA0012 airfoil.  Checks induced drag.
//
// Copyright (C) 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <cmath>
#include <iostream>
#include <fstream>

#include <vortexje/solver.hpp>
#include <vortexje/lifting-surface-builder.hpp>
#include <vortexje/shape-generators/airfoils/naca4-airfoil-generator.hpp>
#include <vortexje/surface-writers/vtk-surface-writer.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

static const double pi = 3.141592653589793238462643383279502884;

#define DELTA_CONVERGENCE 1e-1
#define DELTA_T           1e-2

#define C_D_TEST_TOLERANCE 1e-3

//#define OUTPUT_RESULTS

#ifdef OUTPUT_RESULTS
static ofstream f;
#endif

// Cosine rule:
static double
cosine_rule(int n_points, int i)
{
    return 0.5 * (1 - cos(pi * i / (double) n_points));
}

// Run a test for a single angle of attack:
bool
run_test(double alpha)
{
    // Set up parameters for unsteady simulation:
    Parameters::unsteady_bernoulli = true;
    Parameters::convect_wake       = true;
    
    // Create wing:
    shared_ptr<LiftingSurface> wing(new LiftingSurface("main"));
    
    LiftingSurfaceBuilder surface_builder(*wing);
    
    const double AR = 6.0;
    const double chord = 0.5;
    const double span = AR * chord;
    
    const int n_points_per_airfoil = 40;
    const int n_airfoils = AR * n_points_per_airfoil / 2;
    
    int trailing_edge_point_id;
    vector<int> prev_airfoil_nodes;
    
    vector<vector<int> > node_strips;
    vector<vector<int> > panel_strips;
    
    for (int i = 0; i < n_airfoils; i++) {
        double y        = -span / 2.0 + span * cosine_rule(n_airfoils - 1, i);
        double chord_y  = chord * sqrt(1.0 - pow(2.0 * y / span, 2));
        double offset_x = (chord - chord_y) / 2.0;
        
        vector_aligned<Vector3d> airfoil_points;
        if (i == 0 || i == n_airfoils - 1) {
            Vector3d tip_point(0.0, 0.0, 0.0);
            for (int j = 0; j < n_points_per_airfoil; j++)
                airfoil_points.push_back(tip_point);
                
        } else
            airfoil_points = NACA4AirfoilGenerator::generate(0, 0, 0.12, true, chord_y, n_points_per_airfoil, trailing_edge_point_id);
        
        for (int j = 0; j < (int) airfoil_points.size(); j++) {
            airfoil_points[j](0) += offset_x;
            airfoil_points[j](2) += y;
        }
             
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

    // Create body:
    shared_ptr<Body> body(new Body(string("elliptic-wing")));
    body->add_lifting_surface(wing);
    
    // Set up solver:
    Solver solver("test-elliptic-planform-log");
    solver.add_body(body);
    
    Vector3d freestream_velocity(30, 0, 0);
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

    // Compute lift and drag coefficients:
    Vector3d F_a = solver.force(body);

    double q = 0.5 * fluid_density * freestream_velocity.squaredNorm();

    double S = pi * chord * span / 4.0;

    double C_L = F_a(1) / (q * S);
    double C_D = F_a(0) / (q * S);
    
    double C_D_ref = 1.0 / pi * S / pow(span, 2) * pow(C_L, 2);
    
    // Log lift and drag coefficients:
#ifdef OUTPUT_RESULTS
    f << C_L << ' ' << C_D << ' ' << C_D_ref << endl;
#endif
    
    // Compare induced drag coefficients:
    if (fabs(C_D - C_D_ref) > C_D_TEST_TOLERANCE) {
        cerr << " *** TEST FAILED *** " << endl;
        cerr << " C_D(ref) = " << C_D_ref << endl;
        cerr << " C_D = " << C_D << endl;
        cerr << " ******************* " << endl;
        
        return false;
    }

    // Done.
    return true;
}

int
main (int argc, char **argv)
{
    double alpha_min = 0.0;
    double alpha_max = 8.0;
    double dalpha    = 1.0;
    
#ifdef OUTPUT_RESULTS
    f.open("test-elliptic-planform.txt");
#endif
    
    double alpha = alpha_min;
    while (alpha <= alpha_max) {     
        if (!run_test(alpha / 180.0 * pi))
            exit(1);
            
        alpha += dalpha;
    }

#ifdef OUTPUT_RESULTS    
    f.close();
#endif
    
    return 0;
}
