//
// Vortexje -- NACA0012 airfoil section in harmonic oscillation.
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

int
main (int argc, char **argv)
{
    // Set up parameters for unsteady simulation:
    Parameters::unsteady_bernoulli = true;
    Parameters::convect_wake       = true;
    
    // Create wing:
    shared_ptr<LiftingSurface> wing(new LiftingSurface("main"));
    
    LiftingSurfaceBuilder surface_builder(*wing);

    const double chord = 0.5;
    const double span = 1.0;
    
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
    
    // Create body:
    shared_ptr<Body> body(new Body(string("wing-section")));
    body->add_lifting_surface(wing);
    
    // Set up oscillation:
    double alpha_max = 10.0 / 180.0 * pi;
    
    double omega = 2 * pi / 1.0;
    
    // Set up solver:
    Solver solver("naca0012-section-oscillation-log");
    solver.add_body(body);
    
    Vector3d freestream_velocity(30, 0, 0);
    solver.set_freestream_velocity(freestream_velocity);
    
    double fluid_density = 1.2;
    solver.set_fluid_density(fluid_density);
    
    double dt = 0.01;
    
    // Set up logging:
    VTKSurfaceWriter writer;
    
    ofstream f;
    f.open("naca0012-section-oscillation-log/forces.txt");
    
    // Compute:
    double alpha = 0.0;
    double t = 0.0;
    int i = 0;

    solver.initialize_wakes(dt);

    while (t < 30.0) {
        // Solve:
        solver.solve(dt);
        
        // Log source, doublet, and pressure distributions:
        solver.log(i, writer);
        
        // Log lift and drag coefficients:
        Vector3d F_a = solver.force(body);

        double q = 0.5 * fluid_density * freestream_velocity.squaredNorm();

        double S = chord * span;

        double C_L = F_a(1) / (q * S);
        double C_D = F_a(0) / (q * S);

        f << alpha << ' ' << C_L << ' ' << C_D << endl;
        
        // Rotate wing:
        alpha = alpha_max * sin(omega * t);
        Quaterniond attitude = AngleAxis<double>(alpha, Vector3d::UnitZ()) * Quaterniond(1, 0, 0, 0);
        body->set_attitude(attitude);
        
        // Update rotational velocity:
        double dalphadt = alpha_max * omega * cos(omega * t);
        body->set_rotational_velocity(Vector3d(0, 0, dalphadt));
        
        // Update wake:
        solver.update_wakes(dt);
        
        // Step time:
        t += dt;
        i++;
    }  
    
    f.close();

    // Done.
    return 0;
}
