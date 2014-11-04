//
// Vortexje -- Test velocity interpolation layer.
//
// Copyright (C) 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <iostream>
#include <fstream>
#include <limits>

#include <vortexje/solver.hpp>
#include <vortexje/lifting-surface-builder.hpp>
#include <vortexje/shape-generators/airfoils/naca4-airfoil-generator.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

static const double pi = 3.141592653589793238462643383279502884;

#define TEST_TOLERANCE  1e-12
#define INSIDE_DISTANCE 1e-2

int
main (int argc, char **argv)
{
    // Set parameters:
    Parameters::interpolation_layer_thickness = 2e-2;
    Parameters::convect_wake                  = false;
    Parameters::static_wake_length            = 1e2;
    
    // Create wing:
    shared_ptr<LiftingSurface> wing(new LiftingSurface());
    
    LiftingSurfaceBuilder surface_builder(*wing);

    const double chord = 1.0;
    const double span = 1.0;
    
    const int n_points_per_airfoil = 10;
    const int n_airfoils = 2;
    
    int trailing_edge_point_id;
    vector<int> prev_airfoil_nodes;
    
    vector<vector<int> > node_strips;
    vector<vector<int> > panel_strips;
    
    double alpha = pi * 5.0 / 180.0;
    AngleAxis<double> rotation(-alpha, Vector3d::UnitZ());
    
    vector<Vector3d, Eigen::aligned_allocator<Vector3d> > airfoil_points =
        NACA4AirfoilGenerator::generate(0, 0, 0.12, true, chord, n_points_per_airfoil, trailing_edge_point_id);
    for (int i = 0; i < n_points_per_airfoil; i++)
        airfoil_points[i] = rotation * airfoil_points[i];
        
    for (int i = 0; i < n_airfoils; i++) {
        vector<Vector3d, Eigen::aligned_allocator<Vector3d> > translated_airfoil_points = airfoil_points;
        for (int j = 0; j < (int) airfoil_points.size(); j++)
            translated_airfoil_points[j](2) += i * span / (double) (n_airfoils - 1);
             
        vector<int> airfoil_nodes = surface_builder.create_nodes_for_points(translated_airfoil_points);
        node_strips.push_back(airfoil_nodes);
        
        if (i > 0) {
            vector<int> airfoil_panels = surface_builder.create_panels_between_shapes(airfoil_nodes, prev_airfoil_nodes, trailing_edge_point_id);
            panel_strips.push_back(airfoil_panels);
        }
            
        prev_airfoil_nodes = airfoil_nodes;
    }

    surface_builder.finish(node_strips, panel_strips, trailing_edge_point_id);
    
    // Create wake:
    shared_ptr<Wake> wake(new Wake(wing));
    
    // Create surface body:
    shared_ptr<Body> body(new Body(string("test-wing")));
    body->add_lifting_surface(wing, wake);
    
    // Set up solver:
    Solver solver("test-interpolation-layer-log");
    solver.add_body(body);
    
    Vector3d freestream_velocity(30.0, 0, 0);
    solver.set_freestream_velocity(freestream_velocity);
    
    double fluid_density = 1.2;
    solver.set_fluid_density(fluid_density);
    
    // Run simulation:
    solver.initialize_wakes();
    solver.solve();

    // Check velocities on and above panel collocation points:
    vector<double> perturbations;
    perturbations.push_back(0);
    perturbations.push_back(1e-1);
    
    for (int j = 0; j < (int) perturbations.size(); j++) {
        for (int i = 0; i < n_points_per_airfoil; i++) {  
            Vector3d point, velocity, reference_velocity;
            
            Vector3d collocation_point = wing->panel_collocation_point(i, false);
            
            Vector3d inner_velocity = solver.surface_velocity(wing, i);
            Vector3d outer_velocity = solver.velocity(collocation_point - (Parameters::interpolation_layer_thickness + Parameters::inversion_tolerance) * wing->panel_normal(i));
            
            int next_i;
            if (i == n_points_per_airfoil - 1)
                next_i = 0;
            else
                next_i = i + 1;
                
            Vector3d perturbation_dir = airfoil_points[next_i] - airfoil_points[i];
            
            Vector3d perturbed_collocation_point = collocation_point + perturbations[j] * perturbation_dir;
            
            // Check surface velocity:
            point = perturbed_collocation_point;
            
            velocity           = solver.velocity(point);
            reference_velocity = inner_velocity;
            
            if ((velocity - reference_velocity).norm() > TEST_TOLERANCE) {
                cerr << " *** INTERPOLATION LAYER SURFACE VELOCITY TEST FAILED *** " << endl;
                cerr << " panel = " << i << endl;
                cerr << " perturbation = " << perturbations[j] << endl;
                cerr << " |V_ref| = " << reference_velocity.norm() << endl;
                cerr << " |V| = " << velocity.norm() << endl;
                cerr << " ******************* " << endl;
                
                exit(1);
            }
            
            // Check velocity inside interpolation layer:
            point = perturbed_collocation_point - 0.25 * Parameters::interpolation_layer_thickness * wing->panel_normal(i);
            
            velocity           = solver.velocity(point);
            reference_velocity = 0.75 * inner_velocity + 0.25 * outer_velocity;
            
            if ((velocity - reference_velocity).norm() > TEST_TOLERANCE) {
                cerr << " *** INTERPOLATION LAYER INNER VELOCITY TEST FAILED *** " << endl;
                cerr << " panel = " << i << endl;
                cerr << " perturbation = " << perturbations[j] << endl;
                cerr << " |V_ref| = " << reference_velocity.norm() << endl;
                cerr << " |V| = " << velocity.norm() << endl;
                cerr << " ******************* " << endl;
                
                exit(1);
            }
            
            // Check velocity at edge of interpolation layer:
            point = perturbed_collocation_point - (Parameters::interpolation_layer_thickness - Parameters::inversion_tolerance) * wing->panel_normal(i);
            
            velocity           = solver.velocity(point);
            reference_velocity = outer_velocity;
            
            if ((velocity - reference_velocity).norm() > TEST_TOLERANCE) {
                cerr << " *** INTERPOLATION LAYER EDGE VELOCITY TEST FAILED *** " << endl;
                cerr << " panel = " << i << endl;
                cerr << " perturbation = " << perturbations[j] << endl;
                cerr << " |V_ref| = " << reference_velocity.norm() << endl;
                cerr << " |V| = " << velocity.norm() << endl;
                cerr << " ******************* " << endl;
                
                exit(1);
            }
            
            // Check velocity inside body:
            point = perturbed_collocation_point + INSIDE_DISTANCE * wing->panel_normal(i);
            
            velocity           = solver.velocity(point);
            reference_velocity = freestream_velocity;
            
            if ((velocity - reference_velocity).norm() > TEST_TOLERANCE) {
                cerr << " *** INSIDE BODY VELOCITY TEST FAILED *** " << endl;
                cerr << " panel = " << i << endl;
                cerr << " perturbation = " << perturbations[j] << endl;
                cerr << " |V_ref| = " << reference_velocity.norm() << endl;
                cerr << " |V| = " << velocity.norm() << endl;
                cerr << " ******************* " << endl;
                
                exit(1);
            }
        }
    }
    
    // Check velocities on and above corner points:
    for (int i = 0; i < n_points_per_airfoil; i++) {
        Vector3d point, velocity, reference_velocity;
        
        int prev_i;
        if (i == 0)
            prev_i = n_points_per_airfoil - 1;
        else
            prev_i = i - 1;
        
        Vector3d collocation_point      = wing->panel_collocation_point(i     , false);
        Vector3d prev_collocation_point = wing->panel_collocation_point(prev_i, false);
        
        Vector3d inner_velocity = 0.5 * (solver.surface_velocity(wing, prev_i) + solver.surface_velocity(wing, i));
        Vector3d outer_velocity = 0.5 * (solver.velocity(prev_collocation_point - (Parameters::interpolation_layer_thickness + Parameters::inversion_tolerance) * wing->panel_normal(prev_i)) +
                                         solver.velocity(collocation_point - (Parameters::interpolation_layer_thickness + Parameters::inversion_tolerance) * wing->panel_normal(i)));
                                         
        Vector3d airfoil_point = airfoil_points[i];

        Vector3d direction = wing->panel_normal(prev_i) + wing->panel_normal(i);
        direction.normalize();
            
        // Check surface velocity:
        point = airfoil_point;
        
        velocity           = solver.velocity(point);
        reference_velocity = inner_velocity;
        
        if ((velocity - reference_velocity).norm() > TEST_TOLERANCE) {
            cerr << " *** INTERPOLATION LAYER (CORNER) SURFACE VELOCITY TEST FAILED *** " << endl;
            cerr << " node = " << i << endl;
            cerr << " |V_ref| = " << reference_velocity.norm() << endl;
            cerr << " |V| = " << velocity.norm() << endl;
            cerr << " ******************* " << endl;
            
            exit(1);
        }
        
        // Check velocity inside interpolation layer:
        point = airfoil_point - 0.25 * Parameters::interpolation_layer_thickness * direction;
        
        velocity           = solver.velocity(point);
        reference_velocity = 0.75 * inner_velocity + 0.25 * outer_velocity;
        
        if ((velocity - reference_velocity).norm() > TEST_TOLERANCE) {
            cerr << " *** INTERPOLATION LAYER (CORNER) INNER VELOCITY TEST FAILED *** " << endl;
            cerr << " node = " << i << endl;
            cerr << " |V_ref| = " << reference_velocity.norm() << endl;
            cerr << " |V| = " << velocity.norm() << endl;
            cerr << " ******************* " << endl;
            
            exit(1);
        }
        
        // Check velocity at edge of interpolation layer:
        point = airfoil_point - (Parameters::interpolation_layer_thickness - Parameters::inversion_tolerance) * direction;
        
        velocity           = solver.velocity(point);
        reference_velocity = outer_velocity;
        
        if ((velocity - reference_velocity).norm() > TEST_TOLERANCE) {
            cerr << " *** INTERPOLATION LAYER (CORNER) EDGE VELOCITY TEST FAILED *** " << endl;
            cerr << " node = " << i << endl;
            cerr << " |V_ref| = " << reference_velocity.norm() << endl;
            cerr << " |V| = " << velocity.norm() << endl;
            cerr << " ******************* " << endl;
            
            exit(1);
        }
        
        // Check velocity inside body:
        point = airfoil_point + 1e-4 * direction;

        velocity           = solver.velocity(point);
        reference_velocity = freestream_velocity;
        
        if ((velocity - reference_velocity).norm() > TEST_TOLERANCE) {
            cerr << " *** INSIDE BODY (CORNER) VELOCITY TEST FAILED *** " << endl;
            cerr << " node = " << i << endl;
            cerr << " |V_ref| = " << reference_velocity.norm() << endl;
            cerr << " |V| = " << velocity.norm() << endl;
            cerr << " ******************* " << endl;
            
            exit(1);
        }
    }
    
    // Done:
    return 0;
}
