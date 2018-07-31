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

#define TEST_TOLERANCE 1e-9

int
main (int argc, char **argv)
{
    // Set parameters:
    Parameters::interpolation_layer_thickness = 2e-2;
    Parameters::convect_wake                  = false;
    Parameters::static_wake_length            = 1e2;
    
    // Create wing:
    shared_ptr<LiftingSurface> wing(new LiftingSurface("main"));
    
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
    
    vector_aligned<Vector3d> airfoil_points =
        NACA4AirfoilGenerator::generate(0.09, 0.4, 0.12, true, chord, n_points_per_airfoil, trailing_edge_point_id);
    for (int i = 0; i < n_points_per_airfoil; i++)
        airfoil_points[i] = rotation * airfoil_points[i];
        
    for (int i = 0; i < n_airfoils; i++) {
        vector_aligned<Vector3d> translated_airfoil_points = airfoil_points;
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
    shared_ptr<Body> body(new Body(string("wing-section")));
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

    // Check for a number of perturbations:
    vector<double> perturbations;
    perturbations.push_back(0);
    perturbations.push_back(0.1);
    perturbations.push_back(0.99);
    
    for (int j = 0; j < (int) perturbations.size(); j++) {
        // Check velocities on and above panel collocation points:
        for (int i = 0; i < n_points_per_airfoil; i++) {  
            Vector3d point, velocity, reference_velocity;
            
            Vector3d collocation_point = wing->panel_collocation_point(i, false);
            
            Vector3d inner_velocity = solver.surface_velocity(wing, i);
            
            int next_i;
            if (i == n_points_per_airfoil - 1)
                next_i = 0;
            else
                next_i = i + 1;
                
            Vector3d perturbation_max = 0.5 * (airfoil_points[next_i] - airfoil_points[i]);
            
            Vector3d perturbed_collocation_point = collocation_point + perturbations[j] * perturbation_max;
            
            Vector3d outer_velocity = solver.velocity(perturbed_collocation_point - (Parameters::interpolation_layer_thickness + Parameters::zero_threshold) * wing->panel_normal(i));
            
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
            point = perturbed_collocation_point - (Parameters::interpolation_layer_thickness - Parameters::zero_threshold) * wing->panel_normal(i);
            
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
            
            // Don't run the below tests on the trailing edge:
            if (i == 4 || i == 5)
                continue;
            
            // Don't run the below tests on sharp edges, when perturbed:
            if ((i == 0 || i == 9 || i == 4 || i == 5) && j > 0)
                continue;

            // Check velocity inside body. 
            point = perturbed_collocation_point + 1e-2 * Parameters::interpolation_layer_thickness * wing->panel_normal(i);
            
            velocity           = solver.velocity(point);
            reference_velocity = (1.0 - 1e-2) * inner_velocity + 1e-2 * freestream_velocity;
            
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
    
        // Check velocities on and above corner points:
        for (int i = 0; i < n_points_per_airfoil; i++) {
            Vector3d point, velocity, reference_velocity;
            double w;
            
            int prev_i;
            if (i == 0)
                prev_i = n_points_per_airfoil - 1;
            else
                prev_i = i - 1;
                
            Vector3d inner_point = airfoil_points[i] + Vector3d(0, 0, 0.5);
                
            Vector3d prev_normal = wing->panel_normal(prev_i);
            Vector3d normal      = wing->panel_normal(i);
            
            double lambda_1 = 0.5 + 0.5 * perturbations[j];
            double lambda_2 = 1.0 - lambda_1;
                
            Vector3d direction = lambda_1 * prev_normal + lambda_2 * normal;
            direction.normalize();
            
            Vector3d inner_velocity = 0.5 * solver.surface_velocity(wing, prev_i) + 0.5 * solver.surface_velocity(wing, i);
            
            Vector3d outer_point = inner_point - Parameters::interpolation_layer_thickness * direction;

            Vector3d outer_velocity = solver.velocity(outer_point);
            
            bool convex = (i <= 5 || i == 9);
                                         
            // Check surface velocity:
            point = inner_point;
            
            velocity           = solver.velocity(point);
            reference_velocity = inner_velocity;
            
            if ((velocity - reference_velocity).norm() > TEST_TOLERANCE) {
                cerr << " *** INTERPOLATION LAYER (CORNER) SURFACE VELOCITY TEST FAILED *** " << endl;
                cerr << " node = " << i << endl;
                cerr << " perturbation = " << perturbations[j] << endl;
                cerr << " |V_ref| = " << reference_velocity.norm() << endl;
                cerr << " |V| = " << velocity.norm() << endl;
                cerr << " ******************* " << endl;
                
                exit(1);
            }
            
            // Don't run the below tests on sharp edges, when perturbed:
            if ((i == 0 || i == 9 || i == 4 || i == 5) && j > 0)
                continue;
           
            // Check velocity inside interpolation layer:
            point = inner_point - 0.25 * Parameters::interpolation_layer_thickness * direction;
            
            if (convex) { 
                // Convex corners:
                reference_velocity = (1.0 - 0.25) * inner_velocity + 0.25 * outer_velocity;
                
            } else {
                // Concave corners:
                Vector3d prev_transformed = wing->panel_coordinate_transformation(prev_i) * point;
                Vector3d transformed      = wing->panel_coordinate_transformation(i) * point;
                
                lambda_1 = -prev_transformed(2);
                lambda_2 = -transformed(2);
                
                prev_transformed(2) = -Parameters::interpolation_layer_thickness;
                transformed(2)      = -Parameters::interpolation_layer_thickness;
                
                Vector3d prev_transformed_inv = wing->panel_coordinate_transformation(prev_i).inverse() * prev_transformed;
                Vector3d transformed_inv      = wing->panel_coordinate_transformation(i).inverse() * transformed;
                
                prev_transformed = wing->panel_coordinate_transformation(i) * prev_transformed_inv;
                transformed      = wing->panel_coordinate_transformation(prev_i) * transformed_inv;
                
                double lambda_3 = -prev_transformed(2);
                double lambda_4 = -transformed(2);
                
                prev_transformed(2) = -Parameters::interpolation_layer_thickness;
                transformed(2)      = -Parameters::interpolation_layer_thickness;
                
                prev_transformed_inv = wing->panel_coordinate_transformation(i).inverse() * prev_transformed;
                transformed_inv      = wing->panel_coordinate_transformation(prev_i).inverse() * transformed;
                
                lambda_1 /= Parameters::interpolation_layer_thickness;
                lambda_2 /= Parameters::interpolation_layer_thickness;
                lambda_3 /= Parameters::interpolation_layer_thickness;
                lambda_4 /= Parameters::interpolation_layer_thickness;

                reference_velocity = 0.5 * (((1 - lambda_1) + lambda_2 * (1 - lambda_4)) * solver.surface_velocity(wing, prev_i)
                    + ((1 - lambda_2) + lambda_1 * (1 - lambda_3)) * solver.surface_velocity(wing, i)
                    + lambda_1 * lambda_3 * solver.velocity(prev_transformed_inv)
                    + lambda_2 * lambda_4 * solver.velocity(transformed_inv));
                    
            }
            
            velocity = solver.velocity(point);
            
            if ((velocity - reference_velocity).norm() > TEST_TOLERANCE) {
                cerr << " *** INTERPOLATION LAYER (CORNER) INNER VELOCITY TEST FAILED *** " << endl;
                cerr << " node = " << i << endl;
                cerr << " perturbation = " << perturbations[j] << endl;
                cerr << " |V_ref| = " << reference_velocity.norm() << endl;
                cerr << " |V| = " << velocity.norm() << endl;
                cerr << " ******************* " << endl;
                
                exit(1);
            }
            
            // Check velocity at edge of interpolation layer:
            point = inner_point - (Parameters::interpolation_layer_thickness - Parameters::zero_threshold) * direction;
            
            velocity           = solver.velocity(point);
            reference_velocity = outer_velocity;
            
            if ((velocity - reference_velocity).norm() > TEST_TOLERANCE) {
                cerr << " *** INTERPOLATION LAYER (CORNER) EDGE VELOCITY TEST FAILED *** " << endl;
                cerr << " node = " << i << endl;
                cerr << " perturbation = " << perturbations[j] << endl;
                cerr << " |V_ref| = " << reference_velocity.norm() << endl;
                cerr << " |V| = " << velocity.norm() << endl;
                cerr << " ******************* " << endl;
                
                exit(1);
            }
            
            // Don't run the below tests on sharp edges, or when perturbed:
            if ((i == 0 || i == 9 || i == 4 || i == 5) || j > 0)
                continue;
            
            // Check velocity inside body:
            point = inner_point + 1e-2 * Parameters::interpolation_layer_thickness * direction;

            if (convex) {
                // Convex corners:
                lambda_1 = (wing->panel_coordinate_transformation(prev_i) * point)(2);
                lambda_2 = (wing->panel_coordinate_transformation(i) * point)(2);
            
                w = 0.5 * (lambda_1 + lambda_2) / Parameters::interpolation_layer_thickness;      
            } else {
                // Concave corners:
                w = 1e-2;
            }
                
            velocity           = solver.velocity(point);
            reference_velocity = (1.0 - w) * inner_velocity + w * freestream_velocity;
            
            if ((velocity - reference_velocity).norm() > TEST_TOLERANCE) {
                cerr << " *** INSIDE BODY (CORNER) VELOCITY TEST FAILED *** " << endl;
                cerr << " node = " << i << endl;
                cerr << " perturbation = " << perturbations[j] << endl;
                cerr << " |V_ref| = " << reference_velocity.norm() << endl;
                cerr << " |V| = " << velocity.norm() << endl;
                cerr << " ******************* " << endl;
                
                exit(1);
            }
        }
    }
    
    // Done:
    return 0;
}
