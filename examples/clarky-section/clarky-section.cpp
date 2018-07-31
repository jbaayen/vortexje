//
// Vortexje -- Clark-Y wing section construction example.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <cmath>
#include <iostream>

#include <vortexje/solver.hpp>
#include <vortexje/parameters.hpp>
#include <vortexje/lifting-surface-builder.hpp>
#include <vortexje/surface-writers/vtk-surface-writer.hpp>
#include <vortexje/field-writers/vtk-field-writer.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

static const double pi = 3.141592653589793238462643383279502884;

// Load airfoil data from file:
static vector_aligned<Vector3d>
read_airfoil(const std::string &filename, int &trailing_edge_point_id)
{
    vector_aligned<Vector3d> upper_points;
    vector_aligned<Vector3d> lower_points;
    
    // Parse file:
    ifstream f;
    f.open(filename.c_str());
    
    int side = -1;
    
    while (f.good()) {
        string line;
        getline(f, line);
        
        // Empty line signifies new section:
        bool empty_line = false;
        if (line.find_first_not_of("\r\n\t ") == string::npos)
            empty_line = true;
            
        if (empty_line) {
            side++;
            continue;
        }
            
        // Read x, y coordinates:
        if (side >= 0) {
            istringstream tokens(line);
            
            double x, y;
            tokens >> x >> y;
            
            Vector3d point(x, y, 0.0);
            if (side == 0)
                upper_points.push_back(point);
            else
                lower_points.push_back(point);
        }
    }
    
    // Close file:
    f.close();
    
    // Assemble entire airfoil:
    vector_aligned<Vector3d> points;
    
    for (int i = 0; i < (int) upper_points.size() - 1; i++)
        points.push_back(upper_points[i]);
    
    // Use a thin trailing edge:
    Vector3d trailing_edge_point = 0.5 * (upper_points[(int) upper_points.size() - 1] + lower_points[(int) lower_points.size() - 1]);
    points.push_back(trailing_edge_point);
    trailing_edge_point_id = (int) points.size() - 1;
    
    for (int i = (int) lower_points.size() - 2; i > 0; i--)
        points.push_back(lower_points[i]);

    // Done:
    return points;
}

// Main:
int
main(int argc, char **argv)
{
    // Enable wake convection:
    Parameters::convect_wake = true;
    
    // Load airfoil data:
    int trailing_edge_point_id;
    vector_aligned<Vector3d> clarky_airfoil = read_airfoil("clarky.dat", trailing_edge_point_id);
    
    // Create lifting surface object:
    shared_ptr<LiftingSurface> wing(new LiftingSurface("main"));

    // Construct wing section:
    LiftingSurfaceBuilder surface_builder(*wing);

    const int n_airfoils = 21;
    
    const double chord = 1.0;
    const double span = 5.0;
    
    vector<int> prev_airfoil_nodes;
    
    vector<vector<int> > node_strips;
    vector<vector<int> > panel_strips;
    
    for (int i = 0; i < n_airfoils; i++) {
        vector_aligned<Vector3d> airfoil_points;
        for (int j = 0; j < (int) clarky_airfoil.size(); j++)
            airfoil_points.push_back(Vector3d(chord * clarky_airfoil[j](0), chord * clarky_airfoil[j](1), i * span / (double) (n_airfoils - 1)));
             
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
    wing->translate(translation);
    
    // Prescribe angle of attack:
    double alpha = 5.0 / 180.0 * pi;
    wing->rotate(Vector3d::UnitZ(), -alpha);
    
    // Create surface body:
    shared_ptr<Body> body(new Body(string("wing-section")));
    body->add_lifting_surface(wing);
    
    // Set up solver:
    Solver solver("clarky-section-log");
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
        //field_writer.write_velocity_field(solver, "velocity-field.vtk", -0.5, 2.0, -3.0, 3.0, -1.0, 1.0, 0.2, 0.2, 0.2);
        
        // Update wake:
        solver.update_wakes(dt);
        
        // Step time:
        t += dt;
        step_number++;
    }
    
    // Done:
    return 0;
}
