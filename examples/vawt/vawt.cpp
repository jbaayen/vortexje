//
// Vortexje -- Simple VAWT example.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <cmath>
#include <iostream>
#include <fstream>

#include <vortexje/solver.hpp>
#include <vortexje/lifting-surface-builder.hpp>
#include <vortexje/shape-generators/airfoils/naca4-airfoil-generator.hpp>
#include <vortexje/shape-generators/ellipse-generator.hpp>
#include <vortexje/surface-writers/vtk-surface-writer.hpp>
#include <vortexje/boundary-layers/dummy-boundary-layer.hpp>
#include <vortexje/empirical-wakes/ramasamy-leishman-wake.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

static const double pi = 3.141592653589793238462643383279502884;

#define N_BLADES        2
#define MILL_RADIUS     2.5
#define TIP_SPEED_RATIO 5
#define WIND_VELOCITY   6
//#define INCLUDE_TOWER

class Blade : public LiftingSurface
{
public:
    // Constructor:
    Blade(const string &id)
        : LiftingSurface(id)
    {
        // Create blade:
        LiftingSurfaceBuilder surface_builder(*this);
        
        const int n_points_per_airfoil = 32;
        const int n_airfoils = 21;
        
        const double chord = 0.75;
        const double span = 4.5;
        
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
        
        // Translate and rotate into the canonical coordinate system:
        Vector3d translation(-chord / 3.0, 0.0, -span / 2.0);
        translate(translation);
        
        rotate(Vector3d::UnitZ(), -pi / 2.0);
    }
};

class Tower : public Surface
{
public:
    // Constructor:
    Tower()
        : Surface("tower")
    {
        // Create cylinder:      
        SurfaceBuilder surface_builder(*this);
        
        const double r = 0.1;
        const double h = 4.5;
        
        const int n_points = 32;
        const int n_layers = 21;
        
        vector<int> prev_nodes;
        
        for (int i = 0; i < n_layers; i++) {
            vector_aligned<Vector3d> points =
                EllipseGenerator::generate(r, r, n_points);
            for (int j = 0; j < (int) points.size(); j++)
                points[j](2) += i * h / (double) (n_layers - 1);
                 
            vector<int> nodes = surface_builder.create_nodes_for_points(points);
            
            if (i > 0)
                vector<int> airfoil_panels = surface_builder.create_panels_between_shapes(nodes, prev_nodes);
                
            prev_nodes = nodes;
        }

        surface_builder.finish();
        
        // Translate into the canonical coordinate system:
        Vector3d translation(0.0, 0.0, -h / 2.0);
        translate(translation);
    }
};

class VAWT : public Body
{
public:
    double rotor_radius;
    
    // Constructor:
    VAWT(string   id,
         double   rotor_radius,
         int      n_blades,
         Vector3d position,
         double   theta_0,
         double   dthetadt) :
         Body(id), rotor_radius(rotor_radius)
    {
        // Initialize kinematics:
        this->position = position;
        this->velocity = Vector3d(0, 0, 0);
        this->attitude = AngleAxis<double>(theta_0, Vector3d::UnitZ());
        this->rotational_velocity = Vector3d(0, 0, dthetadt);
        
#ifdef INCLUDE_TOWER
        // Initialize tower:
        shared_ptr<Tower> tower(new Tower());
        add_non_lifting_surface(tower);
#endif
        
        // Initialize blades:
        for (int i = 0; i < n_blades; i++) {
            stringstream ss;
            ss << "blade_" << i;
            shared_ptr<Blade> blade(new Blade(ss.str()));
            
            Vector3d translation(rotor_radius, 0, 0);
            blade->translate(translation);
            
            double theta_blade = theta_0 + 2 * pi / n_blades * i;
            blade->rotate(Vector3d::UnitZ(), theta_blade);
            
            blade->translate(position);
            
            shared_ptr<RamasamyLeishmanWake> wake(new RamasamyLeishmanWake(blade));           
            add_lifting_surface(blade, wake);
        }
    }
    
    // Rotate:
    void
    rotate(double dt)
    { 
        // Compute new kinematic state:
        Quaterniond new_attitude = AngleAxis<double>(rotational_velocity(2) * dt, Vector3d::UnitZ()) * attitude;
        set_attitude(new_attitude);
    }
};

int
main (int argc, char **argv)
{    
    // Set simulation parameters:
    Parameters::convect_wake                  = true;
    Parameters::interpolation_layer_thickness = 1e-1;
    Parameters::wake_vortex_core_radius       = 1e-3;
    
    // Set up VAWT:
    Vector3d position(0, 0, 0);
    
    shared_ptr<VAWT> vawt(new VAWT(string("vawt"),
                                   MILL_RADIUS,
                                   N_BLADES,
                                   position,
                                   pi / 6.0,
                                   TIP_SPEED_RATIO * WIND_VELOCITY / MILL_RADIUS));
    
    // Set up solver:
    Solver solver("vawt-log");
    solver.add_body(vawt);
    
    Vector3d freestream_velocity(WIND_VELOCITY, 0, 0);
    solver.set_freestream_velocity(freestream_velocity);
    
    double fluid_density = 1.2;
    solver.set_fluid_density(fluid_density);
    
    // Set up file format for logging:
    VTKSurfaceWriter surface_writer;
    
    // Log shaft moments:
    ofstream f;
    f.open("vawt-log/shaft_moment.txt");
    
    // Run simulation:
    double t = 0.0;
    double dt = 0.0033;
    int step_number = 0;
    
    solver.initialize_wakes(dt);
    while (t < 60) {
        // Solve:
        solver.solve(dt);
        
        // Log coefficients:
        solver.log(step_number, surface_writer);
        
        // Log shaft moment:
        Vector3d M = solver.moment(vawt, position);
        f << M(2) << endl;

        // Rotate blades:
        vawt->rotate(dt);
        
        // Update wakes:
        solver.update_wakes(dt);
        
        // Step time:
        t += dt;
        step_number++;
    }
    
    // Close shaft log file:
    f.close();
    
    // Done:
    return 0;
}
