//
// Vortexje -- Simple VAWT example.
//
// Copyright (C) 2012 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <iostream>

#include <vortexje/solver.hpp>
#include <vortexje/lifting-surface-builder.hpp>
#include <vortexje/airfoils/naca4.hpp>

#include <iostream>
#include <fstream>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

#define N_BLADES        2
#define MILL_RADIUS     2.5
#define TIP_SPEED_RATIO 5
#define WIND_VELOCITY   6

class Blade : public LiftingSurface
{
public:
    // Constructor:
    Blade()
    {
        // Set up local coordinate system:
        chord_direction = Vector3d(1, 0, 0);
        span_direction  = Vector3d(0, 0, 1);
        
        // Create blade:
        LiftingSurfaceBuilder surface_builder(*this);
        
        int trailing_edge_point_id;
        vector<int> prev_airfoil_nodes;
        
        const int n_points_per_airfoil = 32;
        const int n_airfoils = 21;
        
        const double chord = 0.75;
        const double span = 4.5;
        
        for (int j = 0; j < n_airfoils; j++) {
            vector<Vector3d, Eigen::aligned_allocator<Vector3d> > airfoil_points =
                Airfoils::NACA4::generate(0, 0, 0.12, true, chord, n_points_per_airfoil, trailing_edge_point_id);
            for (int k = 0; k < (int) airfoil_points.size(); k++)
                airfoil_points[k](2) += j * span / (double) (n_airfoils - 1);
                
            vector<int> airfoil_nodes = surface_builder.create_nodes(airfoil_points, trailing_edge_point_id);
            
            if (j > 0)
                surface_builder.create_panels_between(airfoil_nodes, prev_airfoil_nodes, trailing_edge_point_id);
                
            prev_airfoil_nodes = airfoil_nodes;
        }

        compute_topology();
        compute_geometry();
        
        sort_strips();
        
        // Translate and rotate into the canonical coordinate system:
        Vector3d translation(-chord / 3.0, 0.0, -span / 2.0);
        translate(translation);
        
        rotate(Vector3d::UnitZ(), -M_PI / 2.0);
    }
};

class VAWT : public Body
{
public:
    // Constructor:
    double rotor_radius;
    
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
        
        // Initialize blades:
        for (int i = 0; i < n_blades; i++) {
            Blade *blade = new Blade();
            
            Vector3d translation(rotor_radius, 0, 0);
            blade->translate(translation);
            
            double theta_blade = theta_0 + 2 * M_PI / n_blades * i;
            blade->rotate(Vector3d::UnitZ(), theta_blade);
            
            blade->translate(position);
            
            lifting_surfaces.push_back(blade);
            non_wake_surfaces.push_back(blade);
            
            Wake *wake = new Wake(*blade);
            wakes.push_back(wake);
        }
    }

    // Destructor:
    ~VAWT()
    {
        for (int i = 0; i < (int) lifting_surfaces.size(); i++) {
            delete lifting_surfaces[i];
            delete wakes[i];
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
    // Set simulation parameters for dynamic wake convection with Ramasamy-Leishman vortex model:
    Parameters::convect_wake                       = true;
    Parameters::use_ramasamy_leishman_vortex_sheet = true;
        
    // Shut up separation warnings:
    Parameters::min_pressure_coefficient = -1e10;
    
    // Set up VAWT:
    Vector3d position(0, 0, 0);
    
    VAWT vawt(string("vawt"),
              MILL_RADIUS,
              N_BLADES,
              position,
              M_PI / 6.0,
              TIP_SPEED_RATIO * WIND_VELOCITY / MILL_RADIUS);
    
    // Set up solver:
    Solver solver("vawt-log");
    solver.add_body(vawt);
    
    Vector3d freestream_velocity(WIND_VELOCITY, 0, 0);
    solver.set_freestream_velocity(freestream_velocity);
    
    double fluid_density = 1.2;
    solver.set_fluid_density(fluid_density);
    
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
        solver.update_coefficients(dt);
        
        // Log coefficients:
        solver.log_coefficients(step_number, Surface::VTK);
        
        // Log shaft moment:
        Vector3d M = solver.moment(vawt, position);
        f << M(2) << endl;

        // Update wakes:
        solver.update_wakes(dt);
        
        // Rotate blades:
        vawt.rotate(dt);
        
        // Step time:
        t += dt;
        step_number++;
    }
    
    // Close shaft log file:
    f.close();
    
    // Done:
    return 0;
}
