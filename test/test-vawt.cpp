//
// Vortexje -- VAWT test case.
//
// Copyright (C) 2012 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <iostream>

#include <vortexje/solver.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

#define N_BLADES        2
#define MILL_RADIUS     2.5
#define TIP_SPEED_RATIO 5
#define WIND_VELOCITY   6

class VAWT : public Collection
{
public:
    // Constructor:
    double rotor_radius;
    
    VAWT(string   id,
         Mesh    &tower_mesh,
         Mesh    &blade_mesh,
         double   rotor_radius,
         int      n_blades,
         Vector3d position,
         double   theta_0,
         double   dthetadt) :
         Collection(id, tower_mesh), rotor_radius(rotor_radius)
    {
        // Initialize kinematics:
        this->position = position;
        this->velocity = Vector3d(0, 0, 0);
        this->attitude = AngleAxis<double>(theta_0, Vector3d::UnitZ());
        this->rotational_velocity = Vector3d(0, 0, dthetadt);
        
        // Initialize blades:
        Vector3d unit_z = Vector3d::UnitZ();
        for (int i = 0; i < n_blades; i++) {
            Vector3d location(0, 0, 0);
            Vector3d chord_direction(1, 0, 0);
            Vector3d top_direction(0, 1, 0);
            Vector3d span_direction(0, 0, 1);
    
            Wing *blade = new Wing(blade_mesh, location, chord_direction, top_direction, span_direction);
            
            blade->rotate(unit_z, -M_PI / 2.0);
            
            Vector3d translation(rotor_radius, 0, 0);
            blade->translate(translation);
            
            double theta_blade = theta_0 + 2 * M_PI / n_blades * i;
            blade->rotate(unit_z, theta_blade);
            
            blade->translate(position);
            
            wings.push_back(blade);
            meshes_without_wakes.push_back(blade);
            
            Wake *wake = new Wake(*blade);
            wakes.push_back(wake);
        }
    }

    // Destructor:
    ~VAWT()
    {
        for (int i = 0; i < (int) wings.size(); i++) {
            delete wings[i];
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
    
    // Set up mill:
    Mesh tower_mesh;
    Mesh blade_mesh(string("naca0012.msh"));
    
    VAWT vawt(string("test-mill"),
              tower_mesh,
              blade_mesh,
              MILL_RADIUS,
              N_BLADES,
              Vector3d(0, 0, 0),
              M_PI / 6.0,
              TIP_SPEED_RATIO * WIND_VELOCITY / MILL_RADIUS);
    
    // Set up solver:
    Solver solver("test-vawt-log");
    solver.add_collection(vawt);
    
    Vector3d freestream_velocity(WIND_VELOCITY, 0, 0);
    solver.set_freestream_velocity(freestream_velocity);
    
    double fluid_density = 1.2;
    solver.set_fluid_density(fluid_density);
    
    // Run simulation:
    double t = 0.0;
    double dt = 0.0033;
    int step_number = 0;
    
    solver.initialize_wakes(dt);
    while (t < 60) {
        solver.update_coefficients(dt);
        solver.log_coefficients(step_number);

        solver.update_wakes(dt);
        
        vawt.rotate(dt);
        
        t += dt;
        step_number++;
    }
    
    // Done:
    return 0;
}
