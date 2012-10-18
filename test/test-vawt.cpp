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
        for (int i = 0; i < n_blades; i++) {
            Vector3d location(0, 0, 0);
            Vector3d chord_direction(1, 0, 0);
            Vector3d top_direction(0, 1, 0);
            Vector3d span_direction(0, 0, 1);
    
            Wing *blade = new Wing(blade_mesh, location, chord_direction, top_direction, span_direction);
            
            blade->rotate(Vector3d::UnitZ(), -M_PI / 2.0);
            
            blade->translate(Vector3d(rotor_radius, 0, 0));
            
            double theta_blade = theta_0 + 2 * M_PI / n_blades * i;
            blade->rotate(Vector3d::UnitZ(), theta_blade);
            
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
        for (int i = 0; i < wings.size(); i++) {
            delete wings[i];
            delete wakes[i];
        }
    }

    // Rotate:
    void
    rotate(double dt)
    { 
        // Compute new kinematic state:
        Quaterniond q_rotational_velocity(0.0, 0.0, 0.0, 0.5 * rotational_velocity(2) * dt);
        Quaterniond dattitude = q_rotational_velocity * this->attitude;
        Quaterniond attitude = Quaterniond(this->attitude.w() + dattitude.w(), 
                                           this->attitude.x() + dattitude.x(),
                                           this->attitude.y() + dattitude.y(),
                                           this->attitude.z() + dattitude.z());
        set_attitude(attitude);
    }
};

int
main (int argc, char **argv)
{
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
    solver.set_wind_velocity(Vector3d(WIND_VELOCITY, 0, 0));
    solver.set_air_density(1.2);
    
    // Run simulation:
    double t = 0.0;
    double dt = 0.0033;
    int step_number = 0;
    
    solver.initialize_wakes(dt);
    while (t < 60) {
        solver.update_coefficients(dt);
        solver.log_coefficients(step_number);
        
        vawt.rotate(dt);
        
        solver.update_wakes(dt);
        
        t += dt;
        step_number++;
    }
    
    // Done:
    return 0;
}
