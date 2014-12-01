//
// Vortexje -- Simple HAWT example.
//
// Copyright (C) 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <vortexje/solver.hpp>
#include <vortexje/lifting-surface-builder.hpp>
#include <vortexje/shape-generators/airfoils/naca4-airfoil-generator.hpp>
#include <vortexje/surface-writers/vtk-surface-writer.hpp>
#include <vortexje/boundary-layers/dummy-boundary-layer.hpp>
#include <vortexje/empirical-wakes/ramasamy-leishman-wake.hpp>

#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

#define N_BLADES        3
#define ROTOR_RADIUS    5.5320
#define TIP_SPEED_RATIO 3
#define WIND_VELOCITY   6.0

class Blade : public LiftingSurface
{
public:
    // Constructor:
    Blade(const string &id)
        : LiftingSurface(id)
    {
        // Read data:
        read_airfoil("hawt-airfoil.dat");
        read_blade_geometry("hawt-blade.dat");
        
        // Create blade:
        LiftingSurfaceBuilder surface_builder(*this);
        
        vector<int> prev_airfoil_nodes;
        
        vector<vector<int> > node_strips;
        vector<vector<int> > panel_strips;
        
        for (int i = 0; i < (int) blade_dr.size(); i++) {
            vector<Vector3d, Eigen::aligned_allocator<Vector3d> > airfoil_points;
            for (int j = 0; j < (int) unscaled_airfoil_points.size(); j++) {
                Vector3d airfoil_point(blade_chord[i] * (unscaled_airfoil_points[j](0) - 0.25),
                                       blade_thickness[i] / unscaled_airfoil_thickness * unscaled_airfoil_points[j](1),
                                       blade_dr[i]);
                                       
                airfoil_point = AngleAxis<double>(blade_twist[i], Vector3d::UnitZ()) * airfoil_point;
                                       
                airfoil_points.push_back(airfoil_point);
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
        
        rotate(Vector3d::UnitZ(), -M_PI / 2.0);
    }
    
private:
    vector<Vector3d, Eigen::aligned_allocator<Vector3d> > unscaled_airfoil_points;
    double unscaled_airfoil_thickness;
    
    int trailing_edge_point_id;
    
    vector<double> blade_dr;
    vector<double> blade_chord;
    vector<double> blade_twist;
    vector<double> blade_thickness;
    
    // Load airfoil data from file:
    void
    read_airfoil(const std::string &filename)
    {
        vector<Vector3d, Eigen::aligned_allocator<Vector3d> > upper_points;
        vector<Vector3d, Eigen::aligned_allocator<Vector3d> > lower_points;
        
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
        for (int i = 0; i < (int) upper_points.size() - 1; i++)
            unscaled_airfoil_points.push_back(upper_points[i]);
        
        // Use a thin trailing edge:
        Vector3d trailing_edge_point = 0.5 * (upper_points[(int) upper_points.size() - 1] + lower_points[(int) lower_points.size() - 1]);
        unscaled_airfoil_points.push_back(trailing_edge_point);
        trailing_edge_point_id = (int) unscaled_airfoil_points.size() - 1;
        
        for (int i = (int) lower_points.size() - 2; i > 0; i--)
            unscaled_airfoil_points.push_back(lower_points[i]);
            
        // Compute thickness:
        double max_y = 0.0;
        for (int i = 0; i < (int) upper_points.size(); i++) {
            if (upper_points[i](1) > max_y)
                max_y = upper_points[i](1);
        }
        
        double min_y = 0.0;
        for (int i = 0; i < (int) lower_points.size(); i++) {
            if (lower_points[i](1) < min_y)
                min_y = lower_points[i](1);
        }
        
        unscaled_airfoil_thickness = max_y - min_y;
    }
    
    // Load blade geometry from file:
    void
    read_blade_geometry(const std::string &filename)
    {
        // Parse file:
        ifstream f;
        f.open(filename.c_str());
        
        while (f.good()) {
            string line;
            getline(f, line);
            
            // Read section geometry:
            if (line.substr(0, 2) != "dr" && line.length() > 0) {
                istringstream tokens(line);
                
                double dr, chord, twist, thick;
                tokens >> dr >> chord >> twist >> thick;
                
                blade_dr.push_back(dr);
                blade_chord.push_back(chord);
                blade_twist.push_back(M_PI * twist / 180.0);
                blade_thickness.push_back(thick);
            }
        }
        
        // Close file:
        f.close();
    }
};

class HAWT : public Body
{
public:
    // Constructor:
    HAWT(string   id,
         int      n_blades,
         Vector3d position,
         double   theta_0,
         double   dthetadt) :
         Body(id)
    {
        // Initialize kinematics:
        this->position = position;
        this->velocity = Vector3d(0, 0, 0);
        this->attitude = AngleAxis<double>(theta_0, Vector3d::UnitX());
        this->rotational_velocity = Vector3d(-dthetadt, 0, 0);
        
        // Initialize blades:
        for (int i = 0; i < n_blades; i++) {
            stringstream ss;
            ss << "blade_" << i;
            shared_ptr<Blade> blade(new Blade(ss.str()));
            
            double theta_blade = theta_0 + 2 * M_PI / n_blades * i;
            blade->rotate(Vector3d::UnitX(), theta_blade);
            
            blade->translate(position);

            add_lifting_surface(blade);
        }
    }
    
    // Rotate:
    void
    rotate(double dt)
    { 
        // Compute new kinematic state:
        Quaterniond new_attitude = AngleAxis<double>(rotational_velocity(0) * dt, Vector3d::UnitX()) * attitude;
        set_attitude(new_attitude);
    }
};

int
main (int argc, char **argv)
{    
    // Set simulation parameters:
    Parameters::unsteady_bernoulli            = true;
    Parameters::convect_wake                  = true;
    Parameters::interpolation_layer_thickness = 1e-1;
    Parameters::wake_vortex_core_radius       = 1e-3;
    
    // Set up HAWT:
    Vector3d position(0, 0, 0);
    
    shared_ptr<HAWT> hawt(new HAWT(string("hawt"),
                              N_BLADES,
                              position,
                              0.0,
                              TIP_SPEED_RATIO * WIND_VELOCITY / ROTOR_RADIUS));
    
    // Set up solver:
    Solver solver("hawt-log");
    solver.add_body(hawt);
    
    Vector3d freestream_velocity(WIND_VELOCITY, 0, 0);
    solver.set_freestream_velocity(freestream_velocity);
    
    double fluid_density = 1.2;
    solver.set_fluid_density(fluid_density);
    
    // Set up file format for logging:
    VTKSurfaceWriter surface_writer;
    
    // Log shaft moments:
    ofstream f;
    f.open("hawt-log/shaft_moment.txt");
    
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
        Vector3d M = solver.moment(hawt, position);
        f << M(0) << endl;
        
        // Rotate blades:
        hawt->rotate(dt);
        
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
