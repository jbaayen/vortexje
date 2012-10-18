//
// Vortexje -- Solver.
//
// Copyright (C) 2012 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __SOLVER_HPP__
#define __SOLVER_HPP__

#include <vector>
#include <string>
#include <fstream>

#include <Eigen/Core>

#include <vortexje/collection.hpp>

namespace Vortexje
{

class Solver
{
public:
    // Constructor:
    std::string log_folder;
    std::vector<Mesh*> meshes_without_wakes;
    int total_n_panels_without_wakes;
    
    Solver(std::string log_folder);
    
    // Destructor:
    ~Solver();

    // Collection management:
    std::vector<Collection*> collections;
    
    void add_collection(Collection &collection);
    
    // Wind management:
    Eigen::Vector3d wind_velocity;
    
    void set_wind_velocity(Eigen::Vector3d value);
    
    // Air density:
    double air_density;
    
    void set_air_density(double value);
    
    // Solver stepping:
    Eigen::VectorXd source_coefficients;
    Eigen::VectorXd doublet_coefficients;
    Eigen::VectorXd pressure_coefficients;
    Eigen::VectorXd potentials;
    
    void initialize_wakes(double dt);
    void update_coefficients(double dt);
    void update_wakes(double dt);
    
    // Aerodynamic force and moment computation:
    Eigen::Vector3d aerodynamic_force(Collection &collection);
    Eigen::Vector3d aerodynamic_moment(Collection &collection, Eigen::Vector3d x);
    
    // Logging:
    void log_coefficients(int step_number);
    
private:
    void doublet_coefficient_matrix_block(Eigen::MatrixXd &A, Eigen::VectorXd &b,
                                          Mesh &mesh_one, int offset_one, Mesh &mesh_two, int offset_two);
                                          
    void wakes_influence(Eigen::MatrixXd &A, Eigen::VectorXd &b, Mesh &mesh, int offset);
                                          
    double source_coefficient(Mesh &mesh, int panel, Eigen::Vector3d &kinematic_velocity, bool include_wake_influence);
    
    Eigen::Vector3d surface_velocity(Mesh &mesh, int panel, Eigen::VectorXd &doublet_coefficient_field, Eigen::Vector3d &kinematic_velocity);

    double pressure_coefficient(Mesh &mesh, int panel, Eigen::Vector3d &kinematic_velocity,
                                Eigen::VectorXd &doublet_coefficient_field, double dpotentialdt, double v_ref);
    
    double potential(Eigen::Vector3d &x);
    
    Eigen::VectorXd surface_potentials();
    
    Eigen::Vector3d potential_gradient(Eigen::Vector3d &x);
    
    Eigen::Vector3d stream_velocity(Eigen::Vector3d &x, Eigen::Vector3d &kinematic_velocity);
};

};

#endif // __SOLVER_HPP__
