//
// Vortexje -- Solver.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __SOLVER_HPP__
#define __SOLVER_HPP__

#include <string>
#include <fstream>

#include <Eigen/Core>
#include <Eigen/StdVector>

#include <vortexje/collection.hpp>

namespace Vortexje
{

/**
   Class for solution of the panel method equations, and their propagation in time.
   
   @brief Panel method solver.
*/
class Solver
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
    Solver(std::string log_folder);
    
    ~Solver();

    /**
       List of mesh collections.
    */
    std::vector<Collection*> collections;
    
    void add_collection(Collection &collection);
    
    /**
       Freestream velocity.
    */
    Eigen::Vector3d freestream_velocity;
    
    void set_freestream_velocity(const Eigen::Vector3d &value);
    
    /**
       Density of the fluid.
    */
    double fluid_density;
    
    void set_fluid_density(double value);
    
    void initialize_wakes(double dt);
    void update_coefficients(double dt);
    void update_wakes(double dt);
    
    double velocity_potential(const Eigen::Vector3d &x);
    
    Eigen::Vector3d velocity(const Eigen::Vector3d &x);
    
    double pressure_coefficient(const Mesh &mesh, int panel) const;
    
    Eigen::Vector3d aerodynamic_force(const Collection &collection) const;
    Eigen::Vector3d aerodynamic_moment(const Collection &collection, const Eigen::Vector3d &x) const;
    
    void log_coefficients(int step_number) const;

private:
    std::string log_folder;
    std::vector<Mesh*> meshes;
    std::vector<Mesh*> meshes_without_wakes;
    int total_n_panels_without_wakes;
    
    std::map<Mesh*, Collection*> mesh_to_collection;
    
    Eigen::VectorXd source_coefficients;   
    Eigen::VectorXd doublet_coefficients;
    Eigen::VectorXd pressure_coefficients;
    Eigen::VectorXd velocity_potentials;
    
    void doublet_coefficient_matrix_block(Eigen::MatrixXd &doublet_influence_coefficients, 
                                          Eigen::MatrixXd &source_influence_coefficients,
                                          Mesh &mesh_one, int offset_one, Mesh &mesh_two, int offset_two);
                                          
    void wakes_influence(Eigen::MatrixXd &A, Mesh &mesh, int offset);
                                          
    double source_coefficient(Mesh &mesh, int panel, const Eigen::Vector3d &kinematic_velocity, bool include_wake_influence);
    
    Eigen::Vector3d surface_velocity(Mesh &mesh, int panel, const Eigen::VectorXd &doublet_coefficient_field);
    
    double reference_velocity(const Collection &collection) const;

    double pressure_coefficient(Mesh &mesh, int panel, const Eigen::VectorXd &doublet_coefficient_field, double dphidt, double v_ref);   
    
    Eigen::VectorXd surface_velocity_potentials();
    
    Eigen::Vector3d disturbance_potential_gradient(const Eigen::Vector3d &x);
    
    double velocity_potential_time_derivative(const Eigen::VectorXd &velocity_potentials, const Eigen::VectorXd &old_velocity_potentials, int offset, int panel, double dt);
};

};

#endif // __SOLVER_HPP__
