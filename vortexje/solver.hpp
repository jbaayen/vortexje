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
        
    Solver(const std::string log_folder);
    
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
    
    double velocity_potential(const Eigen::Vector3d &x) const;
    
    Eigen::Vector3d velocity(const Eigen::Vector3d &x) const;
    
    double pressure_coefficient(const Mesh &mesh, int panel) const;
    
    Eigen::Vector3d force(const Collection &collection) const;
    Eigen::Vector3d moment(const Collection &collection, const Eigen::Vector3d &x) const;
    
    void log_coefficients(int step_number, Mesh::FileFormat format) const;
    
    void log_fields(int step_number, Mesh::FileFormat format, double dx, double dy, double dz, double x_margin, double y_margin, double z_margin) const;

private:
    std::string log_folder;
    std::vector<Mesh*> meshes;
    std::vector<Mesh*> meshes_without_wakes;
    int total_n_panels_without_wakes;
    
    std::map<const Mesh*, Collection*> mesh_to_collection;
    
    Eigen::VectorXd source_coefficients;   
    Eigen::VectorXd doublet_coefficients;
    Eigen::VectorXd pressure_coefficients;
    Eigen::VectorXd velocity_potentials;
    
    void doublet_coefficient_matrix_block(Eigen::MatrixXd &doublet_influence_coefficients, 
                                          Eigen::MatrixXd &source_influence_coefficients,
                                          const Mesh &mesh_one, int offset_one, const Mesh &mesh_two, int offset_two) const;
                                          
    void wakes_influence(Eigen::MatrixXd &A, Mesh &mesh, int offset) const;
                                          
    double source_coefficient(const Mesh &mesh, int panel, const Eigen::Vector3d &kinematic_velocity, bool include_wake_influence) const;
    
    Eigen::Vector3d surface_velocity(const Mesh &mesh, int panel, const Eigen::VectorXd &doublet_coefficient_field) const;
    
    double reference_velocity(const Collection &collection) const;

    double pressure_coefficient(const Mesh &mesh, int panel, const Eigen::VectorXd &doublet_coefficient_field, double dphidt, double v_ref) const;
    
    double surface_velocity_potential(const Mesh &mesh, int offset, int panel) const;
    
    Eigen::VectorXd surface_velocity_potentials() const;
    
    Eigen::Vector3d disturbance_potential_gradient(const Eigen::Vector3d &x) const;
    
    double velocity_potential_time_derivative(const Eigen::VectorXd &velocity_potentials, const Eigen::VectorXd &old_velocity_potentials, int offset, int panel, double dt) const;
    
    void log_fields_vtk(int step_number, double dx, double dy, double dz, double x_margin, double y_margin, double z_margin) const;
    
    void log_fields_gmsh(int step_number, double dx, double dy, double dz, double x_margin, double y_margin, double z_margin) const;
};

};

#endif // __SOLVER_HPP__
