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

#include <vortexje/body.hpp>

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
       List of surface bodies.
    */
    std::vector<Body*> bodies;
    
    void add_body(Body &body);
    
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
    
    double pressure_coefficient(const Surface &surface, int panel) const;
    
    Eigen::Vector3d force(const Body &body) const;
    Eigen::Vector3d moment(const Body &body, const Eigen::Vector3d &x) const;
    
    void log_coefficients(int step_number, Surface::FileFormat format) const;
    
    void log_fields(int step_number, Surface::FileFormat format, double dx, double dy, double dz, double x_margin, double y_margin, double z_margin) const;

private:
    std::string log_folder;
    std::vector<Surface*> surfaces;
    std::vector<Surface*> non_wake_surfaces;
    int n_non_wake_panels;
    
    std::map<const Surface*, Body*> surface_to_body;
    
    Eigen::VectorXd source_coefficients;   
    Eigen::VectorXd doublet_coefficients;
    Eigen::VectorXd pressure_coefficients;
    Eigen::VectorXd velocity_potentials;
    
    void doublet_coefficient_matrix_block(Eigen::MatrixXd &doublet_influence_coefficients, 
                                          Eigen::MatrixXd &source_influence_coefficients,
                                          const Surface &surface_one, int offset_one, const Surface &surface_two, int offset_two) const;
                                          
    void wakes_influence(Eigen::MatrixXd &A, Surface &surface, int offset) const;
                                          
    double source_coefficient(const Surface &surface, int panel, const Eigen::Vector3d &kinematic_velocity, bool include_wake_influence) const;
    
    Eigen::Vector3d surface_velocity(const Surface &surface, int panel, const Eigen::VectorXd &doublet_coefficient_field) const;
    
    double reference_velocity(const Body &body) const;

    double pressure_coefficient(const Surface &surface, int panel, const Eigen::VectorXd &doublet_coefficient_field, double dphidt, double v_ref) const;
    
    double surface_velocity_potential(const Surface &surface, int offset, int panel) const;
    
    Eigen::VectorXd surface_velocity_potentials() const;
    
    Eigen::Vector3d disturbance_potential_gradient(const Eigen::Vector3d &x) const;
    
    double velocity_potential_time_derivative(const Eigen::VectorXd &velocity_potentials, const Eigen::VectorXd &old_velocity_potentials, int offset, int panel, double dt) const;
    
    void log_fields_vtk(int step_number, double dx, double dy, double dz, double x_margin, double y_margin, double z_margin) const;
    
    void log_fields_gmsh(int step_number, double dx, double dy, double dz, double x_margin, double y_margin, double z_margin) const;
};

};

#endif // __SOLVER_HPP__
