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
#include <vortexje/surface-writer.hpp>

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
        
    Solver(const std::string &log_folder);
    
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
    
    void initialize_wakes(double dt = 0.0);
    
    void update_wakes(double dt = 0.0);
    
    bool solve(double dt = 0.0, bool propagate = true);
    
    void propagate();
    
    double velocity_potential(const Eigen::Vector3d &x) const;
    
    Eigen::Vector3d velocity(const Eigen::Vector3d &x) const;
    
    double surface_velocity_potential(const Surface &surface, int panel) const;
    
    Eigen::Vector3d surface_velocity(const Surface &surface, int panel) const;
    
    double pressure_coefficient(const Surface &surface, int panel) const;
    
    Eigen::Vector3d force(const Body &body) const;
    Eigen::Vector3d moment(const Body &body, const Eigen::Vector3d &x) const;
    
    /**
       Data structure bundling a Surface, a panel ID, and a point on the panel.
       
       @brief Surface, panel ID, and point bundle.
    */
    class SurfacePanelPoint {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        /**
           Constructor.
           
           @param[in]   surface   Associated Surface object.
           @param[in]   panel     Panel ID.
           @param[in]   point     Point.
        */
        SurfacePanelPoint(const Surface &surface, int panel, const Eigen::Vector3d &point) : surface(&surface), panel(panel), point(point) { };
        
        /**
           Associated Surface object.
        */
        const Surface *surface;
        
        /**
           Panel ID.
        */
        int panel;
        
        /**
           Point.
        */
        Eigen::Vector3d point;
    };
    
    std::vector<SurfacePanelPoint> trace_streamline(const SurfacePanelPoint &start) const;
    
    void log(int step_number, SurfaceWriter &writer) const;

private:
    std::string log_folder;
    
    std::vector<Body::SurfaceData*> non_wake_surfaces;
    int n_non_wake_panels;
    
    std::map<int, Body*> surface_id_to_body;
    
    Eigen::VectorXd source_coefficients;   
    Eigen::VectorXd doublet_coefficients;
        
    Eigen::VectorXd surface_velocity_potentials;
    Eigen::MatrixXd surface_velocities;
    Eigen::VectorXd pressure_coefficients;  
    
    Eigen::VectorXd previous_surface_velocity_potentials; 
                                          
    double compute_source_coefficient(const Body &body, const Surface &surface, int panel, const BoundaryLayer &boundary_layer, bool include_wake_influence) const;
    
    double compute_surface_velocity_potential(const Surface &surface, int offset, int panel) const;
    
    double compute_surface_velocity_potential_time_derivative(int offset, int panel, double dt) const;
    
    Eigen::Vector3d compute_surface_velocity(const Body &body, const Surface &surface, int panel) const;
    
    double compute_reference_velocity_squared(const Body &body) const;

    double compute_pressure_coefficient(const Eigen::Vector3d &surface_velocity, double dphidt, double v_ref) const;
    
    Eigen::Vector3d compute_disturbance_velocity(const Eigen::Vector3d &x) const;
    
    double compute_disturbance_velocity_potential(const Eigen::Vector3d &x) const;
    
    Eigen::Vector3d compute_trailing_edge_vortex_displacement(const Body &body, const LiftingSurface &lifting_surface, int index, double dt) const;

    Eigen::Vector3d compute_scalar_field_gradient(const Eigen::VectorXd &scalar_field, const Body &body, const Surface &surface, int panel) const;
    
    int compute_index(const Surface &surface, int panel) const;
};

};

#endif // __SOLVER_HPP__
