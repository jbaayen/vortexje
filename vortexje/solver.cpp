//
// Vortexje -- Solver.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

#include <iostream>
#include <limits>

#ifdef _WIN32
#include <direct.h>
#endif

#include <Eigen/Geometry>
#include <Eigen/IterativeLinearSolvers>

#include <vortexje/solver.hpp>
#include <vortexje/parameters.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

// Helper to create folders:
static void
mkdir_helper(const string folder)
{
#ifdef _WIN32
    if (mkdir(folder.c_str()) < 0)
#else
    if (mkdir(folder.c_str(), S_IRWXU) < 0)
#endif
        if (errno != EEXIST)
            cerr << "Could not create log folder " << folder << ": " << strerror(errno) << endl;
}

/**
   Construct a solver, logging its output into the given folder.
   
   @param[in]   log_folder  Logging output folder.
*/
Solver::Solver(const string log_folder) : log_folder(log_folder)
{ 
    // Initialize wind:
    freestream_velocity = Vector3d(0, 0, 0);
    
    // Initialize fluid density:
    fluid_density = 0.0;
    
    // Total number of panels:
    n_non_wake_panels = 0;
        
    // Open log files:
    mkdir_helper(log_folder);
}

/**
   Destructor.
*/
Solver::~Solver()
{
}

/**
   Adds a surface body to this solver.
   
   @param[in]   body   Body to be added.
*/
void
Solver::add_body(Body &body)
{
    bodies.push_back(&body);
    
    vector<Surface*>::iterator si;
    for (si = body.non_lifting_surfaces.begin(); si != body.non_lifting_surfaces.end(); si++) {
        Surface *non_lifting_surface = *si;
        
        surfaces.push_back(non_lifting_surface);
        non_wake_surfaces.push_back(non_lifting_surface);
           
        surface_to_body[non_lifting_surface] = &body;
        
        n_non_wake_panels = n_non_wake_panels + non_lifting_surface->n_panels();
    }
    
    vector<LiftingSurface*>::iterator lsi;
    vector<Wake*>::iterator wi;
    for (lsi = body.lifting_surfaces.begin(), wi = body.wakes.begin(); lsi != body.lifting_surfaces.end(); lsi++, wi++) {
        LiftingSurface *lifting_surface = *lsi;
        Wake *wake = *wi;
        
        surfaces.push_back(lifting_surface);
        surfaces.push_back(wake);
        non_wake_surfaces.push_back(lifting_surface);
           
        surface_to_body[lifting_surface] = &body;
        surface_to_body[wake] = &body;
        
        n_non_wake_panels = n_non_wake_panels + lifting_surface->n_panels();
    }
    
    doublet_coefficients.resize(n_non_wake_panels);
    doublet_coefficients.setZero();
    
    source_coefficients.resize(n_non_wake_panels);
    source_coefficients.setZero();
    
    surface_velocity_potentials.resize(n_non_wake_panels);
    surface_velocity_potentials.setZero();
    
    surface_velocities.resize(n_non_wake_panels, 3);
    surface_velocities.setZero();
    
    pressure_coefficients.resize(n_non_wake_panels);
    pressure_coefficients.setZero();

    // Open logs:
    string body_log_folder = log_folder + "/" + body.id;
    
    mkdir_helper(body_log_folder);
    
    for (int i = 0; i < (int) body.non_lifting_surfaces.size(); i++) {
        stringstream ss;
        ss << body_log_folder << "/non_lifting_surface_" << i;
        
        string s = ss.str();
        mkdir_helper(s);
    }
    
    for (int i = 0; i < (int) body.lifting_surfaces.size(); i++) {
        stringstream ss;
        ss << body_log_folder << "/lifting_surface_" << i;
        
        string s = ss.str();
        mkdir_helper(s);
        
        stringstream ss2;
        ss2 << body_log_folder << "/wake_" << i;
        
        s = ss2.str();
        mkdir_helper(s);      
    }
}

/**
   Sets the freestream velocity.
   
   @param[in]   value   Freestream velocity.
*/
void
Solver::set_freestream_velocity(const Vector3d &value)
{
    freestream_velocity = value;
}

/**
   Sets the fluid density.
   
   @param[in]   value   Fluid density.
*/
void
Solver::set_fluid_density(double value)
{
    fluid_density = value;
}

// Compute source coefficient for given surface and panel:
double
Solver::source_coefficient(const Surface &surface, int panel, const Vector3d &apparent_velocity, bool include_wake_influence) const
{
    // Main velocity:
    Vector3d velocity = -apparent_velocity;
    
    // Wake contribution:
    if (Parameters::convect_wake && include_wake_influence) {
        vector<Body*>::const_iterator bi;
        for (bi = bodies.begin(); bi != bodies.end(); bi++) {
            const Body *body = *bi;
            
            vector<Wake*>::const_iterator wi;
            for (wi = body->wakes.begin(); wi != body->wakes.end(); wi++) {
                const Wake *wake = *wi;
                
                for (int k = 0; k < wake->n_panels(); k++) {
                    // Use doublet panel - vortex ring equivalence.  Any new wake panels have zero doublet
                    // coefficient, and are therefore not accounted for here.
                    if (Parameters::use_ramasamy_leishman_vortex_sheet)
                        velocity += wake->vortex_ring_ramasamy_leishman_velocity
                            (surface, panel, k, wake->vortex_core_radii[k], wake->doublet_coefficients[k]);
                    else
                        velocity += wake->vortex_ring_unit_velocity(surface, panel, k) * wake->doublet_coefficients[k];
                }
            }
        }
    }
    
    // Take normal component:
    Vector3d normal = surface.panel_normal(panel);
    return -velocity.dot(normal);
}

/**
   Computes the surface velocity for the given panel.
   
   @param[in]   surface                     Reference surface.
   @param[in]   panel                       Reference panel.
   @param[in]   doublet_coefficient_field   Doublet coefficient distribution on given surface.
   
   @returns Surface velocity.
*/
Eigen::Vector3d
Solver::surface_velocity(const Surface &surface, int panel, const Eigen::VectorXd &doublet_coefficient_field) const
{
    // Compute disturbance part of surface velocity.
    Vector3d tangential_velocity;
    if (Parameters::marcov_surface_velocity) {
        Vector3d x = surface.panel_collocation_point(panel, false);
        
        // Use N. Marcov's formula for surface velocity, see L. Dragoş, Mathematical Methods in Aerodynamics, Springer, 2003.
        Vector3d tangential_velocity = disturbance_potential_gradient(x);
        tangential_velocity -= 0.5 * surface.scalar_field_gradient(doublet_coefficient_field, panel);
    } else
        tangential_velocity = -surface.scalar_field_gradient(doublet_coefficient_field, panel);

    // Add flow due to kinematic velocity:
    Body *body = surface_to_body.find(&surface)->second;
    Vector3d apparent_velocity = body->panel_kinematic_velocity(surface, panel) - freestream_velocity;
                                          
    tangential_velocity -= apparent_velocity;
    
    // Remove any normal velocity.  This is the (implicit) contribution of the source term.
    Vector3d normal = surface.panel_normal(panel);
    tangential_velocity -= tangential_velocity.dot(normal) * normal;
    
    // Done:
    return tangential_velocity;
}

/**
   Returns the square of the reference velocity for the given body.
   
   @param[in]   body   Body to establish reference velocity for.
   
   @returns Square of the reference velocity.
*/
double
Solver::reference_velocity_squared(const Body &body) const
{
    return (body.velocity - freestream_velocity).squaredNorm();
}

/**
   Computes the pressure coefficient.
   
   @param[in]   surface_velocity   Surface velocity for the reference panel.
   @param[in]   dphidt             Time-derivative of the velocity potential for the reference panel.
   @param[in]   v_ref              Reference velocity.
   
   @returns Pressure coefficient.
*/
double
Solver::pressure_coefficient(const Vector3d &surface_velocity, double dphidt, double v_ref_squared) const
{
    double C_p = 1 - (surface_velocity.squaredNorm() + 2 * dphidt) / v_ref_squared;
    
    return C_p;
}

/**
   Computes velocity potential at the given point.
   
   @param[in]   x   Reference point.
   
   @returns Velocity potential.
*/
double
Solver::velocity_potential(const Vector3d &x) const
{
    double phi = 0.0;
    
    // Iterate all non-wake surfaces:
    int offset = 0;
    
    vector<Surface*>::const_iterator si;
    for (si = non_wake_surfaces.begin(); si != non_wake_surfaces.end(); si++) {
        const Surface *other_surface = *si;

        for (int i = 0; i < other_surface->n_panels(); i++) {
            phi += other_surface->doublet_influence(x, i) * doublet_coefficients(offset + i);
            phi += other_surface->source_influence(x, i) * source_coefficients(offset + i);
        }
        
        offset += other_surface->n_panels();
    }
    
    // Iterate wakes:
    vector<Body*>::const_iterator bi;
    for (bi = bodies.begin(); bi != bodies.end(); bi++) {
        const Body *body = *bi;
        
        vector<Wake*>::const_iterator wi;
        for (wi = body->wakes.begin(); wi != body->wakes.end(); wi++) {
            const Wake *wake = *wi;
            
            for (int i = 0; i < wake->n_panels(); i++)
                phi += wake->doublet_influence(x, i) * wake->doublet_coefficients[i];
        }
    }
                    
    // Sum with freestream velocity potential:
    return phi + freestream_velocity.dot(x);
}

/**
   Returns velocity potential value on the body surface.
   
   @returns Surface potential value.
*/
double
Solver::surface_velocity_potential(const Surface &surface, int offset, int panel) const
{
    if (Parameters::marcov_surface_velocity) {
        // Since we use N. Marcov's formula for surface velocity, we also compute the surface velocity
        // potential directly.
        return velocity_potential(surface.panel_collocation_point(panel, false));
        
    } else {
        double phi = -doublet_coefficients(offset + panel);
        
        // Add flow potential due to kinematic velocity:
        Body *body = surface_to_body.find(&surface)->second;
        Vector3d apparent_velocity = body->panel_kinematic_velocity(surface, panel) - freestream_velocity;
        
        phi -= apparent_velocity.dot(surface.panel_collocation_point(panel, false));
        
        return phi;
    }
}

/**
   Computes velocity potential values on the body surface.
   
   @returns Vector of velocity potential values, ordered by panel number.
*/
Eigen::VectorXd
Solver::surface_surface_velocity_potentials() const
{
    cout << "Solver: Computing surface potential values." << endl;
    
    VectorXd surface_surface_velocity_potentials(n_non_wake_panels);
    
    int offset = 0;  
 
    vector<Body*>::const_iterator bi;
    for (bi = bodies.begin(); bi != bodies.end(); bi++) {
        const Body *body = *bi;
        
        vector<Surface*>::const_iterator si;
        for (si = body->non_lifting_surfaces.begin(); si != body->non_lifting_surfaces.end(); si++) {
            const Surface *non_lifting_surface = *si;
            
            for (int i = 0; i < non_lifting_surface->n_panels(); i++)
                surface_surface_velocity_potentials(offset + i) = surface_velocity_potential(*non_lifting_surface, offset, i);
            
            offset += non_lifting_surface->n_panels();
        }
        
        vector<LiftingSurface*>::const_iterator lsi;
        for (lsi = body->lifting_surfaces.begin(); lsi != body->lifting_surfaces.end(); lsi++) {
            const LiftingSurface *lifting_surface = *lsi;
            
            for (int i = 0; i < lifting_surface->n_panels(); i++)
                surface_surface_velocity_potentials(offset + i) = surface_velocity_potential(*lifting_surface, offset, i);
            
            offset += lifting_surface->n_panels();
        }
    }
    
    return surface_surface_velocity_potentials;
}

/**
   Computes disturbance potential gradient at the given point.
   
   @param[in]  x    Reference point.
   
   @returns Disturbance potential gradient.
*/ 
Eigen::Vector3d
Solver::disturbance_potential_gradient(const Eigen::Vector3d &x) const
{
    Vector3d gradient(0, 0, 0);
    
    // Iterate all non-wake surfaces:
    int offset = 0;
    
    vector<Surface*>::const_iterator si;
    for (si = non_wake_surfaces.begin(); si != non_wake_surfaces.end(); si++) {
        const Surface *other_surface = *si;

        for (int i = 0; i < other_surface->n_panels(); i++) {
            gradient += other_surface->vortex_ring_unit_velocity(x, i) * doublet_coefficients(offset + i);
            gradient += other_surface->source_unit_velocity(x, i) * source_coefficients(offset + i);
        }
        
        offset += other_surface->n_panels();
    }
    
    // Iterate wakes:
    vector<Body*>::const_iterator bi;
    for (bi = bodies.begin(); bi != bodies.end(); bi++) {
        const Body *body = *bi;
        
        vector<LiftingSurface*>::const_iterator lsi;
        vector<Wake*>::const_iterator wi;
        for (lsi = body->lifting_surfaces.begin(), wi = body->wakes.begin(); lsi != body->lifting_surfaces.end(); lsi++, wi++) {
            const LiftingSurface *lifting_surface = *lsi;
            const Wake *wake = *wi;
            
            if (wake->n_panels() >= lifting_surface->n_spanwise_panels()) {
                for (int i = wake->n_panels() - lifting_surface->n_spanwise_panels(); i < wake->n_panels(); i++)
                    gradient += wake->vortex_ring_unit_velocity(x, i) * wake->doublet_coefficients[i];

                for (int i = 0; i < wake->n_panels() - lifting_surface->n_spanwise_panels(); i++) {
                    if (Parameters::use_ramasamy_leishman_vortex_sheet)
                        gradient += wake->vortex_ring_ramasamy_leishman_velocity(x, i, wake->vortex_core_radii[i], wake->doublet_coefficients[i]);
                    else
                        gradient += wake->vortex_ring_unit_velocity(x, i) * wake->doublet_coefficients[i];
                }
            }
        }
    }
               
    // Done:
    return gradient;
}

/**
   Computes velocity potential time derivative at the given panel.
   
   @param[in]  surface_velocity_potentials       Current potential values.
   @param[in]  old_surface_velocity_potentials   Previous potential values
   @param[in]  offset                            Offset to requested Surface
   @param[in]  panel                             Panel number.
   @param[in]  dt                                Time step size.
   
   @returns Velocity potential time derivative.
*/ 
double
Solver::velocity_potential_time_derivative(const Eigen::VectorXd &surface_velocity_potentials, const Eigen::VectorXd &old_surface_velocity_potentials, int offset, int panel, double dt) const
{
    double dphidt;
    
    // Evaluate the time-derivative of the potential in a body-fixed reference frame, as in
    //   J. P. Giesing, Nonlinear Two-Dimensional Unsteady Potential Flow with Lift, Journal of Aircraft, 1968.
    if (Parameters::unsteady_bernoulli && dt > 0.0)
        dphidt = (surface_velocity_potentials(offset + panel) - old_surface_velocity_potentials(offset + panel)) / dt;
    else
        dphidt = 0.0;
        
    return dphidt;
}

/**
   Computes the total stream velocity at the given point.
   
   @param[in]   x   Reference point.
   
   @returns Stream velocity.
*/
Eigen::Vector3d
Solver::velocity(const Eigen::Vector3d &x) const
{
    // Find closest surface and panel:
    double distance = numeric_limits<double>::max();
    const Surface *close_surface = NULL;
    int close_panel = -1;
    int close_offset = - 1;
    bool close_near_trailing_edge = false;
    int offset = 0;
    
    vector<Surface*>::const_iterator si;
    for (si = non_wake_surfaces.begin(); si != non_wake_surfaces.end(); si++) {
        const Surface *surface_candidate = *si;
        
        int panel_candidate;
        double distance_candidate;
        bool near_trailing_edge = surface_candidate->closest_panel(x, panel_candidate, distance_candidate);
        
        if (distance_candidate < distance) {
            distance                 = distance_candidate;
            close_surface            = surface_candidate;
            close_panel              = panel_candidate;
            close_offset             = offset;
            close_near_trailing_edge = near_trailing_edge;
        }
        
        offset += surface_candidate->n_panels();
    }
 
    // Compute velocity potential gradients near the body:
    vector<Vector3d> disturbance_potential_gradients;
    vector<Vector3d> near_exterior_points;
    if (distance < Parameters::interpolation_layer_thickness && !close_near_trailing_edge) {   
        for (int i = 0; i < (int) close_surface->panel_nodes[close_panel].size(); i++) {
            Vector3d near_exterior_point = close_surface->near_exterior_point(close_surface->panel_nodes[close_panel][i]);
            near_exterior_points.push_back(near_exterior_point);
            
            disturbance_potential_gradients.push_back(disturbance_potential_gradient(near_exterior_point));
        }
        
    } else {
        disturbance_potential_gradients.push_back(disturbance_potential_gradient(x));
        
        close_surface = NULL;
        close_panel   = -1;
    }
    
    Vector3d velocity(0, 0, 0);
    if (close_surface != NULL) { 
        // Compute stream velocity near the boundary using interpolation, see
        //   K. Dixon, C. S. Ferreira, C. Hofemann, G. van Brussel, G. van Kuik,
        //   A 3D Unsteady Panel Method for Vertical Axis Wind Turbines, DUWIND, 2008.
        Vector3d x = close_surface->panel_collocation_point(close_panel, false);
        
        VectorXd doublet_coefficient_field(close_surface->n_panels());
        for (int i = 0; i < close_surface->n_panels(); i++)
            doublet_coefficient_field(i) = doublet_coefficients(close_offset + i);
        
        Vector3d normal = close_surface->panel_normal(close_panel);
        
        double total_weight = 0.0;
        
        for (int i = 0; i < (int) disturbance_potential_gradients.size(); i++) {
            Vector3d layer_point_distance = x - near_exterior_points[i];
            layer_point_distance = layer_point_distance - layer_point_distance.dot(normal) * normal;
            
            double weight = distance * (close_surface->panel_diameter(close_panel) - layer_point_distance.norm());
                 
            velocity += weight * (disturbance_potential_gradients[i] + freestream_velocity);
            total_weight += weight;
        }
        
        double weight = Parameters::interpolation_layer_thickness - distance;
        velocity += weight * surface_velocity(*close_surface, close_panel, doublet_coefficient_field);
        total_weight += weight;
        
        velocity /= total_weight;
        
    } else {
        // We are not close to the boundary.  Use sum of disturbance potential and freestream velocity.
        velocity = disturbance_potential_gradients[0] + freestream_velocity;
    }
    
    return velocity;
}

/**
   Initializes the wakes by adding a first layer of vortex ring panels.
   
   @param[in]   dt   Time step size.
*/
void
Solver::initialize_wakes(double dt)
{
    // Add initial wake layers:
    vector<Body*>::iterator bi;
    for (bi = bodies.begin(); bi != bodies.end(); bi++) {
        Body *body = *bi;
        
        Vector3d body_apparent_velocity = body->velocity - freestream_velocity;
        
        vector<Wake*>::const_iterator wi;
        for (wi = body->wakes.begin(); wi != body->wakes.end(); wi++) {
            Wake *wake = *wi;
            
            wake->add_layer();
            for (int i = 0; i < wake->n_nodes(); i++) {
                if (Parameters::convect_wake)    
                    wake->nodes[i] -= body_apparent_velocity * dt;
                else
                    wake->nodes[i] -= Parameters::static_wake_length * body_apparent_velocity / body_apparent_velocity.norm();
            }
            
            wake->add_layer();
        }
    }
}

// Create block of linear system:
void
Solver::doublet_coefficient_matrix_block(MatrixXd &doublet_influence_coefficients,
                                         MatrixXd &source_influence_coefficients,
                                         const Surface &surface_one, int offset_one, const Surface &surface_two, int offset_two) const
{
    for (int i = 0; i < surface_one.n_panels(); i++) {
        for (int j = 0; j < surface_two.n_panels(); j++) {
            if ((&surface_one == &surface_two) && (i == j))
                doublet_influence_coefficients(offset_one + i, offset_two + j) = -0.5;
            else
                doublet_influence_coefficients(offset_one + i, offset_two + j) = surface_two.doublet_influence(surface_one, i, j);
            
            source_influence_coefficients(offset_one + i, offset_two + j) = surface_two.source_influence(surface_one, i, j);
        }
    }
}

// Add doublet influence of wakes to system matrix:
void
Solver::wake_influence(MatrixXd &A, Surface &surface, int offset) const
{
    for (int i = 0; i < surface.n_panels(); i++) {
        int lifting_surface_offset = 0;
        
        vector<Body*>::const_iterator bi;
        for (bi = bodies.begin(); bi != bodies.end(); bi++) {
            const Body *body = *bi;
            
            vector<Surface*>::const_iterator si;
            for (si = body->non_lifting_surfaces.begin(); si != body->non_lifting_surfaces.end(); si++)
                lifting_surface_offset += (*si)->n_panels();
            
            vector<LiftingSurface*>::const_iterator lsi;
            vector<Wake*>::const_iterator wi;
            for (lsi = body->lifting_surfaces.begin(), wi = body->wakes.begin(); lsi != body->lifting_surfaces.end(); lsi++, wi++) {
                const LiftingSurface *lifting_surface = *lsi;
                const Wake *wake = *wi;
                
                int wake_panel_offset = wake->n_panels() - lifting_surface->n_spanwise_panels();
                for (int j = 0; j < lifting_surface->n_spanwise_panels(); j++) {  
                    int pa = lifting_surface->trailing_edge_upper_panel(j);
                    int pb = lifting_surface->trailing_edge_lower_panel(j);
                    
                    // Account for the influence of the new wake panels.  The doublet strength of these panels
                    // is set according to the Kutta condition;  see below.
                    A(offset + i, lifting_surface_offset + pa) += wake->doublet_influence(surface, i, wake_panel_offset + j);
                    A(offset + i, lifting_surface_offset + pb) -= wake->doublet_influence(surface, i, wake_panel_offset + j);
                }
                
                lifting_surface_offset += lifting_surface->n_panels();
            }
        }
    }
}

/**
   Computes new source, doublet, and pressure distributions.
   
   @param[in]   dt   Time step size.
*/
void
Solver::update_coefficients(double dt)
{
    // Compute source distribution:
    cout << "Solver: Computing source distribution with wake influence." << endl;
    
    int idx = 0;
    
    vector<Body*>::iterator bi;
    for (bi = bodies.begin(); bi != bodies.end(); bi++) {
        Body *body = *bi;
        
        vector<Surface*>::iterator si;
        for (si = body->non_lifting_surfaces.begin(); si != body->non_lifting_surfaces.end(); si++) {
            Surface *non_lifting_surface = *si;
            
            for (int i = 0; i < non_lifting_surface->n_panels(); i++) {
                Vector3d apparent_velocity = body->panel_kinematic_velocity(*non_lifting_surface, i) - freestream_velocity;
            
                source_coefficients(idx) = source_coefficient(*non_lifting_surface, i, apparent_velocity, true);
                idx++;
            }
        }
        
        vector<LiftingSurface*>::iterator lsi;
        for (lsi = body->lifting_surfaces.begin(); lsi != body->lifting_surfaces.end(); lsi++) {
            LiftingSurface *lifting_surface = *lsi;
            
            for (int i = 0; i < lifting_surface->n_panels(); i++) {
                Vector3d apparent_velocity = body->panel_kinematic_velocity(*lifting_surface, i) - freestream_velocity;
            
                source_coefficients(idx) = source_coefficient(*lifting_surface, i, apparent_velocity, true);
                idx++;
            }
        }
    }
  
    // Compute doublet distribution:
    MatrixXd A(n_non_wake_panels, n_non_wake_panels);
    MatrixXd source_influence_coefficients(n_non_wake_panels, n_non_wake_panels);
    
    int offset_one = 0, offset_two = 0;
    for (int i = 0; i < (int) non_wake_surfaces.size(); i++) {
        offset_two = 0;
        for (int j = 0; j < (int) non_wake_surfaces.size(); j++) {
            doublet_coefficient_matrix_block(A,
                                             source_influence_coefficients,
                                             *non_wake_surfaces[i], offset_one, *non_wake_surfaces[j], offset_two);
            
            offset_two = offset_two + non_wake_surfaces[j]->n_panels();
        }
        
        wake_influence(A, *non_wake_surfaces[i], offset_one);
        
        offset_one = offset_one + non_wake_surfaces[i]->n_panels();
    }
    
    VectorXd b = source_influence_coefficients * source_coefficients;
    
    BiCGSTAB<MatrixXd> solver(A);
    solver.setMaxIterations(Parameters::linear_solver_max_iterations);
    solver.setTolerance(Parameters::linear_solver_tolerance);

    cout << "Solver: Computing doublet distribution." << endl;
    
    doublet_coefficients = solver.solveWithGuess(b, doublet_coefficients);
    
    if (solver.info() != Success) {
        cerr << "Solver: Computing doublet distribution failed (" << solver.iterations();
        cerr << " iterations with estimated error=" << solver.error() << ")." << endl;
        exit(1);
    }
    
    cout << "Solver: Done in " << solver.iterations() << " iterations with estimated error " << solver.error() << "." << endl;
    
    // Set new wake panel doublet coefficients:
    int offset = 0;
    for (bi = bodies.begin(); bi != bodies.end(); bi++) {
        Body *body = *bi;
        
        vector<Surface*>::iterator si;
        for (si = body->non_lifting_surfaces.begin(); si != body->non_lifting_surfaces.end(); si++)
            offset += (*si)->n_panels();
        
        vector<LiftingSurface*>::iterator lsi;
        vector<Wake*>::iterator wi;
        for (lsi = body->lifting_surfaces.begin(), wi = body->wakes.begin(); lsi != body->lifting_surfaces.end(); lsi++, wi++) {
            LiftingSurface *lifting_surface = *lsi;
            Wake *wake = *wi;
                     
            // Set panel doublet coefficient:
            for (int i = 0; i < lifting_surface->n_spanwise_panels(); i++) {
                double doublet_coefficient_top    = doublet_coefficients(offset + lifting_surface->trailing_edge_upper_panel(i));
                double doublet_coefficient_bottom = doublet_coefficients(offset + lifting_surface->trailing_edge_lower_panel(i));
                
                // Use the trailing-edge Kutta condition to compute the doublet coefficients of the new wake panels.
                double doublet_coefficient = doublet_coefficient_top - doublet_coefficient_bottom;
                
                int idx = wake->n_panels() - lifting_surface->n_spanwise_panels() + i;
                wake->doublet_coefficients[idx] = doublet_coefficient;
            }
            
            // Update offset:
            offset += lifting_surface->n_panels();
        }
    }

    if (Parameters::convect_wake) {
        // Recompute source distribution without wake influence:
        cout << "Solver: Recomputing source distribution without wake influence." << endl;
        
        idx = 0;
        
        for (bi = bodies.begin(); bi != bodies.end(); bi++) {
            Body *body = *bi;
            
            vector<Surface*>::iterator si;
            for (si = body->non_lifting_surfaces.begin(); si != body->non_lifting_surfaces.end(); si++) {
                Surface *non_lifting_surface = *si;
                
                for (int i = 0; i < non_lifting_surface->n_panels(); i++) {
                    Vector3d apparent_velocity = body->panel_kinematic_velocity(*non_lifting_surface, i) - freestream_velocity;
                
                    source_coefficients(idx) = source_coefficient(*non_lifting_surface, i, apparent_velocity, false);
                    idx++;
                }
            }
            
            vector<LiftingSurface*>::iterator lsi;
            for (lsi = body->lifting_surfaces.begin(); lsi != body->lifting_surfaces.end(); lsi++) {
                LiftingSurface *lifting_surface = *lsi;
                
                for (int i = 0; i < lifting_surface->n_panels(); i++) {
                    Vector3d apparent_velocity = body->panel_kinematic_velocity(*lifting_surface, i) - freestream_velocity;
                
                    source_coefficients(idx) = source_coefficient(*lifting_surface, i, apparent_velocity, false);
                    idx++;
                }
            }
        }
    }
    
    // Compute potential values on new body with new coefficients:
    VectorXd old_surface_velocity_potentials = surface_velocity_potentials;
    surface_velocity_potentials = surface_surface_velocity_potentials();

    // Compute pressure distribution:
    cout << "Solver: Computing pressure distribution." << endl;

    idx = 0;
    
    offset = 0;

    for (bi = bodies.begin(); bi != bodies.end(); bi++) {
        Body *body = *bi;
        
        double v_ref_squared = reference_velocity_squared(*body);
        
        vector<Surface*>::iterator si;
        for (si = body->non_lifting_surfaces.begin(); si != body->non_lifting_surfaces.end(); si++) {
            Surface *non_lifting_surface = *si;
            
            VectorXd doublet_coefficient_field(non_lifting_surface->n_panels());
            for (int i = 0; i < non_lifting_surface->n_panels(); i++)
                doublet_coefficient_field(i) = doublet_coefficients(offset + i); 
                
            for (int i = 0; i < non_lifting_surface->n_panels(); i++) {
                double dphidt = velocity_potential_time_derivative(surface_velocity_potentials, old_surface_velocity_potentials, offset, i, dt);
                
                Vector3d V = surface_velocity(*non_lifting_surface, i, doublet_coefficient_field);
                surface_velocities.block<1, 3>(idx, 0) = V;
 
                pressure_coefficients(idx) = pressure_coefficient(V, dphidt, v_ref_squared);
                
                idx++;
            }   
            
            offset += non_lifting_surface->n_panels();      
        }
        
        vector<LiftingSurface*>::iterator lsi;
        for (lsi = body->lifting_surfaces.begin(); lsi != body->lifting_surfaces.end(); lsi++) {
            LiftingSurface *lifting_surface = *lsi;
            
            VectorXd doublet_coefficient_field(lifting_surface->n_panels());
            for (int i = 0; i < lifting_surface->n_panels(); i++)
                doublet_coefficient_field(i) = doublet_coefficients(offset + i); 
                
            for (int i = 0; i < lifting_surface->n_panels(); i++) {
                double dphidt = velocity_potential_time_derivative(surface_velocity_potentials, old_surface_velocity_potentials, offset, i, dt);
                
                Vector3d V = surface_velocity(*lifting_surface, i, doublet_coefficient_field);
                surface_velocities.block<1, 3>(idx, 0) = V;
 
                pressure_coefficients(idx) = pressure_coefficient(V, dphidt, v_ref_squared);
                
                idx++;
            }   
            
            offset += lifting_surface->n_panels();      
        }
    }
}

/**
   Convects existing wake nodes, and emits a new layer of wake panels.
   
   @param[in]   dt   Time step size.
*/
void
Solver::update_wakes(double dt)
{
    // Do we convect wake panels?
    if (Parameters::convect_wake) {
        // Compute velocity values at wake nodes;
        std::vector<std::vector<Vector3d> > wake_velocities;
        
        vector<Body*>::iterator bi;
        for (bi = bodies.begin(); bi != bodies.end(); bi++) {
            Body *body = *bi;
                 
            vector<Wake*>::iterator wi;
            for (wi = body->wakes.begin(); wi != body->wakes.end(); wi++) {
                Wake *wake = *wi;
                
                std::vector<Vector3d> local_wake_velocities;
                
                for (int i = 0; i < wake->n_nodes(); i++)
                    local_wake_velocities.push_back(velocity(wake->nodes[i]));
                
                wake_velocities.push_back(local_wake_velocities);
            }
        }
        
        // Add new wake panels at trailing edges, and convect all vertices:
        int idx = 0;
        
        for (bi = bodies.begin(); bi != bodies.end(); bi++) {
            Body *body = *bi;
            
            vector<LiftingSurface*>::iterator lsi;
            vector<Wake*>::iterator wi;
            for (lsi = body->lifting_surfaces.begin(), wi = body->wakes.begin(); lsi != body->lifting_surfaces.end(); lsi++, wi++) {
                LiftingSurface *lifting_surface = *lsi;
                Wake *wake = *wi;
                
                // Retrieve local wake velocities:
                std::vector<Vector3d> &local_wake_velocities = wake_velocities[idx];
                idx++;
                
                // Convect wake nodes that coincide with the trailing edge nodes with the freestream velocity.
                // Alternative options are discussed in
                //   K. Dixon, The Near Wake Structure of a Vertical Axis Wind Turbine, M.Sc. Thesis, TU Delft, 2008.
                for (int i = wake->n_nodes() - lifting_surface->n_spanwise_nodes(); i < wake->n_nodes(); i++) {
                    Vector3d apparent_velocity = body->node_kinematic_velocity(*wake, i) - freestream_velocity;
                                                               
                    wake->nodes[i] -= apparent_velocity * dt;
                }                
                
                // Convect all other wake nodes according to the local wake velocity:
                for (int i = 0; i < wake->n_nodes() - lifting_surface->n_spanwise_nodes(); i++)
                    wake->nodes[i] += local_wake_velocities[i] * dt;
                    
                // Update vortex core radii:
                for (int i = 0; i < wake->n_panels(); i++)
                    wake->update_ramasamy_leishman_vortex_core_radii(i, dt);

                // Add new vertices:
                // (This call also updates the geometry)
                wake->add_layer();
            }
        }
        
    } else {
        // No wake convection.  Re-position wake:
        vector<Body*>::iterator bi;
        for (bi = bodies.begin(); bi != bodies.end(); bi++) {
            Body *body = *bi;
            
            Vector3d body_apparent_velocity = body->velocity - freestream_velocity;
            
            vector<LiftingSurface*>::iterator lsi;
            vector<Wake*>::iterator wi;
            for (lsi = body->lifting_surfaces.begin(), wi = body->wakes.begin(); lsi != body->lifting_surfaces.end(); lsi++, wi++) {
                LiftingSurface *lifting_surface = *lsi;
                Wake *wake = *wi;
                
                for (int i = 0; i < lifting_surface->n_spanwise_nodes(); i++) {
                    // Connect wake to trailing edge nodes:                             
                    wake->nodes[lifting_surface->n_spanwise_nodes() + i] = lifting_surface->nodes[lifting_surface->trailing_edge_node(i)];
                    
                    // Point wake in direction of body kinematic velocity:
                    wake->nodes[i] = lifting_surface->nodes[lifting_surface->trailing_edge_node(i)]
                                     - Parameters::static_wake_length * body_apparent_velocity / body_apparent_velocity.norm();
                }
                
                // Need to update geometry:
                wake->compute_geometry();
            }
        }
    }
}

/**
   Returns the surface velocity potential for the given panel.
   
   @param[in]   body    Reference body.
   @param[in]   panel   Reference panel.
  
   @returns Surface velocity potential.
*/
double
Solver::surface_velocity_potential(const Surface &surface, int panel) const
{
    int offset = 0;
    
    vector<Surface*>::const_iterator si;
    for (si = non_wake_surfaces.begin(); si != non_wake_surfaces.end(); si++) {
        const Surface *tmp_surface = *si;
        
        if (&surface == tmp_surface)
            return surface_velocity_potentials(offset + panel);
        
        offset += tmp_surface->n_panels();
    }
    
    cerr << "Solver::surface_velocity_potential():  Panel " << panel << " not found on surface " << surface.id << "." << endl;
    
    return 0.0;
}

/**
   Returns the surface velocity for the given panel.
   
   @param[in]   body    Reference body.
   @param[in]   panel   Reference panel.
  
   @returns Surface velocity.
*/
Vector3d
Solver::surface_velocity(const Surface &surface, int panel) const
{
    int offset = 0;
    
    vector<Surface*>::const_iterator si;
    for (si = non_wake_surfaces.begin(); si != non_wake_surfaces.end(); si++) {
        const Surface *tmp_surface = *si;
        
        if (&surface == tmp_surface)
            return surface_velocities.block<1, 3>(offset + panel, 0);
        
        offset += tmp_surface->n_panels();
    }
    
    cerr << "Solver::surface_velocity():  Panel " << panel << " not found on surface " << surface.id << "." << endl;
    
    return Vector3d(0, 0, 0);
}

/**
   Returns the pressure coefficient for the given panel.
   
   @param[in]   body    Reference body.
   @param[in]   panel   Reference panel.
  
   @returns Pressure coefficient.
*/
double
Solver::pressure_coefficient(const Surface &surface, int panel) const
{
    int offset = 0;
    
    vector<Surface*>::const_iterator si;
    for (si = non_wake_surfaces.begin(); si != non_wake_surfaces.end(); si++) {
        const Surface *tmp_surface = *si;
        
        if (&surface == tmp_surface)
            return pressure_coefficients(offset + panel);
        
        offset += tmp_surface->n_panels();
    }
    
    cerr << "Solver::pressure_coefficient():  Panel " << panel << " not found on surface " << surface.id << "." << endl;
    
    return 0.0;
}

/**
   Computes the force caused by the pressure distribution on the given body.
   
   @param[in]   body   Reference body.
  
   @returns Force vector.
*/
Eigen::Vector3d
Solver::force(const Body &body) const
{
    // Dynamic pressure:
    double q = 0.5 * fluid_density * reference_velocity_squared(body);
        
    // Total force on body:
    Vector3d F(0, 0, 0);
    int offset = 0;
    
    vector<Surface*>::const_iterator si;
    for (si = non_wake_surfaces.begin(); si != non_wake_surfaces.end(); si++) {
        const Surface *surface = *si;
        
        for (int i = 0; i < surface->n_panels(); i++) {
            Vector3d normal = surface->panel_normal(i);
            double surface_area = surface->panel_surface_area(i);
            F += q * surface_area * pressure_coefficients(offset + i) * normal;
        }
        
        offset += surface->n_panels();
    }
    
    // Done:
    return F;      
}

/**
   Computes the moment caused by the pressure distribution on the given body, relative to the given point.
   
   @param[in]   body   Reference body.
   @param[in]   x      Reference point.
  
   @returns Moment vector.
*/
Eigen::Vector3d
Solver::moment(const Body &body, const Eigen::Vector3d &x) const
{
    // Dynamic pressure:
    double q = 0.5 * fluid_density * reference_velocity_squared(body);
        
    // Total moment on body:
    Vector3d M(0, 0, 0);
    int offset = 0;
    
    vector<Surface*>::const_iterator si;
    for (si = non_wake_surfaces.begin(); si != non_wake_surfaces.end(); si++) {
        const Surface *surface = *si;
        
        for (int i = 0; i < surface->n_panels(); i++) {                                    
            Vector3d normal = surface->panel_normal(i);
            double surface_area = surface->panel_surface_area(i);
            Vector3d F = q * surface_area * pressure_coefficients(offset + i) * normal;
            Vector3d r = surface->panel_collocation_point(i, false) - x;
            M += r.cross(F);
        }
        
        offset += surface->n_panels();
    }
    
    // Done:
    return M;
}

/**
   Logs source, doublet, and pressure coefficients into files in the logging folder.
   
   @param[in]   step_number   Step number used to name the output files.
   @param[in]   format        File format to log data in.
*/
void
Solver::log_coefficients(int step_number, SurfaceWriter &writer) const
{   
    // Log coefficients: 
    int offset = 0;
    int save_node_offset = 0;
    int save_panel_offset = 0;
    int idx;
    
    vector<Body*>::const_iterator bi;
    for (bi = bodies.begin(); bi != bodies.end(); bi++) {
        const Body *body = *bi;
        
        // Iterate non-lifting surfaces:
        idx = 0;
        vector<Surface*>::const_iterator si;
        for (si = body->non_lifting_surfaces.begin(); si != body->non_lifting_surfaces.end(); si++) {
            const Surface *non_lifting_surface = *si;
            
            // Log non-lifting surface coefficients:
            VectorXd non_lifting_surface_doublet_coefficients(non_lifting_surface->n_panels());
            VectorXd non_lifting_surface_source_coefficients(non_lifting_surface->n_panels());
            VectorXd non_lifting_surface_pressure_coefficients(non_lifting_surface->n_panels());
            for (int i = 0; i < non_lifting_surface->n_panels(); i++) {
                non_lifting_surface_doublet_coefficients(i)  = doublet_coefficients(offset + i);
                non_lifting_surface_source_coefficients(i)   = source_coefficients(offset + i);
                non_lifting_surface_pressure_coefficients(i) = pressure_coefficients(offset + i);
            }
            
            offset += non_lifting_surface->n_panels();
            
            vector<string> view_names;
            vector<VectorXd> view_data;
            
            view_names.push_back("DoubletDistribution");
            view_data.push_back(non_lifting_surface_doublet_coefficients);
            
            view_names.push_back("SourceDistribution");
            view_data.push_back(non_lifting_surface_source_coefficients);
            
            view_names.push_back("PressureDistribution");
            view_data.push_back(non_lifting_surface_pressure_coefficients);
            
            std::stringstream ss;
            ss << log_folder << "/" << body->id << "/non_lifting_surface_" << idx << "/step_" << step_number << writer.file_extension();

            writer.write(*non_lifting_surface, ss.str(), save_node_offset, save_panel_offset, view_names, view_data);
            
            save_node_offset += non_lifting_surface->n_nodes();
            save_panel_offset += non_lifting_surface->n_panels();
            
            idx++;
        }   
        
        // Iterate lifting surfaces:
        idx = 0;
        vector<LiftingSurface*>::const_iterator lsi;
        for (lsi = body->lifting_surfaces.begin(); lsi != body->lifting_surfaces.end(); lsi++) {
            const LiftingSurface *lifting_surface = *lsi;
            
            // Log lifting surface coefficients:
            VectorXd lifting_surface_doublet_coefficients(lifting_surface->n_panels());
            VectorXd lifting_surface_source_coefficients(lifting_surface->n_panels());
            VectorXd lifting_surface_pressure_coefficients(lifting_surface->n_panels());
            for (int i = 0; i < lifting_surface->n_panels(); i++) {
                lifting_surface_doublet_coefficients(i)  = doublet_coefficients(offset + i);
                lifting_surface_source_coefficients(i)   = source_coefficients(offset + i);
                lifting_surface_pressure_coefficients(i) = pressure_coefficients(offset + i);
            }
            
            offset += lifting_surface->n_panels();
            
            vector<string> view_names;
            vector<VectorXd> view_data;
            
            view_names.push_back("DoubletDistribution");
            view_data.push_back(lifting_surface_doublet_coefficients);
            
            view_names.push_back("SourceDistribution");
            view_data.push_back(lifting_surface_source_coefficients);
            
            view_names.push_back("PressureDistribution");
            view_data.push_back(lifting_surface_pressure_coefficients);
            
            std::stringstream ss;
            ss << log_folder << "/" << body->id << "/lifting_surface_" << idx << "/step_" << step_number << writer.file_extension();

            writer.write(*lifting_surface, ss.str(), save_node_offset, save_panel_offset, view_names, view_data);
            
            save_node_offset += lifting_surface->n_nodes();
            save_panel_offset += lifting_surface->n_panels();
            
            idx++;
        }
        
        // Iterate wakes:
        idx = 0;
        vector<Wake*>::const_iterator wi;
        for (wi = body->wakes.begin(); wi != body->wakes.end(); wi++) {
            const Wake *wake = *wi;
    
            // Log wake surface and coefficients:
            VectorXd wake_doublet_coefficients(wake->doublet_coefficients.size());
            for (int i = 0; i < (int) wake->doublet_coefficients.size(); i++)
                wake_doublet_coefficients(i) = wake->doublet_coefficients[i];

            vector<string> view_names;
            vector<VectorXd> view_data;
            
            view_names.push_back("DoubletDistribution");
            view_data.push_back(wake_doublet_coefficients);
            
            std::stringstream ss;
            ss << log_folder << "/" << body->id << "/wake_" << idx << "/step_" << step_number << writer.file_extension();
            
            writer.write(*wake, ss.str(), 0, save_panel_offset, view_names, view_data);
            
            save_node_offset += wake->n_nodes();
            save_panel_offset += wake->n_panels();
            
            idx++;
        }
    }
}
