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
    if (mkdir(folder.c_str(), S_IRWXU) < 0)
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
    total_n_panels_without_wakes = 0;
    
    // Properly size and zero doublet coefficient vector:
    doublet_coefficients.resize(total_n_panels_without_wakes);
    source_coefficients.resize(total_n_panels_without_wakes);
    pressure_coefficients.resize(total_n_panels_without_wakes);
    velocity_potentials.resize(total_n_panels_without_wakes);
    for (int i = 0; i < total_n_panels_without_wakes; i++) {
        doublet_coefficients(i) = 0.0;
        source_coefficients(i) = 0.0;
        pressure_coefficients(i) = 0.0;
        velocity_potentials(i) = 0.0;
    }
        
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
   Adds a mesh collection to this solver.
   
   @param[in]   collection  Collection to be added.
*/
void
Solver::add_collection(Collection &collection)
{
    collections.push_back(&collection);
    
    meshes.push_back(&collection.nolift_mesh);
    meshes_without_wakes.push_back(&collection.nolift_mesh);
    
    mesh_to_collection[&collection.nolift_mesh] = &collection;
 
    total_n_panels_without_wakes = total_n_panels_without_wakes + collection.nolift_mesh.n_panels();
    
    for (int i = 0; i < (int) collection.wings.size(); i++) {
        meshes.push_back(collection.wings[i]);
        meshes.push_back(collection.wakes[i]);
        meshes_without_wakes.push_back(collection.wings[i]);
           
        mesh_to_collection[collection.wings[i]] = &collection;
        mesh_to_collection[collection.wakes[i]] = &collection;
        
        total_n_panels_without_wakes = total_n_panels_without_wakes + collection.wings[i]->n_panels();
    }
    
    doublet_coefficients.resize(total_n_panels_without_wakes);
    source_coefficients.resize(total_n_panels_without_wakes);
    pressure_coefficients.resize(total_n_panels_without_wakes);
    velocity_potentials.resize(total_n_panels_without_wakes);
    for (int i = 0; i < total_n_panels_without_wakes; i++) {
        doublet_coefficients(i) = 0.0;
        source_coefficients(i) = 0.0;
        pressure_coefficients(i) = 0.0;
        velocity_potentials(i) = 0.0;
    }
        
    // Open logs:
    string collection_log_folder = log_folder + "/" + collection.id;
    
    mkdir_helper(collection_log_folder);
    
    stringstream ss;
    ss << collection_log_folder << "/nolift_mesh";
   
    string s = ss.str();
    mkdir_helper(s);
    
    for (int i = 0; i < (int) collection.wings.size(); i++) {
        stringstream ss;
        ss << collection_log_folder << "/wake_" << i;
        
        string s = ss.str();
        mkdir_helper(s);
        
        stringstream ss2;
        ss2 << collection_log_folder << "/wing_" << i;
        
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

// Add doublet influence of wakes to system matrix:
void
Solver::wakes_influence(MatrixXd &A, Mesh &mesh, int offset) const
{
    for (int j = 0; j < mesh.n_panels(); j++) {
        int wing_offset = 0;
        
        for (int k = 0; k < (int) collections.size(); k++) {
            Collection *collection = collections[k];
            
            wing_offset += collection->nolift_mesh.n_panels();
            
            for (int l = 0; l < (int) collection->wakes.size(); l++) {
                Wing *wing = collection->wings[l];
                Wake *wake = collection->wakes[l];
                    
                int idx = 0;
                for (int m = wake->n_panels() - (int) wing->trailing_edge_top_panels.size(); m < (int) wake->n_panels(); m++) {
                    int pa = wing->trailing_edge_top_panels[idx];
                    int pb = wing->trailing_edge_bottom_panels[idx];
                    
                    // Account for the influence of the new wake panels.  The doublet strength of these panels
                    // is set according to the Kutta condition;  see below.
                    A(offset + j, wing_offset + pa) += wake->doublet_influence(mesh, j, m);
                    A(offset + j, wing_offset + pb) -= wake->doublet_influence(mesh, j, m);
                    
                    idx++;
                }
                
                wing_offset += wing->n_panels();
            }
        }
    }
}

// Compute source coefficient for given mesh and panel:
double
Solver::source_coefficient(const Mesh &mesh, int panel, const Vector3d &kinematic_velocity, bool include_wake_influence) const
{
    // Main velocity:
    Vector3d velocity = -kinematic_velocity;
    
    // Wake contribution:
    if (Parameters::convect_wake && include_wake_influence) {
        for (int i = 0; i < (int) collections.size(); i++) {
            Collection *collection = collections[i];
            
            for (int j = 0; j < (int) collection->wakes.size(); j++) {
                Wake *wake = collection->wakes[j];
                
                for (int k = 0; k < wake->n_panels(); k++) {
                    // Use doublet panel - vortex ring equivalence.  Any new wake panels have zero doublet
                    // coefficient, and are therefore not accounted for here.
                    if (Parameters::use_ramasamy_leishman_vortex_sheet)
                        velocity += wake->vortex_ring_ramasamy_leishman_velocity
                            (mesh, panel, k, wake->vortex_core_radii[k], wake->doublet_coefficients[k]);
                    else
                        velocity += wake->vortex_ring_unit_velocity(mesh, panel, k) * wake->doublet_coefficients[k];
                }
            }
        }
    }
    
    // Take normal component:
    Vector3d normal = mesh.panel_normal(panel);
    return -velocity.dot(normal);
}

/**
   Computes the surface velocity for the given panel.
   
   @param[in]   mesh                        Reference mesh.
   @param[in]   panel                       Reference panel.
   @param[in]   doublet_coefficient_field   Doublet coefficient distribution on given mesh.
   
   @returns Surface velocity.
*/
Eigen::Vector3d
Solver::surface_velocity(const Mesh &mesh, int panel, const Eigen::VectorXd &doublet_coefficient_field) const
{
    // Compute disturbance part of surface velocity.
    Vector3d tangential_velocity;
    if (Parameters::marcov_surface_velocity) {
        Vector3d x = mesh.panel_collocation_point(panel, false);
        
        // Use N. Marcov's formula for surface velocity, see L. Dragoş, Mathematical Methods in Aerodynamics, Springer, 2003.
        Vector3d tangential_velocity = disturbance_potential_gradient(x);
        tangential_velocity -= 0.5 * mesh.scalar_field_gradient(doublet_coefficient_field, panel);
    } else
        tangential_velocity = -mesh.scalar_field_gradient(doublet_coefficient_field, panel);

    // Add flow due to kinematic velocity:
    Collection *collection = mesh_to_collection.find(&mesh)->second;
    Vector3d kinematic_velocity = mesh.panel_deformation_velocity(panel)
                                      + collection->panel_kinematic_velocity(mesh, panel)
                                      - freestream_velocity;
                                          
    tangential_velocity -= kinematic_velocity;
    
    // Remove any normal velocity.  This is the (implicit) contribution of the source term.
    Vector3d normal = mesh.panel_normal(panel);
    tangential_velocity -= tangential_velocity.dot(normal) * normal;
    
    // Done:
    return tangential_velocity;
}

/**
   Establishes the reference velocity the given collection.
   
   @param[in]   collection   Collection to establish reference velocity for.
   
   @returns Reference velocity.
*/
double
Solver::reference_velocity(const Collection &collection) const
{
    return (collection.velocity - freestream_velocity).norm();
}

/**
   Computes the pressure coefficient for the given panel.
   
   @param[in]   mesh                        Reference Mesh.
   @param[in]   panel                       Reference panel.
   @param[in]   doublet_coefficient_field   Doublet coefficient distribution on given mesh.
   @param[in]   dphidt                      Time-derivative of the velocity potential of the reference panel.
   @param[in]   v_ref                       Reference velocity for unit removal.
   
   @returns Pressure coefficient.
*/
double
Solver::pressure_coefficient(const Mesh &mesh, int panel, const Eigen::VectorXd &doublet_coefficient_field, double dphidt, double v_ref) const
{
    double C_p = 1 - (surface_velocity(mesh, panel, doublet_coefficient_field).squaredNorm() + 2 * dphidt) / pow(v_ref, 2);
    if (C_p < Parameters::min_pressure_coefficient)
        cerr << "Solver: Pressure coefficient on mesh " << mesh.id << ", panel " << panel << " is less than minimum." << endl;
    
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
    
    // Iterate all non-wake meshes:
    int offset = 0;
    for (int i = 0; i < (int) meshes_without_wakes.size(); i++) {
        Mesh *other_mesh = meshes_without_wakes[i];

        for (int j = 0; j < other_mesh->n_panels(); j++) {
            phi += other_mesh->doublet_influence(x, j) * doublet_coefficients(offset + j);
            phi += other_mesh->source_influence(x, j) * source_coefficients(offset + j);
        }
        
        offset += other_mesh->n_panels();
    }
    
    // Iterate wakes:
    for (int i = 0; i < (int) collections.size(); i++) {
        Collection *collection = collections[i];
        
        for (int j = 0; j < (int) collection->wings.size(); j++) {
            Wake *wake = collection->wakes[j];
            
            for (int k = 0; k < wake->n_panels(); k++)
                phi += wake->doublet_influence(x, k) * wake->doublet_coefficients[k];
        }
    }
                    
    // Sum with freestream velocity potential:
    return phi + freestream_velocity.dot(x);
}

/**
   Computes velocity potential values on the body surface.
   
   @returns Vector of velocity potential values, ordered by panel number.
*/
Eigen::VectorXd
Solver::surface_velocity_potentials() const
{
    cout << "Solver: Computing surface potential values." << endl;
    
    VectorXd surface_velocity_potentials(total_n_panels_without_wakes);
    
    int surface_velocity_potentials_idx = 0;
    
    // The potential functions are not singular on the panel collocation points.  We can
    // therefoce compute the surface potential values directly.  
    for (int i = 0; i < (int) collections.size(); i++) {
        Collection *collection = collections[i];
        
        for (int j = 0; j < collection->nolift_mesh.n_panels(); j++) {
            surface_velocity_potentials(surface_velocity_potentials_idx) = velocity_potential(collection->nolift_mesh.panel_collocation_point(j, false));
            surface_velocity_potentials_idx++;
        }
        
        for (int j = 0; j < (int) collection->wings.size(); j++) {
            Wing *wing = collection->wings[j];
            
            for (int k = 0; k < wing->n_panels(); k++) {
                surface_velocity_potentials(surface_velocity_potentials_idx) = velocity_potential(wing->panel_collocation_point(k, false));
                surface_velocity_potentials_idx++;
            }
        }
    }
    
    return surface_velocity_potentials;
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
    
    // Iterate all non-wake meshes:
    int offset = 0;
    for (int i = 0; i < (int) meshes_without_wakes.size(); i++) {
        Mesh *other_mesh = meshes_without_wakes[i];

        for (int j = 0; j < other_mesh->n_panels(); j++) {
            gradient += other_mesh->vortex_ring_unit_velocity(x, j) * doublet_coefficients(offset + j);
            gradient += other_mesh->source_unit_velocity(x, j) * source_coefficients(offset + j);
        }
        
        offset += other_mesh->n_panels();
    }
    
    // Iterate wakes:
    for (int i = 0; i < (int) collections.size(); i++) {
        Collection *collection = collections[i];
        
        for (int j = 0; j < (int) collection->wings.size(); j++) {
            Wake *wake = collection->wakes[j];
            Wing *wing = collection->wings[j];
            
            if (wake->n_panels() >= (int) wing->trailing_edge_top_panels.size()) {
                for (int k = wake->n_panels() - (int) wing->trailing_edge_top_panels.size(); k < wake->n_panels(); k++)
                    gradient += wake->vortex_ring_unit_velocity(x, k) * wake->doublet_coefficients[k];

                for (int k = 0; k < wake->n_panels() - (int) wing->trailing_edge_top_panels.size(); k++) {
                    if (Parameters::use_ramasamy_leishman_vortex_sheet)
                        gradient += wake->vortex_ring_ramasamy_leishman_velocity(x, k, wake->vortex_core_radii[k], wake->doublet_coefficients[k]);
                    else
                        gradient += wake->vortex_ring_unit_velocity(x, k) * wake->doublet_coefficients[k];
                }
            }
        }
    }
               
    // Done:
    return gradient;
}

/**
   Computes velocity potential time derivative at the given panel.
   
   @param[in]  velocity_potentials         Current potential values.
   @param[in]  old_velocity_potentials     Previous potential values
   @param[in]  offset                      Offset to requested Mesh
   @param[in]  panel                       Panel number.
   @param[in]  dt                          Time step.
   
   @returns Velocity potential time derivative.
*/ 
double
Solver::velocity_potential_time_derivative(const Eigen::VectorXd &velocity_potentials, const Eigen::VectorXd &old_velocity_potentials, int offset, int panel, double dt) const
{
    double dphidt;
    
    // Evaluate the time-derivative of the potential in a body-fixed reference frame, as in
    //   J. P. Giesing, Nonlinear Two-Dimensional Unsteady Potential Flow with Lift, Journal of Aircraft, 1968.
    if (Parameters::unsteady_bernoulli && dt > 0.0)
        dphidt = (velocity_potentials(offset + panel) - old_velocity_potentials(offset + panel)) / dt;
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
    // Find closest mesh and panel:
    double distance = numeric_limits<double>::max();
    Mesh *close_mesh = NULL;
    int close_panel = -1;
    int close_offset = - 1;
    bool close_near_sharp_edge = false;
    int offset = 0;
    
    for (int i = 0; i < (int) meshes_without_wakes.size(); i++) {
        Mesh *mesh_candidate = meshes_without_wakes[i];
        
        int panel_candidate;
        double distance_candidate;
        bool near_sharp_edge = mesh_candidate->closest_panel(x, panel_candidate, distance_candidate);
        
        if (distance_candidate < distance) {
            distance              = distance_candidate;
            close_mesh            = mesh_candidate;
            close_panel           = panel_candidate;
            close_offset          = offset;
            close_near_sharp_edge = near_sharp_edge;
        }
        
        offset += meshes_without_wakes[i]->n_panels();
    }
 
    // Compute velocity potential gradients near the body:
    vector<Vector3d> disturbance_potential_gradients;
    vector<Vector3d> close_to_body_points;
    if (distance < Parameters::interpolation_layer_thickness && !close_near_sharp_edge) {   
        for (int i = 0; i < (int) close_mesh->panel_nodes[close_panel].size(); i++) {
            Vector3d close_to_body_point = close_mesh->close_to_body_point(close_mesh->panel_nodes[close_panel][i]);
            close_to_body_points.push_back(close_to_body_point);
            
            disturbance_potential_gradients.push_back(disturbance_potential_gradient(close_to_body_point));
        }
        
    } else {
        disturbance_potential_gradients.push_back(disturbance_potential_gradient(x));
        
        close_mesh  = NULL;
        close_panel = -1;
    }
    
    Vector3d velocity(0, 0, 0);
    if (close_mesh != NULL) { 
        // Compute stream velocity near the boundary using interpolation, see
        //   K. Dixon, C. S. Ferreira, C. Hofemann, G. van Brussel, G. van Kuik,
        //   A 3D Unsteady Panel Method for Vertical Axis Wind Turbines, DUWIND, 2008.
        Vector3d x = close_mesh->panel_collocation_point(close_panel, false);
        
        VectorXd doublet_coefficient_field(close_mesh->n_panels());
        for (int i = 0; i < close_mesh->n_panels(); i++)
            doublet_coefficient_field(i) = doublet_coefficients(close_offset + i);
        
        Vector3d normal = close_mesh->panel_normal(close_panel);
        
        double total_weight = 0.0;
        
        for (int i = 0; i < (int) disturbance_potential_gradients.size(); i++) {
            Vector3d layer_point_distance = x - close_to_body_points[i];
            layer_point_distance = layer_point_distance - layer_point_distance.dot(normal) * normal;
            
            double weight = distance * (close_mesh->panel_diameter(close_panel) - layer_point_distance.norm());
                 
            velocity += weight * (disturbance_potential_gradients[i] + freestream_velocity);
            total_weight += weight;
        }
        
        double weight = Parameters::interpolation_layer_thickness - distance;
        velocity += weight * surface_velocity(*close_mesh, close_panel, doublet_coefficient_field);
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
   
   @param[in]   dt  Time step size.
*/
void
Solver::initialize_wakes(double dt)
{
    // Add initial wake layers:
    for (int i = 0; i < (int) collections.size(); i++) {
        Collection *collection = collections[i];
        
        Vector3d collection_kinematic_velocity = collection->velocity - freestream_velocity;
        
        for (int j = 0; j < (int) collection->wings.size(); j++) {
            Wake *wake = collection->wakes[j];
            
            wake->add_layer();
            for (int k = 0; k < wake->n_nodes(); k++) {                                  
                if (Parameters::convect_wake)    
                    wake->nodes[k] -= collection_kinematic_velocity * dt;
                else
                    wake->nodes[k] -= Parameters::static_wake_length * collection_kinematic_velocity / collection_kinematic_velocity.norm();
            }
            
            wake->add_layer();
        }
    }
}

// Create block of linear system:
void
Solver::doublet_coefficient_matrix_block(MatrixXd &doublet_influence_coefficients,
                                         MatrixXd &source_influence_coefficients,
                                         const Mesh &mesh_one, int offset_one, const Mesh &mesh_two, int offset_two) const
{
    for (int i = 0; i < mesh_one.n_panels(); i++) {
        for (int j = 0; j < mesh_two.n_panels(); j++) {
            if ((&mesh_one == &mesh_two) && (i == j))
                doublet_influence_coefficients(offset_one + i, offset_two + j) = -0.5;
            else
                doublet_influence_coefficients(offset_one + i, offset_two + j) = mesh_two.doublet_influence(mesh_one, i, j);
            
            source_influence_coefficients(offset_one + i, offset_two + j) = mesh_two.source_influence(mesh_one, i, j);
        }
    }
}

/**
   Computes new source, doublet, and pressure distributions.
   
   @param[in]   dt  Time step size.
*/
void
Solver::update_coefficients(double dt)
{
    // Compute source distribution:
    cout << "Solver: Computing source distribution with wake influence." << endl;
    
    int source_coefficients_idx = 0;
    
    for (int i = 0; i < (int) collections.size(); i++) {
        Collection *collection = collections[i];
        
        for (int j = 0; j < collection->nolift_mesh.n_panels(); j++) {
            Vector3d kinematic_velocity = collection->nolift_mesh.panel_deformation_velocity(j)
                                          + collection->panel_kinematic_velocity(collection->nolift_mesh, j)
                                          - freestream_velocity;
            
            source_coefficients(source_coefficients_idx) = source_coefficient(collection->nolift_mesh, j, kinematic_velocity, true);
            source_coefficients_idx++;
        }
        
        for (int j = 0; j < (int) collection->wings.size(); j++) {
            Wing *wing = collection->wings[j];
            
            for (int k = 0; k < wing->n_panels(); k++) {
                Vector3d kinematic_velocity = wing->panel_deformation_velocity(k)
                                              + collection->panel_kinematic_velocity(*wing, k) 
                                              - freestream_velocity;
            
                source_coefficients(source_coefficients_idx) = source_coefficient(*wing, k, kinematic_velocity, true);
                source_coefficients_idx++;
            }
        }
    }
  
    // Compute doublet distribution:
    MatrixXd A(total_n_panels_without_wakes, total_n_panels_without_wakes);
    MatrixXd source_influence_coefficients(total_n_panels_without_wakes, total_n_panels_without_wakes);
    
    int offset_one = 0, offset_two = 0;
    for (int i = 0; i < (int) meshes_without_wakes.size(); i++) {
        offset_two = 0;
        for (int j = 0; j < (int) meshes_without_wakes.size(); j++) {
            doublet_coefficient_matrix_block(A,
                                             source_influence_coefficients,
                                             *meshes_without_wakes[i], offset_one, *meshes_without_wakes[j], offset_two);
            
            offset_two = offset_two + meshes_without_wakes[j]->n_panels();
        }
        
        wakes_influence(A, *meshes_without_wakes[i], offset_one);
        
        offset_one = offset_one + meshes_without_wakes[i]->n_panels();
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
    for (int i = 0; i < (int) collections.size(); i++) {
        Collection *collection = collections[i];
        
        offset += collection->nolift_mesh.n_panels();
        
        for (int j = 0; j < (int) collection->wings.size(); j++) {
            Wing *wing = collection->wings[j];
            Wake *wake = collection->wakes[j];
                     
            // Set panel doublet coefficient:
            int trailing_edge_n_nodes = wing->trailing_edge_nodes.size();
            int trailing_edge_n_panels = trailing_edge_n_nodes - 1;
            for (int k = 0; k < trailing_edge_n_panels; k++) {
                double doublet_coefficient_top    = doublet_coefficients(offset + wing->trailing_edge_top_panels[k]);
                double doublet_coefficient_bottom = doublet_coefficients(offset + wing->trailing_edge_bottom_panels[k]);
                
                // Use the trailing-edge Kutta condition to compute the doublet coefficients of the new wake panels.
                double doublet_coefficient = doublet_coefficient_top - doublet_coefficient_bottom;
                
                int idx = wake->n_panels() - trailing_edge_n_panels + k;
                wake->doublet_coefficients[idx] = doublet_coefficient;
            }
            
            // Update offset:
            offset += wing->n_panels();
        }
    }

    if (Parameters::convect_wake) {
        // Recompute source distribution without wake influence:
        cout << "Solver: Recomputing source distribution without wake influence." << endl;
        
        source_coefficients_idx = 0;
        
        for (int i = 0; i < (int) collections.size(); i++) {
            Collection *collection = collections[i];
            
            for (int j = 0; j < collection->nolift_mesh.n_panels(); j++) {
                Vector3d kinematic_velocity = collection->nolift_mesh.panel_deformation_velocity(j)
                                              + collection->panel_kinematic_velocity(collection->nolift_mesh, j) 
                                              - freestream_velocity;
                
                source_coefficients(source_coefficients_idx) = source_coefficient(collection->nolift_mesh, j, kinematic_velocity, false);
                source_coefficients_idx++;
            }
            
            for (int j = 0; j < (int) collection->wings.size(); j++) {
                Wing *wing = collection->wings[j];
                
                for (int k = 0; k < wing->n_panels(); k++) {
                    Vector3d kinematic_velocity = wing->panel_deformation_velocity(k)
                                                  + collection->panel_kinematic_velocity(*wing, k) 
                                                  - freestream_velocity;
                
                    source_coefficients(source_coefficients_idx) = source_coefficient(*wing, k, kinematic_velocity, false);
                    source_coefficients_idx++;
                }
            }
        }
    }
    
    // Compute potential values on new body with new coefficients:
    VectorXd old_velocity_potentials = velocity_potentials;
    if (Parameters::unsteady_bernoulli)
        velocity_potentials = surface_velocity_potentials();

    // Compute pressure distribution:
    cout << "Solver: Computing pressure distribution." << endl;

    int pressure_coefficients_idx = 0;
    
    offset = 0;

    for (int i = 0; i < (int) collections.size(); i++) {
        const Collection *collection = collections[i];
        
        double v_ref = reference_velocity(*collection);
        
        VectorXd doublet_coefficient_field(collection->nolift_mesh.n_panels());
        for (int j = 0; j < collection->nolift_mesh.n_panels(); j++)
            doublet_coefficient_field(j) = doublet_coefficients(offset + j); 
            
        for (int j = 0; j < collection->nolift_mesh.n_panels(); j++) {                          
            double dphidt = velocity_potential_time_derivative(velocity_potentials, old_velocity_potentials, offset, j, dt);
                                                 
            pressure_coefficients(pressure_coefficients_idx) = pressure_coefficient(collection->nolift_mesh, j,
                                                                                    doublet_coefficient_field,
                                                                                    dphidt, v_ref);
            pressure_coefficients_idx++;
        }
        
        offset += collection->nolift_mesh.n_panels();
        
        for (int j = 0; j < (int) collection->wings.size(); j++) {
            Wing *wing = collection->wings[j];
            
            VectorXd doublet_coefficient_field(wing->n_panels());
            for (int k = 0; k < wing->n_panels(); k++)
                doublet_coefficient_field(k) = doublet_coefficients(offset + k); 
                
            Vector3d wing_force(0, 0, 0);
                
            for (int k = 0; k < wing->n_panels(); k++) {
                double dphidt = velocity_potential_time_derivative(velocity_potentials, old_velocity_potentials, offset, k, dt);
 
                pressure_coefficients(pressure_coefficients_idx) = pressure_coefficient(*wing, k,
                                                                                        doublet_coefficient_field,
                                                                                        dphidt, v_ref);
                
                pressure_coefficients_idx++;
            }   
            
            offset += wing->n_panels();      
        }
    }
}

/**
   Convects existing wake nodes, and emits a new layer of wake panels.
   
   @param[in]   dt  Time step size.
*/
void
Solver::update_wakes(double dt)
{
    // Do we convect wake panels?
    if (Parameters::convect_wake) {
        // Compute velocity values at wake nodes;
        std::vector<std::vector<Vector3d> > wake_velocities;
        
        for (int i = 0; i < (int) collections.size(); i++) {
            Collection *collection = collections[i];
                 
            for (int j = 0; j < (int) collection->wings.size(); j++) {
                Wake *wake = collection->wakes[j];
                
                std::vector<Vector3d> local_wake_velocities;
                
                for (int k = 0; k < wake->n_nodes(); k++)
                    local_wake_velocities.push_back(velocity(wake->nodes[k]));
                
                wake_velocities.push_back(local_wake_velocities);
            }
        }
        
        // Add new wake panels at trailing edges, and convect all vertices:
        int idx = 0;
        
        for (int i = 0; i < (int) collections.size(); i++) {
            Collection *collection = collections[i];
            
            for (int j = 0; j < (int) collection->wings.size(); j++) {
                Wing *wing = collection->wings[j];
                Wake *wake = collection->wakes[j];
                
                // Retrieve local wing velocities:
                std::vector<Vector3d> &local_wake_velocities = wake_velocities[idx];
                idx++;
                
                // Convect wake nodes that coincide with the trailing edge nodes with the freestream velocity.
                // Alternative options are discussed in
                //   K. Dixon, The Near Wake Structure of a Vertical Axis Wind Turbine, M.Sc. Thesis, TU Delft, 2008.
                for (int k = wake->n_nodes() - (int) wing->trailing_edge_nodes.size(); k < wake->n_nodes(); k++) {
                     Vector3d kinematic_velocity = wake->node_deformation_velocities[k]
                                                   + collection->node_kinematic_velocity(*wake, k)
                                                   - freestream_velocity;
                                                               
                     wake->nodes[k] -= kinematic_velocity * dt;
                }                
                
                // Convect all other wake nodes according to the local wake velocity:
                for (int k = 0; k < wake->n_nodes() - (int) wing->trailing_edge_nodes.size(); k++)         
                    wake->nodes[k] += local_wake_velocities[k] * dt;
                    
                // Update vortex core radii:
                for (int k = 0; k < wake->n_panels(); k++)
                    wake->update_ramasamy_leishman_vortex_core_radii(k, dt);

                // Add new vertices:
                // (This call also updates the geometry)
                wake->add_layer();
            }
        }
        
    } else {
        // No wake convection.  Re-position wake:
        for (int i = 0; i < (int) collections.size(); i++) {
            Collection *collection = collections[i];
            
            Vector3d collection_kinematic_velocity = collection->velocity - freestream_velocity;
            
            for (int j = 0; j < (int) collection->wings.size(); j++) {
                Wing *wing = collection->wings[j];
                Wake *wake = collection->wakes[j];
                
                for (int k = 0; k < (int) wing->trailing_edge_nodes.size(); k++) {
                    // Connect wake to trailing edge nodes:                             
                    wake->nodes[wing->trailing_edge_nodes.size() + k] = wing->nodes[wing->trailing_edge_nodes[k]];
                    
                    // Point wake in direction of collection kinematic velocity:
                    wake->nodes[k] = wing->nodes[wing->trailing_edge_nodes[k]]
                                     - Parameters::static_wake_length * collection_kinematic_velocity / collection_kinematic_velocity.norm();
                }
                
                // Need to update geometry:
                wake->compute_geometry();
            }
        }
    }
}

/**
   Computes the aerodynamic force on a given collection.
   
   @param[in]   collection  Reference collection.
  
   @returns Aerodynamic force.
*/
double
Solver::pressure_coefficient(const Mesh &mesh, int panel) const
{
    int offset = 0;
    
    for (int i = 0; i < (int) meshes_without_wakes.size(); i++) {
        Mesh *tmp_mesh = meshes_without_wakes[i];
        
        if (&mesh == tmp_mesh)
            return pressure_coefficients(offset + panel);
        
        offset += tmp_mesh->n_panels();
    }
    
    cerr << "Solver::pressure_coefficient():  Panel " << panel << " not found on mesh " << mesh.id << "." << endl;
    
    return 0.0;
}

/**
   Computes the aerodynamic force on a given collection.
   
   @param[in]   collection  Reference collection.
  
   @returns Aerodynamic force.
*/
Eigen::Vector3d
Solver::aerodynamic_force(const Collection &collection) const
{
    double v_ref = reference_velocity(collection);
        
    Vector3d F(0, 0, 0);
    int offset = 0;
    
    for (int i = 0; i < (int) meshes_without_wakes.size(); i++) {
        Mesh *mesh = meshes_without_wakes[i];
        
        for (int k = 0; k < mesh->n_panels(); k++) {                                    
            Vector3d normal = mesh->panel_normal(k);            
            double surface_area = mesh->panel_surface_area(k);
            F += 0.5 * fluid_density * pow(v_ref, 2) * surface_area * pressure_coefficients(offset + k) * normal;
        }
        
        offset += mesh->n_panels();
    }
    
    return F;      
}

/**
   Computes the aerodynamic moment on a given point.
   
   @param[in]   collection  Reference collection.
   @param[in]   x           Reference point.
  
   @returns Aerodynamic moment.
*/
Eigen::Vector3d
Solver::aerodynamic_moment(const Collection &collection, const Eigen::Vector3d &x) const
{
    double v_ref = reference_velocity(collection);
        
    Vector3d M(0, 0, 0);
    int offset = 0;
    
    for (int i = 0; i < (int) meshes_without_wakes.size(); i++) {
        Mesh *mesh = meshes_without_wakes[i];
        
        for (int k = 0; k < mesh->n_panels(); k++) {                                    
            Vector3d normal = mesh->panel_normal(k);            
            double surface_area = mesh->panel_surface_area(k);
            Vector3d F = 0.5 * fluid_density * pow(v_ref, 2) * surface_area * pressure_coefficients(offset + k) * normal;
            Vector3d r = mesh->panel_collocation_point(k, false) - x;
            M += r.cross(F);
        }
        
        offset += mesh->n_panels();
    }
    
    return M;
}

/**
   Logs source, doublet, and pressure coefficients into files into the logging folder.
   
   @param[in]   step_number     Step number used to name the output files.
*/
void
Solver::log_coefficients(int step_number) const
{   
    // Log coefficients: 
    int offset = 0;
    int save_node_offset = 0;
    int save_panel_offset = 0;
    
    for (int i = 0; i < (int) collections.size(); i++) {
        Collection *collection = collections[i];
        
        // Log no-lift mesh:
        VectorXd nolift_mesh_doublet_coefficients(collection->nolift_mesh.n_panels());
        VectorXd nolift_mesh_source_coefficients(collection->nolift_mesh.n_panels());
        VectorXd nolift_mesh_pressure_coefficients(collection->nolift_mesh.n_panels());
        for (int k = 0; k < collection->nolift_mesh.n_panels(); k++) {
            nolift_mesh_doublet_coefficients(k)  = doublet_coefficients(offset + k);
            nolift_mesh_source_coefficients(k)   = source_coefficients(offset + k);
            nolift_mesh_pressure_coefficients(k) = pressure_coefficients(offset + k);
        }
        
        offset += collection->nolift_mesh.n_panels();
        
        vector<string> view_names;
        vector<VectorXd> view_data;
        
        view_names.push_back("Doublet distribution");
        view_data.push_back(nolift_mesh_doublet_coefficients);
        
        view_names.push_back("Source distribution");
        view_data.push_back(nolift_mesh_source_coefficients);
        
        view_names.push_back("Pressure distribution");
        view_data.push_back(nolift_mesh_pressure_coefficients);
        
        std::stringstream ss;
        ss << log_folder << "/" << collection->id << "/nolift_mesh/step_" << step_number << ".msh";

        collection->nolift_mesh.save(ss.str(), view_names, view_data, save_node_offset, save_panel_offset);
        save_node_offset += collection->nolift_mesh.n_nodes();
        save_panel_offset += collection->nolift_mesh.n_panels();
        
        // Iterate wings:
        for (int j = 0; j < (int) collection->wings.size(); j++) {
            Wing *wing = collection->wings[j];
            Wake *wake = collection->wakes[j];
            
            // Log wing coefficients:
            VectorXd wing_doublet_coefficients(wing->n_panels());
            VectorXd wing_source_coefficients(wing->n_panels());
            VectorXd wing_pressure_coefficients(wing->n_panels());
            for (int k = 0; k < wing->n_panels(); k++) {
                wing_doublet_coefficients(k)  = doublet_coefficients(offset + k);
                wing_source_coefficients(k)   = source_coefficients(offset + k);
                wing_pressure_coefficients(k) = pressure_coefficients(offset + k);
            }
            
            offset += wing->n_panels();
            
            vector<string> view_names;
            vector<VectorXd> view_data;
            
            view_names.push_back("Doublet distribution");
            view_data.push_back(wing_doublet_coefficients);
            
            view_names.push_back("Source distribution");
            view_data.push_back(wing_source_coefficients);
            
            view_names.push_back("Pressure distribution");
            view_data.push_back(wing_pressure_coefficients);
            
            std::stringstream ss;
            ss << log_folder << "/" << collection->id << "/wing_" << j << "/step_" << step_number << ".msh";

            wing->save(ss.str(), view_names, view_data, save_node_offset, save_panel_offset);
            save_node_offset += wing->n_nodes();
            save_panel_offset += wing->n_panels();
            
            // Log wake mesh and coefficients:
            VectorXd wake_doublet_coefficients(wake->doublet_coefficients.size());
            for (int k = 0; k < (int) wake->doublet_coefficients.size(); k++)
                wake_doublet_coefficients(k) = wake->doublet_coefficients[k];

            view_names.clear();
            view_data.clear();
            
            view_names.push_back("Doublet distribution");
            view_data.push_back(wake_doublet_coefficients);
            
            std::stringstream ss2;
            ss2 << log_folder << "/" << collection->id << "/wake_" << j << "/step_" << step_number << ".msh";
            wake->save(ss2.str(), view_names, view_data, 0, save_panel_offset);
            save_node_offset += wake->n_nodes();
            save_panel_offset += wake->n_panels();
        }
    }
}
