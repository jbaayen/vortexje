//
// Vortexje -- Solver.
//
// Copyright (C) 2012 Baayen & Heinz GmbH.
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
mkdir_helper(string folder)
{
    if (mkdir(folder.c_str(), S_IRWXU) < 0)
        if (errno != EEXIST)
            cerr << "Could not create log folder " << folder << ": " << strerror(errno) << endl;
}

/**
   Construct a solver, logging its output into the given folder.
   
   @param[in]   log_folder  Logging output folder.
*/
Solver::Solver(string log_folder) : log_folder(log_folder)
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
    potentials.resize(total_n_panels_without_wakes);
    for (int i = 0; i < total_n_panels_without_wakes; i++) {
        doublet_coefficients(i) = 0.0;
        source_coefficients(i) = 0.0;
        pressure_coefficients(i) = 0.0;
        potentials(i) = 0.0;
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
    total_n_panels_without_wakes = total_n_panels_without_wakes + collection.nolift_mesh.n_panels();
    
    for (int i = 0; i < (int) collection.wings.size(); i++) {
        meshes.push_back(collection.wings[i]);
        meshes.push_back(collection.wakes[i]);
        meshes_without_wakes.push_back(collection.wings[i]);      
        total_n_panels_without_wakes = total_n_panels_without_wakes + collection.wings[i]->n_panels();
    }
    
    doublet_coefficients.resize(total_n_panels_without_wakes);
    source_coefficients.resize(total_n_panels_without_wakes);
    pressure_coefficients.resize(total_n_panels_without_wakes);
    potentials.resize(total_n_panels_without_wakes);
    for (int i = 0; i < total_n_panels_without_wakes; i++) {
        doublet_coefficients(i) = 0.0;
        source_coefficients(i) = 0.0;
        pressure_coefficients(i) = 0.0;
        potentials(i) = 0.0;
    }
        
    // Open logs:
    string collection_log_folder = log_folder + "/" + collection.id;
    
    mkdir_helper(collection_log_folder);
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

// Solver stepping:
void
Solver::doublet_coefficient_matrix_block(MatrixXd &A, VectorXd &b, 
                                         Mesh &mesh_one, int offset_one, Mesh &mesh_two, int offset_two)
{
    for (int i = 0; i < mesh_one.n_panels(); i++) {
        for (int j = 0; j < mesh_two.n_panels(); j++) {
            A(offset_one + i, offset_two + j) = -mesh_two.doublet_influence(mesh_one, i, j);
            
            b(offset_one + i) = b(offset_one + i) - mesh_two.source_influence(mesh_one, i, j) * source_coefficients(offset_two + j);         
        }
    }
}

// Add doublet influence of wakes to system matrix:
void
Solver::wakes_influence(MatrixXd &A, VectorXd &b, Mesh &mesh, int offset)
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
                    
                    A(offset + j, wing_offset + pa) -= wake->doublet_influence(mesh, j, m);
                    A(offset + j, wing_offset + pb) += wake->doublet_influence(mesh, j, m);
                    
                    idx++;
                }
                
                wing_offset += wing->n_panels();
            }
        }
    }
}

// Compute source coefficient for given mesh and panel:
double
Solver::source_coefficient(Mesh &mesh, int panel, const Vector3d &kinematic_velocity, bool include_wake_influence)
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
                    // Use doublet panel - vortex ring equivalence:
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
    return velocity.dot(normal);
}

/**
   Computes the surface velocity for the given panel.
   
   @param[in]   mesh                        Reference mesh.
   @param[in]   panel                       Reference panel.
   @param[in]   doublet_coefficient_field   Doublet coefficient distribution on given mesh.
   @param[in]   kinematic_velocity          Reference kinematic velocity.
   
   @returns Surface velocity.
*/
Eigen::Vector3d
Solver::surface_velocity(Mesh &mesh, int panel, const Eigen::VectorXd &doublet_coefficient_field, const Eigen::Vector3d &kinematic_velocity)
{
    Vector3d x = mesh.panel_collocation_point(panel, false);
    
    // Use N. Marcov's formula for surface velocity, see L. Dragoş, Mathematical Methods in Aerodynamics, Springer, 2003.
    Vector3d tangential_velocity = potential_gradient(x);
    tangential_velocity -= 0.5 * mesh.scalar_field_gradient(doublet_coefficient_field, panel);

    // Add flow due to kinematic velocity:
    tangential_velocity -= kinematic_velocity;
    
    // Remove any normal velocity:
    Vector3d normal = mesh.panel_normal(panel);
    tangential_velocity -= tangential_velocity.dot(normal) * normal;
    
    // Done:
    return tangential_velocity;
}

/**
   Computes the pressure coefficient for the given panel.
   
   @param[in]   mesh                        Reference Mesh.
   @param[in]   panel                       Reference panel.
   @param[in]   kinematic_velocity          Reference kinematic velocity.
   @param[in]   doublet_coefficient_field   Doublet coefficient distribution on given mesh.
   @param[in]   dpotentialdt                Time-derivative of the disturbance potential of the reference panel.
   @param[in]   v_ref                       Reference velocity for unit removal.
   
   @returns Pressure coefficient.
*/
double
Solver::pressure_coefficient(Mesh &mesh, int panel, const Eigen::Vector3d &kinematic_velocity,
                             const Eigen::VectorXd &doublet_coefficient_field, double dpotentialdt, double v_ref)
{
    double C_p = 1 - (surface_velocity(mesh, panel, doublet_coefficient_field, kinematic_velocity).squaredNorm() + 2 * dpotentialdt) / pow(v_ref, 2);
    if (C_p < Parameters::min_pressure_coefficient)
        cerr << "Solver: Pressure coefficient on mesh " << mesh.id << ", panel " << panel << " is less than minimum." << endl;
    
    return C_p;
}

/**
   Computes disturbance potential at the given point.
   
   @param[in]   x   Reference point.
   
   @returns Disturbance potential.
*/
double
Solver::potential(const Vector3d &x)
{
    double potential = 0.0;
    
    // Iterate all non-wake meshes:
    int offset = 0;
    for (int i = 0; i < (int) meshes_without_wakes.size(); i++) {
        Mesh *other_mesh = meshes_without_wakes[i];

        for (int j = 0; j < other_mesh->n_panels(); j++) {
            potential += other_mesh->doublet_influence(x, j) * doublet_coefficients(offset + j);
            potential -= other_mesh->source_influence(x, j) * source_coefficients(offset + j);
        }
        
        offset += other_mesh->n_panels();
    }
    
    // Iterate wakes:
    for (int i = 0; i < (int) collections.size(); i++) {
        Collection *collection = collections[i];
        
        for (int j = 0; j < (int) collection->wings.size(); j++) {
            Wake *wake = collection->wakes[j];
            
            for (int k = 0; k < wake->n_panels(); k++)
                potential += wake->doublet_influence(x, k) * wake->doublet_coefficients[k];
        }
    }
                    
    return potential;
}

/**
   Computes disturbance potential values on the body surface.
   
   @returns Vector of potential values, ordered by panel number.
*/
Eigen::VectorXd
Solver::surface_potentials()
{
    cout << "Solver: Computing surface potential values." << endl;
    
    VectorXd surface_potentials(total_n_panels_without_wakes);
    
    int surface_potentials_idx = 0;
    
    for (int i = 0; i < (int) collections.size(); i++) {
        Collection *collection = collections[i];
        
        for (int j = 0; j < collection->nolift_mesh.n_panels(); j++) {
            surface_potentials(surface_potentials_idx) = -doublet_coefficients(surface_potentials_idx);
            surface_potentials_idx++;
        }
        
        for (int j = 0; j < (int) collection->wings.size(); j++) {
            Wing *wing = collection->wings[j];
            
            for (int k = 0; k < wing->n_panels(); k++) {
                surface_potentials(surface_potentials_idx) = -doublet_coefficients(surface_potentials_idx);
                surface_potentials_idx++;
            }
        }
    }
    
    return surface_potentials;
}

/**
   Computes disturbance potential gradient at the given point.
   
   @param[in]  x    Reference point.
   
   @returns Disturbance potential gradient.
*/ 
Eigen::Vector3d
Solver::potential_gradient(const Eigen::Vector3d &x)
{
    Vector3d gradient(0, 0, 0);
    
    // Iterate all non-wake meshes:
    int offset = 0;
    for (int i = 0; i < (int) meshes_without_wakes.size(); i++) {
        Mesh *other_mesh = meshes_without_wakes[i];

        for (int j = 0; j < other_mesh->n_panels(); j++) {
            gradient += other_mesh->vortex_ring_unit_velocity(x, j) * doublet_coefficients(offset + j);
            gradient -= other_mesh->source_unit_velocity(x, j) * source_coefficients(offset + j);
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
                    
    return gradient;
}

/**
   Compute stream velocity at a given point.
   
   @param[in]   x                   Reference point.
   @param[in]   kinematic_velocity  Reference kinematic velocity.
   
   @returns Stream velocity.
*/
Eigen::Vector3d
Solver::stream_velocity(const Eigen::Vector3d &x, const Eigen::Vector3d &kinematic_velocity)
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
 
    // Compute disturbance potential gradients near the body:
    vector<Vector3d> potential_gradients;
    vector<Vector3d> close_to_body_points;
    if (distance < Parameters::interpolation_layer_thickness && !close_near_sharp_edge) {   
        for (int i = 0; i < (int) close_mesh->panel_nodes[close_panel].size(); i++) {
            Vector3d close_to_body_point = close_mesh->close_to_body_point(close_mesh->panel_nodes[close_panel][i]);
            close_to_body_points.push_back(close_to_body_point);
            
            potential_gradients.push_back(potential_gradient(close_to_body_point));
        }
        
    } else {
        potential_gradients.push_back(potential_gradient(x));
        
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
        
        for (int i = 0; i < (int) potential_gradients.size(); i++) {
            Vector3d layer_point_distance = x - close_to_body_points[i];
            layer_point_distance = layer_point_distance - layer_point_distance.dot(normal) * normal;
            
            double weight = distance * (close_mesh->panel_diameter(close_panel) - layer_point_distance.norm());
                 
            velocity += weight * (potential_gradients[i] - kinematic_velocity);
            total_weight += weight;
        }
        
        double weight = Parameters::interpolation_layer_thickness - distance;
        velocity += weight * surface_velocity(*close_mesh, close_panel, doublet_coefficient_field, kinematic_velocity);
        total_weight += weight;
        
        velocity /= total_weight;
        
    } else {
        // We are not close to the boundary.  Use sum of disturbance potential and flow due to kinematic velocity.
        velocity = potential_gradients[0] - kinematic_velocity;
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
            
            wake->add_layer(meshes);
            for (int k = 0; k < wake->n_nodes(); k++) {                                  
                if (Parameters::convect_wake)    
                    wake->nodes[k] -= collection_kinematic_velocity * dt;
                else
                    wake->nodes[k] -= Parameters::static_wake_length * collection_kinematic_velocity / collection_kinematic_velocity.norm();
            }
            
            wake->add_layer(meshes);
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
    VectorXd b(total_n_panels_without_wakes);
    for (int i = 0; i < total_n_panels_without_wakes; i++)
        b(i) = 0.0;
    
    int offset_one = 0, offset_two = 0;
    for (int i = 0; i < (int) meshes_without_wakes.size(); i++) {
        offset_two = 0;
        for (int j = 0; j < (int) meshes_without_wakes.size(); j++) {
            doublet_coefficient_matrix_block(A, b, *meshes_without_wakes[i], offset_one, *meshes_without_wakes[j], offset_two);
            
            offset_two = offset_two + meshes_without_wakes[j]->n_panels();
        }
        
        wakes_influence(A, b, *meshes_without_wakes[i], offset_one);
        
        offset_one = offset_one + meshes_without_wakes[i]->n_panels();
    }
    
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
    VectorXd old_potentials = potentials;
    if (Parameters::unsteady_bernoulli)
        potentials = surface_potentials();

    // Compute pressure distribution:
    cout << "Solver: Computing pressure distribution." << endl;

    int pressure_coefficients_idx = 0;
    
    offset = 0;

    for (int i = 0; i < (int) collections.size(); i++) {
        Collection *collection = collections[i];
        
        double v_ref = (collection->velocity - freestream_velocity).norm();
        
        VectorXd doublet_coefficient_field(collection->nolift_mesh.n_panels());
        for (int j = 0; j < collection->nolift_mesh.n_panels(); j++)
            doublet_coefficient_field(j) = doublet_coefficients(offset + j); 
            
        for (int j = 0; j < collection->nolift_mesh.n_panels(); j++) {
            Vector3d kinematic_velocity = collection->nolift_mesh.panel_deformation_velocity(j)
                                          + collection->panel_kinematic_velocity(collection->nolift_mesh, j) 
                                          - freestream_velocity;

            double dpotentialdt;
            if (Parameters::unsteady_bernoulli && dt > 0.0)
                dpotentialdt = (potentials(offset + j) - old_potentials(offset + j)) / dt;
            else
                dpotentialdt = 0.0;
                                                 
            pressure_coefficients(pressure_coefficients_idx) = pressure_coefficient(collection->nolift_mesh, j,
                                                                                    kinematic_velocity,
                                                                                    doublet_coefficient_field,
                                                                                    dpotentialdt, v_ref);
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
                Vector3d kinematic_velocity = wing->panel_deformation_velocity(k) 
                                              + collection->panel_kinematic_velocity(*wing, k) 
                                              - freestream_velocity;
               
                double dpotentialdt;
                if (Parameters::unsteady_bernoulli && dt > 0.0)
                    dpotentialdt = (potentials(offset + k) - old_potentials(offset + k)) / dt;
                else
                    dpotentialdt = 0.0;
 
                pressure_coefficients(pressure_coefficients_idx) = pressure_coefficient(*wing, k,
                                                                                        kinematic_velocity,
                                                                                        doublet_coefficient_field,
                                                                                        dpotentialdt, v_ref);
                
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
                
                for (int k = 0; k < wake->n_nodes(); k++) {  
                    Vector3d kinematic_velocity = - freestream_velocity;
                    
                    Vector3d x = wake->nodes[k];

                    Vector3d velocity = stream_velocity(x, kinematic_velocity);
                    
                    local_wake_velocities.push_back(velocity);
                }
                
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
                
                // Convect wake nodes that coincide with the trailing edge nodes with the freestream velocity,
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
                wake->add_layer(meshes);
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
                
                // Need to update influence coefficients:
                wake->invalidate_cache();
            }
        }
    }
}

/**
   Computes the aerodynamic force on a given collection.
   
   @param[in]   collection  Reference collection.
  
   @returns Aerodynamic force.
*/
Eigen::Vector3d
Solver::aerodynamic_force(Collection &collection)
{
    double v_ref = (collection.velocity - freestream_velocity).norm();
        
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
Solver::aerodynamic_moment(Collection &collection, const Eigen::Vector3d &x)
{
    double v_ref = (collection.velocity - freestream_velocity).norm();
        
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
Solver::log_coefficients(int step_number)
{   
    // Log coefficients: 
    int offset = 0;
    int save_node_offset = 0;
    int save_panel_offset = 0;
    
    for (int i = 0; i < (int) collections.size(); i++) {
        Collection *collection = collections[i];
        
        offset += collection->nolift_mesh.n_panels();
        
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

#if 0
    // Make a snapshot of the velocity field:
    stringstream str_x;
    str_x << "/home/jorn/snapshot_x_" << step_number << ".txt";
    ofstream snapshot_x(str_x.str().c_str());
    stringstream str_y;
    str_y << "/home/jorn/snapshot_y_" << step_number << ".txt";
    ofstream snapshot_y(str_y.str().c_str());
    stringstream str_p;
    str_p << "/home/jorn/snapshot_p_" << step_number << ".txt";
    ofstream snapshot_p(str_p.str().c_str());
    stringstream str_a;
    str_a << "/home/jorn/snapshot_a_" << step_number << ".txt";
    ofstream snapshot_a(str_a.str().c_str());
    
    double x_min = -0.5;
    double x_max = 1.5;
    double y_min = -0.2;
    double y_max = 0.2;
    double z = -0.1875;
    for (int i = 0; i <= 20; i++) {
        for (int j = 0; j <= 20; j++) {
            if (j > 0) {
                snapshot_x << ' ';
                snapshot_y << ' ';
                snapshot_p << ' ';
            }
            
            Vector3d kinematic_velocity = - freestream_velocity;
                
            Vector3d x = Vector3d(x_min + (x_max - x_min) / 20.0 * i, y_min + (y_max - y_min) / 20.0 * j, z);
            
            Vector3d velocity = stream_velocity(x, kinematic_velocity);
            
            snapshot_x << velocity(0);
            snapshot_y << velocity(1);
            snapshot_p << potential(x);
        }
        
        snapshot_x << endl;
        snapshot_y << endl;
        snapshot_p << endl;
    }
    
    for (int i = 0; i < meshes_without_wakes[2]->nodes.size(); i++) {
        Vector3d n = meshes_without_wakes[2]->nodes[i];
        
        if (meshes_without_wakes[2]->node_panel_neighbors[i].size() == 4)
            snapshot_a << n(0) << ' ' << n(1) << endl;
    }
    
    snapshot_x.close();
    snapshot_y.close();
    snapshot_p.close();
    snapshot_a.close();
#endif
}
