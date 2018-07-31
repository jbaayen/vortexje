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
#include <typeinfo>

#ifdef _WIN32
#include <direct.h>
#endif

#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <Eigen/IterativeLinearSolvers>

#include <vortexje/solver.hpp>
#include <vortexje/parameters.hpp>
#include <vortexje/boundary-layers/dummy-boundary-layer.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

// String constants:
#define VIEW_NAME_SOURCE_DISTRIBUTION   "sigma"
#define VIEW_NAME_DOUBLET_DISTRIBUTION  "mu"
#define VIEW_NAME_PRESSURE_DISTRIBUTION "Cp"
#define VIEW_NAME_VELOCITY_DISTRIBUTION "V"

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
Solver::Solver(const string &log_folder) : log_folder(log_folder)
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
   Adds a body to this solver.
   
   @param[in]   body   Body to be added.
*/
void
Solver::add_body(std::shared_ptr<Body> body)
{
    static shared_ptr<DummyBoundaryLayer> dummy_boundary_layer(new DummyBoundaryLayer());
    
    add_body(body, dummy_boundary_layer);
}

/**
   Adds a body, with boundary layer model, to this solver.
   
   @param[in]   body             Body to be added.
   @param[in]   boundary_layer   Boundary layer model.
*/
void
Solver::add_body(std::shared_ptr<Body> body, std::shared_ptr<BoundaryLayer> boundary_layer)
{
    shared_ptr<BodyData> bd(new BodyData(body, boundary_layer));
    
    bodies.push_back(bd);
    
    vector<shared_ptr<Body::SurfaceData> >::iterator si;
    for (si = body->non_lifting_surfaces.begin(); si != body->non_lifting_surfaces.end(); si++) {
        shared_ptr<Body::SurfaceData> d = *si;
        
        non_wake_surfaces.push_back(d);
           
        surface_to_body[d->surface] = bd;
        
        n_non_wake_panels += d->surface->n_panels();
    }
    
    vector<shared_ptr<Body::LiftingSurfaceData> >::iterator lsi;
    for (lsi = body->lifting_surfaces.begin(); lsi != body->lifting_surfaces.end(); lsi++) {
        shared_ptr<Body::LiftingSurfaceData> d = *lsi;
        
        non_wake_surfaces.push_back(d);
           
        surface_to_body[d->surface] = bd;
        surface_to_body[d->wake]    = bd;
        
        n_non_wake_panels += d->lifting_surface->n_panels();
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
    
    previous_surface_velocity_potentials.resize(n_non_wake_panels);
    previous_surface_velocity_potentials.setZero();

    // Open logs:
    string body_log_folder = log_folder + "/" + body->id;
    
    mkdir_helper(body_log_folder);
    
    for (si = bd->body->non_lifting_surfaces.begin(); si != bd->body->non_lifting_surfaces.end(); si++) {
        shared_ptr<Body::SurfaceData> d = *si;
        
        stringstream ss;
        ss << body_log_folder << "/" << d->surface->id;
        
        mkdir_helper(ss.str());
    }
          
    for (lsi = bd->body->lifting_surfaces.begin(); lsi != bd->body->lifting_surfaces.end(); lsi++) {
        shared_ptr<Body::LiftingSurfaceData> d = *lsi;
        
        stringstream ss;
        ss << body_log_folder << "/" << d->lifting_surface->id;
        
        mkdir_helper(ss.str());
        
        ss.str(string());
        ss.clear();
        ss << body_log_folder << "/" << d->wake->id;
        
        mkdir_helper(ss.str());      
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

/**
   Computes the velocity potential at the given point.
   
   @param[in]   x   Reference point.
   
   @returns Velocity potential.
*/
double
Solver::velocity_potential(const Vector3d &x) const
{
    return compute_velocity_potential(x);
}

/**
   Computes the total stream velocity at the given point.
   
   @param[in]   x   Reference point.
   
   @returns Stream velocity.
*/
Eigen::Vector3d
Solver::velocity(const Eigen::Vector3d &x) const
{
    std::set<int> ignore_set;
    
    return compute_velocity_interpolated(x, ignore_set);
}

/**
   Returns the surface velocity potential for the given panel.
   
   @param[in]   surface   Reference surface.
   @param[in]   panel     Reference panel.
  
   @returns Surface velocity potential.
*/
double
Solver::surface_velocity_potential(const shared_ptr<Surface> &surface, int panel) const
{
    int index = compute_index(surface, panel);
    if (index >= 0)
        return surface_velocity_potentials(index);

    cerr << "Solver::surface_velocity_potential():  Panel " << panel << " not found on surface " << surface->id << "." << endl;
    
    return 0.0;
}

/**
   Returns the surface velocity for the given panel.
   
   @param[in]   surface   Reference surface.
   @param[in]   panel     Reference panel.
  
   @returns Surface velocity.
*/
Vector3d
Solver::surface_velocity(const shared_ptr<Surface> &surface, int panel) const
{
    int index = compute_index(surface, panel);
    if (index >= 0)
        return surface_velocities.row(index);
    
    cerr << "Solver::surface_velocity():  Panel " << panel << " not found on surface " << surface->id << "." << endl;
    
    return Vector3d(0, 0, 0);
}

/**
   Returns the pressure coefficient for the given panel.
   
   @param[in]   surface   Reference surface.
   @param[in]   panel     Reference panel.
  
   @returns Pressure coefficient.
*/
double
Solver::pressure_coefficient(const shared_ptr<Surface> &surface, int panel) const
{
    int index = compute_index(surface, panel);
    if (index >= 0)
        return pressure_coefficients(index);
    
    cerr << "Solver::pressure_coefficient():  Panel " << panel << " not found on surface " << surface->id << "." << endl;
    
    return 0.0;
}

/**
   Computes the force caused by the pressure distribution on the given body.
   
   @param[in]   body   Reference body.
  
   @returns Force vector.
*/
Eigen::Vector3d
Solver::force(const std::shared_ptr<Body> &body) const
{
    // Dynamic pressure:
    double q = 0.5 * fluid_density * compute_reference_velocity_squared(body);
        
    // Total force on body:
    Vector3d F(0, 0, 0);
    int offset = 0;
    
    vector<shared_ptr<Body::SurfaceData> >::const_iterator si;
    for (si = non_wake_surfaces.begin(); si != non_wake_surfaces.end(); si++) {
        const shared_ptr<Body::SurfaceData> &d = *si;
        
        const shared_ptr<BodyData> &bd = surface_to_body.find(d->surface)->second;
        if (body == bd->body) {        
            for (int i = 0; i < d->surface->n_panels(); i++) {
                const Vector3d &normal = d->surface->panel_normal(i);
                double surface_area = d->surface->panel_surface_area(i);
                F += q * surface_area * pressure_coefficients(offset + i) * normal;
                
                F += bd->boundary_layer->friction(d->surface, i);
            }
        }
        
        offset += d->surface->n_panels();
    }
                    
    // Done:
    return F;      
}

/**
   Computes the force caused by the pressure distribution on the given surface.
   
   @param[in]   surface   Reference surface.
  
   @returns Force vector.
*/
Eigen::Vector3d
Solver::force(const std::shared_ptr<Surface> &surface) const
{
    // Total force on surface:
    Vector3d F(0, 0, 0);
    int offset = 0;
    
    vector<shared_ptr<Body::SurfaceData> >::const_iterator si;
    for (si = non_wake_surfaces.begin(); si != non_wake_surfaces.end(); si++) {
        const shared_ptr<Body::SurfaceData> &d = *si;
        
        if (d->surface == surface) {   
            const shared_ptr<BodyData> &bd = surface_to_body.find(d->surface)->second;
            
            // Dynamic pressure:
            double q = 0.5 * fluid_density * compute_reference_velocity_squared(bd->body);
                 
            for (int i = 0; i < d->surface->n_panels(); i++) {
                const Vector3d &normal = d->surface->panel_normal(i);
                double surface_area = d->surface->panel_surface_area(i);
                F += q * surface_area * pressure_coefficients(offset + i) * normal;
                
                F += bd->boundary_layer->friction(d->surface, i);
            }
            
            break;
        }
        
        offset += d->surface->n_panels();
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
Solver::moment(const std::shared_ptr<Body> &body, const Eigen::Vector3d &x) const
{
    // Dynamic pressure:
    double q = 0.5 * fluid_density * compute_reference_velocity_squared(body);
        
    // Total moment on body:
    Vector3d M(0, 0, 0);
    int offset = 0;
    
    vector<shared_ptr<Body::SurfaceData> >::const_iterator si;
    for (si = non_wake_surfaces.begin(); si != non_wake_surfaces.end(); si++) {
        const shared_ptr<Body::SurfaceData> &d = *si;
        
        const shared_ptr<BodyData> &bd = surface_to_body.find(d->surface)->second;
        if (body == bd->body) { 
            for (int i = 0; i < d->surface->n_panels(); i++) {                                    
                const Vector3d &normal = d->surface->panel_normal(i);
                double surface_area = d->surface->panel_surface_area(i);
                Vector3d F = q * surface_area * pressure_coefficients(offset + i) * normal;
                
                F += bd->boundary_layer->friction(d->surface, i);
                    
                Vector3d r = d->surface->panel_collocation_point(i, false) - x;
                M += r.cross(F);
            }
        }
        
        offset += d->surface->n_panels();
    }
    
    // Done:
    return M;
}

/**
   Computes the moment caused by the pressure distribution on the given surface, relative to the given point.
   
   @param[in]   surface   Reference surface.
   @param[in]   x         Reference point.
  
   @returns Moment vector.
*/
Eigen::Vector3d
Solver::moment(const std::shared_ptr<Surface> &surface, const Eigen::Vector3d &x) const
{
    // Total moment on surface:
    Vector3d M(0, 0, 0);
    int offset = 0;
    
    vector<shared_ptr<Body::SurfaceData> >::const_iterator si;
    for (si = non_wake_surfaces.begin(); si != non_wake_surfaces.end(); si++) {
        const shared_ptr<Body::SurfaceData> &d = *si;

        if (d->surface == surface) {
            const shared_ptr<BodyData> &bd = surface_to_body.find(d->surface)->second;
            
            // Dynamic pressure:
            double q = 0.5 * fluid_density * compute_reference_velocity_squared(bd->body);
            
            for (int i = 0; i < d->surface->n_panels(); i++) {                                    
                const Vector3d &normal = d->surface->panel_normal(i);
                double surface_area = d->surface->panel_surface_area(i);
                Vector3d F = q * surface_area * pressure_coefficients(offset + i) * normal;
                
                F += bd->boundary_layer->friction(d->surface, i);
                    
                Vector3d r = d->surface->panel_collocation_point(i, false) - x;
                M += r.cross(F);
            }
            
            break;
        }
        
        offset += d->surface->n_panels();
    }
    
    // Done:
    return M;
}

/**
   Traces a streamline, starting from the given starting point.
   
   @param[in]   start   Starting point for streamline.
   
   @returns A list of points tracing the streamline.
*/
vector_aligned<Solver::SurfacePanelPoint>
Solver::trace_streamline(const SurfacePanelPoint &start) const
{
    vector_aligned<SurfacePanelPoint> streamline;
    
    SurfacePanelPoint cur(start.surface, start.panel, start.point);
    
    Vector3d prev_intersection = start.point;
    int originating_edge = -1;
    
    // Trace until we hit the end of a surface, or until we hit a stagnation point:
    while (true) {
        // Look up panel velocity:
        Vector3d velocity = surface_velocity(cur.surface, cur.panel);
        
        // Stop following the streamline at stagnation points:
        if (velocity.norm() < Parameters::zero_threshold)
            break;
        
        // Transform into a panel frame:
        const Transform<double, 3, Affine> &transformation = cur.surface->panel_coordinate_transformation(cur.panel);
        
        Vector3d transformed_velocity = transformation.linear() * velocity;        
        Vector3d transformed_point = transformation * cur.point;
        
        // Project point onto panel for improved numerical stability:
        transformed_point(2) = 0.0;
        
        // Intersect with one of the panel edges.
        int edge_id = -1;
        double t = numeric_limits<double>::max();
        for (int i = 0; i < (int) cur.surface->panel_nodes[cur.panel].size(); i++) {
            // Do not try to intersect with the edge we are already on:
            if (i == originating_edge)
                continue;
                
            // Compute next node index:
            int next_idx;
            if (i == (int) cur.surface->panel_nodes[cur.panel].size() - 1)
                next_idx = 0;
            else
                next_idx = i + 1;
              
            // Retrieve nodes in panel-local coordinates:  
            const Vector3d &node_a = cur.surface->panel_transformed_points[cur.panel][i];
            const Vector3d &node_b = cur.surface->panel_transformed_points[cur.panel][next_idx];
            
            // Compute edge:
            Vector3d edge = node_b - node_a;
            
            // Find intersection, if any:
            Matrix2d A;
            Vector2d b;
            for (int j = 0; j < 2; j++) {
                A(j, 0) = transformed_velocity(j);
                A(j, 1) = -edge(j);
                b(j)    = -transformed_point(j) + node_a(j);
            }
            
            ColPivHouseholderQR<Matrix2d> solver(A);
            solver.setThreshold(Parameters::zero_threshold);
            if (!solver.isInvertible())
                continue;
            
            Vector2d x = solver.solve(b);
  
            // Only accept positive quadrant:
            if (x(0) < 0 || x(1) < 0)
                continue;
                
            // Do not accept infinitesimally small solutions:
            if (x(0) < Parameters::zero_threshold)
                continue;
                
            // Is this the smallest positive 't' (velocity coefficient)?   
            if (x(0) < t) {
                t = x(0);
                
                edge_id = i;
            }
        }
        
        // Dead end?
        if (edge_id < 0)
            break;
        
        // Compute intersection vector:
        Vector3d transformed_intersection = transformed_point + t * transformed_velocity;
        Vector3d intersection = transformation.inverse() * transformed_intersection;
        
        // Compute mean between intersection points:
        Vector3d mean_point = 0.5 * (intersection + prev_intersection);
        prev_intersection = intersection;
        
        // Add to streamline:
        SurfacePanelPoint n(cur.surface, cur.panel, mean_point);
        streamline.push_back(n);
        
        // Find neighbor across edge: 
        // N.B.:  This code assumes that every panel has at most one neighbor across an edge.  
        //        The rest of Vortexje supports more general geometries, however.  This code needs work.
        const shared_ptr<BodyData> &bd = surface_to_body.find(cur.surface)->second;
        vector<Body::SurfacePanelEdge> neighbors = bd->body->panel_neighbors(cur.surface, cur.panel, edge_id);
        
        // No neighbor?
        if (neighbors.size() == 0)
            break;

        // Verify the direction of the neighboring velocity vector:
        Vector3d neighbor_velocity = surface_velocity(neighbors[0].surface, neighbors[0].panel);
        
        const Vector3d &normal          = cur.surface->panel_normal(cur.panel);
        const Vector3d &neighbor_normal = neighbors[0].surface->panel_normal(neighbors[0].panel);
        
        Quaterniond unfold = Quaterniond::FromTwoVectors(neighbor_normal, normal);
        Vector3d unfolded_neighbor_velocity = unfold * neighbor_velocity;
       
        if (velocity.dot(unfolded_neighbor_velocity) < 0) {
            // Velocity vectors point in opposite directions.
            break;
        }
            
        // Proceed to neighboring panel:
        cur.surface = neighbors[0].surface;
        cur.panel   = neighbors[0].panel;
        cur.point   = intersection;
        
        originating_edge = neighbors[0].edge;
    }
    
    // Done:
    return streamline;
}

/**
   Initializes the wakes by adding a first layer of vortex ring panels.
   
   @param[in]   dt   Time step size.
*/
void
Solver::initialize_wakes(double dt)
{
    // Add initial wake layers:
    vector<shared_ptr<BodyData> >::iterator bdi;
    for (bdi = bodies.begin(); bdi != bodies.end(); bdi++) {
        shared_ptr<BodyData> bd = *bdi;
    
        vector<shared_ptr<Body::LiftingSurfaceData> >::iterator lsi;
        for (lsi = bd->body->lifting_surfaces.begin(); lsi != bd->body->lifting_surfaces.end(); lsi++) {
            shared_ptr<Body::LiftingSurfaceData> d = *lsi;
            
            d->wake->add_layer();
            for (int i = 0; i < d->lifting_surface->n_spanwise_nodes(); i++) {
                if (Parameters::convect_wake) {
                    // Convect wake nodes that coincide with the trailing edge.
                    d->wake->nodes[i] += compute_trailing_edge_vortex_displacement(bd->body, d->lifting_surface, i, dt);
                    
                } else {
                    // Initialize static wake->
                    Vector3d body_apparent_velocity = bd->body->velocity - freestream_velocity;
                    
                    d->wake->nodes[i] -= Parameters::static_wake_length * body_apparent_velocity / body_apparent_velocity.norm();
                }
            }
            
            d->wake->add_layer();
        }
    }
}

/**
   Computes new source, doublet, and pressure distributions.
   
   @param[in]   dt          Time step size.
   @param[in]   propagate   Propagate solution forward in time.
   
   @returns true on success.
*/
bool
Solver::solve(double dt, bool propagate)
{
    int offset;
    
    // Iterate inviscid and boundary layer solutions until convergence.
    VectorXd previous_source_coefficients;
    VectorXd previous_doublet_coefficients;
    
    int boundary_layer_iteration = 0;
    
    while (true) {
        // Copy state:
        previous_source_coefficients  = source_coefficients;
        previous_doublet_coefficients = doublet_coefficients;
        
        // Compute new source distribution:
        cout << "Solver: Computing source distribution with wake influence." << endl;
        
        offset = 0;
        
        vector<shared_ptr<Body::SurfaceData> >::const_iterator si;
        for (si = non_wake_surfaces.begin(); si != non_wake_surfaces.end(); si++) {
            const shared_ptr<Body::SurfaceData> &d = *si;
            int i;
            
            const shared_ptr<BodyData> &bd = surface_to_body.find(d->surface)->second;
            
            #pragma omp parallel
            {
                #pragma omp for schedule(dynamic, 1)
                for (i = 0; i < d->surface->n_panels(); i++)
                    source_coefficients(offset + i) = compute_source_coefficient(bd->body, d->surface, i, bd->boundary_layer, true);
            }
            
            offset += d->surface->n_panels();
        }
      
        // Populate the matrices of influence coefficients:
        cout << "Solver: Computing matrices of influence coefficients." << endl;
        
        MatrixXd A(n_non_wake_panels, n_non_wake_panels);
        MatrixXd source_influence_coefficients(n_non_wake_panels, n_non_wake_panels);
        
        int offset_row = 0, offset_col = 0;
        
        vector<shared_ptr<Body::SurfaceData> >::const_iterator si_row;
        for (si_row = non_wake_surfaces.begin(); si_row != non_wake_surfaces.end(); si_row++) {
            const shared_ptr<Body::SurfaceData> &d_row = *si_row;

            offset_col = 0;
     
            // Influence coefficients between all non-wake surfaces:
            vector<shared_ptr<Body::SurfaceData> >::const_iterator si_col;
            for (si_col = non_wake_surfaces.begin(); si_col != non_wake_surfaces.end(); si_col++) {
                shared_ptr<Body::SurfaceData> d_col = *si_col;
                int i, j;
                
                #pragma omp parallel private(j) 
                {
                    #pragma omp for schedule(dynamic, 1)
                    for (i = 0; i < d_row->surface->n_panels(); i++) {
                        for (j = 0; j < d_col->surface->n_panels(); j++) {
                            d_col->surface->source_and_doublet_influence(d_row->surface, i, j,
                                                                         source_influence_coefficients(offset_row + i, offset_col + j), 
                                                                         A(offset_row + i, offset_col + j));
                        }
                    }
                }
                
                offset_col = offset_col + d_col->surface->n_panels();
            }
            
            // The influence of the new wake panels:
            int i, j, lifting_surface_offset, wake_panel_offset, pa, pb;
            vector<shared_ptr<BodyData> >::const_iterator bdi;
            vector<shared_ptr<Body::SurfaceData> >::const_iterator si;
            vector<shared_ptr<Body::LiftingSurfaceData> >::const_iterator lsi;
            shared_ptr<BodyData> bd;
            shared_ptr<Body::LiftingSurfaceData> d;
            
            #pragma omp parallel private(bdi, si, lsi, lifting_surface_offset, j, wake_panel_offset, pa, pb, bd, d)
            {
                #pragma omp for schedule(dynamic, 1)
                for (i = 0; i < d_row->surface->n_panels(); i++) {
                    lifting_surface_offset = 0;      
                    
                    for (bdi = bodies.begin(); bdi != bodies.end(); bdi++) {
                        bd = *bdi;
                        
                        for (si = bd->body->non_lifting_surfaces.begin(); si != bd->body->non_lifting_surfaces.end(); si++)
                            lifting_surface_offset += (*si)->surface->n_panels();
                                      
                        for (lsi = bd->body->lifting_surfaces.begin(); lsi != bd->body->lifting_surfaces.end(); lsi++) {
                            d = *lsi;
                            
                            wake_panel_offset = d->wake->n_panels() - d->lifting_surface->n_spanwise_panels();
                            for (j = 0; j < d->lifting_surface->n_spanwise_panels(); j++) {  
                                pa = d->lifting_surface->trailing_edge_upper_panel(j);
                                pb = d->lifting_surface->trailing_edge_lower_panel(j);
                                
                                // Account for the influence of the new wake panels.  The doublet strength of these panels
                                // is set according to the Kutta condition.
                                A(offset_row + i, lifting_surface_offset + pa) += d->wake->doublet_influence(d_row->surface, i, wake_panel_offset + j);
                                A(offset_row + i, lifting_surface_offset + pb) -= d->wake->doublet_influence(d_row->surface, i, wake_panel_offset + j);
                            }
                            
                            lifting_surface_offset += d->lifting_surface->n_panels();
                        }
                    }
                }
            }
                
            offset_row = offset_row + d_row->surface->n_panels();
        }
        
        // Compute new doublet distribution:
        cout << "Solver: Computing doublet distribution." << endl;
        
        VectorXd b = source_influence_coefficients * source_coefficients;
        
        BiCGSTAB<MatrixXd, DiagonalPreconditioner<double> > solver(A);
        solver.setMaxIterations(Parameters::linear_solver_max_iterations);
        solver.setTolerance(Parameters::linear_solver_tolerance);

        doublet_coefficients = solver.solveWithGuess(b, previous_doublet_coefficients);
        
        if (solver.info() != Success) {
            cerr << "Solver: Computing doublet distribution failed (" << solver.iterations();
            cerr << " iterations with estimated error=" << solver.error() << ")." << endl;
           
            return false;
        }
        
        cout << "Solver: Done computing doublet distribution in " << solver.iterations() << " iterations with estimated error " << solver.error() << "." << endl;

        // Check for convergence from second iteration onwards.
        bool converged = false;
        if (boundary_layer_iteration > 0) {
            double delta = (source_coefficients - previous_source_coefficients).norm();
            
            cout << "Solver: Boundary layer convergence delta = " << delta << endl;
            
            if (delta < Parameters::boundary_layer_iteration_tolerance)
                converged = true;
        }
        
        // Set new wake panel doublet coefficients:
        cout << "Solver: Updating wake doublet distribution." << endl;
        
        offset = 0;
        
        vector<shared_ptr<BodyData> >::iterator bdi;
        for (bdi = bodies.begin(); bdi != bodies.end(); bdi++) {
            shared_ptr<BodyData> bd = *bdi;
            
            vector<shared_ptr<Body::SurfaceData> >::iterator si;
            for (si = bd->body->non_lifting_surfaces.begin(); si != bd->body->non_lifting_surfaces.end(); si++)
                offset += (*si)->surface->n_panels();
            
            vector<shared_ptr<Body::LiftingSurfaceData> >::iterator lsi;
            for (lsi = bd->body->lifting_surfaces.begin(); lsi != bd->body->lifting_surfaces.end(); lsi++) {
                shared_ptr<Body::LiftingSurfaceData> d = *lsi;
                         
                // Set panel doublet coefficient:
                for (int i = 0; i < d->lifting_surface->n_spanwise_panels(); i++) {
                    double doublet_coefficient_top    = doublet_coefficients(offset + d->lifting_surface->trailing_edge_upper_panel(i));
                    double doublet_coefficient_bottom = doublet_coefficients(offset + d->lifting_surface->trailing_edge_lower_panel(i));
                    
                    // Use the trailing-edge Kutta condition to compute the doublet coefficients of the new wake panels.
                    double doublet_coefficient = doublet_coefficient_top - doublet_coefficient_bottom;
                    
                    int idx = d->wake->n_panels() - d->lifting_surface->n_spanwise_panels() + i;
                    d->wake->doublet_coefficients[idx] = doublet_coefficient;
                }
                
                // Update offset:
                offset += d->lifting_surface->n_panels();
            }
        }
        
        // Compute surface velocity distribution:
        cout << "Solver: Computing surface velocity distribution." << endl;
        
        offset = 0;

        for (bdi = bodies.begin(); bdi != bodies.end(); bdi++) {
            shared_ptr<BodyData> bd = *bdi;
            
            vector<shared_ptr<Body::SurfaceData> >::const_iterator si;
            for (si = bd->body->non_lifting_surfaces.begin(); si != bd->body->non_lifting_surfaces.end(); si++) {
                const shared_ptr<Body::SurfaceData> &d = *si;
                int i;
                
                #pragma omp parallel
                {
                    #pragma omp for schedule(dynamic, 1)
                    for (i = 0; i < d->surface->n_panels(); i++)
                        surface_velocities.row(offset + i) = compute_surface_velocity(bd->body, d->surface, i);
                }
                
                offset += d->surface->n_panels();   
            }
            
            vector<shared_ptr<Body::LiftingSurfaceData> >::const_iterator lsi;
            for (lsi = bd->body->lifting_surfaces.begin(); lsi != bd->body->lifting_surfaces.end(); lsi++) {
                const shared_ptr<Body::LiftingSurfaceData> &d = *lsi;
                int i;
                
                #pragma omp parallel
                {
                    #pragma omp for schedule(dynamic, 1)
                    for (i = 0; i < d->surface->n_panels(); i++)
                        surface_velocities.row(offset + i) = compute_surface_velocity(bd->body, d->surface, i);
                }  
                
                offset += d->surface->n_panels();                           
            }
        }

        // If we converged, then this is the time to break out of the loop.
        if (converged) {
            cout << "Solver: Boundary layer iteration converged in " << boundary_layer_iteration << " steps." << endl;
            
            break;
        }
        
        if (boundary_layer_iteration > Parameters::max_boundary_layer_iterations) {
            cout << "Solver: Maximum number of boundary layer iterations ranged.  Aborting iteration." << endl;
            
            break;
        }
        
        // Recompute the boundary layers.
        offset = 0;
        
        bool have_boundary_layer = false;
        for (bdi = bodies.begin(); bdi != bodies.end(); bdi++) {
            shared_ptr<BodyData> bd = *bdi;
            
            // Count panels on body:
            int body_n_panels = 0;
                       
            vector<shared_ptr<Body::SurfaceData> >::const_iterator si;
            for (si = bd->body->non_lifting_surfaces.begin(); si != bd->body->non_lifting_surfaces.end(); si++) {
                const shared_ptr<Body::SurfaceData> &d = *si;

                body_n_panels += d->surface->n_panels();   
            }
            
            vector<shared_ptr<Body::LiftingSurfaceData> >::const_iterator lsi;
            for (lsi = bd->body->lifting_surfaces.begin(); lsi != bd->body->lifting_surfaces.end(); lsi++) {
                const shared_ptr<Body::LiftingSurfaceData> &d = *lsi;

                body_n_panels += d->surface->n_panels();                           
            }
            
            // Recompute boundary layer:
            if (typeid(*bd->boundary_layer.get()) != typeid(DummyBoundaryLayer)) {
                have_boundary_layer = true;
                
                if (!bd->boundary_layer->recalculate(freestream_velocity, surface_velocities.block(offset, 0, body_n_panels, 3)))
                    return false;
            }
            
            offset += body_n_panels;
        }
        
        // Did we did not find any boundary layers, then there is no need to iterate.
        if (!have_boundary_layer)
            break;
        
        // Increase iteration counter:
        boundary_layer_iteration++;
    }

    if (Parameters::convect_wake) {
        // Recompute source distribution without wake influence:
        cout << "Solver: Recomputing source distribution without wake influence." << endl;
        
        offset = 0;
        
        vector<shared_ptr<Body::SurfaceData> >::const_iterator si;
        for (si = non_wake_surfaces.begin(); si != non_wake_surfaces.end(); si++) {
            const shared_ptr<Body::SurfaceData> &d = *si;
            int i;
            
            const shared_ptr<BodyData> &bd = surface_to_body.find(d->surface)->second;
            
            #pragma omp parallel
            {
                #pragma omp for schedule(dynamic, 1)
                for (i = 0; i < d->surface->n_panels(); i++)
                    source_coefficients(offset + i) = compute_source_coefficient(bd->body, d->surface, i, bd->boundary_layer, false);
            }
            
            offset += d->surface->n_panels();
        }
    }

    // Compute pressure distribution:
    cout << "Solver: Computing pressure distribution." << endl;
    
    offset = 0;

    vector<shared_ptr<Body::SurfaceData> >::const_iterator si;   
    for (si = non_wake_surfaces.begin(); si != non_wake_surfaces.end(); si++) {
        const shared_ptr<Body::SurfaceData> &d = *si;
        int i;
        
        const shared_ptr<BodyData> &bd = surface_to_body.find(d->surface)->second;
        double v_ref_squared = compute_reference_velocity_squared(bd->body);
        
        double dphidt;
            
        #pragma omp parallel private(dphidt)
        {
            #pragma omp for schedule(dynamic, 1)
            for (i = 0; i < d->surface->n_panels(); i++) {
                // Velocity potential:
                surface_velocity_potentials(offset + i) = compute_surface_velocity_potential(d->surface, offset, i);
                
                // Pressure coefficient:
                dphidt = compute_surface_velocity_potential_time_derivative(offset, i, dt);
                pressure_coefficients(offset + i) = compute_pressure_coefficient(surface_velocities.row(offset + i), dphidt, v_ref_squared);
            }
        }
        
        offset += d->surface->n_panels();      
    }
    
    // Propagate solution forward in time, if requested.
    if (propagate)
        this->propagate();
    
    // Done:
    return true;
}

/**
   Propagates solution forward in time.  Relevant in unsteady mode only.
*/
void
Solver::propagate()
{
    // Store previous values of the surface velocity potentials:
    previous_surface_velocity_potentials = surface_velocity_potentials;
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
        cout << "Solver: Convecting wakes." << endl;
        
        // Compute velocity values at wake nodes, with the wakes in their original state:
        vector_aligned<vector_aligned<Vector3d> > wake_velocities;
        
        vector<shared_ptr<BodyData> >::const_iterator bdi;
        for (bdi = bodies.begin(); bdi != bodies.end(); bdi++) {
            const shared_ptr<BodyData> &bd = *bdi;
                 
            vector<shared_ptr<Body::LiftingSurfaceData> >::const_iterator lsi;
            for (lsi = bd->body->lifting_surfaces.begin(); lsi != bd->body->lifting_surfaces.end(); lsi++) {
                const shared_ptr<Body::LiftingSurfaceData> &d = *lsi;
                
                vector_aligned<Vector3d> local_wake_velocities;
                local_wake_velocities.resize(d->wake->n_nodes());
                
                int i;
                
                #pragma omp parallel
                {
                    #pragma omp for schedule(dynamic, 1)
                    for (i = 0; i < d->wake->n_nodes(); i++)
                        local_wake_velocities[i] = velocity(d->wake->nodes[i]);
                }
                
                wake_velocities.push_back(local_wake_velocities);
            }
        }
        
        // Add new wake panels at trailing edges, and convect all vertices:
        int idx = 0;
        
        for (bdi = bodies.begin(); bdi != bodies.end(); bdi++) {
            shared_ptr<BodyData> bd = *bdi;
            
            vector<shared_ptr<Body::LiftingSurfaceData> >::iterator lsi;
            for (lsi = bd->body->lifting_surfaces.begin(); lsi != bd->body->lifting_surfaces.end(); lsi++) {
                shared_ptr<Body::LiftingSurfaceData> d = *lsi;
                
                // Retrieve local wake velocities:
                vector_aligned<Vector3d> &local_wake_velocities = wake_velocities[idx];
                idx++;
                
                // Convect wake nodes that coincide with the trailing edge.
                for (int i = 0; i < d->lifting_surface->n_spanwise_nodes(); i++) {                                                  
                    d->wake->nodes[d->wake->n_nodes() - d->lifting_surface->n_spanwise_nodes() + i]
                        += compute_trailing_edge_vortex_displacement(bd->body, d->lifting_surface, i, dt);
                }                
                
                // Convect all other wake nodes according to the local wake velocity:
                int i;
                
                #pragma omp parallel
                {
                    #pragma omp for schedule(dynamic, 1)
                    for (i = 0; i < d->wake->n_nodes() - d->lifting_surface->n_spanwise_nodes(); i++)
                        d->wake->nodes[i] += local_wake_velocities[i] * dt;
                }
                    
                // Run internal wake update:
                d->wake->update_properties(dt);

                // Add new vertices:
                // (This call also updates the geometry)
                d->wake->add_layer();
            }
        }
        
    } else {
        cout << "Solver: Re-positioning wakes." << endl;
        
        // No wake convection.  Re-position wake:
        vector<shared_ptr<BodyData> >::iterator bdi;
        for (bdi = bodies.begin(); bdi != bodies.end(); bdi++) {
            shared_ptr<BodyData> bd = *bdi;
            
            Vector3d body_apparent_velocity = bd->body->velocity - freestream_velocity;
            
            vector<shared_ptr<Body::LiftingSurfaceData> >::iterator lsi;
            for (lsi = bd->body->lifting_surfaces.begin(); lsi != bd->body->lifting_surfaces.end(); lsi++) {
                shared_ptr<Body::LiftingSurfaceData> d = *lsi;
            
                for (int i = 0; i < d->lifting_surface->n_spanwise_nodes(); i++) {
                    // Connect wake to trailing edge nodes:                             
                    d->wake->nodes[d->lifting_surface->n_spanwise_nodes() + i] = d->lifting_surface->nodes[d->lifting_surface->trailing_edge_node(i)];
                    
                    // Point wake in direction of body kinematic velocity:
                    d->wake->nodes[i] = d->lifting_surface->nodes[d->lifting_surface->trailing_edge_node(i)]
                                     - Parameters::static_wake_length * body_apparent_velocity / body_apparent_velocity.norm();
                }
                
                // Need to update geometry:
                d->wake->compute_geometry();
            }
        }
    }
}

/**
   Logs source and doublet distributions, as well as the pressure coefficients, into files in the logging folder
   tagged with the specified step number.
   
   @param[in]   step_number   Step number used to name the output files.
   @param[in]   writer        SurfaceWriter object to use.
*/
void
Solver::log(int step_number, SurfaceWriter &writer) const
{   
    // Log coefficients: 
    int offset = 0;
    int save_node_offset = 0;
    int save_panel_offset = 0;
    int idx;
    
    vector<shared_ptr<BodyData> >::const_iterator bdi;
    for (bdi = bodies.begin(); bdi != bodies.end(); bdi++) {
        const shared_ptr<BodyData> &bd = *bdi;
        
        // Iterate non-lifting surfaces:
        idx = 0;
        
        vector<shared_ptr<Body::SurfaceData> >::const_iterator si;
        for (si = bd->body->non_lifting_surfaces.begin(); si != bd->body->non_lifting_surfaces.end(); si++) {
            const shared_ptr<Body::SurfaceData> &d = *si;
            
            // Log non-lifting surface coefficients:
            MatrixXd non_lifting_surface_doublet_coefficients(d->surface->n_panels(), 1);
            MatrixXd non_lifting_surface_source_coefficients(d->surface->n_panels(), 1);
            MatrixXd non_lifting_surface_pressure_coefficients(d->surface->n_panels(), 1);
            MatrixXd non_lifting_surface_velocity_vectors(d->surface->n_panels(), 3);
            for (int i = 0; i < d->surface->n_panels(); i++) {
                non_lifting_surface_doublet_coefficients(i, 0)  = doublet_coefficients(offset + i);
                non_lifting_surface_source_coefficients(i, 0)   = source_coefficients(offset + i);
                non_lifting_surface_pressure_coefficients(i, 0) = pressure_coefficients(offset + i);
                non_lifting_surface_velocity_vectors.row(i)     = surface_velocities.row(offset + i);
            }
            
            offset += d->surface->n_panels();
            
            vector<string> view_names;
            vector_aligned<MatrixXd> view_data;
            
            view_names.push_back(VIEW_NAME_DOUBLET_DISTRIBUTION);
            view_data.push_back(non_lifting_surface_doublet_coefficients);
            
            view_names.push_back(VIEW_NAME_SOURCE_DISTRIBUTION);
            view_data.push_back(non_lifting_surface_source_coefficients);
            
            view_names.push_back(VIEW_NAME_PRESSURE_DISTRIBUTION);
            view_data.push_back(non_lifting_surface_pressure_coefficients);
            
            view_names.push_back(VIEW_NAME_VELOCITY_DISTRIBUTION);
            view_data.push_back(non_lifting_surface_velocity_vectors);
            
            stringstream ss;
            ss << log_folder << "/" << bd->body->id << "/" << d->surface->id << "/step_" << step_number << writer.file_extension();

            writer.write(d->surface, ss.str(), save_node_offset, save_panel_offset, view_names, view_data);
            
            save_node_offset += d->surface->n_nodes();
            save_panel_offset += d->surface->n_panels();
            
            idx++;
        }   
        
        // Iterate lifting surfaces:
        idx = 0;
        
        vector<shared_ptr<Body::LiftingSurfaceData> >::const_iterator lsi;
        for (lsi = bd->body->lifting_surfaces.begin(); lsi != bd->body->lifting_surfaces.end(); lsi++) {
            const shared_ptr<Body::LiftingSurfaceData> &d = *lsi;
            
            // Log lifting surface coefficients:
            MatrixXd lifting_surface_doublet_coefficients(d->lifting_surface->n_panels(), 1);
            MatrixXd lifting_surface_source_coefficients(d->lifting_surface->n_panels(), 1);
            MatrixXd lifting_surface_pressure_coefficients(d->lifting_surface->n_panels(), 1);
            MatrixXd lifting_surface_velocity_vectors(d->surface->n_panels(), 3);
            for (int i = 0; i < d->lifting_surface->n_panels(); i++) {
                lifting_surface_doublet_coefficients(i, 0)  = doublet_coefficients(offset + i);
                lifting_surface_source_coefficients(i, 0)   = source_coefficients(offset + i);
                lifting_surface_pressure_coefficients(i, 0) = pressure_coefficients(offset + i);
                lifting_surface_velocity_vectors.row(i)     = surface_velocities.row(offset + i);
            }
            
            offset += d->lifting_surface->n_panels();
            
            vector<string> view_names;
            vector_aligned<MatrixXd> view_data;
            
            view_names.push_back(VIEW_NAME_DOUBLET_DISTRIBUTION);
            view_data.push_back(lifting_surface_doublet_coefficients);
            
            view_names.push_back(VIEW_NAME_SOURCE_DISTRIBUTION);
            view_data.push_back(lifting_surface_source_coefficients);
            
            view_names.push_back(VIEW_NAME_PRESSURE_DISTRIBUTION);
            view_data.push_back(lifting_surface_pressure_coefficients);
            
            view_names.push_back(VIEW_NAME_VELOCITY_DISTRIBUTION);
            view_data.push_back(lifting_surface_velocity_vectors);
            
            stringstream ss;
            ss << log_folder << "/" << bd->body->id << "/" << d->lifting_surface->id << "/step_" << step_number << writer.file_extension();

            writer.write(d->lifting_surface, ss.str(), save_node_offset, save_panel_offset, view_names, view_data);
            
            save_node_offset += d->lifting_surface->n_nodes();
            save_panel_offset += d->lifting_surface->n_panels();
    
            // Log wake surface and coefficients:
            MatrixXd wake_doublet_coefficients(d->wake->doublet_coefficients.size(), 1);
            for (int i = 0; i < (int) d->wake->doublet_coefficients.size(); i++)
                wake_doublet_coefficients(i, 0) = d->wake->doublet_coefficients[i];

            view_names.clear();
            view_data.clear();
            
            view_names.push_back(VIEW_NAME_DOUBLET_DISTRIBUTION);
            view_data.push_back(wake_doublet_coefficients);
            
            stringstream ssw;
            ssw << log_folder << "/" << bd->body->id << "/" << d->wake->id << "/step_" << step_number << writer.file_extension();
            
            writer.write(d->wake, ssw.str(), 0, save_panel_offset, view_names, view_data);
            
            save_node_offset += d->wake->n_nodes();
            save_panel_offset += d->wake->n_panels();
            
            idx++;
        }
    }
}
 
// Compute source coefficient for given surface and panel:
double
Solver::compute_source_coefficient(const std::shared_ptr<Body> &body, const std::shared_ptr<Surface> &surface, int panel, const std::shared_ptr<BoundaryLayer> &boundary_layer, bool include_wake_influence) const
{
    // Start with apparent velocity:
    Vector3d velocity = body->panel_kinematic_velocity(surface, panel) - freestream_velocity;
    
    // Wake contribution:
    if (Parameters::convect_wake && include_wake_influence) {
        vector<shared_ptr<BodyData> >::const_iterator bdi;
        for (bdi = bodies.begin(); bdi != bodies.end(); bdi++) {
            const shared_ptr<BodyData> &bd = *bdi;
            
            vector<shared_ptr<Body::LiftingSurfaceData> >::const_iterator lsi;
            for (lsi = bd->body->lifting_surfaces.begin(); lsi != bd->body->lifting_surfaces.end(); lsi++) {
                const shared_ptr<Body::LiftingSurfaceData> &d = *lsi;
                
                // Add influence of old wake panels.  That is, those wake panels which already have a doublet
                // strength assigned to them.
                for (int k = 0; k < d->wake->n_panels() - d->lifting_surface->n_spanwise_panels(); k++) {
                    // Use doublet panel - vortex ring equivalence.
                    velocity -= d->wake->vortex_ring_unit_velocity(surface->panel_collocation_point(panel, true), k)
                        * d->wake->doublet_coefficients[k];
                }
            }
        }
    }
    
    // Take normal component, and subtract blowing velocity:
    const Vector3d &normal = surface->panel_normal(panel);
    
    double blowing_velocity = boundary_layer->blowing_velocity(surface, panel);
        
    return velocity.dot(normal) - blowing_velocity;
}

/**
   Returns velocity potential value on the body surface.
   
   @returns Surface potential value.
*/
double
Solver::compute_surface_velocity_potential(const std::shared_ptr<Surface> &surface, int offset, int panel) const
{
    double phi = -doublet_coefficients(offset + panel);
    
    // Add flow potential due to kinematic velocity:
    const shared_ptr<BodyData> &bd = surface_to_body.find(surface)->second;
    Vector3d apparent_velocity = bd->body->panel_kinematic_velocity(surface, panel) - freestream_velocity;
    
    phi -= apparent_velocity.dot(surface->panel_collocation_point(panel, false));
    
    return phi;
}

/**
   Computes velocity potential time derivative at the given panel.
   
   @param[in]  surface_velocity_potentials            Current potential values.
   @param[in]  previous_surface_velocity_potentials   Previous potential values
   @param[in]  offset                                 Offset to requested Surface
   @param[in]  panel                                  Panel number.
   @param[in]  dt                                     Time step size.
   
   @returns Velocity potential time derivative.
*/ 
double
Solver::compute_surface_velocity_potential_time_derivative(int offset, int panel, double dt) const
{
    double dphidt;
    
    // Evaluate the time-derivative of the potential in a body-fixed reference frame, as in
    //   J. P. Giesing, Nonlinear Two-Dimensional Unsteady Potential Flow with Lift, Journal of Aircraft, 1968.
    if (Parameters::unsteady_bernoulli && dt > 0.0)
        dphidt = (surface_velocity_potentials(offset + panel) - previous_surface_velocity_potentials(offset + panel)) / dt;
    else
        dphidt = 0.0;
        
    return dphidt;
}

/**
   Computes the on-body gradient of a scalar field.
   
   @param[in]   scalar_field   Scalar field, ordered by panel number.
   @param[in]   surface        Surface to which the panel belongs.
   @param[in]   panel          Panel on which the on-body gradient is evaluated.
   
   @returns On-body gradient vector.
*/
Vector3d
Solver::compute_scalar_field_gradient(const Eigen::VectorXd &scalar_field, const std::shared_ptr<Body> &body, const std::shared_ptr<Surface> &surface, int panel) const
{
    // We compute the scalar field gradient by fitting a linear model.
    
    // Retrieve panel neighbors.
    vector<Body::SurfacePanelEdge> neighbors = body->panel_neighbors(surface, panel);

    // Set up a transformation such that panel normal becomes unit Z vector:
    Transform<double, 3, Affine> transformation = surface->panel_coordinate_transformation(panel);
    
    // Set up model equations:
    MatrixXd A(neighbors.size(), 2);
    VectorXd b(neighbors.size());
    
    // The model is centered on panel:
    double panel_value = scalar_field(compute_index(surface, panel));
    
    for (int i = 0; i < (int) neighbors.size(); i++) {
        Body::SurfacePanelEdge neighbor_panel = neighbors[i];
        
        // Add neighbor relative to panel:
        Vector3d neighbor_vector_normalized = transformation * neighbor_panel.surface->panel_collocation_point(neighbor_panel.panel, false);
    
        A(i, 0) = neighbor_vector_normalized(0);
        A(i, 1) = neighbor_vector_normalized(1);
    
        b(i) = scalar_field(compute_index(neighbor_panel.surface, neighbor_panel.panel)) - panel_value;
    }
    
    // Solve model equations:
    JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);
    svd.setThreshold(Parameters::zero_threshold);
    
    VectorXd model_coefficients = svd.solve(b);
    
    // Extract gradient in local frame:
    Vector3d gradient_normalized = Vector3d(model_coefficients(0), model_coefficients(1), 0.0);
    
    // Transform gradient to global frame:
    return transformation.linear().transpose() * gradient_normalized;
}

/**
   Computes the surface velocity for the given panel.
   
   @param[in]   surface   Reference surface.
   @param[it]   offset    Doublet coefficien vector offset.
   @param[in]   panel     Reference panel.
   
   @returns Surface velocity.
*/
Eigen::Vector3d
Solver::compute_surface_velocity(const std::shared_ptr<Body> &body, const std::shared_ptr<Surface> &surface, int panel) const
{
    // Compute doublet surface gradient:
    Vector3d tangential_velocity = -compute_scalar_field_gradient(doublet_coefficients, body, surface, panel);

    // Add flow due to kinematic velocity:
    Vector3d apparent_velocity = body->panel_kinematic_velocity(surface, panel) - freestream_velocity;
                                          
    tangential_velocity -= apparent_velocity;
    
    // Remove any normal velocity.  This is the (implicit) contribution of the source term.
    const Vector3d &normal = surface->panel_normal(panel);
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
Solver::compute_reference_velocity_squared(const std::shared_ptr<Body> &body) const
{
    return (body->velocity - freestream_velocity).squaredNorm();
}

/**
   Computes the pressure coefficient.
   
   @param[in]   surface_velocity   Surface velocity for the reference panel.
   @param[in]   dphidt             Time-derivative of the velocity potential for the reference panel.
   @param[in]   v_ref              Reference velocity.
   
   @returns Pressure coefficient.
*/
double
Solver::compute_pressure_coefficient(const Vector3d &surface_velocity, double dphidt, double v_ref_squared) const
{
    double C_p = 1 - (surface_velocity.squaredNorm() + 2 * dphidt) / v_ref_squared;
    
    return C_p;
}

/**
   Computes the velocity potential at the given point.
   
   @param[in]   x   Reference point.
   
   @returns Velocity potential.
*/
double
Solver::compute_velocity_potential(const Vector3d &x) const
{
    double phi = 0.0;
    
    // Iterate all non-wake surfaces:
    int offset = 0;
    
    vector<shared_ptr<Body::SurfaceData> >::const_iterator si;
    for (si = non_wake_surfaces.begin(); si != non_wake_surfaces.end(); si++) {
        const shared_ptr<Body::SurfaceData> &d = *si;

        for (int i = 0; i < d->surface->n_panels(); i++) {
            double source_influence, doublet_influence;
            
            d->surface->source_and_doublet_influence(x, i, source_influence, doublet_influence);
            
            phi += doublet_influence * doublet_coefficients(offset + i);
            phi += source_influence * source_coefficients(offset + i);
        }
        
        offset += d->surface->n_panels();
    }
    
    // Iterate wakes:
    vector<shared_ptr<BodyData> >::const_iterator bdi;
    for (bdi = bodies.begin(); bdi != bodies.end(); bdi++) {
        const shared_ptr<BodyData> &bd = *bdi;
        
        vector<shared_ptr<Body::LiftingSurfaceData> >::const_iterator lsi;
        for (lsi = bd->body->lifting_surfaces.begin(); lsi != bd->body->lifting_surfaces.end(); lsi++) {
            const shared_ptr<Body::LiftingSurfaceData> &d = *lsi;

            for (int i = 0; i < d->wake->n_panels(); i++)
                phi += d->wake->doublet_influence(x, i) * d->wake->doublet_coefficients[i];
        }
    }
                    
    // Done:
    return phi + freestream_velocity.dot(x);
}

/**
   Computes velocity at the given point, interpolating if close to the body. 
   
   The interpolation code assumes the outer angles between panels to be over 90 degrees.
   
   @param[in]   x            Reference point.
   @param[in]   ignore_set   Set of panel IDs not to interpolate for.
   
   @returns Velocity vector.
*/ 
Eigen::Vector3d
Solver::compute_velocity_interpolated(const Eigen::Vector3d &x, std::set<int> &ignore_set) const
{
    // Lists of close velocities, ordered by primacy:
    vector_aligned<Vector3d> close_panel_velocities;
    vector_aligned<Vector3d> close_panel_edge_velocities;
    vector_aligned<Vector3d> interior_close_panel_velocities;
    vector_aligned<Vector3d> interior_close_panel_edge_velocities;
    
    vector<vector_aligned<Vector3d>* > velocity_lists;
    velocity_lists.push_back(&close_panel_velocities);
    velocity_lists.push_back(&close_panel_edge_velocities);
    velocity_lists.push_back(&interior_close_panel_velocities);
    velocity_lists.push_back(&interior_close_panel_edge_velocities);
    
    // Iterate bodies:
    int offset = 0;
    
    vector<shared_ptr<BodyData> >::const_iterator bdi;
    for (bdi = bodies.begin(); bdi != bodies.end(); bdi++) {
        const shared_ptr<BodyData> &bd = *bdi;
        
        // Iterate surfaces:
        vector<shared_ptr<Body::SurfaceData> > surfaces;
        vector<shared_ptr<Body::LiftingSurfaceData> >::const_iterator lsi;
        vector<shared_ptr<Body::SurfaceData> >::const_iterator si; 
        for (si = bd->body->non_lifting_surfaces.begin(); si != bd->body->non_lifting_surfaces.end(); si++)
            surfaces.push_back(*si);
        for (lsi = bd->body->lifting_surfaces.begin(); lsi != bd->body->lifting_surfaces.end(); lsi++)
            surfaces.push_back(*lsi);
        
        for (si = surfaces.begin(); si != surfaces.end(); si++) {
            const shared_ptr<Body::SurfaceData> &d = *si;

            for (int i = 0; i < d->surface->n_panels(); i++) {
                // Ignore the given set of panels.
                if (ignore_set.find(i) != ignore_set.end())
                    continue;
                    
                // Transform the point 'x' into the panel coordinate system:
                Vector3d x_transformed = d->surface->panel_coordinate_transformation(i) * x;
                
                // Are we in the exterior, relative to the panel?
                bool in_exterior;
                if (x_transformed(2) < Parameters::zero_threshold)
                    in_exterior = true;
                else
                    in_exterior = false;
   
                // Compute normal distance of the point 'x' from panel:
                double normal_distance = fabs(x_transformed(2));
                
                // We have three zones: 
                // The boundary layer, followed by the interpolation layer, followed by the rest of the control volume.
                double boundary_layer_thickness      = bd->boundary_layer->thickness(d->surface, i);
                double interpolation_layer_thickness = Parameters::interpolation_layer_thickness;
                double total_thickness               = boundary_layer_thickness + interpolation_layer_thickness;
                
                // Are we inside one of the first two layers?
                if (normal_distance < total_thickness) {
                    // Yes.  Check whether A) the point projection lies inside the panel, and whether B) we are close to one of the panel's edges:
                    bool projection_in_panel = true;
                    double panel_edge_distance = total_thickness;
                    Vector3d panel_to_point_direction(0, 0, 0);
                    for (int l = 0; l < (int) d->surface->panel_nodes[i].size(); l++) {
                        int next_l;
                        if (l == (int) d->surface->panel_nodes[i].size() - 1)
                            next_l = 0;
                        else
                            next_l = l + 1;
                            
                        Vector3d point_a = d->surface->panel_transformed_points[i][l];
                        Vector3d point_b = d->surface->panel_transformed_points[i][next_l];
                            
                        Vector3d edge = point_b - point_a;
                        Vector3d normal(-edge(1), edge(0), 0.0);
                        normal.normalize();

                        // We are above the panel if the projection lies inside all four panel edges:
                        double normal_component = (x_transformed - point_a).dot(normal);
                        if (normal_component <= 0)
                            projection_in_panel = false;
                            
                        double edge_distance = sqrt(pow(normal_component, 2) + pow(x_transformed(2), 2));
                        if (edge_distance < panel_edge_distance) {
                            // Does the point lie beside the panel edge?
                            if (edge.dot(x_transformed - point_a) * edge.dot(x_transformed - point_b) < 0) {
                                panel_edge_distance = edge_distance;
                                if (edge_distance > 0)
                                    panel_to_point_direction = (normal * normal_component + Vector3d(0, 0, x_transformed(2))).normalized();
                            }
                                
                            // Is the point close to the panel vertex?
                            Vector3d delta = x_transformed - point_a;
                            double node_distance = delta.norm();
                            if (node_distance < panel_edge_distance) {
                                panel_edge_distance = node_distance;
                                if (node_distance > 0)
                                    panel_to_point_direction = delta.normalized();
                            }
                        }
                    }
                    
                    // Compute distance to panel:
                    double panel_distance;
                    Vector3d to_point_direction;
                    if (projection_in_panel) {
                        panel_distance = normal_distance;
                        to_point_direction = Vector3d(0, 0, -1);
                    } else {
                        panel_distance = panel_edge_distance;
                        to_point_direction = panel_to_point_direction;
                    }
                    
                    // Are we close to the panel?
                    if (panel_distance < total_thickness) {                            
                        // Yes. 
                        Vector3d velocity;
                        
                        if (panel_distance < boundary_layer_thickness) {
                            // We are in the boundary layer:
                            velocity = bd->boundary_layer->velocity(d->surface, i, panel_distance);
                                
                        } else if (panel_distance > 0) {
                            // We are in the interpolation layer.
                            // Interpolate between the surface velocity, and the velocity away from the body:
                            Vector3d lower_velocity = surface_velocity(d->surface, i);
   
                            // This point lies in the control volume only, if A) no other body lies in the way, and B) the exterior angles\ are more than 90 degrees each.
                            Vector3d upper_point_transformed = x_transformed + (total_thickness - panel_distance) * to_point_direction;
                            Vector3d upper_point = d->surface->panel_coordinate_transformation(i).inverse() * upper_point_transformed;
                            
                            Vector3d upper_velocity;
                            if (in_exterior) {
                                // Compute the upper velocity again using interpolation, in case we are now close to another panel.  This can happen in concave corners.
                                // We must take care, however, to avoid the possibility of an infinite loop.
                                set<int> ignore_set_copy(ignore_set);
                                ignore_set_copy.insert(i);
                                upper_velocity = compute_velocity_interpolated(upper_point, ignore_set_copy);
                                
                            } else {
                                // In the interior, we have the undisturbed freestream velocity.
                                upper_velocity = freestream_velocity;
                                
                            }
                                
                            // Interpolate:
                            double interpolation_distance = panel_distance - boundary_layer_thickness;
                            velocity = (interpolation_distance * upper_velocity + (interpolation_layer_thickness - interpolation_distance) * lower_velocity) 
                                        / interpolation_layer_thickness;
                        } else {
                            // We are on the panel.  Use surface velocity:
                            velocity = surface_velocity(d->surface, i);
                            
                        }
                        
                        // Store interpolated velocity.  We cannot return here.  In concave corners, a point may be close to more than one panel.
                        if (in_exterior) {
                            if (projection_in_panel)
                                close_panel_velocities.push_back(velocity);
                            else
                                close_panel_edge_velocities.push_back(velocity);
                        } else {
                            if (projection_in_panel)
                                interior_close_panel_velocities.push_back(velocity);
                            else
                                interior_close_panel_edge_velocities.push_back(velocity);
                        }
                    }
                }
            }
            
            offset += d->surface->n_panels();
        }
    }
    
    // Are we close to any panels?   
    for (int i = 0; i < (int) velocity_lists.size(); i++) {
        // Yes, at primacy level i:
        if (velocity_lists[i]->size() > 0) {
            // Average:
            Vector3d velocity = Vector3d(0, 0, 0);
            vector_aligned<Vector3d>::iterator it;
            for (it = velocity_lists[i]->begin(); it != velocity_lists[i]->end(); it++)
                velocity += *it;
                
            // Normalize velocity.  The weights sum up to (n - 1) times the sum of the distances.
            velocity /= velocity_lists[i]->size();
            
            return velocity;
        }
    }

    // No close panels.  Compute potential velocity:
    return compute_velocity(x);
}

/**
   Computes the potential velocity at the given point.
   
   @param[in]   x   Reference point.
   
   @returns Potential velocity vector.
*/ 
Eigen::Vector3d
Solver::compute_velocity(const Eigen::Vector3d &x) const
{
    Vector3d velocity = Vector3d(0, 0, 0);
    
    int offset = 0;

    // Iterate bodies:
    vector<shared_ptr<BodyData> >::const_iterator bdi;
    for (bdi = bodies.begin(); bdi != bodies.end(); bdi++) {
        const shared_ptr<BodyData> &bd = *bdi;
        
        // Iterate surfaces:
        vector<shared_ptr<Body::SurfaceData> > surfaces;
        vector<shared_ptr<Body::LiftingSurfaceData> >::const_iterator lsi;
        vector<shared_ptr<Body::SurfaceData> >::const_iterator si; 
        for (si = bd->body->non_lifting_surfaces.begin(); si != bd->body->non_lifting_surfaces.end(); si++)
            surfaces.push_back(*si);
        for (lsi = bd->body->lifting_surfaces.begin(); lsi != bd->body->lifting_surfaces.end(); lsi++)
            surfaces.push_back(*lsi);
        
        for (si = surfaces.begin(); si != surfaces.end(); si++) {
            const shared_ptr<Body::SurfaceData> &d = *si;

            for (int i = 0; i < d->surface->n_panels(); i++) {
                // If no close panels were detected so far, add the influence of this panel:
                velocity += d->surface->vortex_ring_unit_velocity(x, i) * doublet_coefficients(offset + i);
                velocity += d->surface->source_unit_velocity(x, i) * source_coefficients(offset + i);
            }
            
            offset += d->surface->n_panels();
        }
        
        // If no close panels were detected so far, add the influence of the wakes:
        for (lsi = bd->body->lifting_surfaces.begin(); lsi != bd->body->lifting_surfaces.end(); lsi++) {
            const shared_ptr<Body::LiftingSurfaceData> &d = *lsi;
            
            if (d->wake->n_panels() >= d->lifting_surface->n_spanwise_panels()) {
                for (int i = 0; i < d->wake->n_panels(); i++)
                    velocity += d->wake->vortex_ring_unit_velocity(x, i) * d->wake->doublet_coefficients[i];
            }
        }
    }
    
    // Done:
    return velocity + freestream_velocity;
}

/**
   Computes the vector by which the first wake vortex is offset from the trailing edge.
   
   @param[in]   body              Reference body.
   @param[in]   lifting_surface   Reference lifting surface.
   @param[in]   index             Trailing edge index.
   @param[in]   dt                Time step size.
   
   @returns The trailing edge vortex displacement.
*/
Eigen::Vector3d
Solver::compute_trailing_edge_vortex_displacement(const std::shared_ptr<Body> &body, const std::shared_ptr<LiftingSurface> &lifting_surface, int index, double dt) const
{
    Vector3d apparent_velocity = body->node_kinematic_velocity(lifting_surface, lifting_surface->trailing_edge_node(index)) - freestream_velocity;
                    
    Vector3d wake_emission_velocity = lifting_surface->wake_emission_velocity(apparent_velocity, index);
    
    return Parameters::wake_emission_distance_factor * wake_emission_velocity * dt;   
}

/**
   Computes the result index for a given (surface, panel)-pair.
   
   @param[in]   surface   Reference surface.
   @param[in]   panel     Reference panel.
   
   @returns The index.
*/
int
Solver::compute_index(const std::shared_ptr<Surface> &surface, int panel) const
{
    int offset = 0;
    vector<shared_ptr<Body::SurfaceData> >::const_iterator si;
    for (si = non_wake_surfaces.begin(); si != non_wake_surfaces.end(); si++) {
        const shared_ptr<Body::SurfaceData> &d = *si;
        
        if (surface == d->surface)
            return offset + panel;
        
        offset += d->surface->n_panels();
    }
    
    return -1;
}
