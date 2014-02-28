//
// Vortexje -- Surface.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>
#include <cmath>

#include <Eigen/Geometry>
#include <Eigen/SVD>

#include <vortexje/surface.hpp>
#include <vortexje/parameters.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

// Static counter to give every surface a unique ID.
static int id_counter = 0;

// Avoid having to divide by 4 pi all the time:
static const double one_over_4pi = 1.0 / (4 * M_PI);

/**
   Constructs an empty surface.
*/
Surface::Surface()
{
    // Set ID:
    id = ++id_counter;
}

/**
   Clears the node-panel neighbor data structure, and frees up its memory.
*/
void
Surface::clear_node_panel_neighbors()
{
    vector<vector<int> *> unique;
    
    for (int i = 0; i < (int) node_panel_neighbors.size(); i++) {
        bool found = false;
        for (int j = 0; j < (int) unique.size(); j++) {
            if (node_panel_neighbors[i] == unique[j]) {
                found = true;
                break;
            }
        }
        
        if (!found)
            unique.push_back(node_panel_neighbors[i]);
    }
    
    for (int i = 0; i < (int) unique.size(); i++)
        delete unique[i];
        
    node_panel_neighbors.clear();
}

/**
   Destructor.
*/
Surface::~Surface()
{
    clear_node_panel_neighbors();
}

/**
   Adds a triangle to this surface, following the Gmsh orientation convention.
    
   @param[in]   node_a   Node number to form vertex A.
   @param[in]   node_b   Node number to form vertex B.
   @param[in]   node_c   Node number to form vertex C.
    
   @returns New panel number.
    
   @note This method does not update the panel-panel neighbor data structure.  One must call
   compute_topology() as well as compute_geometry() when done adding panels.
*/
int
Surface::add_triangle(int node_a, int node_b, int node_c)
{
    vector<int> single_panel_nodes;
    single_panel_nodes.push_back(node_a);
    single_panel_nodes.push_back(node_b);
    single_panel_nodes.push_back(node_c);
    
    int panel_id = panel_nodes.size();
    panel_nodes.push_back(single_panel_nodes);
    
    node_panel_neighbors[node_a]->push_back(panel_id);
    node_panel_neighbors[node_b]->push_back(panel_id);
    node_panel_neighbors[node_c]->push_back(panel_id);
    
    return panel_id;
}

/**
   Adds a quadrangle to this surface, following the gmsh orientation convention.
    
   @param[in]   node_a   Node number to form vertex A.
   @param[in]   node_b   Node number to form vertex B.
   @param[in]   node_c   Node number to form vertex C.
   @param[in]   node_d   Node number to form vertex D.
    
   @returns New panel number.
    
   @note This method does not update the panel-panel neighbor data structure.  One must call
   compute_panel_neighbors() when done adding panels.
*/
int
Surface::add_quadrangle(int node_a, int node_b, int node_c, int node_d)
{
    vector<int> single_panel_nodes;
    single_panel_nodes.push_back(node_a);
    single_panel_nodes.push_back(node_b);
    single_panel_nodes.push_back(node_c);       
    single_panel_nodes.push_back(node_d);
    
    int panel_id = panel_nodes.size();
    panel_nodes.push_back(single_panel_nodes);
    
    node_panel_neighbors[node_a]->push_back(panel_id);
    node_panel_neighbors[node_b]->push_back(panel_id);
    node_panel_neighbors[node_c]->push_back(panel_id);
    node_panel_neighbors[node_d]->push_back(panel_id);
    
    return panel_id;
}

/**
   Computes neighboring panels of panels, based on existing node-panel data structures.
*/
void
Surface::compute_topology()
{   
    // Compute panel neighbors:
    for (int i = 0; i < (int) panel_nodes.size(); i++) {
        vector<int> single_panel_potential_neighbors;
        vector<int> single_panel_neighbors;
        
        for (int j = 0; j < (int) panel_nodes[i].size(); j++) {
            int node = panel_nodes[i][j];
            for (int k = 0; k < (int) node_panel_neighbors[node]->size(); k++) {
                int potential_neighbor = (*node_panel_neighbors[node])[k];
                if (potential_neighbor == i)
                    continue;
                
                bool found = false;
                for (int l = 0; l < (int) single_panel_potential_neighbors.size(); l++) {
                    if (potential_neighbor == single_panel_potential_neighbors[l]) {
                        found = true;
                        break;
                    }
                }
                
                // Must have two nodes in common.
                if (found)
                    single_panel_neighbors.push_back(potential_neighbor);
                else
                    single_panel_potential_neighbors.push_back(potential_neighbor);
            }
        }
        
        panel_neighbors.push_back(single_panel_neighbors);
    }
}


/**
   Computes the normals, collocation points, surface areas, and diameters of all panels.
*/
void
Surface::compute_geometry()
{
    // Normals:
    cout << "Surface " << id << ": Generating panel normals." << endl;
        
    panel_normals.clear();
    panel_normals.reserve(n_panels());
    
    for (int i = 0; i < n_panels(); i++) {
        vector<int> &single_panel_nodes = panel_nodes[i];
        
        Vector3d normal;
        if (single_panel_nodes.size() == 3) {
            Vector3d AB = nodes[single_panel_nodes[1]] - nodes[single_panel_nodes[0]];
            Vector3d AC = nodes[single_panel_nodes[2]] - nodes[single_panel_nodes[0]];
            
            normal = AB.cross(AC);
            
        } else { // 4 sides
            Vector3d AC = nodes[single_panel_nodes[2]] - nodes[single_panel_nodes[0]];
            Vector3d BD = nodes[single_panel_nodes[3]] - nodes[single_panel_nodes[1]];
            
            normal = AC.cross(BD);
        }

        normal.normalize();

        panel_normals.push_back(normal);
    }
    
    // Collocation points: 
    cout << "Surface " << id << ": Generating panel collocation points." << endl;
    
    panel_collocation_points[0].clear();
    panel_collocation_points[0].reserve(n_panels());
    
    panel_collocation_points[1].clear();
    panel_collocation_points[1].reserve(n_panels());

    for (int i = 0; i < n_panels(); i++) {
        vector<int> &single_panel_nodes = panel_nodes[i];
        
        Vector3d collocation_point(0, 0, 0);
        for (int j = 0; j < (int) single_panel_nodes.size(); j++)
            collocation_point = collocation_point + nodes[single_panel_nodes[j]];

        collocation_point = collocation_point / single_panel_nodes.size();
            
        panel_collocation_points[0].push_back(collocation_point);
        
        Vector3d below_surface_collocation_point = collocation_point + Parameters::collocation_point_delta * panel_normal(i);
        panel_collocation_points[1].push_back(below_surface_collocation_point);
    }
    
    // Coordinate transformations:
    cout << "Surface " << id << ": Generating panel coordinate transformations." << endl;
        
    panel_coordinate_transformations.clear();
    panel_coordinate_transformations.reserve(n_panels());
    
    panel_transformed_points.clear();
    panel_transformed_points.reserve(n_panels());
    
    for (int i = 0; i < n_panels(); i++) {
        vector<int> &single_panel_nodes = panel_nodes[i];
        
        Vector3d AB = nodes[single_panel_nodes[1]] - nodes[single_panel_nodes[0]];
        AB.normalize();
        
        const Vector3d &normal = panel_normal(i);
        
        Matrix3d rotation;
        rotation.block<1, 3>(0, 0) = AB;
        rotation.block<1, 3>(1, 0) = normal.cross(AB);
        rotation.block<1, 3>(2, 0) = normal;
        
        Transform<double, 3, Affine> transformation = rotation * Translation<double, 3>(-panel_collocation_point(i, false));

        panel_coordinate_transformations.push_back(transformation);
        
        // Create transformed points.
        vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > single_panel_transformed_points;
        single_panel_transformed_points.reserve(single_panel_nodes.size());
        for (int j = 0; j < (int) single_panel_nodes.size(); j++)
            single_panel_transformed_points.push_back(transformation * nodes[single_panel_nodes[j]]);
        panel_transformed_points.push_back(single_panel_transformed_points);
    }
    
    // Surface areas:
    cout << "Surface " << id << ": Generating panel surface area cache." << endl;
    
    panel_surface_areas.clear();
    panel_surface_areas.reserve(n_panels());
    
    for (int i = 0; i < n_panels(); i++) {
        vector<int> &single_panel_nodes = panel_nodes[i];
        
        double surface_area = 0.0;
        if (single_panel_nodes.size() == 3) {
            Vector3d AB = nodes[single_panel_nodes[1]] - nodes[single_panel_nodes[0]];
            Vector3d AC = nodes[single_panel_nodes[2]] - nodes[single_panel_nodes[0]];
            
            surface_area = 0.5 * AB.cross(AC).norm();
            
        } else { // 4 sides
            Vector3d AC = nodes[single_panel_nodes[2]] - nodes[single_panel_nodes[0]];
            Vector3d BD = nodes[single_panel_nodes[3]] - nodes[single_panel_nodes[1]];
            
            surface_area = 0.5 * AC.cross(BD).norm();
        }
        
        panel_surface_areas.push_back(surface_area);
    }
    
    // Diameters:
    cout << "Surface " << id << ": Generating panel diameter cache." << endl;
    
    panel_diameters.clear();
    panel_diameters.reserve(n_panels());
    
    for (int i = 0; i < n_panels(); i++) {
        double diameter = numeric_limits<double>::min();
        
        for (int j = 0; j < (int) panel_nodes[i].size(); j++) {
            Vector3d a = nodes[panel_nodes[i][j]];
            
            for (int k = 0; k < j; k++) {
                Vector3d b = nodes[panel_nodes[i][k]];
                
                double diameter_candidate = (b - a).norm();
                if (diameter_candidate > diameter)
                    diameter = diameter_candidate;
            }
        }
        
        panel_diameters.push_back(diameter);
    }
}

/**
   Returns the number of nodes contained in this surface.
   
   @returns Number of nodes.
*/
int
Surface::n_nodes() const
{
    return nodes.size();
}

/**
   Returns the number of panels contained in this surface.
   
   @returns Number of panels.
*/
int
Surface::n_panels() const
{
    return panel_nodes.size();
}

/**
   Rotates this surface.
   
   @param[in]   axis    Axis of rotation.
   @param[in]   angle   Angle of rotation.
*/
void
Surface::rotate(const Eigen::Vector3d &axis, double angle)
{
    transform(AngleAxis<double>(angle, axis).toRotationMatrix());
}

/**
   Transforms this surface.
   
   @param[in]   transformation   Transformation matrix.
*/
void
Surface::transform(const Eigen::Matrix3d &transformation)
{
    transform(Transform<double, 3, Eigen::Affine>(transformation));
}

/**
   Transforms this surface.
   
   @param[in]   transformation   Affine transformation.
*/
void
Surface::transform(const Eigen::Transform<double, 3, Eigen::Affine> &transformation)
{
    for (int i = 0; i < n_nodes(); i++)
        nodes[i] = transformation * nodes[i];   
            
    for (int j = 0; j < 2; j++)
        for (int i = 0; i < n_panels(); i++)
            panel_collocation_points[j][i] = transformation * panel_collocation_points[j][i];
    
    for (int i = 0; i < n_panels(); i++)
        panel_normals[i] = transformation * panel_normals[i];
        
    for (int i = 0; i < n_panels(); i++)
        panel_coordinate_transformations[i] = panel_coordinate_transformations[i] * transformation.inverse();
}

/**
   Translates this surface.
   
   @param[in]   translation   Translation vector.
*/
void
Surface::translate(const Eigen::Vector3d &translation)
{
    for (int i = 0; i < n_nodes(); i++)
        nodes[i] = nodes[i] + translation;
            
    for (int j = 0; j < 2; j++)
        for (int i = 0; i < n_panels(); i++)
            panel_collocation_points[j][i] = panel_collocation_points[j][i] + translation;
            
    for (int i = 0; i < n_panels(); i++)
        panel_coordinate_transformations[i] = panel_coordinate_transformations[i] * Translation<double, 3>(-translation);
}

/**
   Returns the collocation point of the given panel.
   
   @param[in]   panel           Panel of which the collocation point is returned.
   @param[in]   below_surface   true to request the collocation point lying underneath the surface.
   
   @returns Collocation point.
*/
const Vector3d &
Surface::panel_collocation_point(int panel, bool below_surface) const
{
    return panel_collocation_points[below_surface][panel];
}

/**
   Returns the inward-pointing normal of the given panel.
   
   @param[in]   panel   Panel of which the inward-pointing normal is returned.
   
   @returns Inward-pointing normal.
*/
const Vector3d &
Surface::panel_normal(int panel) const
{
    return panel_normals[panel];
}

/**
   Returns the panel coordinate transformation for the given panel.
   
   @param[in]   panel   Panel of which the coordinate transformation is returned.
   
   @returns Panel coordinate transformation.
*/
const Transform<double, 3, Affine> &
Surface::panel_coordinate_transformation(int panel) const
{
    return panel_coordinate_transformations[panel];
}

/**
   Returns the surface area of the given panel.
   
   @param[in]   panel   Panel of which the surface area is returned.
   
   @returns Panel surface area.
*/
double
Surface::panel_surface_area(int panel) const
{  
    return panel_surface_areas[panel];
}

/**
   Returns the diameter of the given panel.
   
   @param[in]   panel   Panel of which the diameter is returned.
   
   @returns Panel diameter.
*/
double
Surface::panel_diameter(int panel) const
{
    return panel_diameters[panel];
}

/**
   Computes the on-body gradient of a scalar field.
   
   @param[in]   scalar_field   Scalar field, ordered by panel number.
   @param[in]   offset         Scalar field offset
   @param[in]   this_panel     Panel on which the on-body gradient is evaluated.
   
   @returns On-body gradient.
*/
Vector3d
Surface::scalar_field_gradient(const Eigen::VectorXd &scalar_field, int offset, int this_panel) const
{
    // We compute the scalar field gradient by fitting a linear model.
    const Vector3d &this_normal = panel_normal(this_panel);

    // Set up a transformation such that panel normal becomes unit Z vector:
    Transform<double, 3, Affine> transformation = panel_coordinate_transformation(this_panel);
    
    // Set up model equations:
    MatrixXd A(panel_neighbors[this_panel].size() + 1, 3);
    VectorXd b(panel_neighbors[this_panel].size() + 1);
    
    // The model is centered on this_panel:
    A(0, 0) = 0.0;
    A(0, 1) = 0.0;
    A(0, 2) = 1.0;
    b(0) = scalar_field(offset + this_panel);
    
    for (int i = 0; i < (int) panel_neighbors[this_panel].size(); i++) {
        int neighbor_panel = panel_neighbors[this_panel][i];
        
        const Vector3d &neighbor_normal = panel_normal(neighbor_panel);
        if (this_normal.dot(neighbor_normal) >= 0) {
            // Add neighbor relative to this_panel:
            Vector3d neighbor_vector_normalized = transformation * panel_collocation_point(neighbor_panel, false);
        
            A(i + 1, 0) = neighbor_vector_normalized(0);
            A(i + 1, 1) = neighbor_vector_normalized(1);
            A(i + 1, 2) = 1.0;
        
            b(i + 1) = scalar_field(offset + neighbor_panel);
        } else {
            // Don't differentiate along sharp edges, such as the trailing edge.
            A(i + 1, 0) = A(i + 1, 1) = A(i + 1, 2) = 0.0;
        }
    }
    
    // Solve model equations:
    JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);
    
    VectorXd model_coefficients = svd.solve(b);
    
    // Extract gradient in local frame:
    Vector3d gradient_normalized = Vector3d(model_coefficients(0), model_coefficients(1), 0.0);
    
    // Transform gradient to global frame:
    return transformation.linear().transpose() * gradient_normalized;
}

// Simultaneously compute influence of source and doublet panel edges on given point.
static void
source_and_doublet_edge_influence(const Vector3d &x, const Vector3d &node_a, const Vector3d &node_b, double *source_edge_influence, double *doublet_edge_influence)
{
    double d = sqrt(pow(node_b(0) - node_a(0), 2) + pow(node_b(1) - node_a(1), 2));

    if (d < Parameters::inversion_tolerance) {
        if (source_edge_influence != NULL)
            *source_edge_influence = 0.0;
        if (doublet_edge_influence != NULL)
            *doublet_edge_influence = 0.0;
            
        return;
    }
        
    double z = x(2);
    
    double m = (node_b(1) - node_a(1)) / (node_b(0) - node_a(0));
    
    double e1 = pow(x(0) - node_a(0), 2) + pow(z, 2);
    double e2 = pow(x(0) - node_b(0), 2) + pow(z, 2);
    
    double r1 = sqrt(e1 + pow(x(1) - node_a(1), 2));
    double r2 = sqrt(e2 + pow(x(1) - node_b(1), 2));
    
    double h1 = (x(0) - node_a(0)) * (x(1) - node_a(1));
    double h2 = (x(0) - node_b(0)) * (x(1) - node_b(1));
    
    // IEEE-754 floating point division by zero results in +/- inf, and atan(inf) = pi / 2.
    double u = (m * e1 - h1) / (z * r1);
    double v = (m * e2 - h2) / (z * r2);
    
    double delta_theta;
    if (u == v)
        delta_theta = 0.0;
    else
        delta_theta = atan2(u - v, 1 + u * v);
    
    if (source_edge_influence != NULL)
        *source_edge_influence  = ((x(0) - node_a(0)) * (node_b(1) - node_a(1)) - (x(1) - node_a(1)) * (node_b(0) - node_a(0))) / d * log((r1 + r2 + d) / (r1 + r2 - d)) - fabs(z) * delta_theta; 
    if (doublet_edge_influence != NULL)
        *doublet_edge_influence = delta_theta;
}

/**
   Simultaneously computes the potential influences induced by source and doublet panels of unit strength.  
   
   @param[in]   x                   Point at which the influence coefficient is evaluated.
   @param[in]   this_panel          Panel on which the doublet panel is located.
   @param[out]  source_influence    Source influence value.
   @param[out]  doublet_influence   Doublet influence value.
*/
void
Surface::source_and_doublet_influence(const Eigen::Vector3d &x, int this_panel, double &source_influence, double &doublet_influence) const
{
    // Transform such that panel normal becomes unit Z vector:    
    const Transform<double, 3, Affine> &transformation = panel_coordinate_transformation(this_panel);
    
    Vector3d x_normalized = transformation * x;
    
    // Compute influence coefficient according to Hess:
    source_influence  = 0.0;
    doublet_influence = 0.0;
    
    for (int i = 0; i < (int) panel_nodes[this_panel].size(); i++) {
        int prev_idx;
        if (i == 0)
            prev_idx = panel_nodes[this_panel].size() - 1;
        else
            prev_idx = i - 1;
            
        const Vector3d &node_a = panel_transformed_points[this_panel][prev_idx];
        const Vector3d &node_b = panel_transformed_points[this_panel][i];
        
        double source_edge_influence, doublet_edge_influence;
        
        source_and_doublet_edge_influence(x_normalized, node_a, node_b, &source_edge_influence, &doublet_edge_influence);
        
        source_influence  += source_edge_influence;
        doublet_influence += doublet_edge_influence;
    }   
    
    source_influence  *= -one_over_4pi;
    doublet_influence *=  one_over_4pi;
}

/**
   Computes the potential influence induced by a source panel of unit strength.  
   
   @param[in]   x            Point at which the influence coefficient is evaluated.
   @param[in]   this_panel   Panel on which the source panel is located.
   
   @returns Influence coefficient.
*/
double
Surface::source_influence(const Eigen::Vector3d &x, int this_panel) const
{
    // Transform such that panel normal becomes unit Z vector:    
    const Transform<double, 3, Affine> &transformation = panel_coordinate_transformation(this_panel);
    
    Vector3d x_normalized = transformation * x;
    
    // Compute influence coefficient according to Hess:
    double influence = 0.0;
    for (int i = 0; i < (int) panel_nodes[this_panel].size(); i++) {
        int prev_idx;
        if (i == 0)
            prev_idx = panel_nodes[this_panel].size() - 1;
        else
            prev_idx = i - 1;
            
        const Vector3d &node_a = panel_transformed_points[this_panel][prev_idx];
        const Vector3d &node_b = panel_transformed_points[this_panel][i];
        
        double edge_influence;
        source_and_doublet_edge_influence(x_normalized, node_a, node_b, &edge_influence, NULL);
        
        influence += edge_influence;
    }   
    
    return -one_over_4pi * influence;
}

/**
   Computes the potential influence induced by a doublet panel of unit strength.  
   
   @param[in]   x            Point at which the influence coefficient is evaluated.
   @param[in]   this_panel   Panel on which the doublet panel is located.
   
   @returns Influence coefficient.
*/
double
Surface::doublet_influence(const Eigen::Vector3d &x, int this_panel) const
{
    // Transform such that panel normal becomes unit Z vector:
    const Transform<double, 3, Affine> &transformation = panel_coordinate_transformation(this_panel);
    
    Vector3d x_normalized = transformation * x;
    
    // Compute influence coefficient according to Hess:
    double influence = 0.0;
    for (int i = 0; i < (int) panel_nodes[this_panel].size(); i++) {
        int prev_idx;
        if (i == 0)
            prev_idx = panel_nodes[this_panel].size() - 1;
        else
            prev_idx = i - 1;
            
        const Vector3d &node_a = panel_transformed_points[this_panel][prev_idx];
        const Vector3d &node_b = panel_transformed_points[this_panel][i];
        
        double edge_influence;
        source_and_doublet_edge_influence(x_normalized, node_a, node_b, NULL, &edge_influence);
        
        influence += edge_influence;
    }
    
    return one_over_4pi * influence;
}

// Compute velocity induced by an edge of a source panel:
static Vector3d
source_edge_unit_velocity(const Vector3d &x, const Vector3d &node_a, const Vector3d &node_b)
{   
    double d = sqrt(pow(node_b(0) - node_a(0), 2) + pow(node_b(1) - node_a(1), 2));
    
    if (d < Parameters::inversion_tolerance)
        return Vector3d(0, 0, 0);
        
    double z = x(2);
    
    double m = (node_b(1) - node_a(1)) / (node_b(0) - node_a(0));
    
    double e1 = pow(x(0) - node_a(0), 2) + pow(z, 2);
    double e2 = pow(x(0) - node_b(0), 2) + pow(z, 2);
    
    double r1 = sqrt(e1 + pow(x(1) - node_a(1), 2));
    double r2 = sqrt(e2 + pow(x(1) - node_b(1), 2));
    
    double h1 = (x(0) - node_a(0)) * (x(1) - node_a(1));
    double h2 = (x(0) - node_b(0)) * (x(1) - node_b(1));
    
    // IEEE-754 floating point division by zero results in +/- inf, and atan(inf) = pi / 2.
    double u = (m * e1 - h1) / (z * r1);
    double v = (m * e2 - h2) / (z * r2);
    
    double delta_theta;
    if (u == v)
        delta_theta = 0.0;
    else
        delta_theta = atan2(u - v, 1 + u * v);
        
    double l = log((r1 + r2 - d) / (r1 + r2 + d));
    
    return Vector3d((node_b(1) - node_a(1)) / d * l,
                    (node_a(0) - node_b(0)) / d * l,
                    delta_theta);
}

/**
   Computes the velocity induced by a source panel of unit strength.  
   
   @param[in]   x            Point at which the velocity is evaluated.
   @param[in]   this_panel   The panel on which the vortex ring is located.
   
   @returns Velocity induced by the source panel.
*/
Vector3d
Surface::source_unit_velocity(const Eigen::Vector3d &x, int this_panel) const
{   
    // Transform such that panel normal becomes unit Z vector:
    const Transform<double, 3, Affine> &transformation = panel_coordinate_transformation(this_panel);
    
    Vector3d x_normalized = transformation * x;
    
    // Compute influence coefficient according to Hess:
    Vector3d velocity(0, 0, 0);
    for (int i = 0; i < (int) panel_nodes[this_panel].size(); i++) {
        int prev_idx;
        if (i == 0)
            prev_idx = panel_nodes[this_panel].size() - 1;
        else
            prev_idx = i - 1;
            
        const Vector3d &node_a = panel_transformed_points[this_panel][prev_idx];
        const Vector3d &node_b = panel_transformed_points[this_panel][i];
        
        velocity += source_edge_unit_velocity(x_normalized, node_a, node_b);
    }   
    
    // Transform back:
    velocity = transformation.linear().transpose() * velocity;
    
    // Done:
    return one_over_4pi * velocity;
}

/**
   Computes the velocity induced by a vortex ring of unit strength.
   
   @param[in]   x            Point at which the velocity is evaluated.
   @param[in]   this_panel   Panel on which the vortex ring is located.
   
   @returns Velocity induced by the vortex ring.
*/
Vector3d
Surface::vortex_ring_unit_velocity(const Eigen::Vector3d &x, int this_panel) const
{    
    Vector3d velocity(0, 0, 0);
    
    for (int i = 0; i < (int) panel_nodes[this_panel].size(); i++) {
        int previous_idx;
        if (i == 0)
            previous_idx = panel_nodes[this_panel].size() - 1;
        else
            previous_idx = i - 1;
            
        const Vector3d &node_a = nodes[panel_nodes[this_panel][previous_idx]];
        const Vector3d &node_b = nodes[panel_nodes[this_panel][i]];
        
        Vector3d r_0 = node_b - node_a;
        Vector3d r_1 = node_a - x;
        Vector3d r_2 = node_b - x;
        
        double r_1_norm = r_1.norm();
        double r_2_norm = r_2.norm();
        
        Vector3d r_1xr_2 = r_1.cross(r_2);
        double r_1xr_2_sqnorm = r_1xr_2.squaredNorm();
        
        if (r_1_norm < Parameters::inversion_tolerance ||
            r_2_norm < Parameters::inversion_tolerance ||
            r_1xr_2_sqnorm < Parameters::inversion_tolerance)
            continue;

        velocity += r_1xr_2 / r_1xr_2_sqnorm * r_0.dot(r_1 / r_1_norm - r_2 / r_2_norm);
    }

    return one_over_4pi * velocity;
}

/**
   Computes the potential influence induced by a doublet panel of unit strength.  If the influence
   has been computed before, the cached influence coefficient is returned.  If not, the coefficient
   is computed and cached.
   
   @param[in]   other         Surface on which the influence coefficient is evaluated.
   @param[in]   other_panel   Panel on which the influence coefficient is evaluated.
   @param[in]   this_panel    Panel on which the doublet panel is located.
   
   @returns Influence coefficient.
*/
double
Surface::doublet_influence(const Surface &other, int other_panel, int this_panel) const
{ 
    if ((this == &other) && (this_panel == other_panel))
        return -0.5;
    else
        return doublet_influence(other.panel_collocation_point(other_panel, true), this_panel);
}

/**
   Computes the potential influence induced by a source panel of unit strength.  If the influence
   has been computed before, the cached influence coefficient is returned.  If not, the coefficient
   is computed and cached.
   
   @param[in]   other         Surface on which the influence coefficient is evaluated.
   @param[in]   other_panel   Panel on which the influence coefficient is evaluated.
   @param[in]   this_panel    Panel on which the source panel is located.
   
   @returns Influence coefficient.
*/
double
Surface::source_influence(const Surface &other, int other_panel, int this_panel) const
{
    return source_influence(other.panel_collocation_point(other_panel, true), this_panel);
}

/**
   Simultaneously computes the potential influences induced by source and doublet panels of unit strength.  
   
   @param[in]   other               Surface on which the velocity is evaluated.
   @param[in]   other_panel         Panel on which the velocity is evaluated.
   @param[in]   this_panel          Panel on which the doublet panel is located.
   @param[out]  source_influence    Source influence value.
   @param[out]  doublet_influence   Doublet influence value.
*/
void
Surface::source_and_doublet_influence(const Surface &other, int other_panel, int this_panel, double &source_influence, double &doublet_influence) const
{
    if ((this == &other) && (this_panel == other_panel)) {
        doublet_influence = -0.5;
        
        source_influence = this->source_influence(other, other_panel, this_panel);
        
    } else
        source_and_doublet_influence(other.panel_collocation_point(other_panel, true), this_panel, source_influence, doublet_influence);
}

/**
   Computes the velocity induced by a source panel of unit strength.  If the velocity
   has been computed before, the cached velocity is returned.  If not, the quantity
   is computed and cached.
   
   @param[in]   other         Surface on which the velocity is evaluated.
   @param[in]   other_panel   Panel on which the velocity is evaluated.
   @param[in]   this_panel    Panel on which the vortex ring is located.
   
   @returns Velocity induced by the source panel.
*/
Vector3d
Surface::source_unit_velocity(const Surface &other, int other_panel, int this_panel) const
{ 
    return source_unit_velocity(other.panel_collocation_point(other_panel, true), this_panel);
}

/**
   Computes the velocity induced by a vortex ring of unit strength.  If the velocity
   has been computed before, the cached velocity is returned.  If not, the quantity
   is computed and cached.
   
   @param[in]   other         Surface on which the velocity is evaluated.
   @param[in]   other_panel   Panel on which the velocity is evaluated.
   @param[in]   this_panel    Panel on which the vortex ring is located.
   
   @returns Velocity induced by the vortex ring.
*/
Vector3d
Surface::vortex_ring_unit_velocity(const Surface &other, int other_panel, int this_panel) const
{
    return vortex_ring_unit_velocity(other.panel_collocation_point(other_panel, true), this_panel);
}
