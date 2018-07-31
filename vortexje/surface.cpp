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
#include <algorithm>

#include <Eigen/Geometry>

#include <vortexje/surface.hpp>
#include <vortexje/parameters.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

// Avoid having to divide by 4 pi all the time:
static const double pi = 3.141592653589793238462643383279502884;
static const double one_over_4pi = 1.0 / (4 * pi);

/**
   Constructs an empty surface.
*/
Surface::Surface(const string &id) : id(id)
{
}

/**
   Destructor.
*/
Surface::~Surface()
{
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
   compute_topology() when done adding panels.
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
        vector<vector<pair<int, int> > > single_panel_neighbors;
        single_panel_neighbors.resize(panel_nodes[i].size());
        
        // Every node gives rise to one edge, which gives rise to at most one neighbor.
        for (int j = 0; j < (int) panel_nodes[i].size(); j++) {          
            // Compute index for next node.
            int next_j;
            if (j == (int) panel_nodes[i].size() - 1)
                next_j = 0;
            else
                next_j = j + 1;
                
            int node      = panel_nodes[i][j];
            int next_node = panel_nodes[i][next_j];
            
            for (int k = 0; k < (int) node_panel_neighbors[node]->size(); k++) {
                int potential_neighbor = (*node_panel_neighbors[node])[k];
                if (potential_neighbor == i)
                    continue;
                
                // Is this neighbor shared with the next node?   
                if (find(node_panel_neighbors[next_node]->begin(),
                         node_panel_neighbors[next_node]->end(),
                         potential_neighbor) != node_panel_neighbors[next_node]->end()) {
                    // Yes: Establish neighbor relationship.
                    int potential_neighbor_edge = -1;
                    for (int l = 0; l < (int) panel_nodes[potential_neighbor].size(); l++) {
                        int next_l;
                        if (l == (int) panel_nodes[potential_neighbor].size() - 1)
                            next_l = 0;
                        else
                            next_l = l + 1;
                            
                        int potential_neighbor_node      = panel_nodes[potential_neighbor][l];
                        int potential_neighbor_next_node = panel_nodes[potential_neighbor][next_l];
                            
                        if (find(node_panel_neighbors[potential_neighbor_node]->begin(),
                                 node_panel_neighbors[potential_neighbor_node]->end(),
                                 i) != node_panel_neighbors[potential_neighbor_node]->end() &&
                            find(node_panel_neighbors[potential_neighbor_next_node]->begin(),
                                 node_panel_neighbors[potential_neighbor_next_node]->end(),
                                 i) != node_panel_neighbors[potential_neighbor_next_node]->end()) {
                            potential_neighbor_edge = l;
                            break;
                        }
                    }
                    
                    assert(potential_neighbor_edge >= 0);
                    
                    single_panel_neighbors[j].push_back(make_pair(potential_neighbor, potential_neighbor_edge));
                }
            }
        }
        
        panel_neighbors.push_back(single_panel_neighbors);
    }
}

/**
   Cuts the panel neighbor relationship between two given panels.
   
   @param[in]   panel_a   First panel.
   @param[in]   panel_b   Second panel.
*/
void
Surface::cut_panels(int panel_a, int panel_b)
{
    for (int i = 0; i < 2; i++) {
        int panel_ids[2];
        if (i == 0) {
            panel_ids[0] = panel_a;
            panel_ids[1] = panel_b;
        } else {
            panel_ids[0] = panel_b;
            panel_ids[1] = panel_a;
        }
        
        vector<vector<pair<int, int> > > &single_panel_neighbors = panel_neighbors[panel_ids[0]];
        vector<vector<pair<int, int> > >::iterator it;
        for (it = single_panel_neighbors.begin(); it != single_panel_neighbors.end(); it++) {
            vector<pair<int, int> > &edge_neighbors = *it;
            vector<pair<int, int> >::iterator sit;
            for (sit = edge_neighbors.begin(); sit != edge_neighbors.end(); ) {
                if (sit->first == panel_ids[1])
                    sit = edge_neighbors.erase(sit);
                else  
                    sit++;
            }
        }
    }  
}

/**
   Computes the normal, collocation point, and surface area of the specified panel.
   
   @param[in]   panel   Reference panel.
*/
void
Surface::compute_geometry(int panel)
{
    // Resize arrays, if necessary:
    panel_normals.resize(n_panels());
    panel_collocation_points[0].resize(n_panels());
    panel_collocation_points[1].resize(n_panels());
    panel_coordinate_transformations.resize(n_panels());
    panel_transformed_points.resize(n_panels());
    panel_surface_areas.resize(n_panels());
    
    // Get panel nodes:
    vector<int> &single_panel_nodes = panel_nodes[panel];
    
    // Normal:
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

    panel_normals[panel] = normal;
    
    // Collocation point:
    Vector3d collocation_point(0, 0, 0);
    for (int j = 0; j < (int) single_panel_nodes.size(); j++)
        collocation_point = collocation_point + nodes[single_panel_nodes[j]];

    collocation_point = collocation_point / single_panel_nodes.size();
        
    panel_collocation_points[0][panel] = collocation_point;
    
    Vector3d below_surface_collocation_point = collocation_point + Parameters::collocation_point_delta * normal;
    panel_collocation_points[1][panel] = below_surface_collocation_point;
    
    // Coordinate transformation:
    Vector3d AB = nodes[single_panel_nodes[1]] - nodes[single_panel_nodes[0]];
    AB.normalize();
    
    Matrix3d rotation;
    rotation.row(0) = AB;
    rotation.row(1) = normal.cross(AB).normalized(); // Should be normalized already.
    rotation.row(2) = normal;
    
    Transform<double, 3, Affine> transformation = rotation * Translation<double, 3>(-collocation_point);

    panel_coordinate_transformations[panel] = transformation;
    
    // Create transformed points.
    vector_aligned<Eigen::Vector3d> single_panel_transformed_points;
    single_panel_transformed_points.reserve(single_panel_nodes.size());
    for (int j = 0; j < (int) single_panel_nodes.size(); j++)
        single_panel_transformed_points.push_back(transformation * nodes[single_panel_nodes[j]]);
        
    panel_transformed_points[panel] = single_panel_transformed_points;
    
    // Surface area: 
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
    
    panel_surface_areas[panel] = surface_area;
}

/**
   Computes the normals, collocation points, and surface areas of all panels.
*/
void
Surface::compute_geometry()
{
    cout << "Surface " << id << ": Computing geometry." << endl;
    
    for (int i = 0; i < n_panels(); i++)
        compute_geometry(i);
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
        panel_normals[i] = transformation.linear() * panel_normals[i];
        
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

// Simultaneously compute influence of source and doublet panel edges on given point.
static void
source_and_doublet_edge_influence(const Vector3d &x, const Vector3d &node_a, const Vector3d &node_b, double *source_edge_influence, double *doublet_edge_influence)
{
    double d = sqrt(pow(node_b(0) - node_a(0), 2) + pow(node_b(1) - node_a(1), 2));

    if (d < Parameters::zero_threshold) {
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
        int next_idx;
        if (i == (int) panel_nodes[this_panel].size() - 1)
            next_idx = 0;
        else
            next_idx = i + 1;
            
        const Vector3d &node_a = panel_transformed_points[this_panel][i];
        const Vector3d &node_b = panel_transformed_points[this_panel][next_idx];
        
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
        int next_idx;
        if (i == (int) panel_nodes[this_panel].size() - 1)
            next_idx = 0;
        else
            next_idx = i + 1;
            
        const Vector3d &node_a = panel_transformed_points[this_panel][i];
        const Vector3d &node_b = panel_transformed_points[this_panel][next_idx];
        
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
        int next_idx;
        if (i == (int) panel_nodes[this_panel].size() - 1)
            next_idx = 0;
        else
            next_idx = i + 1;
            
        const Vector3d &node_a = panel_transformed_points[this_panel][i];
        const Vector3d &node_b = panel_transformed_points[this_panel][next_idx];
        
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
    
    if (d < Parameters::zero_threshold)
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
        int next_idx;
        if (i == (int) panel_nodes[this_panel].size() - 1)
            next_idx = 0;
        else
            next_idx = i + 1;
            
        const Vector3d &node_a = panel_transformed_points[this_panel][i];
        const Vector3d &node_b = panel_transformed_points[this_panel][next_idx];
        
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
        
        if (r_1_norm < Parameters::zero_threshold ||
            r_2_norm < Parameters::zero_threshold ||
            r_1xr_2_sqnorm < Parameters::zero_threshold)
            continue;

        velocity += r_1xr_2 / r_1xr_2_sqnorm * r_0.dot(r_1 / r_1_norm - r_2 / r_2_norm);
    }

    return one_over_4pi * velocity;
}

/**
   Computes the potential influence induced by a doublet panel of unit strength. 
   
   @param[in]   other         Surface on which the influence coefficient is evaluated.
   @param[in]   other_panel   Panel on which the influence coefficient is evaluated.
   @param[in]   this_panel    Panel on which the doublet panel is located.
   
   @returns Influence coefficient.
*/
double
Surface::doublet_influence(const std::shared_ptr<Surface> &other, int other_panel, int this_panel) const
{ 
    if ((this == other.get()) && (this_panel == other_panel))
        return -0.5;
    else
        return doublet_influence(other->panel_collocation_point(other_panel, true), this_panel);
}

/**
   Computes the potential influence induced by a source panel of unit strength.  
   
   @param[in]   other         Surface on which the influence coefficient is evaluated.
   @param[in]   other_panel   Panel on which the influence coefficient is evaluated.
   @param[in]   this_panel    Panel on which the source panel is located.
   
   @returns Influence coefficient.
*/
double
Surface::source_influence(const std::shared_ptr<Surface> &other, int other_panel, int this_panel) const
{
    return source_influence(other->panel_collocation_point(other_panel, true), this_panel);
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
Surface::source_and_doublet_influence(const std::shared_ptr<Surface> &other, int other_panel, int this_panel, double &source_influence, double &doublet_influence) const
{
    if ((this == other.get()) && (this_panel == other_panel)) {
        doublet_influence = -0.5;
        
        source_influence = this->source_influence(other, other_panel, this_panel);
        
    } else
        source_and_doublet_influence(other->panel_collocation_point(other_panel, true), this_panel, source_influence, doublet_influence);
}

/**
   Computes the velocity induced by a source panel of unit strength.  
   
   @param[in]   other         Surface on which the velocity is evaluated.
   @param[in]   other_panel   Panel on which the velocity is evaluated.
   @param[in]   this_panel    Panel on which the vortex ring is located.
   
   @returns Velocity induced by the source panel.
*/
Vector3d
Surface::source_unit_velocity(const std::shared_ptr<Surface> &other, int other_panel, int this_panel) const
{ 
    return source_unit_velocity(other->panel_collocation_point(other_panel, true), this_panel);
}

/**
   Computes the velocity induced by a vortex ring of unit strength.  
   
   @param[in]   other         Surface on which the velocity is evaluated.
   @param[in]   other_panel   Panel on which the velocity is evaluated.
   @param[in]   this_panel    Panel on which the vortex ring is located.
   
   @returns Velocity induced by the vortex ring.
*/
Vector3d
Surface::vortex_ring_unit_velocity(const std::shared_ptr<Surface> &other, int other_panel, int this_panel) const
{
    return vortex_ring_unit_velocity(other->panel_collocation_point(other_panel, true), this_panel);
}
