//
// Vortexje -- Wing.
//
// Copyright (C) 2012 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <iostream>
#include <limits>
#include <algorithm>

#include <Eigen/Geometry>

#include <vortexje/wing.hpp>
#include <vortexje/parameters.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

/**
   Constructs an empty wing.
*/
Wing::Wing() : Mesh()
{
    // Do nothing.  Let subclass construct wing.
}

/**
   Constructs a wing from a given mesh and aerodynamic interpretation.
  
   @param[in]   mesh                Source mesh containing the wing geometry.
   @param[in]   location            Reference location for new wing 
   @param[in]   chord_direction     Unit vector pointing along the wing chord.
   @param[in]   top_direction       Unit vector pointing upwards, out of, and perpendicular to the wing.
   @param[in]   span_direction      Unit vector pointing along the wing span.
*/
Wing::Wing(Mesh &mesh, const Eigen::Vector3d &location, const Eigen::Vector3d &chord_direction, const Eigen::Vector3d &top_direction, const Eigen::Vector3d &span_direction)
    : Mesh(), location(location), chord_direction(chord_direction), top_direction(top_direction), span_direction(span_direction)
{   
    // Load wing data from existing mesh:
    nodes                       = mesh.nodes;
    for (int i = 0; i < (int) mesh.node_panel_neighbors.size(); i++) {
        vector<int> *single_node_panel_neighbors = new vector<int>;
        for (int j = 0; j < (int) mesh.node_panel_neighbors[i]->size(); j++)
            single_node_panel_neighbors->push_back((*mesh.node_panel_neighbors[i])[j]);
        node_panel_neighbors.push_back(single_node_panel_neighbors);
    }
    node_deformation_velocities = mesh.node_deformation_velocities;
    panel_nodes                 = mesh.panel_nodes;
    panel_neighbors             = mesh.panel_neighbors;
    
    // Autodetect trailing edge nodes:
    double max_x = numeric_limits<double>::min(), min_x = numeric_limits<double>::max();
    double max_y = numeric_limits<double>::min(), min_y = numeric_limits<double>::max();
    double max_z = numeric_limits<double>::min(), min_z = numeric_limits<double>::max();
    for (int i = 0; i < n_nodes(); i++) {
        if (nodes[i](0) > max_x)
            max_x = nodes[i](0);
        if (nodes[i](0) < min_x)
            min_x = nodes[i](0);
        if (nodes[i](1) > max_y)
            max_y = nodes[i](1);
        if (nodes[i](1) < min_y)
            min_y = nodes[i](1);
        if (nodes[i](2) > max_z)
            max_z = nodes[i](2);
        if (nodes[i](2) < min_z)
            min_z = nodes[i](2);
    }
    
    for (int i = 0; i < n_nodes(); i++)
        if (nodes[i](0) == max_x)
            trailing_edge_nodes.push_back(i);
    
    // Initialize trailing edge panels:
    for (int i = 0; i < (int) trailing_edge_nodes.size(); i++) {
        int node = trailing_edge_nodes[i];
        for (int j = 0; j < (int) node_panel_neighbors[node]->size(); j++) {
            int neighbor_panel = (*node_panel_neighbors[node])[j];
            
            double direction_coefficient = panel_normal(neighbor_panel).dot(Vector3d::UnitY());
            if (direction_coefficient < -Parameters::sharp_edge_threshold) {
                double found = false;
                for (int k = 0; k < (int) trailing_edge_top_panels.size(); k++) {
                    if (trailing_edge_top_panels[k] == neighbor_panel) {
                        found = true;
                        break;
                    }
                }
                
                if (!found)
                    trailing_edge_top_panels.push_back(neighbor_panel);
                    
            } else if (direction_coefficient > Parameters::sharp_edge_threshold) {
                double found = false;
                for (int k = 0; k < (int) trailing_edge_bottom_panels.size(); k++) {
                    if (trailing_edge_bottom_panels[k] == neighbor_panel) {
                        found = true;
                        break;
                    }
                }
                
                if (!found)
                    trailing_edge_bottom_panels.push_back(neighbor_panel);
            }
        }
    }
    
    // Sort trailing edge:
    sort_trailing_edge();
}

/**
   Sorts the trailing edge node and panel lists.
   
   @note Sorted trailing edge node and panel lists are a prerequisite for wake emission.
*/
void
Wing::sort_trailing_edge()
{
    // Sort trailing edge nodes:
    bool swapped = true;
    int n = trailing_edge_nodes.size();
    while (swapped) {
       swapped = false;
       for (int i = 1; i < n; i++) {
            if (nodes[trailing_edge_nodes[i - 1]].dot(span_direction) > nodes[trailing_edge_nodes[i]].dot(span_direction)) {
                int temp = trailing_edge_nodes[i];
                trailing_edge_nodes[i] = trailing_edge_nodes[i - 1];
                trailing_edge_nodes[i - 1] = temp;
                swapped = true;
            }
       }
       n--;
    }
    
    // Sort trailing edge panels:
    swapped = true;
    n = trailing_edge_top_panels.size();
    while (swapped) {
        swapped = false;
        for (int i = 1; i < n; i++) {
            if (panel_collocation_point(trailing_edge_top_panels[i - 1], false).dot(span_direction) > panel_collocation_point(trailing_edge_top_panels[i], false).dot(span_direction)) {
                int temp = trailing_edge_top_panels[i];
                trailing_edge_top_panels[i] = trailing_edge_top_panels[i - 1];
                trailing_edge_top_panels[i - 1] = temp;
                swapped = true;
            }
        }
        n--;
    }
    
    swapped = true;
    n = trailing_edge_bottom_panels.size();
    while (swapped) {
        swapped = false;
        for (int i = 1; i < n; i++) {
            if (panel_collocation_point(trailing_edge_bottom_panels[i - 1], false).dot(span_direction) > panel_collocation_point(trailing_edge_bottom_panels[i], false).dot(span_direction)) {
                int temp = trailing_edge_bottom_panels[i];
                trailing_edge_bottom_panels[i] = trailing_edge_bottom_panels[i - 1];
                trailing_edge_bottom_panels[i - 1] = temp;
                swapped = true;
            }
        }
        n--;
    }
}

// Translate mesh.
void
Wing::translate(const Eigen::Vector3d &translation)
{
    vector<Mesh*> empty;
    translate(translation, empty);
}

void
Wing::translate(const Eigen::Vector3d &translation, std::vector<Mesh*> &corotating_meshes)
{
    this->Mesh::translate(translation, corotating_meshes);
    
    location += translation;
}

// Rotate mesh.
void
Wing::transform(const Eigen::Matrix3d &transformation)
{
    vector<Mesh*> empty;
    transform(transformation, empty);
}

void
Wing::transform(const Eigen::Matrix3d &transformation, std::vector<Mesh*> &corotating_meshes)
{
    this->Mesh::transform(transformation, corotating_meshes);
    
    chord_direction = transformation * chord_direction;
    top_direction   = transformation * top_direction;
    span_direction  = transformation * span_direction;
}

// Find closest panel to given point, returning 'true' if this panel borders the trailing edge:
bool
Wing::closest_panel(const Eigen::Vector3d &x, int &panel, double &distance)
{
    this->Mesh::closest_panel(x, panel, distance);
    
    bool trailing_edge = false;
    for (int i = 0; i < (int) trailing_edge_nodes.size(); i++) {
        Vector3d &trailing_edge_node = nodes[trailing_edge_nodes[i]];
        
        Vector3d x_direction = x - trailing_edge_node;
        x_direction -= x_direction.dot(span_direction) * span_direction;
        if (x_direction.norm() < Parameters::inversion_tolerance) {
            trailing_edge = true;
            break;
        }
        
        x_direction /= x_direction.norm();
        
        double dot_product = chord_direction.dot(x_direction);  
        double angle = acos(dot_product);
        if (angle < Parameters::interpolation_layer_notch_angle) {
            trailing_edge = true;
            break;
        }
    }
    
    return trailing_edge;
}

// Compute out-of-body point close to given node:
Vector3d
Wing::close_to_body_point(int node)
{
    // For the trailing edge, the close to body point is set equal to the trailing edge point. 
    // This is to guarantee wake emission at the correct angle; see
    //   K. Dixon, C. S. Ferreira, C. Hofemann, G. van Brussel, G. van Kuik,
    //   A 3D Unsteady Panel Method for Vertical Axis Wind Turbines, DUWIND, 2008.
    bool trailing_edge = false;
    for (int i = 0; i < (int) trailing_edge_nodes.size(); i++) {
        if (node == trailing_edge_nodes[i]) {
            trailing_edge = true;
            break;
        }
    }
    
    if (trailing_edge)
        return nodes[node];
    else
        return this->Mesh::close_to_body_point(node);
}
