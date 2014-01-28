//
// Vortexje -- Lifting surface.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <algorithm>

#include <vortexje/lifting-surface.hpp>
#include <vortexje/parameters.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

/**
   Constructs an empty LiftingSurface.
*/
LiftingSurface::LiftingSurface() : Surface()
{
    // Do nothing.
}

/**
   Returns the number of chordwise nodes.

   @returns The number of chordwise nodes.
*/
int 
LiftingSurface::n_chordwise_nodes() const
{
    return (int) upper_nodes.rows();
}

/**
   Returns the number of chordwise panels.
   
   @returns The number of chordwise panels.
*/
int 
LiftingSurface::n_chordwise_panels() const
{
    return (int) upper_panels.rows();
}

/**
   Returns the number of spanwise nodes.

   @returns The number of spanwise nodes.
*/
int 
LiftingSurface::n_spanwise_nodes() const
{
    return (int) upper_nodes.cols();
}

/**
   Returns the number of spanwise panels.
   
   @returns The number of spanwise panels.
*/
int 
LiftingSurface::n_spanwise_panels() const
{
    return (int) upper_panels.cols();
}

/**
   Returns the index'th trailing edge node.
  
   @param[in]   index   Trailing edge node index.
   
   @returns The node number of the index'th trailing edge node.
*/
int 
LiftingSurface::trailing_edge_node(int index) const
{
    return upper_nodes(upper_nodes.rows() - 1, index);
}

/**
   Returns the index'th upper trailing edge panel.
  
   @param[in]   index   Trailing edge panel index.
   
   @returns The panel number of the index'th upper trailing edge panel.
*/
int
LiftingSurface::trailing_edge_upper_panel(int index) const
{
    return upper_panels(upper_panels.rows() - 1, index);
}

/**
   Returns the index'th lower trailing edge panel.
  
   @param[in]   index   Trailing edge panel index.
   
   @returns The panel number of the index'th lower trailing edge panel.
*/
int
LiftingSurface::trailing_edge_lower_panel(int index) const
{
    return lower_panels(lower_panels.rows() - 1, index);
}

// Find panel closest to given point, returning 'true' if the given point lies with the notch
// of the interpolation layer close to the trailing edge. See
//   K. Dixon, C. S. Ferreira, C. Hofemann, G. van Brussel, G. van Kuik,
//   A 3D Unsteady Panel Method for Vertical Axis Wind Turbines, DUWIND, 2008.
bool
LiftingSurface::closest_panel(const Eigen::Vector3d &x, int &panel, double &distance) const
{
    this->Surface::closest_panel(x, panel, distance);
    
    bool trailing_edge = false;
    for (int i = 0; i < n_spanwise_nodes(); i++) {
        const Vector3d &trailing_edge_point = nodes[trailing_edge_node(i)];
        
        // Compute local coordinate system:
        const Vector3d &upper_chordwise_neighbor = nodes[upper_nodes(n_chordwise_nodes() - 2, i)];
        const Vector3d &lower_chordwise_neighbor = nodes[lower_nodes(n_chordwise_nodes() - 2, i)];
        
        Vector3d upper_chord_direction = trailing_edge_point - upper_chordwise_neighbor;
        Vector3d lower_chord_direction = trailing_edge_point - lower_chordwise_neighbor;
        upper_chord_direction.normalize();
        lower_chord_direction.normalize();
        
        Vector3d chord_direction = upper_chord_direction + lower_chord_direction;
        chord_direction.normalize();
        
        Vector3d span_direction(0, 0, 0);
        
        if (i > 0) {
            const Vector3d &left_spanwise_neighbor = nodes[trailing_edge_node(i - 1)];
            
            Vector3d left_span_direction = trailing_edge_point - left_spanwise_neighbor;
            left_span_direction.normalize();
            
            span_direction += left_span_direction;
        }
        
        if (i < n_spanwise_nodes() - 1) {
            const Vector3d &right_spanwise_neighbor = nodes[trailing_edge_node(i + 1)];
            
            Vector3d right_span_direction = trailing_edge_point - right_spanwise_neighbor;
            right_span_direction.normalize();
            
            span_direction -= right_span_direction;
        }
        
        span_direction.normalize(); 
        
        // Test for interpolation layer notch:
        Vector3d x_direction = x - trailing_edge_point;
        x_direction -= x_direction.dot(span_direction) * span_direction;
        if (x_direction.norm() < Parameters::inversion_tolerance) {
            trailing_edge = true;
            break;
        }
        
        x_direction.normalize();
        
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
LiftingSurface::near_exterior_point(int node) const
{
    // For the trailing edge, the close to body point is set equal to the trailing edge point. 
    // This is to guarantee wake emission at the correct angle; see
    //   K. Dixon, C. S. Ferreira, C. Hofemann, G. van Brussel, G. van Kuik,
    //   A 3D Unsteady Panel Method for Vertical Axis Wind Turbines, DUWIND, 2008.
    bool trailing_edge = false;
    for (int i = 0; i < n_spanwise_nodes(); i++) {
        if (node == trailing_edge_node(i)) {
            trailing_edge = true;
            break;
        }
    }
    
    if (trailing_edge)
        return nodes[node];
    else
        return this->Surface::near_exterior_point(node);
}
