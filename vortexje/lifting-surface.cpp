//
// Vortexje -- Lifting surface.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

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

/**
   Finishes the set up of the trailing edge.  
   
   This function computes the trailing edge bisectors, the initial wake strip normals, and terminates the neigbor relationships
   along the trailing edge.
*/
void
LiftingSurface::finish_trailing_edge() 
{
    // Compute trailing edge bisectors and normals to the initial wake strip surface:
    trailing_edge_bisectors.resize(n_spanwise_nodes(), 3);
    wake_normals.resize(n_spanwise_nodes(), 3);
    
    if (n_chordwise_nodes() > 1) {
        for (int i = 0; i < n_spanwise_nodes(); i++) {
            // Compute bisector: 
            Vector3d upper = nodes[upper_nodes(upper_nodes.rows() - 1, i)] - nodes[upper_nodes(upper_nodes.rows() - 2, i)];
            Vector3d lower = nodes[lower_nodes(lower_nodes.rows() - 1, i)] - nodes[lower_nodes(lower_nodes.rows() - 2, i)];
            
            upper.normalize();
            lower.normalize();
            
            Vector3d trailing_edge_bisector = upper + lower;
            trailing_edge_bisector.normalize();
            
            trailing_edge_bisectors.row(i) = trailing_edge_bisector;
            
            // Compute normal to the initial wake strip surface, spanned by the bisector and by the span direction:
            int prev_node, next_node;
            
            if (i > 0)
                prev_node = trailing_edge_node(i - 1);
            else
                prev_node = trailing_edge_node(i);
            
            if (i < n_spanwise_nodes() - 1)
                next_node = trailing_edge_node(i + 1);
            else
                next_node = trailing_edge_node(i);
                
            Vector3d wake_normal(0, 0, 0);
            
            if (prev_node != next_node) {
                Vector3d span_direction = nodes[next_node] - nodes[prev_node];
                
                wake_normal = span_direction.cross(trailing_edge_bisector);
                wake_normal.normalize();
            }
            
            wake_normals.row(i) = wake_normal;
        }
        
    } else {
        // No bisector information available:
        trailing_edge_bisectors.setZero();
        wake_normals.setZero();
    }
    
    // Terminate neighbor relationships across trailing edge.
    for (int i = 0; i < n_spanwise_panels(); i++)
        cut_panels(trailing_edge_upper_panel(i), trailing_edge_lower_panel(i));
}

/**
   Transforms this lifting surface.
   
   @param[in]   transformation   Affine transformation.
*/
void
LiftingSurface::transform(const Eigen::Transform<double, 3, Eigen::Affine> &transformation)
{
    // Call super:
    this->Surface::transform(transformation);
    
    // Transform bisectors and wake normals:
    for (int i = 0; i < n_spanwise_nodes(); i++) {
        Vector3d trailing_edge_bisector = trailing_edge_bisectors.row(i);
        trailing_edge_bisectors.row(i) = transformation.linear() * trailing_edge_bisector;
        
        Vector3d wake_normal = wake_normals.row(i);
        wake_normals.row(i) = transformation.linear() * wake_normal;
    }
}

/**
   Returns the wake emission velocity, at the node_index'th trailing edge node.
  
   @param[in]   apparent_velocity   Apparent velocity.
   @param[in]   node_index          Trailing edge node index.
   
   @returns The wake emission velocity, at the node_index'th trailinge edge node.
*/
Eigen::Vector3d
LiftingSurface::wake_emission_velocity(const Eigen::Vector3d &apparent_velocity, int node_index) const
{
    Vector3d wake_emission_velocity;
    
    if (Parameters::wake_emission_follow_bisector) {
        Vector3d wake_normal = wake_normals.row(node_index);
        
        // Project apparent velocity onto wake emission plane:
        wake_emission_velocity = -(apparent_velocity - apparent_velocity.dot(wake_normal) * wake_normal);
        
    } else {
        // Emit wake in direction of apparent velocity:
        wake_emission_velocity = -apparent_velocity;
        
    }
    
    // Done:
    return wake_emission_velocity;
}
