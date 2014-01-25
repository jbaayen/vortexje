//
// Vortexje -- Lifting surface builder.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <iostream>
#include <algorithm>

#include <Eigen/Geometry>

#include <vortexje/lifting-surface-builder.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

/**
   Constructs a new LiftingSurfaceBuilder object for the given LiftingSurface.
   
   @param[in]   lifting_surface    LiftingSurface object to construct.
*/
LiftingSurfaceBuilder::LiftingSurfaceBuilder(LiftingSurface &lifting_surface) : SurfaceBuilder(lifting_surface), lifting_surface(lifting_surface)
{
}

/**
   Creates new nodes for the given points, tracing an airfoil, adding the designated one to the trailing edge.
   
   @param[in]   airfoil_points          List of points tracing an airfoil.
   @param[in]   trailing_edge_point_id  Index of trailing edge point in the above list.
   
   @returns A list of new node numbers.
*/
vector<int>
LiftingSurfaceBuilder::create_nodes(const vector<Vector3d, Eigen::aligned_allocator<Vector3d> > &airfoil_points, int trailing_edge_point_id)
{
    // Add points:
    vector<int> new_nodes = this->SurfaceBuilder::create_nodes(airfoil_points);
    
    // Add trailing edge node to trailing edge node list:
    lifting_surface.trailing_edge_nodes.push_back(new_nodes[trailing_edge_point_id]);
    
    // Done:
    return new_nodes;
}

/**
   Connects two lists of nodes, tracing two airfoils, with new panels.
  
   @param[in]   first_airfoil_nodes                    First list of node numbers.
   @param[in]   second_airfoil_nodes                   Second list of node numbers.
   @param[in]   cyclic                         True if the last nodes in the lists are adjacent to the first nodes in the lists.
   @param[in]   trailing_edge_point_id         Index of trailing edge node in lists of node numbers.
   
   @returns A list of new panel numbers.
*/
vector<int>
LiftingSurfaceBuilder::create_panels_between(const vector<int> &first_airfoil_nodes, const vector<int> &second_airfoil_nodes, int trailing_edge_point_id)
{
    // Connect nodes:
    vector<int> new_panels = this->SurfaceBuilder::create_panels_between(first_airfoil_nodes, second_airfoil_nodes, true);
    
    // Construct airfoil strips:
    vector<int> upper_panel_strip;
    vector<int> lower_panel_strip;
    
    for (int i = 0; i < (int) new_panels.size(); i++) {
        int panel_id = new_panels[i];
        
        // Add to relevant strip:
        if (i < trailing_edge_point_id)
            upper_panel_strip.push_back(panel_id);
        else
            lower_panel_strip.insert(lower_panel_strip.begin(), panel_id);
        
        // Mark as trailing edge panel if bordering trailing edge node:
        if (i == trailing_edge_point_id - 1)
            lifting_surface.upper_trailing_edge_panels.push_back(panel_id);
        else if (i == trailing_edge_point_id)
            lifting_surface.lower_trailing_edge_panels.push_back(panel_id);
    }
    
    lifting_surface.upper_panel_strips.push_back(upper_panel_strip);
    lifting_surface.lower_panel_strips.push_back(lower_panel_strip);
    
    // Done:
    return new_panels;
}

/**
   Fills an airfoil with panels.
  
   @param[in]   airfoil_nodes           Node numbers tracing an airfoil.
   @param[in]   trailing_edge_point_id  Index of trailing edge node in list of airfoil nodes.
   @param[in]   z_sign                  Handedness of the created panels.
   
   @returns List of new panels numbers.
*/
vector<int>
LiftingSurfaceBuilder::create_panels_inside(const vector<int> &airfoil_nodes, int trailing_edge_point_id, int z_sign)
{
    // Add middle nodes:
    vector<int> upper_nodes;
    vector<int> lower_nodes;
    vector<int> middle_nodes;
    for (int i = 1; i < trailing_edge_point_id; i++) {
        int upper_node_id    = airfoil_nodes[i];
        int lower_node_id = airfoil_nodes[airfoil_nodes.size() - i];
        
        Vector3d upper_point    = lifting_surface.nodes[upper_node_id];
        Vector3d lower_point = lifting_surface.nodes[lower_node_id];
        
        Vector3d upper_deformation_velocity    = lifting_surface.node_deformation_velocities[upper_node_id];
        Vector3d lower_deformation_velocity = lifting_surface.node_deformation_velocities[lower_node_id];
        
        Vector3d middle_point = 0.5 * (upper_point + lower_point);
        Vector3d middle_deformation_velocity = 0.5 * (upper_deformation_velocity + lower_deformation_velocity);
        
        int middle_node_id = lifting_surface.nodes.size();

        lifting_surface.nodes.push_back(middle_point);
        
        lifting_surface.node_deformation_velocities.push_back(middle_deformation_velocity);
        
        vector<int> *empty_vector = new vector<int>;
        lifting_surface.node_panel_neighbors.push_back(empty_vector);
        
        upper_nodes.push_back(upper_node_id);
        lower_nodes.push_back(lower_node_id);
        middle_nodes.push_back(middle_node_id);
    }
    
    // Close middle part with panels:
    if (z_sign == 1) {
        create_panels_between(middle_nodes, upper_nodes, false);
        create_panels_between(lower_nodes, middle_nodes, false);
    } else { 
        create_panels_between(upper_nodes, middle_nodes, false);
        create_panels_between(middle_nodes, lower_nodes, false);
    }
    
    // Create triangle for leading and trailing edges:
    if (z_sign == 1) { 
        lifting_surface.add_triangle(airfoil_nodes[0], airfoil_nodes[1], middle_nodes[0]);
        lifting_surface.add_triangle(airfoil_nodes[0], middle_nodes[0], airfoil_nodes[airfoil_nodes.size() - 1]);
        
        lifting_surface.add_triangle(middle_nodes[trailing_edge_point_id - 2], airfoil_nodes[trailing_edge_point_id - 1], airfoil_nodes[trailing_edge_point_id]);
        lifting_surface.add_triangle(middle_nodes[trailing_edge_point_id - 2], airfoil_nodes[trailing_edge_point_id], airfoil_nodes[trailing_edge_point_id + 1]);
    } else {
        lifting_surface.add_triangle(airfoil_nodes[0], middle_nodes[0], airfoil_nodes[1]);
        lifting_surface.add_triangle(airfoil_nodes[0], airfoil_nodes[airfoil_nodes.size() - 1], middle_nodes[0]);
        
        lifting_surface.add_triangle(middle_nodes[trailing_edge_point_id - 2], airfoil_nodes[trailing_edge_point_id], airfoil_nodes[trailing_edge_point_id - 1]);
        lifting_surface.add_triangle(middle_nodes[trailing_edge_point_id - 2], airfoil_nodes[trailing_edge_point_id + 1], airfoil_nodes[trailing_edge_point_id]);
    }
    
    // Done:
    return middle_nodes;
}
