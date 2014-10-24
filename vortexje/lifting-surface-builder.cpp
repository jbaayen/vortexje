//
// Vortexje -- Lifting surface builder.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <vortexje/lifting-surface-builder.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

/**
   Constructs a new LiftingSurfaceBuilder object for the given LiftingSurface.
   
   @param[in]   lifting_surface   LiftingSurface object to construct.
*/
LiftingSurfaceBuilder::LiftingSurfaceBuilder(LiftingSurface &lifting_surface) : SurfaceBuilder(lifting_surface), lifting_surface(lifting_surface)
{
}

/**
   Fills an airfoil with panels.
  
   @param[in]   airfoil_nodes            Node numbers tracing an airfoil, starting with the leading edge.
   @param[in]   trailing_edge_point_id   Index of trailing edge node in list of airfoil nodes.
   @param[in]   z_sign                   Handedness of the created panels.
   
   @returns List of new panel numbers.
*/
vector<int>
LiftingSurfaceBuilder::create_panels_inside_airfoil(const vector<int> &airfoil_nodes, int trailing_edge_point_id, int z_sign)
{
    vector<int> new_panels;
    
    // Add middle nodes:
    vector<int> upper_nodes;
    vector<int> lower_nodes;
    vector<int> middle_nodes;
    for (int i = 1; i < trailing_edge_point_id; i++) {
        int upper_node_id = airfoil_nodes[i];
        int lower_node_id = airfoil_nodes[airfoil_nodes.size() - i];
        
        Vector3d upper_point = lifting_surface.nodes[upper_node_id];
        Vector3d lower_point = lifting_surface.nodes[lower_node_id];
        
        Vector3d middle_point = 0.5 * (upper_point + lower_point);
        
        int middle_node_id = lifting_surface.nodes.size();

        lifting_surface.nodes.push_back(middle_point);
        
        shared_ptr<vector<int> > empty_vector = make_shared<vector<int> >();
        lifting_surface.node_panel_neighbors.push_back(empty_vector);
        
        upper_nodes.push_back(upper_node_id);
        lower_nodes.push_back(lower_node_id);
        middle_nodes.push_back(middle_node_id);
    }
    
    // Close middle part with panels:
    vector<int> new_between_panels;
    
    if (z_sign == 1) {
        new_between_panels = create_panels_between_shapes(middle_nodes, upper_nodes, false);
        new_panels.insert(new_panels.end(), new_between_panels.begin(), new_between_panels.end());
        
        new_between_panels = create_panels_between_shapes(lower_nodes, middle_nodes, false);
        new_panels.insert(new_panels.end(), new_between_panels.begin(), new_between_panels.end());
        
    } else { 
        new_between_panels = create_panels_between_shapes(upper_nodes, middle_nodes, false);
        new_panels.insert(new_panels.end(), new_between_panels.begin(), new_between_panels.end());
        
        new_between_panels = create_panels_between_shapes(middle_nodes, lower_nodes, false);
        new_panels.insert(new_panels.end(), new_between_panels.begin(), new_between_panels.end());
        
    }
    
    // Create triangle for leading and trailing edges:
    int new_panel;
    
    if (z_sign == 1) { 
        new_panel = lifting_surface.add_triangle(airfoil_nodes[0], airfoil_nodes[1], middle_nodes[0]);
        new_panels.push_back(new_panel);
        
        new_panel = lifting_surface.add_triangle(airfoil_nodes[0], middle_nodes[0], airfoil_nodes[airfoil_nodes.size() - 1]);
        new_panels.push_back(new_panel);
        
        new_panel = lifting_surface.add_triangle(middle_nodes[trailing_edge_point_id - 2], airfoil_nodes[trailing_edge_point_id - 1], airfoil_nodes[trailing_edge_point_id]);
        new_panels.push_back(new_panel);
        
        new_panel = lifting_surface.add_triangle(middle_nodes[trailing_edge_point_id - 2], airfoil_nodes[trailing_edge_point_id], airfoil_nodes[trailing_edge_point_id + 1]);
        new_panels.push_back(new_panel);
        
    } else {
        new_panel = lifting_surface.add_triangle(airfoil_nodes[0], middle_nodes[0], airfoil_nodes[1]);
        new_panels.push_back(new_panel);
        
        new_panel = lifting_surface.add_triangle(airfoil_nodes[0], airfoil_nodes[airfoil_nodes.size() - 1], middle_nodes[0]);
        new_panels.push_back(new_panel);
        
        new_panel = lifting_surface.add_triangle(middle_nodes[trailing_edge_point_id - 2], airfoil_nodes[trailing_edge_point_id], airfoil_nodes[trailing_edge_point_id - 1]);
        new_panels.push_back(new_panel);
        
        new_panel = lifting_surface.add_triangle(middle_nodes[trailing_edge_point_id - 2], airfoil_nodes[trailing_edge_point_id + 1], airfoil_nodes[trailing_edge_point_id]);
        new_panels.push_back(new_panel);
        
    }
    
    // Done:
    return new_panels;
}

/**
   Finishes the surface construction process by computing the topology as well as various geometrical properties.
  
   @param[in]   node_strips              Chordwise strips of nodes, each strip starting with the leading edge.
   @param[in]   panel_strips             Chordwise strips of panels, each strip starting with the leading edge.
   @param[in]   trailing_edge_point_id   Index of trailing edge node in list of airfoil nodes.
*/
void
LiftingSurfaceBuilder::finish(const std::vector<std::vector<int> > &node_strips, const std::vector<std::vector<int> > &panel_strips, int trailing_edge_point_id)
{
    // Perform basic surface finishing.
    this->SurfaceBuilder::finish();
    
    // Put together the node matrices.
    // We need to apply some trickery to make sure that the leading and trailing edge nodes appear in both matrices.
    lifting_surface.upper_nodes.resize(trailing_edge_point_id + 1, node_strips.size());
    lifting_surface.lower_nodes.resize(node_strips[0].size() - trailing_edge_point_id + 1, node_strips.size());
    
    for (int i = 0; i <= trailing_edge_point_id; i++) {
        for (int j = 0; j < (int) node_strips.size(); j++) {
            lifting_surface.upper_nodes(i, j) = node_strips[j][i];
            
            if (i == 0)
                lifting_surface.lower_nodes(i, j) = node_strips[j][i];
        }
    }
    
    for (int i = trailing_edge_point_id; i < (int) node_strips[0].size(); i++)
        for (int j = 0; j < (int) node_strips.size(); j++)
            lifting_surface.lower_nodes(node_strips[0].size() - i, j) = node_strips[j][i];
    
    // Put together the panel matrices.
    lifting_surface.upper_panels.resize(lifting_surface.upper_nodes.rows() - 1, lifting_surface.upper_nodes.cols() - 1);
    lifting_surface.lower_panels.resize(lifting_surface.lower_nodes.rows() - 1, lifting_surface.lower_nodes.cols() - 1);
    
    for (int i = 0; i < (int) lifting_surface.upper_panels.rows(); i++)
        for (int j = 0; j < (int) lifting_surface.upper_panels.cols(); j++)
            lifting_surface.upper_panels(i, j) = panel_strips[j][i];
    
    for (int i = 0; i < (int) lifting_surface.lower_panels.rows(); i++)
        for (int j = 0; j < (int) lifting_surface.lower_panels.cols(); j++)
            lifting_surface.lower_panels(lifting_surface.lower_panels.rows() - 1 - i, j) = panel_strips[j][lifting_surface.upper_panels.rows() + i];
            
    // Finish trailing edge setup.
    lifting_surface.finish_trailing_edge();
}
