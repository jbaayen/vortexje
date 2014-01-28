//
// Vortexje -- Surface builder.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <iostream>

#include <Eigen/Geometry>

#include <vortexje/surface-builder.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

/**
   Constructs a new SurfaceBuilder object for the given Surface.
   
   @param[in]   surface    Surface object to construct.
*/
SurfaceBuilder::SurfaceBuilder(Surface &surface) : surface(surface)
{
}

/**
   Creates new nodes for the given list of points.
   
   @param[in]   points   List of points.
   
   @returns A list of new node numbers.
*/
vector<int>
SurfaceBuilder::create_nodes(const vector<Vector3d, Eigen::aligned_allocator<Vector3d> > &points)
{
    vector<int> new_nodes;
    
    for (int i = 0; i < (int) points.size(); i++) {
        int node_id = surface.nodes.size();
        
        surface.nodes.push_back(points[i]);
        
        vector<int> *empty_vector = new vector<int>;
        surface.node_panel_neighbors.push_back(empty_vector);
        
        surface.node_deformation_velocities.push_back(Vector3d(0, 0, 0));
            
        new_nodes.push_back(node_id);
    }
    
    return new_nodes;
}

/**
   Connects two lists of nodes with new panels.
  
   @param[in]   first_nodes    First list of node numbers.
   @param[in]   second_nodes   Second list of node numbers.
   @param[in]   cyclic         True if the last nodes in the lists are adjacent to the first nodes in the lists.
   
   @returns A list of new panel numbers.
*/
vector<int>
SurfaceBuilder::create_panels_between(const vector<int> &first_nodes, const vector<int> &second_nodes, bool cyclic)
{
    vector<int> new_panels;
    
    for (int i = 0; i < (int) first_nodes.size(); i++) {
        int next_i;
        if (i == (int) first_nodes.size() - 1) {
            if (cyclic)
                next_i = 0;
            else
                break;
        } else
            next_i = i + 1;

        // Add a planar trapezoidal panel, follosurface
        //    T. Cebeci, An Engineering Approach to the Calculation of Aerodynamic Flows, Springer, 1999.
        vector<int> original_nodes;
        original_nodes.push_back(first_nodes[i]);
        original_nodes.push_back(first_nodes[next_i]);
        original_nodes.push_back(second_nodes[i]);
        original_nodes.push_back(second_nodes[next_i]);
        
        Vector3d first_line  = surface.nodes[first_nodes[next_i]] - surface.nodes[first_nodes[i]];
        Vector3d second_line = surface.nodes[second_nodes[next_i]] - surface.nodes[second_nodes[i]];
        
        Vector3d first_line_deformation_velocity = surface.node_deformation_velocities[first_nodes[next_i]] - surface.node_deformation_velocities[first_nodes[i]];
        Vector3d second_line_deformation_velocity = surface.node_deformation_velocities[second_nodes[next_i]] - surface.node_deformation_velocities[second_nodes[i]];
        
        Vector3d line_direction = first_line + second_line;
        line_direction.normalize();
        
        Vector3d line_direction_deformation_velocity = (first_line_deformation_velocity + second_line_deformation_velocity) / (first_line + second_line).norm() - (first_line + second_line) * (first_line + second_line).dot(first_line_deformation_velocity + second_line_deformation_velocity) / pow((first_line + second_line).norm(), 3);
        
        Vector3d first_mid  = 0.5 * (surface.nodes[first_nodes[i]] + surface.nodes[first_nodes[next_i]]);
        Vector3d second_mid = 0.5 * (surface.nodes[second_nodes[i]] + surface.nodes[second_nodes[next_i]]);
        
        Vector3d first_mid_deformation_velocity = 0.5 * (surface.node_deformation_velocities[first_nodes[i]] + surface.node_deformation_velocities[first_nodes[next_i]]);
        Vector3d second_mid_deformation_velocity = 0.5 * (surface.node_deformation_velocities[second_nodes[i]] + surface.node_deformation_velocities[second_nodes[next_i]]);
        
        Vector3d vertices[4];
        vertices[0] = first_mid - 0.5 * first_line.norm() * line_direction;
        vertices[1] = first_mid + 0.5 * first_line.norm() * line_direction;
        vertices[2] = second_mid - 0.5 * second_line.norm() * line_direction;
        vertices[3] = second_mid + 0.5 * second_line.norm() * line_direction;
        
        Vector3d vertex_deformation_velocities[4];
        vertex_deformation_velocities[0] = first_mid_deformation_velocity - 0.5 * (first_line.dot(first_line_deformation_velocity) / first_line.norm() * line_direction + first_line.norm() * line_direction_deformation_velocity);
        vertex_deformation_velocities[1] = first_mid_deformation_velocity + 0.5 * (first_line.dot(first_line_deformation_velocity) / first_line.norm() * line_direction + first_line.norm() * line_direction_deformation_velocity);
        vertex_deformation_velocities[2] = second_mid_deformation_velocity - 0.5 * (second_line.dot(second_line_deformation_velocity) / second_line.norm() * line_direction + second_line.norm() * line_direction_deformation_velocity);
        vertex_deformation_velocities[3] = second_mid_deformation_velocity + 0.5 * (second_line.dot(second_line_deformation_velocity) / second_line.norm() * line_direction + second_line.norm() * line_direction_deformation_velocity);
        
        int new_nodes[4];
        for (int j = 0; j < 4; j++) {
            // If the new points don't match the original ones, create new nodes:
            if ((vertices[j] - surface.nodes[original_nodes[j]]).norm() < Parameters::inversion_tolerance) {
                new_nodes[j] = original_nodes[j];
            } else {         
                new_nodes[j] = surface.nodes.size();
                
                surface.nodes.push_back(vertices[j]);
                
                surface.node_deformation_velocities.push_back(vertex_deformation_velocities[j]);
                
                surface.node_panel_neighbors.push_back(surface.node_panel_neighbors[original_nodes[j]]);
            }
        }
        
        int panel_id = surface.add_quadrangle(new_nodes[0], new_nodes[2], new_nodes[3], new_nodes[1]);
        
        new_panels.push_back(panel_id);
    }
    
    return new_panels;
}

/**
   Fills a shape with triangular panels, all meeting on the specified tip point.
  
   @param[in]   nodes       Node numbers tracing a shape.
   @param[in]   tip_point   Point where all triangles meet.
   @param[in]   z_sign      Handedness of the created panels.
   
   @returns List of new panel numbers.
*/
vector<int>
SurfaceBuilder::create_panels_inside(const vector<int> &nodes, const Vector3d &tip_point, int z_sign)
{
    vector<int> new_panels;
    
    // Add tip node:
    int tip_node = surface.nodes.size();

    surface.nodes.push_back(tip_point);
    
    surface.node_deformation_velocities.push_back(Vector3d(0, 0, 0));
    
    vector<int> *empty_vector = new vector<int>;
    surface.node_panel_neighbors.push_back(empty_vector);
    
    // Create triangle for leading and trailing edges:
    for (int i = 0; i < (int) nodes.size(); i++) {
        int triangle[3];
        
        if (z_sign == 1) {
            triangle[0] = nodes[i];
            if (i == (int) nodes.size() - 1)
                triangle[1] = nodes[0];
            else
                triangle[1] = nodes[i + 1];
            triangle[2] = tip_node;
                
        } else {
            triangle[0] = nodes[i];
            triangle[1] = tip_node;
            if (i == (int) nodes.size() - 1)
                triangle[2] = nodes[0];
            else
                triangle[2] = nodes[i + 1];
            
        }
            
        int new_panel = surface.add_triangle(triangle[0], triangle[1], triangle[2]);
        new_panels.push_back(new_panel);
    }
    
    // Done:
    return new_panels;
}

/**
   Finishes the surface construction process by computing the topology as well as various geometrical properties.
*/
void
SurfaceBuilder::finish()
{
    surface.compute_topology();
    surface.compute_geometry();
}
