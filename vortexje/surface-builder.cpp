//
// Vortexje -- Surface builder.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <iostream>
#include <limits>

#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>

#include <vortexje/surface-builder.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

/**
   Constructs a new SurfaceBuilder object for the given surface.
   
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
SurfaceBuilder::create_nodes_for_points(const vector_aligned<Vector3d> &points)
{
    vector<int> new_nodes;
    
    for (int i = 0; i < (int) points.size(); i++) {
        int node_id = surface.nodes.size();
        
        surface.nodes.push_back(points[i]);
        
        shared_ptr<vector<int> > empty_vector = make_shared<vector<int> >();
        surface.node_panel_neighbors.push_back(empty_vector);
            
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
SurfaceBuilder::create_panels_between_shapes(const vector<int> &first_nodes, const vector<int> &second_nodes, bool cyclic)
{
    vector<int> new_panels;
    
    for (int i = 0; i < (int) first_nodes.size(); i++) {
        // Bundle panel nodes in appropriate order:
        int next_i;
        if (i == (int) first_nodes.size() - 1) {
            if (cyclic)
                next_i = 0;
            else
                break;
        } else
            next_i = i + 1;
            
        vector<int> original_nodes;
        original_nodes.push_back(first_nodes[i]);
        original_nodes.push_back(second_nodes[i]);
        original_nodes.push_back(second_nodes[next_i]);
        original_nodes.push_back(first_nodes[next_i]); 
        
        // Filter out duplicate nodes while preserving order:
        vector<int> unique_nodes;

        for (int j = 0; j < 4; j++) {
            bool duplicate = false;
            for (int k = j + 1; k < 4; k++) {
                if (original_nodes[j] == original_nodes[k]) {
                    duplicate = true;
                    
                    break;
                }
            }
            
            if (!duplicate)
                unique_nodes.push_back(original_nodes[j]);
        }
        
        // Add panel appropriate for number of nodes:
        int panel_id;
        switch (unique_nodes.size()) {
        case 3:
        {
            // Add triangle:
            panel_id = surface.add_triangle(unique_nodes[0], unique_nodes[1], unique_nodes[2]);

            break;
        }
            
        case 4:
        {
            // Construct a planar quadrangle:
            
            // Center points around their mean:
            Vector3d mean(0, 0, 0);
            for (int i = 0; i < 4; i++)
                mean += surface.nodes[unique_nodes[i]];
            mean /= 4.0;
            
            MatrixXd X(4, 3);
            for (int i = 0; i < 4; i++)
                X.row(i) = surface.nodes[unique_nodes[i]] - mean;
                
            // Perform PCA to find dominant directions:
            SelfAdjointEigenSolver<Matrix3d> solver(X.transpose() * X);
            
            Vector3d eigenvalues = solver.eigenvalues();
            Matrix3d eigenvectors = solver.eigenvectors();

            double min_eigenvalue = numeric_limits<double>::max();
            int min_eigenvalue_index = -1;
            for (int i = 0; i < 3; i++) {
                if (eigenvalues(i) < min_eigenvalue) {
                    min_eigenvalue = eigenvalues(i);
                    min_eigenvalue_index = i;
                }
            }
                
            Vector3d normal = eigenvectors.col(min_eigenvalue_index);
                
            // Create new points by projecting onto surface spanned by dominant directions:
            Vector3d vertices[4];
            for (int i = 0; i < 4; i++)
                vertices[i] = surface.nodes[unique_nodes[i]] - (normal * X.row(i)) * normal; 

            // Add points to surface:
            int new_nodes[4];
            for (int j = 0; j < 4; j++) {
                // If the new points don't match the original ones, create new nodes:
                if ((vertices[j] - surface.nodes[unique_nodes[j]]).norm() < Parameters::zero_threshold) {
                    new_nodes[j] = unique_nodes[j];
                } else {         
                    new_nodes[j] = surface.nodes.size();
                    
                    surface.nodes.push_back(vertices[j]);
                    
                    surface.node_panel_neighbors.push_back(surface.node_panel_neighbors[unique_nodes[j]]);
                }
            }
            
            // Add planar quadrangle:
            panel_id = surface.add_quadrangle(new_nodes[0], new_nodes[1], new_nodes[2], new_nodes[3]);
            
            break;
        }
            
        default:
            // Unknown panel type:
            cerr << "SurfaceBuilder::create_panels_between_shapes: Cannot create panel with " << unique_nodes.size() << " vertices." << endl;
            exit(1);
        }
        
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
SurfaceBuilder::create_panels_inside_shape(const vector<int> &nodes, const Vector3d &tip_point, int z_sign)
{
    vector<int> new_panels;
    
    // Add tip node:
    int tip_node = surface.nodes.size();

    surface.nodes.push_back(tip_point);
    
    shared_ptr<vector<int> > empty_vector = make_shared<vector<int> >();
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
