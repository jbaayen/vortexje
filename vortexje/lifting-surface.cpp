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
   Sorts the strips, as well as the trailing edge node and panel lists, by the spanwise coordinate.
   
   @note Sorted trailing edge node and panel lists are a prerequisite for wake emission.
*/
void
LiftingSurface::sort_strips()
{
    // Note:  Eventually, we want to use std::sort and C++11 lambda functions here.
    bool swapped;
    int n;
        
    // Sort strips:
    {
        vector<vector<int> > *lists[2];
        lists[0] = &upper_panel_strips;
        lists[1] = &lower_panel_strips;
        for (int k = 0; k < 2; k++) {
            vector<vector<int> > &list = *lists[k];
            
            swapped = true;
            n = list.size();
            while (swapped) {
                swapped = false;
                for (int i = 1; i < n; i++) {
                    if (panel_collocation_point(list[i - 1][0], false).dot(span_direction) > panel_collocation_point(list[i][0], false).dot(span_direction)) {
                        vector<int> temp = list[i];
                        list[i] = list[i - 1];
                        list[i - 1] = temp;
                        swapped = true;
                    }
                }
                n--;
            }
        }
    }
    
    // Sort trailing edge nodes:
    {
        swapped = true;
        n = trailing_edge_nodes.size();
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
    }
    
    // Sort trailing edge panels:
    {
        vector<int> *lists[2];
        lists[0] = &upper_trailing_edge_panels;
        lists[1] = &lower_trailing_edge_panels;
        for (int k = 0; k < 2; k++) {
            vector<int> &list = *lists[k];
            
            swapped = true;
            n = list.size();
            while (swapped) {
                swapped = false;
                for (int i = 1; i < n; i++) {
                    if (panel_collocation_point(list[i - 1], false).dot(span_direction) > panel_collocation_point(list[i], false).dot(span_direction)) {
                        int temp = list[i];
                        list[i] = list[i - 1];
                        list[i - 1] = temp;
                        swapped = true;
                    }
                }
                n--;
            }
        }
    }
}

// Transform surface.
void
LiftingSurface::transform(const Eigen::Matrix3d &transformation)
{
    this->Surface::transform(transformation);
    
    chord_direction = transformation * chord_direction;
    span_direction  = transformation * span_direction;
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
    for (int i = 0; i < (int) trailing_edge_nodes.size(); i++) {
        const Vector3d &trailing_edge_node = nodes[trailing_edge_nodes[i]];
        
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
LiftingSurface::close_to_body_point(int node) const
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
        return this->Surface::close_to_body_point(node);
}
