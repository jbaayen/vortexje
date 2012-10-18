//
// Vortexje -- Wake
//
// Copyright (C) 2012 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <iostream>
#include <algorithm>

#include <Eigen/Geometry>

#include <vortexje/wake.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

// Constructor:
Wake::Wake(Wing &wing): wing(wing)
{
}

// Add new layer of wake panels:
void
Wake::add_layer(vector<Mesh*> &meshes_without_wakes)
{
    // Add layer of nodes at trailing edge, and add panels if necessary:
    int trailing_edge_n_nodes = wing.trailing_edge_nodes.size();
    int trailing_edge_n_panels = trailing_edge_n_nodes - 1;
    for (int k = 0; k < trailing_edge_n_nodes; k++) {
        Vector3d new_point = wing.nodes[wing.trailing_edge_nodes[k]];
        
        int node = n_nodes();
        nodes.push_back(new_point);
        
        node_deformation_velocities.push_back(Vector3d(0, 0, 0));
        
        if (k > 0 && nodes.size() > trailing_edge_n_nodes + 1) {
            vector<int> vertices;
            vertices.push_back(node - 1);
            vertices.push_back(node - 1 - trailing_edge_n_nodes);
            vertices.push_back(node - trailing_edge_n_nodes);
            vertices.push_back(node);
            
            int panel = n_panels();
            panel_nodes.push_back(vertices);
        
            vector<int> local_panel_neighbors;
            if (panel - trailing_edge_n_panels >= 0) {
                local_panel_neighbors.push_back(panel - trailing_edge_n_panels);
                panel_neighbors[panel - trailing_edge_n_panels].push_back(panel);
            }
            if (k > 1) {
                local_panel_neighbors.push_back(panel - 1);
                panel_neighbors[panel - 1].push_back(panel);
            }
            
            panel_neighbors.push_back(local_panel_neighbors);
            
            vector<int> *local_node_panel_neighbors = new vector<int>;
            local_node_panel_neighbors->push_back(panel);
            node_panel_neighbors.push_back(local_node_panel_neighbors);
            
            node_panel_neighbors[node - 1]->push_back(panel);
            if (node - 1 - trailing_edge_n_nodes >= 0) {
                node_panel_neighbors[node - trailing_edge_n_nodes]->push_back(panel); 
                node_panel_neighbors[node - 1 - trailing_edge_n_nodes]->push_back(panel);
            }
            
            doublet_coefficients.push_back(0);
            
            vector<double> panel_vortex_core_radii;
            for (int i = 0; i < 4; i++)
                panel_vortex_core_radii.push_back(Parameters::initial_vortex_core_radius);
            vortex_core_radii.push_back(panel_vortex_core_radii);
            
            vector<double> edge_lengths;
            for (int i = 0; i < 4; i++) {
                int prev_idx;
                if (i == 0)
                    prev_idx = 3;
                else
                    prev_idx = i - 1;
                    
                Vector3d &node_a = nodes[vertices[prev_idx]];
                Vector3d &node_b = nodes[vertices[i]];
                
                Vector3d edge = node_b - node_a;
                edge_lengths.push_back(edge.norm());
            }
            base_edge_lengths.push_back(edge_lengths);

        } else {
            vector<int> *empty = new vector<int>;
            node_panel_neighbors.push_back(empty);
        }
    }
        
    // Invalidate caches:
    invalidate_cache();
    
    for (int k = 0; k < meshes_without_wakes.size(); k++) {
        meshes_without_wakes[k]->doublet_influence_cache.erase(id);
        meshes_without_wakes[k]->source_influence_cache.erase(id);
        meshes_without_wakes[k]->source_unit_velocity_cache.erase(id);
        meshes_without_wakes[k]->vortex_ring_unit_velocity_cache.erase(id);
        meshes_without_wakes[k]->vortex_ring_ramasamy_leishman_velocity_cache.erase(id);
    }
}

// Geometry of trailing edge nodes:
void
Wake::translate_trailing_edge(Eigen::Vector3d translation)
{
    if (n_nodes() < wing.trailing_edge_nodes.size())
        return;
        
    int k0;
    if (Parameters::convect_wake)
        k0 = n_nodes() - wing.trailing_edge_nodes.size();
    else
        k0 = 0;
        
    for (int k = k0; k < n_nodes(); k++)                
        nodes[k] += translation;
        
    invalidate_cache();
}

void
Wake::transform_trailing_edge(Eigen::Matrix3d transformation)
{
    if (n_nodes() < wing.trailing_edge_nodes.size())
        return;
        
    int k0;
    if (Parameters::convect_wake)
        k0 = n_nodes() - wing.trailing_edge_nodes.size();
    else
        k0 = 0;
        
    for (int k = k0; k < n_nodes(); k++)                
        nodes[k] = transformation * nodes[k];
        
    invalidate_cache();
}

// Update Ramasamy-Leishman vortex core radii according to 
//    M. Ramasamy and J. G. Leishman, Reynolds Number Based Blade Tip Vortex Model, University of Maryland, 2005.
void
Wake::update_ramasamy_leishman_vortex_core_radii(int panel, double dt)
{
    for (int i = 0; i < 4; i++) {
        int prev_idx;
        if (i == 0)
            prev_idx = 3;
        else
            prev_idx = i - 1;
            
        double vortex_reynolds_number = fabs(doublet_coefficients[panel]) / Parameters::fluid_kinematic_viscosity;
        
        double t_multiplier = 4 * Parameters::lambs_constant *
            (1 + vortex_reynolds_number * Parameters::a_prime) * Parameters::fluid_kinematic_viscosity;
        
        double t = (pow(vortex_core_radii[panel][i], 2) - pow(Parameters::initial_vortex_core_radius, 2)) / t_multiplier;
        
        double vortex_core_size_0 = sqrt(pow(Parameters::initial_vortex_core_radius, 2) + t_multiplier * (t + dt));
        
        Vector3d node_a = nodes[panel_nodes[panel][prev_idx]];
        Vector3d node_b = nodes[panel_nodes[panel][i]];
        
        Vector3d edge = node_b - node_a;
        
        double strain = (edge.norm() - base_edge_lengths[panel][i]) / base_edge_lengths[panel][i];
        
        vortex_core_radii[panel][i] = fmax(Parameters::min_vortex_core_radius, vortex_core_size_0 / sqrt(1 + strain));
    }
}
