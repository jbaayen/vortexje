//
// Vortexje -- Wing.
//
// Copyright (C) 2012 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <iostream>
#include <algorithm>

#include <Eigen/Geometry>

#include <vortexje/wing-builder.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

// Constructor:
WingBuilder::WingBuilder(Wing &wing) : wing(wing)
{
}

// Return list of vectors representing a NACA airfoil in the (x, y)-plane.  See:
//    J. Moran, An Introduction to Theoretical and Computational Aerodynamics, Dover, 2010.
static double
cosine_rule(int n_points, int i)
{
    return 0.5 * (1 - cos(M_PI * i / (double) n_points));
}

// NACA airfoil generation:
static double
compute_y_c(double x, double max_camber, double max_camber_dist, double max_thickness, double chord)
{
    if (x < max_camber_dist * chord)
        return max_camber * x / pow(max_camber_dist, 2) * (2 * max_camber_dist - x / chord);
    else
        return max_camber * (chord - x) / pow(1 - max_camber_dist, 2) * (1 + x / chord - 2 * max_camber_dist);
}

static double
compute_theta(double x, double max_camber, double max_camber_dist, double max_thickness, double chord)
{
    if (x < max_camber_dist * chord)
        return atan(max_camber / pow(max_camber_dist, 2) * (2 * max_camber_dist - x / chord) 
                    - max_camber * x / pow(max_camber_dist, 2) / chord);
    else
        return atan(-max_camber / pow(1 - max_camber_dist, 2) * (1 + x / chord - 2 * max_camber_dist) 
                    + max_camber * (chord - x) / pow(1 - max_camber_dist, 2) / chord);
}

static double
compute_y_t(double x, double max_camber, double max_camber_dist, double max_thickness, double chord)
{
    return max_thickness / 0.2 * chord * (0.2969 * sqrt(x / chord) - 0.1260 * (x / chord) 
                                          - 0.3516 * pow(x / chord, 2) + 0.2843 * pow(x / chord, 3) - 0.1036 * pow(x / chord, 4));
}

vector<Vector3d>
WingBuilder::generate_naca_airfoil(double max_camber, double max_camber_dist, double max_thickness, double chord, int n_points, int &trailing_edge_point_id)
{
    if (n_points % 2 == 1) {
        cerr << "Wing::add_naca_airfoil(): n_nodes must be even." << endl;
        exit(1);
    }
    
    vector<Vector3d> airfoil_points;
    
    // Add upper nodes:
    for (int i = 0; i < n_points / 2; i++) {
        double x     = chord * cosine_rule(n_points / 2, i);
        
        double y_c   = compute_y_c  (x, max_camber, max_camber_dist, max_thickness, chord);
        double theta = compute_theta(x, max_camber, max_camber_dist, max_thickness, chord);
        double y_t   = compute_y_t  (x, max_camber, max_camber_dist, max_thickness, chord);
        
        Vector3d upper_point(x - y_t * sin(theta), y_c + y_t * cos(theta), 0.0);
        airfoil_points.push_back(upper_point);
    }
    
    // Add lower nodes:
    for (int i = 0; i < n_points / 2; i++) {
        double x     = chord * (1 - cosine_rule(n_points / 2, i));
        
        double y_c   = compute_y_c  (x, max_camber, max_camber_dist, max_thickness, chord);
        double theta = compute_theta(x, max_camber, max_camber_dist, max_thickness, chord);
        double y_t   = compute_y_t  (x, max_camber, max_camber_dist, max_thickness, chord);
        
        Vector3d lower_point(x + y_t * sin(theta), y_c - y_t * cos(theta), 0.0);
        airfoil_points.push_back(lower_point);
    }
    
    // Done.
    trailing_edge_point_id = n_points / 2;
    
    return airfoil_points;
}

// Clark-Y airfoil generation:
#define CLARKY_DATA_SIZE 61

static struct {
    double x;
    double y;
} clarky_upper_data[CLARKY_DATA_SIZE] = {
    {0.0000000, 0.0000000},
    {0.0005000, 0.0023390},
    {0.0010000, 0.0037271},
    {0.0020000, 0.0058025},
    {0.0040000, 0.0089238},
    {0.0080000, 0.0137350},
    {0.0120000, 0.0178581},
    {0.0200000, 0.0253735},
    {0.0300000, 0.0330215},
    {0.0400000, 0.0391283},
    {0.0500000, 0.0442753},
    {0.0600000, 0.0487571},
    {0.0800000, 0.0564308},
    {0.1000000, 0.0629981},
    {0.1200000, 0.0686204},
    {0.1400000, 0.0734360},
    {0.1600000, 0.0775707},
    {0.1800000, 0.0810687},
    {0.2000000, 0.0839202},
    {0.2200000, 0.0861433},
    {0.2400000, 0.0878308},
    {0.2600000, 0.0890840},
    {0.2800000, 0.0900016},
    {0.3000000, 0.0906804},
    {0.3200000, 0.0911857},
    {0.3400000, 0.0915079},
    {0.3600000, 0.0916266},
    {0.3800000, 0.0915212},
    {0.4000000, 0.0911712},
    {0.4200000, 0.0905657},
    {0.4400000, 0.0897175},
    {0.4600000, 0.0886427},
    {0.4800000, 0.0873572},
    {0.5000000, 0.0858772},
    {0.5200000, 0.0842145},
    {0.5400000, 0.0823712},
    {0.5600000, 0.0803480},
    {0.5800000, 0.0781451},
    {0.6000000, 0.0757633},
    {0.6200000, 0.0732055},
    {0.6400000, 0.0704822},
    {0.6600000, 0.0676046},
    {0.6800000, 0.0645843},
    {0.7000000, 0.0614329},
    {0.7200000, 0.0581599},
    {0.7400000, 0.0547675},
    {0.7600000, 0.0512565},
    {0.7800000, 0.0476281},
    {0.8000000, 0.0438836},
    {0.8200000, 0.0400245},
    {0.8400000, 0.0360536},
    {0.8600000, 0.0319740},
    {0.8800000, 0.0277891},
    {0.9000000, 0.0235025},
    {0.9200000, 0.0191156},
    {0.9400000, 0.0146239},
    {0.9600000, 0.0100232},
    {0.9700000, 0.0076868},
    {0.9800000, 0.0053335},
    {0.9900000, 0.0029690},
    {1.0000000, 0.0005993}
};

static struct {
    double x;
    double y;
} clarky_lower_data[CLARKY_DATA_SIZE] = {
    {0.0000000, 0.0000000},
    {0.0005000, -.0046700},
    {0.0010000, -.0059418},
    {0.0020000, -.0078113},
    {0.0040000, -.0105126},
    {0.0080000, -.0142862},
    {0.0120000, -.0169733},
    {0.0200000, -.0202723},
    {0.0300000, -.0226056},
    {0.0400000, -.0245211},
    {0.0500000, -.0260452},
    {0.0600000, -.0271277},
    {0.0800000, -.0284595},
    {0.1000000, -.0293786},
    {0.1200000, -.0299633},
    {0.1400000, -.0302404},
    {0.1600000, -.0302546},
    {0.1800000, -.0300490},
    {0.2000000, -.0296656},
    {0.2200000, -.0291445},
    {0.2400000, -.0285181},
    {0.2600000, -.0278164},
    {0.2800000, -.0270696},
    {0.3000000, -.0263079},
    {0.3200000, -.0255565},
    {0.3400000, -.0248176},
    {0.3600000, -.0240870},
    {0.3800000, -.0233606},
    {0.4000000, -.0226341},
    {0.4200000, -.0219042},
    {0.4400000, -.0211708},
    {0.4600000, -.0204353},
    {0.4800000, -.0196986},
    {0.5000000, -.0189619},
    {0.5200000, -.0182262},
    {0.5400000, -.0174914},
    {0.5600000, -.0167572},
    {0.5800000, -.0160232},
    {0.6000000, -.0152893},
    {0.6200000, -.0145551},
    {0.6400000, -.0138207},
    {0.6600000, -.0130862},
    {0.6800000, -.0123515},
    {0.7000000, -.0116169},
    {0.7200000, -.0108823},
    {0.7400000, -.0101478},
    {0.7600000, -.0094133},
    {0.7800000, -.0086788},
    {0.8000000, -.0079443},
    {0.8200000, -.0072098},
    {0.8400000, -.0064753},
    {0.8600000, -.0057408},
    {0.8800000, -.0050063},
    {0.9000000, -.0042718},
    {0.9200000, -.0035373},
    {0.9400000, -.0028028},
    {0.9600000, -.0020683},
    {0.9700000, -.0017011},
    {0.9800000, -.0013339},
    {0.9900000, -.0009666},
    {1.0000000, -.0005993}
};

vector<Vector3d>
WingBuilder::generate_clarky_airfoil(double chord, int n_points, int &trailing_edge_point_id)
{
    if (n_points % 2 == 1) {
        cerr << "Wing::add_naca_airfoil(): n_nodes must be even." << endl;
        exit(1);
    }
    
    vector<Vector3d> airfoil_points;
    
    // Add upper nodes:
    for (int i = 0; i < n_points / 2; i++) {
        double x = cosine_rule(n_points / 2, i);
        
        int clarky_upper_idx;
        for (clarky_upper_idx = 0; clarky_upper_idx < CLARKY_DATA_SIZE; clarky_upper_idx++) {
            if (x <= clarky_upper_data[clarky_upper_idx].x)
                break;
        }
        
        double y;
        if (clarky_upper_idx == 0)
            y = clarky_upper_data[0].y;
        else
            y = clarky_upper_data[clarky_upper_idx - 1].y + (clarky_upper_data[clarky_upper_idx].y - clarky_upper_data[clarky_upper_idx - 1].y) / (clarky_upper_data[clarky_upper_idx].x - clarky_upper_data[clarky_upper_idx - 1].x) * (x - clarky_upper_data[clarky_upper_idx - 1].x);
        
        Vector3d upper_point(chord * x, chord * y, 0.0);
        airfoil_points.push_back(upper_point);
    }
    
    // Add lower nodes:
    for (int i = 0; i < n_points / 2; i++) {
        double x = 1 - cosine_rule(n_points / 2, i);
        
        int clarky_lower_idx;
        for (clarky_lower_idx = 0; clarky_lower_idx < CLARKY_DATA_SIZE; clarky_lower_idx++) {
            if (x <= clarky_lower_data[clarky_lower_idx].x)
                break;
        }
        
        double y;
        if (clarky_lower_idx == 0)
            y = clarky_lower_data[0].y;
        else
            y = clarky_lower_data[clarky_lower_idx - 1].y + (clarky_lower_data[clarky_lower_idx].y - clarky_lower_data[clarky_lower_idx - 1].y) / (clarky_lower_data[clarky_lower_idx].x - clarky_lower_data[clarky_lower_idx - 1].x) * (x - clarky_lower_data[clarky_lower_idx - 1].x);
        
        Vector3d lower_point(chord * x, chord * y, 0.0);
        airfoil_points.push_back(lower_point);
    }
    
    // Done.
    trailing_edge_point_id = n_points / 2;
    
    return airfoil_points;
}

// Add points to wing:
vector<int>
WingBuilder::add_points(vector<Vector3d> &points, int trailing_edge_point_id)
{
    vector<int> added_nodes;
    
    for (int i = 0; i < points.size(); i++) {
        int node_id = wing.nodes.size();
        
        wing.nodes.push_back(points[i]);
        
        vector<int> *empty_vector = new vector<int>;
        wing.node_panel_neighbors.push_back(empty_vector);
        
        wing.node_deformation_velocities.push_back(Vector3d(0, 0, 0));
        
        if (i == trailing_edge_point_id)
            wing.trailing_edge_nodes.push_back(node_id);
            
        added_nodes.push_back(node_id);
    }
    
    return added_nodes;
}

// Connect two lists of nodes with panels.  If trianges = true, then triangles are used instead of quadrangles.
void
WingBuilder::connect_nodes(vector<int> &first_nodes, vector<int> &second_nodes,
                           int trailing_edge_point_id, int &trailing_edge_top_panel_id, int &trailing_edge_bottom_panel_id,
                           bool cyclic, ConnectNodesMode mode)
{
    for (int i = 0; i < first_nodes.size(); i++) {
        int next_i;
        if (i == first_nodes.size() - 1) {
            if (cyclic)
                next_i = 0;
            else
                break;
        } else
            next_i = i + 1;
            
        // Create panel(s)
        int panel_id, middle_node_id;
        Vector3d middle_point;
        Vector3d middle_deformation_velocity;
        vector<int> original_nodes;
        vector<int> *empty_vector;
        Vector3d vertices[4];
        Vector3d vertex_deformation_velocities[4];
        int new_nodes[4];
        switch (mode) {
        case TRIANGLES_A:
            panel_id;
            
            panel_id = wing.add_triangle(first_nodes[i], second_nodes[i], second_nodes[next_i]);
            
            // Mark as trailing edge panel if bordering trailing edge node.
            if (i == trailing_edge_point_id)
                trailing_edge_bottom_panel_id = panel_id;
                
            panel_id = wing.add_triangle(first_nodes[i], second_nodes[next_i], first_nodes[next_i]);
            
            // Mark as trailing edge panel if bordering trailing edge node.
            if (i == trailing_edge_point_id - 1)
                trailing_edge_top_panel_id = panel_id;

            break;
            
        case TRIANGLES_B:       
            panel_id = wing.add_triangle(first_nodes[i], second_nodes[i], first_nodes[next_i]);
            
            // Mark as trailing edge panel if bordering trailing edge node.
            if (i == trailing_edge_point_id)
                trailing_edge_bottom_panel_id = panel_id;
                
            panel_id = wing.add_triangle(second_nodes[i], second_nodes[next_i], first_nodes[next_i]);
            
            // Mark as trailing edge panel if bordering trailing edge node.
            if (i == trailing_edge_point_id - 1)
                trailing_edge_top_panel_id = panel_id;   
        
            break;
            
        case TRIANGLES_X:
            // Add node in the middle of the would-be quadrangle:
            middle_node_id = wing.nodes.size();

            middle_point = (wing.nodes[first_nodes[i]] + wing.nodes[first_nodes[next_i]] + wing.nodes[second_nodes[i]] + wing.nodes[second_nodes[next_i]]) / 4.0;
            wing.nodes.push_back(middle_point);
            
            middle_deformation_velocity = (wing.node_deformation_velocities[first_nodes[i]] + wing.node_deformation_velocities[first_nodes[next_i]] + wing.node_deformation_velocities[second_nodes[i]] + wing.node_deformation_velocities[second_nodes[next_i]]) / 4.0;
            wing.node_deformation_velocities.push_back(middle_deformation_velocity);
            
            empty_vector = new vector<int>;
            wing.node_panel_neighbors.push_back(empty_vector);
            
            // Add four triangles:
            panel_id = wing.add_triangle(first_nodes[i], second_nodes[i], middle_node_id);
            
            // Mark as trailing edge panel if bordering trailing edge node.
            if (i == trailing_edge_point_id)
                trailing_edge_bottom_panel_id = panel_id;
                
            panel_id = wing.add_triangle(second_nodes[next_i], first_nodes[next_i], middle_node_id);
            
            // Mark as trailing edge panel if bordering trailing edge node.
            if (i == trailing_edge_point_id - 1)
                trailing_edge_top_panel_id = panel_id;
                
            panel_id = wing.add_triangle(first_nodes[i], middle_node_id, first_nodes[next_i]);
            panel_id = wing.add_triangle(second_nodes[i], second_nodes[next_i], middle_node_id);
            
            break;
            
        case QUADRANGLES:
            // Add a planar trapezoidal panel, following
            //    T. Cebeci, An Engineering Approach to the Calculation of Aerodynamic Flows, Springer, 1999.
            original_nodes.push_back(first_nodes[i]);
            original_nodes.push_back(first_nodes[next_i]);
            original_nodes.push_back(second_nodes[i]);
            original_nodes.push_back(second_nodes[next_i]);
            
            Vector3d first_line  = wing.nodes[first_nodes[next_i]] - wing.nodes[first_nodes[i]];
            Vector3d second_line = wing.nodes[second_nodes[next_i]] - wing.nodes[second_nodes[i]];
            
            Vector3d first_line_deformation_velocity = wing.node_deformation_velocities[first_nodes[next_i]] - wing.node_deformation_velocities[first_nodes[i]];
            Vector3d second_line_deformation_velocity = wing.node_deformation_velocities[second_nodes[next_i]] - wing.node_deformation_velocities[second_nodes[i]];
            
            Vector3d line_direction = first_line + second_line;
            line_direction.normalize();
            
            Vector3d line_direction_deformation_velocity = (first_line_deformation_velocity + second_line_deformation_velocity) / (first_line + second_line).norm() - (first_line + second_line) * (first_line + second_line).dot(first_line_deformation_velocity + second_line_deformation_velocity) / pow((first_line + second_line).norm(), 3);
            
            Vector3d first_mid  = 0.5 * (wing.nodes[first_nodes[i]] + wing.nodes[first_nodes[next_i]]);
            Vector3d second_mid = 0.5 * (wing.nodes[second_nodes[i]] + wing.nodes[second_nodes[next_i]]);
            
            Vector3d first_mid_deformation_velocity = 0.5 * (wing.node_deformation_velocities[first_nodes[i]] + wing.node_deformation_velocities[first_nodes[next_i]]);
            Vector3d second_mid_deformation_velocity = 0.5 * (wing.node_deformation_velocities[second_nodes[i]] + wing.node_deformation_velocities[second_nodes[next_i]]);
            
            vertices[0] = first_mid - 0.5 * first_line.norm() * line_direction;
            vertices[1] = first_mid + 0.5 * first_line.norm() * line_direction;
            vertices[2] = second_mid - 0.5 * second_line.norm() * line_direction;
            vertices[3] = second_mid + 0.5 * second_line.norm() * line_direction;
            
            vertex_deformation_velocities[0] = first_mid_deformation_velocity - 0.5 * (first_line.dot(first_line_deformation_velocity) / first_line.norm() * line_direction + first_line.norm() * line_direction_deformation_velocity);
            vertex_deformation_velocities[1] = first_mid_deformation_velocity + 0.5 * (first_line.dot(first_line_deformation_velocity) / first_line.norm() * line_direction + first_line.norm() * line_direction_deformation_velocity);
            vertex_deformation_velocities[2] = second_mid_deformation_velocity - 0.5 * (second_line.dot(second_line_deformation_velocity) / second_line.norm() * line_direction + second_line.norm() * line_direction_deformation_velocity);
            vertex_deformation_velocities[3] = second_mid_deformation_velocity + 0.5 * (second_line.dot(second_line_deformation_velocity) / second_line.norm() * line_direction + second_line.norm() * line_direction_deformation_velocity);
            
            for (int j = 0; j < 4; j++) {
                // If the new points don't match the original ones, create new nodes:
                if ((vertices[j] - wing.nodes[original_nodes[j]]).norm() < Parameters::inversion_tolerance) {
                    new_nodes[j] = original_nodes[j];
                } else {         
                    new_nodes[j] = wing.nodes.size();
                    
                    wing.nodes.push_back(vertices[j]);
                    
                    wing.node_deformation_velocities.push_back(vertex_deformation_velocities[j]);
                    
                    wing.node_panel_neighbors.push_back(wing.node_panel_neighbors[original_nodes[j]]);
                }
            }
            
            panel_id = wing.add_quadrangle(new_nodes[0], new_nodes[2], new_nodes[3], new_nodes[1]);
            
            // Mark as trailing edge panel if bordering trailing edge node.
            if (i == trailing_edge_point_id - 1)
                trailing_edge_top_panel_id = panel_id;
            else if (i == trailing_edge_point_id)
                trailing_edge_bottom_panel_id = panel_id;
                
            break;
        }
    }
}

// Fill airfoil with panels.
vector<int>
WingBuilder::fill_airfoil(vector<int> airfoil_nodes, int trailing_edge_point_id, int z_sign)
{
    // Add middle nodes:
    vector<int> top_nodes;
    vector<int> bottom_nodes;
    vector<int> middle_nodes;
    for (int i = 1; i < trailing_edge_point_id; i++) {
        int top_node_id    = airfoil_nodes[i];
        int bottom_node_id = airfoil_nodes[airfoil_nodes.size() - i];
        
        Vector3d top_point    = wing.nodes[top_node_id];
        Vector3d bottom_point = wing.nodes[bottom_node_id];
        
        Vector3d top_deformation_velocity    = wing.node_deformation_velocities[top_node_id];
        Vector3d bottom_deformation_velocity = wing.node_deformation_velocities[bottom_node_id];
        
        Vector3d middle_point = 0.5 * (top_point + bottom_point);
        Vector3d middle_deformation_velocity = 0.5 * (top_deformation_velocity + bottom_deformation_velocity);
        
        int middle_node_id = wing.nodes.size();

        wing.nodes.push_back(middle_point);
        
        wing.node_deformation_velocities.push_back(middle_deformation_velocity);
        
        vector<int> *empty_vector = new vector<int>;
        wing.node_panel_neighbors.push_back(empty_vector);
        
        top_nodes.push_back(top_node_id);
        bottom_nodes.push_back(bottom_node_id);
        middle_nodes.push_back(middle_node_id);
    }
    
    // Close middle part with panels:
    int empty;
    if (z_sign == 1) {
        connect_nodes(middle_nodes, top_nodes, -1, empty, empty, false, QUADRANGLES);
        connect_nodes(bottom_nodes, middle_nodes, -1, empty, empty, false, QUADRANGLES);
    } else { 
        connect_nodes(top_nodes, middle_nodes, -1, empty, empty, false, QUADRANGLES);
        connect_nodes(middle_nodes, bottom_nodes, -1, empty, empty, false, QUADRANGLES);
    }
    
    // Create triangle for leading and trailing edges:
    if (z_sign == 1) { 
        wing.add_triangle(airfoil_nodes[0], airfoil_nodes[1], middle_nodes[0]);
        wing.add_triangle(airfoil_nodes[0], middle_nodes[0], airfoil_nodes[airfoil_nodes.size() - 1]);
        
        wing.add_triangle(middle_nodes[trailing_edge_point_id - 2], airfoil_nodes[trailing_edge_point_id - 1], airfoil_nodes[trailing_edge_point_id]);
        wing.add_triangle(middle_nodes[trailing_edge_point_id - 2], airfoil_nodes[trailing_edge_point_id], airfoil_nodes[trailing_edge_point_id + 1]);
    } else {
        wing.add_triangle(airfoil_nodes[0], middle_nodes[0], airfoil_nodes[1]);
        wing.add_triangle(airfoil_nodes[0], airfoil_nodes[airfoil_nodes.size() - 1], middle_nodes[0]);
        
        wing.add_triangle(middle_nodes[trailing_edge_point_id - 2], airfoil_nodes[trailing_edge_point_id], airfoil_nodes[trailing_edge_point_id - 1]);
        wing.add_triangle(middle_nodes[trailing_edge_point_id - 2], airfoil_nodes[trailing_edge_point_id + 1], airfoil_nodes[trailing_edge_point_id]);
    }
    
    // Done:
    return middle_nodes;
}
