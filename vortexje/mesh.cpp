//
// Vortexje -- Mesh.
//
// Copyright (C) 2012 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>
#include <cmath>

#include <Eigen/Geometry>

#include <vortexje/mesh.hpp>
#include <vortexje/parameters.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

// Static counter to give every mesh a unique ID.
static int id_counter = 0;

/**
   Constructs an empty mesh.
*/
Mesh::Mesh()
{
    // Set ID:
    id = ++id_counter;
}

/**
   Constructs a mesh from the given gmsh MSH file.
   
   @param[in]   file    Filename pointing to a gmsh MSH file.
*/
Mesh::Mesh(std::string file)
{
    // Set ID:
    id = ++id_counter;
    
    // Load mesh from file:
    load(file);
}

/**
   Clears the node-panel neighbor data structure, and frees up its memory.
*/
void
Mesh::clear_node_panel_neighbors()
{
    vector<vector<int> *> unique;
    
    for (int i = 0; i < (int) node_panel_neighbors.size(); i++) {
        bool found = false;
        for (int j = 0; j < (int) unique.size(); j++) {
            if (node_panel_neighbors[i] == unique[j]) {
                found = true;
                break;
            }
        }
        
        if (!found)
            unique.push_back(node_panel_neighbors[i]);
    }
    
    for (int i = 0; i < (int) unique.size(); i++)
        delete unique[i];
        
    node_panel_neighbors.clear();
}

/**
   Destructor.
*/
Mesh::~Mesh()
{
    clear_node_panel_neighbors();
}

/**
   Adds a triangle to this mesh, following the gmsh orientation convention.
    
   @param[in]  node_a  Node number to form vertex A.
   @param[in]  node_b  Node number to form vertex B.
   @param[in]  node_c  Node number to form vertex C.
    
   @returns New panel number.
    
   @note This method does not update the panel-panel neighbor data structure.  One must call
   compute_panel_neighbors() when done adding panels.
*/
int
Mesh::add_triangle(int node_a, int node_b, int node_c)
{
    vector<int> single_panel_nodes;
    single_panel_nodes.push_back(node_a);
    single_panel_nodes.push_back(node_b);
    single_panel_nodes.push_back(node_c);
    
    int panel_id = panel_nodes.size();
    panel_nodes.push_back(single_panel_nodes);
    
    node_panel_neighbors[node_a]->push_back(panel_id);
    node_panel_neighbors[node_b]->push_back(panel_id);
    node_panel_neighbors[node_c]->push_back(panel_id);
    
    return panel_id;
}

/**
   Adds a quadrangle to this mesh, following the gmsh orientation convention.
    
   @param[in]  node_a  Node number to form vertex A.
   @param[in]  node_b  Node number to form vertex B.
   @param[in]  node_c  Node number to form vertex C.
   @param[in]  node_d  Node number to form vertex D.
    
   @returns New panel number.
    
   @note This method does not update the panel-panel neighbor data structure.  One must call
   compute_panel_neighbors() when done adding panels.
*/
int
Mesh::add_quadrangle(int node_a, int node_b, int node_c, int node_d)
{
    vector<int> single_panel_nodes;
    single_panel_nodes.push_back(node_a);
    single_panel_nodes.push_back(node_b);
    single_panel_nodes.push_back(node_c);       
    single_panel_nodes.push_back(node_d);
    
    int panel_id = panel_nodes.size();
    panel_nodes.push_back(single_panel_nodes);
    
    node_panel_neighbors[node_a]->push_back(panel_id);
    node_panel_neighbors[node_b]->push_back(panel_id);
    node_panel_neighbors[node_c]->push_back(panel_id);
    node_panel_neighbors[node_d]->push_back(panel_id);
    
    return panel_id;
}

/**
   Reconstructs a mesh from a gmsh MSH file.
   
   @param[in]   file   gmsh MSH file filename.
   
   @returns true on success.
*/
bool
Mesh::load(string file)
{
    cout << "Mesh " << id << ": Loading from " << file << "." << endl;
        
    // Clear cache:
    invalidate_cache();
    
    // Load mesh from gmsh MSH file:
    ifstream f;
    f.open(file.c_str());
    
    bool in_nodes = false, in_elements = false;
    int current_panel = 0;
    while (f.good()) {
        string line;
        getline(f, line);
        
        if (line[0] == '$') {
            if        (line == "$MeshFormat"   ) {
                getline(f, line);
                
                istringstream tokens(line);
                
                string version;
                int file_type, data_size;
                tokens >> version >> file_type >> data_size;
                
                if (version != "2.2" || file_type != 0 || data_size != 8) {
                    cerr << "Mesh " << id << ": Unknown data format in " << f << "." << endl;
                    
                    f.close();
                    
                    return false;
                }
                
                getline(f, line);
                
            } else if (line == "$Nodes"        ) {
                getline(f, line);
                
                in_nodes = true;
              
            } else if (line == "$EndNodes"     ) {
                in_nodes = false;
                
            } else if (line == "$Elements"     ) {
                getline(f, line);
                
                in_elements = true;
                
            } else if (line == "$EndElements"  ) {
                in_elements = false;
                
            }
            
        } else if (in_nodes) {
            istringstream tokens(line);
            
            int node_number;
            double x, y, z;
            tokens >> node_number >> x >> y >> z;
            
            nodes.push_back(Vector3d(x, y, z));            
            
            vector<int> *neighbor_list = new vector<int>;
            node_panel_neighbors.push_back(neighbor_list);
            
            node_deformation_velocities.push_back(Vector3d(0, 0, 0));
            
        } else if (in_elements) {
            istringstream tokens(line);
            
            int single_number, single_type, number_of_tags, number_of_nodes;
            tokens >> single_number >> single_type >> number_of_tags;
            
            switch (single_type) {
            case 2:
                number_of_nodes = 3;
                break;
            case 3:
                number_of_nodes = 4;
                break;
            default:
                continue; // Only read triangles and quatrangles.
            }
            
            for (int i = 0; i < number_of_tags; i++) {
                int tag;
                tokens >> tag;
            }
            
            vector<int> single_panel_nodes;
            for (int i = 0; i < number_of_nodes; i++) {
                int node;
                tokens >> node;
                
                single_panel_nodes.push_back(node - 1);
                
                node_panel_neighbors[node - 1]->push_back(current_panel);
            }

            panel_nodes.push_back(single_panel_nodes);
            
            current_panel++;
        }  
    }
    
    f.close();
    
    compute_panel_neighbors();
    
    return true;
}

/**
   Computes neighboring panels of panels, based on existing node-panel data structures.
*/
void
Mesh::compute_panel_neighbors()
{   
    // Compute panel neighbors:
    for (int i = 0; i < (int) panel_nodes.size(); i++) {
        vector<int> single_panel_potential_neighbors;
        vector<int> single_panel_neighbors;
        
        for (int j = 0; j < (int) panel_nodes[i].size(); j++) {
            int node = panel_nodes[i][j];
            for (int k = 0; k < (int) node_panel_neighbors[node]->size(); k++) {
                int potential_neighbor = (*node_panel_neighbors[node])[k];
                if (potential_neighbor == i)
                    continue;
                
                bool found = false;
                for (int l = 0; l < (int) single_panel_potential_neighbors.size(); l++) {
                    if (potential_neighbor == single_panel_potential_neighbors[l]) {
                        found = true;
                        break;
                    }
                }
                
                if (found)
                    single_panel_neighbors.push_back(potential_neighbor);
                else
                    single_panel_potential_neighbors.push_back(potential_neighbor);
            }
        }
        
        panel_neighbors.push_back(single_panel_neighbors);
    }
}

/**
   Saves this mesh to a gmsh MSH file, including data vectors associating numerical values to each panel -- views,
   in gmsh terminology.
  
   @param[in]   file            Destination filename.
   @param[in]   view_names      List of names of data vectors to be stored.
   @param[in]   view_data       List of data vectors to be stored.
   @param[in]   node_offset     Node numbering offset in output file.
   @param[in]   panel_offset    Panel numbering offset in output file.
*/
void
Mesh::save(std::string file, std::vector<std::string> &view_names, std::vector<Eigen::VectorXd> &view_data, int node_offset, int panel_offset)
{
    cout << "Mesh " << id << ": Saving to " << file << "." << endl;
    
    // Save mesh to gmsh file:
    ofstream f;
    f.open(file.c_str());
    
    f << "$MeshFormat" << endl;
    f << "2.2 0 8" << endl; 
    f << "$EndMeshFormat" << endl;
    
    f << "$Nodes" << endl;
    f << n_nodes() << endl;
    
    for (int i = 0; i < n_nodes(); i++) {
        f << i + node_offset + 1;
        
        for (int j = 0; j < 3; j++) {
            f << ' ';
            f << nodes[i](j);
        }
        
        f << endl;
    }
    
    f << "$EndNodes" << endl;
    f << "$Elements" << endl;
    f << n_panels() << endl;
    
    for (int i = 0; i < n_panels(); i++) {
        int element_type;
        switch (panel_nodes[i].size()) {
        case 3:
            element_type = 2;
            break;
        case 4:
            element_type = 3;
            break;
        default:
            cerr << "Mesh " << id << ": Unknown polygon at panel " << i << "." << endl;
            continue;
        }
        
        f << i + panel_offset + 1;
        
        f << ' ';
        
        f << element_type;
        
        f << ' ';
        
        f << 0;
        
        for (int j = 0; j < (int) panel_nodes[i].size(); j++) {
            f << ' ';
            f << panel_nodes[i][j] + node_offset + 1;
        }
        
        f << endl;
    }
    
    f << "$EndElements" << endl;
    
    for (int k = 0; k < (int) view_names.size(); k++) {
        f << "$ElementData" << endl;
        f << "1" << endl;
        f << '"' << view_names[k] << '"' << endl;
        f << "1" << endl;
        f << 0.0 << endl;
        f << "3" << endl;
        f << 0 << endl;
        f << 1 << endl;
        f << n_panels() << endl;
        
        for (int i = 0; i < n_panels(); i++) {            
            f << i + panel_offset + 1;
            
            f << ' ';
            
            f << view_data[k][i];
            
            f << endl;       
        }
        
        f << "$EndElementData" << endl;
    }
    
    f.close();
}

/**
   Saves this mesh to a gmsh MSH file.
  
   @param[in]   file            Destination filename.
   @param[in]   node_offset     Node numbering offset in output file.
   @param[in]   panel_offset    Panel numbering offset in output file.
*/
void
Mesh::save(std::string file, int node_offset, int panel_offset)
{
    vector<string> empty_names;
    vector<VectorXd> empty_data;
    
    save(file, empty_names, empty_data, node_offset, panel_offset);
}

/**
   Returns the number of nodes contained in this mesh.
   
   @returns Number of nodes.
*/
int
Mesh::n_nodes()
{
    return nodes.size();
}

/**
   Returns the number of panels contained in this mesh.
   
   @returns Number of panels.
*/
int
Mesh::n_panels()
{
    return panel_nodes.size();
}

/**
   Rotates this mesh.
   
   @param[in]   axis                    Axis of rotation.
   @param[in]   angle                   Angle of rotation.
*/
void
Mesh::rotate(Eigen::Vector3d &axis, double angle)
{
    Eigen::Matrix3d transformation = AngleAxis<double>(angle, axis).toRotationMatrix();
    transform(transformation);
}

/**
   Transforms this mesh.
   
   @param[in]   transformation             Transformation matrix.
*/
void
Mesh::transform(Eigen::Matrix3d &transformation)
{
    vector<Mesh*> empty;
    transform(transformation, empty);
}

/**
   Translates this mesh.
   
   @param[in]   translation             Translation vector.
*/
void
Mesh::translate(Eigen::Vector3d &translation)
{
    vector<Mesh*> empty;
    translate(translation, empty);
}

/**
   Rotates this mesh, under the assumption that the list of co-rotating meshes experiences
   the same rotation.  The caches of the co-rotating meshes are updated accordingly.
   
   @param[in]   axis                    Axis of rotation.
   @param[in]   angle                   Angle of rotation.
   @param[in]   corotating_meshes       List of co-rotating meshes.
*/
void
Mesh::rotate(Eigen::Vector3d &axis, double angle, std::vector<Mesh*> &corotating_meshes)
{
    Eigen::Matrix3d transformation = AngleAxis<double>(angle, axis).toRotationMatrix();
    transform(transformation, corotating_meshes);
}

/**
   Transforms this mesh, under the assumption that the list of co-transforming meshes experiences
   the same transformation.  The caches of the co-transforming meshes are updated accordingly.
   
   @param[in]   transformation          Transformation matrix.
   @param[in]   cotransforming_meshes   List of co-transforming meshes.
*/
void
Mesh::transform(Eigen::Matrix3d &transformation, std::vector<Mesh*> &cotransforming_meshes)
{
    for (int i = 0; i < n_nodes(); i++)
        nodes[i] = transformation * nodes[i];
            
    for (int j = 0; j < 2; j++) {
        if (!panel_collocation_point_cache[j].empty()) {
            for (int i = 0; i < n_panels(); i++)
                panel_collocation_point_cache[j][i] = transformation * panel_collocation_point_cache[j][i];
        }
    }
    
    if (!panel_normal_cache.empty()) {
        for (int i = 0; i < n_panels(); i++)
            panel_normal_cache[i] = transformation * panel_normal_cache[i];
    }
    
    std::vector<std::vector<double> > *doublet_influence_backup =
        new std::vector<std::vector<double> >[cotransforming_meshes.size()];
    std::vector<std::vector<double> > *source_influence_backup  =
        new std::vector<std::vector<double> >[cotransforming_meshes.size()];
    std::vector<std::vector<Vector3d, Eigen::aligned_allocator<Vector3d> > > *vortex_ring_unit_velocity_backup =
        new std::vector<std::vector<Vector3d, Eigen::aligned_allocator<Vector3d> > >[cotransforming_meshes.size()];
    std::vector<std::vector<Vector3d, Eigen::aligned_allocator<Vector3d> > > *source_unit_velocity_backup      =
        new std::vector<std::vector<Vector3d, Eigen::aligned_allocator<Vector3d> > >[cotransforming_meshes.size()];
    
    for (int i = 0; i < (int) cotransforming_meshes.size(); i++) {
        // Ensure the caches exist.
        if (n_panels() > 0 && cotransforming_meshes[i]->n_panels() > 0) {
            doublet_influence(*cotransforming_meshes[i], 0, 0);
            source_influence(*cotransforming_meshes[i], 0, 0);
            vortex_ring_unit_velocity(*cotransforming_meshes[i], 0, 0);
            source_unit_velocity(*cotransforming_meshes[i], 0, 0);
        }
        
        doublet_influence_backup[i]         = doublet_influence_cache[cotransforming_meshes[i]->id];
        source_influence_backup[i]          = source_influence_cache[cotransforming_meshes[i]->id];
        vortex_ring_unit_velocity_backup[i] = vortex_ring_unit_velocity_cache[cotransforming_meshes[i]->id];
        source_unit_velocity_backup[i]      = source_unit_velocity_cache[cotransforming_meshes[i]->id];
    }
    
    doublet_influence_cache.clear();
    source_influence_cache.clear();
    vortex_ring_unit_velocity_cache.clear();
    source_unit_velocity_cache.clear();
    vortex_ring_ramasamy_leishman_velocity_cache.clear();
    
    for (int i = 0; i < (int) cotransforming_meshes.size(); i++) {
        doublet_influence_cache[cotransforming_meshes[i]->id]         = doublet_influence_backup[i];
        source_influence_cache[cotransforming_meshes[i]->id]          = source_influence_backup[i];
        vortex_ring_unit_velocity_cache[cotransforming_meshes[i]->id] = vortex_ring_unit_velocity_backup[i];
        source_unit_velocity_cache[cotransforming_meshes[i]->id]      = source_unit_velocity_backup[i];
    }
    
    delete[] doublet_influence_backup;
    delete[] source_influence_backup;
    delete[] vortex_ring_unit_velocity_backup;
    delete[] source_unit_velocity_backup;
}

/**
   Translates this mesh, under the assumption that the list of co-translating meshes experiences
   the same translation.  The caches of the co-translating meshes are updated accordingly.
   
   @param[in]   translation             Translation vector.
   @param[in]   cotranslating_meshes    List of co-translating meshes.
*/
void
Mesh::translate(Eigen::Vector3d &translation, std::vector<Mesh*> &cotranslating_meshes)
{
    for (int i = 0; i < n_nodes(); i++)
        nodes[i] = nodes[i] + translation;
            
    for (int j = 0; j < 2; j++) {
        if (!panel_collocation_point_cache[j].empty()) {
            for (int i = 0; i < n_panels(); i++)
                panel_collocation_point_cache[j][i] = panel_collocation_point_cache[j][i] + translation;
        }
    }
    
    std::vector<std::vector<double> > *doublet_influence_backup =
        new std::vector<std::vector<double> >[cotranslating_meshes.size()];
    std::vector<std::vector<double> > *source_influence_backup  =
        new std::vector<std::vector<double> >[cotranslating_meshes.size()];
    std::vector<std::vector<Vector3d, Eigen::aligned_allocator<Vector3d> > > *vortex_ring_unit_velocity_backup =
        new std::vector<std::vector<Vector3d, Eigen::aligned_allocator<Vector3d> > >[cotranslating_meshes.size()];
    std::vector<std::vector<Vector3d, Eigen::aligned_allocator<Vector3d> > > *source_unit_velocity_backup      =
        new std::vector<std::vector<Vector3d, Eigen::aligned_allocator<Vector3d> > >[cotranslating_meshes.size()];
    
    for (int i = 0; i < (int) cotranslating_meshes.size(); i++) {
        // Ensure the caches exist.
        if (n_panels() > 0 && cotranslating_meshes[i]->n_panels() > 0) {
            doublet_influence(*cotranslating_meshes[i], 0, 0);
            source_influence(*cotranslating_meshes[i], 0, 0);
            vortex_ring_unit_velocity(*cotranslating_meshes[i], 0, 0);
            source_unit_velocity(*cotranslating_meshes[i], 0, 0);
        }
        
        doublet_influence_backup[i]         = doublet_influence_cache[cotranslating_meshes[i]->id];
        source_influence_backup[i]          = source_influence_cache[cotranslating_meshes[i]->id];
        vortex_ring_unit_velocity_backup[i] = vortex_ring_unit_velocity_cache[cotranslating_meshes[i]->id];
        source_unit_velocity_backup[i]      = source_unit_velocity_cache[cotranslating_meshes[i]->id];
    }
    
    doublet_influence_cache.clear();
    source_influence_cache.clear();
    vortex_ring_unit_velocity_cache.clear();
    source_unit_velocity_cache.clear();
    vortex_ring_ramasamy_leishman_velocity_cache.clear();
    
    for (int i = 0; i < (int) cotranslating_meshes.size(); i++) {
        doublet_influence_cache[cotranslating_meshes[i]->id]         = doublet_influence_backup[i];
        source_influence_cache[cotranslating_meshes[i]->id]          = source_influence_backup[i];
        vortex_ring_unit_velocity_cache[cotranslating_meshes[i]->id] = vortex_ring_unit_velocity_backup[i];
        source_unit_velocity_cache[cotranslating_meshes[i]->id]      = source_unit_velocity_backup[i];
    }
    
    delete[] doublet_influence_backup;
    delete[] source_influence_backup;
    delete[] vortex_ring_unit_velocity_backup;
    delete[] source_unit_velocity_backup;
}

/**
   Computes the distance between a point and a given panel.
   
   @param[in]   x       Reference point.
   @param[in]   panel   Panel to measure distance to.
   
   @returns Distance between reference point and panel.
*/
double
Mesh::distance_to_panel(Eigen::Vector3d &x, int panel)
{
    double distance = numeric_limits<double>::max();
    
    vector<int> &single_panel_nodes = panel_nodes[panel];
    
    // Distance to vertices:
    Vector3d mid(0, 0, 0);
    for (int i = 0; i < (int) single_panel_nodes.size(); i++) {           
        Vector3d &vertex = nodes[single_panel_nodes[i]];
        
        double vertex_distance = (x - vertex).norm();
        if (vertex_distance < distance)
            distance = vertex_distance;
            
        mid += vertex;
    }
    mid /= single_panel_nodes.size();
    
    // Distance to edges:
    for (int i = 0; i < (int) single_panel_nodes.size(); i++) {
        int prev_idx;
        if (i == 0)
            prev_idx = single_panel_nodes.size() - 1;
        else
            prev_idx = i - 1;
            
        Vector3d &vertex = nodes[single_panel_nodes[i]];
        Vector3d &prev_vertex = nodes[single_panel_nodes[prev_idx]];
        
        Vector3d edge = vertex - prev_vertex;
        Vector3d edge_direction = edge / edge.norm();
            
        Vector3d vector_to_vertex = vertex - x;
        Vector3d vector_to_edge = vector_to_vertex - vector_to_vertex.dot(edge_direction) * edge_direction;
        
        Vector3d projection_on_edge = x + vector_to_edge;
        
        Vector3d delta_vertex = projection_on_edge - vertex;
        Vector3d delta_prev_vertex = projection_on_edge - prev_vertex;
        if (delta_vertex.dot(delta_prev_vertex) <= 0) {
            double edge_distance = vector_to_edge.norm();
            if (edge_distance < distance)
                distance = edge_distance;
        }
    }
    
    // Distance to interior:
    Vector3d normal = panel_normal(panel);
    Vector3d vector_to_plane = normal.dot(mid - x) * normal;
    Vector3d projection_on_plane = x + vector_to_plane;
    
    int n_right_of_side = 0;
    for (int i = 0; i < (int) single_panel_nodes.size(); i++) {
        int prev_idx;
        if (i == 0)
            prev_idx = single_panel_nodes.size() - 1;
        else
            prev_idx = i - 1;
            
        Vector3d &vertex = nodes[single_panel_nodes[i]];
        Vector3d &prev_vertex = nodes[single_panel_nodes[prev_idx]];
        
        Vector3d edge = vertex - prev_vertex;
        Vector3d edge_direction = edge / edge.norm();
        Vector3d edge_normal = normal.cross(edge_direction);
        
        Vector3d vector_to_edge = vertex - projection_on_plane;
        if (vector_to_edge.dot(edge_normal) >= 0)
            n_right_of_side++;
    }
    
    if (n_right_of_side == (int) single_panel_nodes.size()) {
        double plane_distance = vector_to_plane.norm();
        if (plane_distance < distance)
            distance = plane_distance;
    }
    
    // Done:
    return distance;
}

/**
   Finds the panel closest to the given point, and reports the distance.
   
   @param[in]   x           Reference point.
   @param[out]  panel       Closest panel number.
   @param[out]  distance    Distance to closest panel.
   
   @returns true if the closest panel borders a trailing edge.
*/
bool
Mesh::closest_panel(Eigen::Vector3d &x, int &panel, double &distance)
{
    distance = numeric_limits<double>::max();
    panel = -1;
    
    for (int i = 0; i < n_panels(); i++) {
        double distance_candidate = distance_to_panel(x, i);
        if (distance_candidate < distance) {
            panel    = i;
            distance = distance_candidate;
        }
    }

    return false;
}

/**
   Computes the collocation point of a given panel.  If the point has been computed before,
   the cached collocation point is returned.  If not, the collocation point is computed and cached.
   
   @param[in]   panel           Panel of which the collocation point is evaluated.
   @param[in]   below_surface   true to request the collocation point lying underneath the surface.
   
   @returns Collocation point.
*/
Vector3d
Mesh::panel_collocation_point(int panel, bool below_surface)
{
    if (panel_collocation_point_cache[below_surface].empty()) {
        cout << "Mesh " << id << ": Generating panel collocation point cache " << below_surface << "." << endl;
        
        panel_collocation_point_cache[below_surface].reserve(n_panels());
        
        for (int i = 0; i < n_panels(); i++) {
            vector<int> &single_panel_nodes = panel_nodes[i];
            
            Vector3d collocation_point(0, 0, 0);
            double max_edge = numeric_limits<double>::min();
            double min_edge = numeric_limits<double>::max();
            for (int j = 0; j < (int) single_panel_nodes.size(); j++) {
                collocation_point = collocation_point + nodes[single_panel_nodes[j]];
                
                int prev_idx;
                if (j == 0)
                    prev_idx = single_panel_nodes.size() - 1;
                else
                    prev_idx = j - 1;
                    
                Vector3d edge = nodes[single_panel_nodes[j]] - nodes[single_panel_nodes[prev_idx]];
                if (edge.norm() > max_edge)
                    max_edge = edge.norm();
                if (edge.norm() < min_edge)
                    min_edge = edge.norm();
            }
            collocation_point = collocation_point / single_panel_nodes.size();
            
            if (below_surface) {
                double collocation_point_delta = Parameters::collocation_point_delta_factor * min_edge;
                collocation_point = collocation_point + collocation_point_delta * panel_normal(i);
            }
                
            panel_collocation_point_cache[below_surface].push_back(collocation_point);
        }
    }
    
    return panel_collocation_point_cache[below_surface][panel];
}

/**
   Computes the inward-pointing normal of a given panel.  If the normal has been computed before,
   the cached normal is returned.  If not, the normal is computed and cached.
   
   @param[in]   panel   Panel of which the inward-pointing normal is evaluated.
   
   @returns Inward-pointing normal.
*/
Vector3d
Mesh::panel_normal(int panel)
{
    if (panel_normal_cache.empty()) {
        cout << "Mesh " << id << ": Generating panel normal cache." << endl;
        
        panel_normal_cache.reserve(n_panels());
        
        for (int i = 0; i < n_panels(); i++) {
            vector<int> &single_panel_nodes = panel_nodes[i];
            
            Vector3d normal;
            if (single_panel_nodes.size() == 3) {
                Vector3d AB = nodes[single_panel_nodes[1]] - nodes[single_panel_nodes[0]];
                Vector3d AC = nodes[single_panel_nodes[2]] - nodes[single_panel_nodes[0]];
                
                normal = AB.cross(AC);
                
            } else { // 4 sides
                Vector3d AC = nodes[single_panel_nodes[2]] - nodes[single_panel_nodes[0]];
                Vector3d BD = nodes[single_panel_nodes[3]] - nodes[single_panel_nodes[1]];
                
                normal = AC.cross(BD);
            }

            normal = normal / normal.norm();

            panel_normal_cache.push_back(normal);
        }
    }
    
    return panel_normal_cache[panel];
}

/**
   Computes the surface area of a given panel.  If the surface area has been computed before,
   the cached surface area is returned.  If not, the surface area is computed and cached.
   
   @param[in]   panel   Panel of which the surface area is evaluated.
   
   @returns Panel surface area.
*/
double
Mesh::panel_surface_area(int panel)
{
    if (panel_surface_area_cache.empty()) {
        cout << "Mesh " << id << ": Generating panel surface area cache." << endl;
        
        panel_surface_area_cache.reserve(n_panels());
        
        for (int i = 0; i < n_panels(); i++) {
            vector<int> &single_panel_nodes = panel_nodes[i];
            
            double surface_area = 0.0;
            if (single_panel_nodes.size() == 3) {
                Vector3d AB = nodes[single_panel_nodes[1]] - nodes[single_panel_nodes[0]];
                Vector3d AC = nodes[single_panel_nodes[2]] - nodes[single_panel_nodes[0]];
                
                surface_area = 0.5 * AB.cross(AC).norm();
                
            } else { // 4 sides
                Vector3d AC = nodes[single_panel_nodes[2]] - nodes[single_panel_nodes[0]];
                Vector3d BD = nodes[single_panel_nodes[3]] - nodes[single_panel_nodes[1]];
                
                surface_area = 0.5 * AC.cross(BD).norm();
            }
            
            panel_surface_area_cache.push_back(surface_area);
        }
    }
    
    return panel_surface_area_cache[panel];
}

/**
   Computes the diameter of a given panel.  If the diameter has been computed before,
   the cached diameter is returned.  If not, the diameter is computed and cached.
   
   @param[in]   panel   Panel of which the diameter is evaluated.
   
   @returns Panel diameter.
*/
double
Mesh::panel_diameter(int panel)
{
    if (panel_diameter_cache.empty()) {
        cout << "Mesh " << id << ": Generating panel diameter cache." << endl;
        
        panel_diameter_cache.reserve(n_panels());
        
        for (int i = 0; i < n_panels(); i++) {
            double diameter = numeric_limits<double>::min();
            
            for (int j = 0; j < (int) panel_nodes[i].size(); j++) {
                Vector3d a = nodes[panel_nodes[i][j]];
                
                for (int k = 0; k < j; k++) {
                    Vector3d b = nodes[panel_nodes[i][k]];
                    
                    double diameter_candidate = (b - a).norm();
                    if (diameter_candidate > diameter)
                        diameter = diameter_candidate;
                }
            }
            
            panel_diameter_cache.push_back(diameter);
        }
    }
    
    return panel_diameter_cache[panel];
}

/**
   Computes the deformation velocity of the given panel.
   
   @param[in]   panel   Panel number.
   
   @returns Deformation velocity.
*/
Vector3d
Mesh::panel_deformation_velocity(int panel)
{
    Vector3d deformation_velocity(0, 0, 0);
    for (int i = 0; i < (int) panel_nodes[panel].size(); i++)
        deformation_velocity += node_deformation_velocities[panel_nodes[panel][i]];
    return deformation_velocity / (double) panel_nodes[panel].size();
}

/**
   Computes a point, located outside of the body, that lies close to a given node.
   
   This method is used by the wake-wake and wake-body interaction code to evaluate velocities and
   influence coefficients close to the body, where the solutions would otherwise become numerically singular.
   
   @param[in]   node    Node number.
   
   @returns Point, located outside of the body, close to the given node.
*/
Vector3d
Mesh::close_to_body_point(int node)
{
    Vector3d layer_direction(0, 0, 0);
    for (int i = 0; i < (int) node_panel_neighbors[node]->size(); i++)
        layer_direction += panel_normal((*node_panel_neighbors[node])[i]);
    layer_direction /= layer_direction.norm();
    
    return nodes[node] - Parameters::interpolation_layer_thickness * layer_direction;
}

/**
   Computes the on-body gradient of a scalar field.
   
   @param[in]   scalar_field    Scalar field, ordered by panel number.
   @param[in]   panel           Panel on which the on-body gradient is evaluated.
   
   @returns On-body gradient.
*/
Vector3d
Mesh::scalar_field_gradient(Eigen::VectorXd &scalar_field, int panel)
{
    Vector3d x = panel_collocation_point(panel, false);
    double v = scalar_field(panel);
    
    Vector3d normal = panel_normal(panel);
    
    Vector3d gradient(0, 0, 0);
        
    for (int i = 0; i < 3; i++) {
        Vector3d derivative_direction(0, 0, 0);
        derivative_direction(i) = 1.0;
        
        double a = 0.0, b = 0.0;
        
        for (int n = 0; n < (int) panel_neighbors[panel].size(); n++) {
            Vector3d neighbor_normal = panel_normal(panel_neighbors[panel][n]);
            if (normal.dot(neighbor_normal) < 0)
                continue; // Don't differentiate along sharp corners.
            
            Vector3d neighbor_vector = panel_collocation_point(panel_neighbors[panel][n], false) - x;
            double neighbor_distance = neighbor_vector.norm();
            if (neighbor_distance < Parameters::inversion_tolerance)
                continue; // Don't differentiate along too small distances.
                
            Vector3d neighbor_direction = neighbor_vector / neighbor_distance;

            double neighbor_gradient = (scalar_field(panel_neighbors[panel][n]) - v) / neighbor_distance;
            double direction_coefficient = neighbor_direction.dot(derivative_direction);
            
            a += direction_coefficient * neighbor_gradient;
            b += fabs(direction_coefficient);
        }
        
        if (fabs(b) > Parameters::inversion_tolerance)
            gradient += derivative_direction * a / b;
    }
    
    gradient -= gradient.dot(normal) * normal;

    return gradient;
}

// Compute a matrix that rotates x to y.
static Matrix3d
x_to_y_rotation(const Vector3d &unit_x, const Vector3d &unit_y)
{   
    Quaterniond q;
    q.setFromTwoVectors(unit_x, unit_y);
    return q.toRotationMatrix();
}

// Compute influence of doublet panel edge on given point.
static double
doublet_edge_influence(Vector3d &x, Vector3d &node_a, Vector3d &node_b)
{   
    double r1 = (x - node_a).norm();
    double r2 = (x - node_b).norm();
    
    if (fabs(node_b(0) - node_a(0)) < Parameters::inversion_tolerance)
        return 0.0;
    
    double m = (node_b(1) - node_a(1)) / (node_b(0) - node_a(0));
    
    double e1 = pow(x(0) - node_a(0), 2) + pow(x(2), 2);
    double e2 = pow(x(0) - node_b(0), 2) + pow(x(2), 2);
    
    double h1 = (x(0) - node_a(0)) * (x(1) - node_a(1));
    double h2 = (x(0) - node_b(0)) * (x(1) - node_b(1));

    double delta_theta = atan((m * e1 - h1) / (x(2) * r1)) - atan((m * e2 - h2) / (x(2) * r2));
    
    return delta_theta;
}

/**
   Computes the potential influence induced by a doublet panel of unit strength.  
   
   @param[in]   x           Point at which the influence coefficient is evaluated.
   @param[in]   this_panel  Panel on which the doublet panel is located.
   
   @returns Influence coefficient.
*/
double
Mesh::doublet_influence(Eigen::Vector3d &x, int this_panel)
{
    // Transform such that panel normal becomes unit Z vector:
    Matrix3d rotation = x_to_y_rotation(panel_normal(this_panel), Vector3d::UnitZ());
    
    Vector3d x_normalized = rotation * x;
    
    vector<Vector3d, Eigen::aligned_allocator<Vector3d> > normalized_panel_nodes;
    for (int i = 0; i < (int) panel_nodes[this_panel].size(); i++)
        normalized_panel_nodes.push_back(rotation * nodes[panel_nodes[this_panel][i]]);
    
    Vector3d translation(0, 0, -normalized_panel_nodes[0](2));
    
    x_normalized += translation;
    for (int i = 0; i < (int) panel_nodes[this_panel].size(); i++)
        normalized_panel_nodes[i] += translation;
    
    // Compute influence coefficient according to Hess:
    double influence = 0.0;
    for (int i = 0; i < (int) normalized_panel_nodes.size(); i++) {
        int prev_idx;
        if (i == 0)
            prev_idx = normalized_panel_nodes.size() - 1;
        else
            prev_idx = i - 1;
            
        Vector3d &node_a = normalized_panel_nodes[prev_idx];
        Vector3d &node_b = normalized_panel_nodes[i];
        
        influence += doublet_edge_influence(x_normalized, node_a, node_b);
    }
    
    return influence / (4 * M_PI);
}

// Compute influence of source panel edge on given point.
static double
source_edge_influence(Vector3d &x, Vector3d &node_a, Vector3d &node_b)
{
    Vector3d edge = node_b - node_a;
    Vector3d Jedge(-edge(1), edge(0), 0.0);
    
    double d = edge.norm();
    
    double r1 = (x - node_a).norm();
    double r2 = (x - node_b).norm();
    
    if (fabs(node_b(0) - node_a(0)) < Parameters::inversion_tolerance ||
        d                           < Parameters::inversion_tolerance ||
        fabs(r1 + r2 - d)           < Parameters::inversion_tolerance)
        return 0.0;
    
    double m = (node_b(1) - node_a(1)) / (node_b(0) - node_a(0));
    
    double e1 = pow(x(0) - node_a(0), 2) + pow(x(2), 2);
    double e2 = pow(x(0) - node_b(0), 2) + pow(x(2), 2);
    
    double h1 = (x(0) - node_a(0)) * (x(1) - node_a(1));
    double h2 = (x(0) - node_b(0)) * (x(1) - node_b(1));
    
    double delta_theta = atan((m * e1 - h1) / (x(2) * r1)) - atan((m * e2 - h2) / (x(2) * r2));
    
    return (x - node_a).dot(-Jedge) / d * log((r1 + r2 + d) / (r1 + r2 - d)) + fabs(x(2)) * delta_theta;
}

/**
   Computes the potential influence induced by a source panel of unit strength.  
   
   @param[in]   x           Point at which the influence coefficient is evaluated.
   @param[in]   this_panel  Panel on which the source panel is located.
   
   @returns Influence coefficient.
*/
double
Mesh::source_influence(Eigen::Vector3d &x, int this_panel)
{
    // Transform such that panel normal becomes unit Z vector:
    Matrix3d rotation = x_to_y_rotation(panel_normal(this_panel), Vector3d::UnitZ());
    
    Vector3d x_normalized = rotation * x;
    
    vector<Vector3d, Eigen::aligned_allocator<Vector3d> > normalized_panel_nodes;
    for (int i = 0; i < (int) panel_nodes[this_panel].size(); i++)
        normalized_panel_nodes.push_back(rotation * nodes[panel_nodes[this_panel][i]]);
    
    Vector3d translation(0, 0, -normalized_panel_nodes[0](2));
    
    x_normalized += translation;
    for (int i = 0; i < (int) panel_nodes[this_panel].size(); i++)
        normalized_panel_nodes[i] += translation;
    
    // Compute influence coefficient according to Hess:
    double influence = 0.0;
    for (int i = 0; i < (int) normalized_panel_nodes.size(); i++) {
        int prev_idx;
        if (i == 0)
            prev_idx = normalized_panel_nodes.size() - 1;
        else
            prev_idx = i - 1;
            
        Vector3d &node_a = normalized_panel_nodes[prev_idx];
        Vector3d &node_b = normalized_panel_nodes[i];
        
        influence += source_edge_influence(x_normalized, node_a, node_b);
    }   
    
    return -influence / (4 * M_PI);
}

// Compute velocity induced by an edge of a source panel:
static Vector3d
source_edge_unit_velocity(Vector3d &x, Vector3d &node_a, Vector3d &node_b)
{
    Vector3d edge = node_b - node_a;
    
    double d = edge.norm();
    
    double r1 = (x - node_a).norm();
    double r2 = (x - node_b).norm();
    
    if (fabs(node_b(0) - node_a(0)) < Parameters::inversion_tolerance ||
        d                           < Parameters::inversion_tolerance ||
        fabs(r1 + r2 + d)           < Parameters::inversion_tolerance)
        return Vector3d(0, 0, 0);
    
    double m = (node_b(1) - node_a(1)) / (node_b(0) - node_a(0));
    
    double e1 = pow(x(0) - node_a(0), 2) + pow(x(2), 2);
    double e2 = pow(x(0) - node_b(0), 2) + pow(x(2), 2);
    
    double h1 = (x(0) - node_a(0)) * (x(1) - node_a(1));
    double h2 = (x(0) - node_b(0)) * (x(1) - node_b(1));
    
    double delta_theta = atan((m * e1 - h1) / (x(2) * r1)) - atan((m * e2 - h2) / (x(2) * r2));
    
    double u = (node_b(1) - node_a(1)) / d * log((r1 + r2 - d) / (r1 + r2 + d));
    double v = (node_a(0) - node_b(0)) / d * log((r1 + r2 - d) / (r1 + r2 + d));
    double w = delta_theta;
    
    return Vector3d(u, v, w);
}

/**
   Computes the velocity induced by a source panel of unit strength.  
   
   @param[in]   x           Point at which the velocity is evaluated.
   @param[in]   this_panel  The panel on which the vortex ring is located.
   
   @returns Velocity induced by the source panel.
*/
Vector3d
Mesh::source_unit_velocity(Eigen::Vector3d &x, int this_panel)
{   
    // Transform such that panel normal becomes unit Z vector:
    Matrix3d rotation = x_to_y_rotation(panel_normal(this_panel), Vector3d::UnitZ());
    
    Vector3d x_normalized = rotation * x;
    
    vector<Vector3d, Eigen::aligned_allocator<Vector3d> > normalized_panel_nodes;
    for (int i = 0; i < (int) panel_nodes[this_panel].size(); i++)
        normalized_panel_nodes.push_back(rotation * nodes[panel_nodes[this_panel][i]]);
    
    Vector3d translation(0, 0, -normalized_panel_nodes[0](2));
    
    x_normalized += translation;
    for (int i = 0; i < (int) panel_nodes[this_panel].size(); i++)
        normalized_panel_nodes[i] += translation;
    
    // Compute influence coefficient according to Hess:
    Vector3d velocity(0, 0, 0);
    for (int i = 0; i < (int) normalized_panel_nodes.size(); i++) {
        int prev_idx;
        if (i == 0)
            prev_idx = normalized_panel_nodes.size() - 1;
        else
            prev_idx = i - 1;
            
        Vector3d &node_a = normalized_panel_nodes[prev_idx];
        Vector3d &node_b = normalized_panel_nodes[i];
        
        velocity += source_edge_unit_velocity(x_normalized, node_a, node_b);
    }   
    
    // Transform back:
    velocity = rotation.transpose() * velocity;
    
    // Done:
    return velocity / (4 * M_PI);
}

/**
   Computes the velocity induced by a vortex ring of unit strength.
   
   @param[in]   x           Point at which the velocity is evaluated.
   @param[in]   this_panel  Panel on which the vortex ring is located.
   
   @returns Velocity induced by the vortex ring.
*/
Vector3d
Mesh::vortex_ring_unit_velocity(Eigen::Vector3d &x, int this_panel)
{    
    Vector3d velocity(0, 0, 0);
    
    for (int i = 0; i < (int) panel_nodes[this_panel].size(); i++) {
        int previous_idx;
        if (i == 0)
            previous_idx = panel_nodes[this_panel].size() - 1;
        else
            previous_idx = i - 1;
            
        Vector3d &node_a = nodes[panel_nodes[this_panel][previous_idx]];
        Vector3d &node_b = nodes[panel_nodes[this_panel][i]];
        
        Vector3d r_0 = node_b - node_a;
        Vector3d r_1 = node_a - x;
        Vector3d r_2 = node_b - x;
        
        double r_1_norm = r_1.norm();
        double r_2_norm = r_2.norm();
        
        Vector3d r_1xr_2 = r_1.cross(r_2);
        double r_1xr_2_sqnorm = r_1xr_2.squaredNorm();
        
        if (r_1_norm < Parameters::inversion_tolerance ||
            r_2_norm < Parameters::inversion_tolerance ||
            r_1xr_2_sqnorm < Parameters::inversion_tolerance)
            continue;

        velocity += r_1xr_2 / r_1xr_2_sqnorm * r_0.dot(r_1 / r_1_norm - r_2 / r_2_norm);
    }

    return velocity / (4 * M_PI);
}

/**
   Computes the velocity induced by a Ramasamy-Leishman vortex ring.
   
   @param[in]   x           Point at which the velocity is evaluated.
   @param[in]   this_panel  Panel on which the vortex ring is located.
   @param[in]   core_radii  Radii of the filaments forming the vortex ring.
   @param[in]   vorticity   Strength of vortex ring.
   
   @returns Velocity induced by the Ramasamy-Leishman vortex ring.
   
   @note See M. Ramasamy and J. G. Leishman, Reynolds Number Based Blade Tip Vortex Model, University of Maryland, 2005.
*/
Vector3d
Mesh::vortex_ring_ramasamy_leishman_velocity(Eigen::Vector3d &x, int this_panel, std::vector<double> core_radii, double vorticity)
{
    // Ramasamy-Leishman series data:
    typedef struct {
        double vortex_reynolds_number;
        double a_1;
        double b_1;
        double a_2;
        double b_2;
        double b_3;
    } ramasamy_leishman_data_row;
    
    static ramasamy_leishman_data_row ramasamy_leishman_data[12] = { {    1, 1.0000, 1.2560, 0.0000, 0.00000, 0.0000},
                                                                     {  100, 1.0000, 1.2515, 0.0000, 0.00000, 0.0000},
                                                                     { 1000, 1.0000, 1.2328, 0.0000, 0.00000, 0.0000},
                                                                     {10000, 0.8247, 1.2073, 0.1753, 0.02630, 0.0000},
                                                                     {2.5e4, 0.5933, 1.3480, 0.2678, 0.01870, 0.2070},
                                                                     {4.8e4, 0.4602, 1.3660, 0.3800, 0.01380, 0.1674},
                                                                     {7.5e4, 0.3574, 1.3995, 0.4840, 0.01300, 0.1636},
                                                                     {  1e5, 0.3021, 1.4219, 0.5448, 0.01220, 0.1624},
                                                                     {2.5e5, 0.1838, 1.4563, 0.6854, 0.00830, 0.1412},
                                                                     {  5e5, 0.1386, 1.4285, 0.7432, 0.00580, 0.1144},
                                                                     {7.5e5, 0.1011, 1.4462, 0.7995, 0.00480, 0.1078},
                                                                     {  1e6, 0.0792, 1.4716, 0.8352, 0.00420, 0.1077} };
     
    // Compute vortex Reynolds number:                          
    double vortex_reynolds_number = vorticity / Parameters::fluid_kinematic_viscosity;
    
    // Interpolate Ramasamy-Leishman series values piecewise-linearly:
    int less_than_idx;
    for (less_than_idx = 0; less_than_idx < 12; less_than_idx++) {
        ramasamy_leishman_data_row &row = ramasamy_leishman_data[less_than_idx];
        
        if (vortex_reynolds_number < row.vortex_reynolds_number)
            break;
    }
    
    double a[3];
    double b[3];
    if (less_than_idx == 0) {
        a[0] = ramasamy_leishman_data[0].a_1;
        a[1] = ramasamy_leishman_data[0].a_2;
        
        b[0] = ramasamy_leishman_data[0].b_1;
        b[1] = ramasamy_leishman_data[0].b_2;
        b[2] = ramasamy_leishman_data[0].b_3;
    } else if (less_than_idx == 12) {
        a[0] = ramasamy_leishman_data[11].a_1;
        a[1] = ramasamy_leishman_data[11].a_2;
        
        b[0] = ramasamy_leishman_data[11].b_1;
        b[1] = ramasamy_leishman_data[11].b_2;
        b[2] = ramasamy_leishman_data[11].b_3;
    } else {
        double delta_vortex_reynolds_number =
            ramasamy_leishman_data[less_than_idx].vortex_reynolds_number - ramasamy_leishman_data[less_than_idx - 1].vortex_reynolds_number;
        double x = vortex_reynolds_number - ramasamy_leishman_data[less_than_idx - 1].vortex_reynolds_number;
        double slope;
        
        slope = (ramasamy_leishman_data[less_than_idx].a_1 - ramasamy_leishman_data[less_than_idx - 1].a_1) / delta_vortex_reynolds_number;
        a[0] = ramasamy_leishman_data[less_than_idx - 1].a_1 + slope * x;
        
        slope = (ramasamy_leishman_data[less_than_idx].a_2 - ramasamy_leishman_data[less_than_idx - 1].a_2) / delta_vortex_reynolds_number;
        a[1] = ramasamy_leishman_data[less_than_idx - 1].a_2 + slope * x;
        
        slope = (ramasamy_leishman_data[less_than_idx].b_1 - ramasamy_leishman_data[less_than_idx - 1].b_1) / delta_vortex_reynolds_number;
        b[0] = ramasamy_leishman_data[less_than_idx - 1].b_1 + slope * x;
        
        slope = (ramasamy_leishman_data[less_than_idx].b_2 - ramasamy_leishman_data[less_than_idx - 1].b_2) / delta_vortex_reynolds_number;
        b[1] = ramasamy_leishman_data[less_than_idx - 1].b_2 + slope * x;
        
        slope = (ramasamy_leishman_data[less_than_idx].b_3 - ramasamy_leishman_data[less_than_idx - 1].b_3) / delta_vortex_reynolds_number;
        b[2] = ramasamy_leishman_data[less_than_idx - 1].b_3 + slope * x;
    }
    
    a[2] = 1 - a[0] - a[1];
    
    // Compute velocity:
    Vector3d velocity(0, 0, 0);
    
    for (int i = 0; i < (int) panel_nodes[this_panel].size(); i++) {
        int previous_idx;
        if (i == 0)
            previous_idx = panel_nodes[this_panel].size() - 1;
        else
            previous_idx = i - 1;
            
        Vector3d &node_a = nodes[panel_nodes[this_panel][previous_idx]];
        Vector3d &node_b = nodes[panel_nodes[this_panel][i]];
        
        Vector3d r_0 = node_b - node_a;
        Vector3d r_1 = node_a - x;
        Vector3d r_2 = node_b - x;
        
        double r_0_norm = r_0.norm();
        double r_1_norm = r_1.norm();
        double r_2_norm = r_2.norm();
        
        Vector3d r_1xr_2 = r_1.cross(r_2);
        double r_1xr_2_sqnorm = r_1xr_2.squaredNorm();
        double r_1xr_2_norm = sqrt(r_1xr_2_sqnorm);
        
        if (r_0_norm < Parameters::inversion_tolerance ||
            r_1_norm < Parameters::inversion_tolerance ||
            r_2_norm < Parameters::inversion_tolerance ||
            r_1xr_2_sqnorm < Parameters::inversion_tolerance)
            continue;
            
        double d = r_1xr_2_norm / r_0_norm;
            
        double sum = 0;
        for (int j = 0; j < 3; j++)
            sum += a[j] * exp(-b[j] * pow(d / core_radii[i], 2));

        velocity += (1 - sum) * r_1xr_2 / r_1xr_2_sqnorm * r_0.dot(r_1 / r_1_norm - r_2 / r_2_norm);
    }

    return vorticity * velocity / (4 * M_PI);    
}

/**
   Computes the potential influence induced by a doublet panel of unit strength.  If the influence
   has been computed before, the cached influence coefficient is returned.  If not, the coefficient
   is computed and cached.
   
   @param[in]   other       Mesh on which the influence coefficient is evaluated.
   @param[in]   other_panel Panel on which the influence coefficient is evaluated.
   @param[in]   this_panel  Panel on which the doublet panel is located.
   
   @returns Influence coefficient.
*/
double
Mesh::doublet_influence(Mesh &other, int other_panel, int this_panel)
{
    // Retrieve inner matrix:
    map<int, std::vector<std::vector<double> > >::iterator it = doublet_influence_cache.find(other.id);
    if (it == doublet_influence_cache.end()) {
        cout << "Mesh " << id << ": Generating doublet influence expression cache for mesh " << other.id << "." << endl;
        
        std::vector<std::vector<double> > inner;
        inner.reserve(other.n_panels());
        for (int i = 0; i < other.n_panels(); i++) {
            Vector3d x = other.panel_collocation_point(i, true);
            
            std::vector<double> inner_inner;
            inner_inner.reserve(n_panels());
            for (int j = 0; j < n_panels(); j++)
                inner_inner.push_back(doublet_influence(x, j));
            inner.push_back(inner_inner);
        }
        
        doublet_influence_cache[other.id] = inner;
        
        return inner[other_panel][this_panel];
    } else {
        std::vector<std::vector<double> > &inner = it->second;
 
        return inner[other_panel][this_panel];
    }
}

/**
   Computes the potential influence induced by a source panel of unit strength.  If the influence
   has been computed before, the cached influence coefficient is returned.  If not, the coefficient
   is computed and cached.
   
   @param[in]   other       Mesh on which the influence coefficient is evaluated.
   @param[in]   other_panel Panel on which the influence coefficient is evaluated.
   @param[in]   this_panel  Panel on which the source panel is located.
   
   @returns Influence coefficient.
*/
double
Mesh::source_influence(Mesh &other, int other_panel, int this_panel)
{
    // Retrieve inner matrix:
    map<int, std::vector<std::vector<double> > >::iterator it = source_influence_cache.find(other.id);
    if (it == source_influence_cache.end()) {
        cout << "Mesh " << id << ": Generating source influence expression cache for mesh " << other.id << "." << endl;
        
        std::vector<std::vector<double> > inner;
        inner.reserve(other.n_panels());
        for (int i = 0; i < other.n_panels(); i++) {
            Vector3d x = other.panel_collocation_point(i, true);
            
            std::vector<double> inner_inner;
            inner_inner.reserve(n_panels());
            for (int j = 0; j < n_panels(); j++)        
                inner_inner.push_back(source_influence(x, j));
            inner.push_back(inner_inner);
        }
        
        source_influence_cache[other.id] = inner;
        
        return inner[other_panel][this_panel];
    } else {
        std::vector<std::vector<double> > &inner = it->second;
 
        return inner[other_panel][this_panel];
    }
}

/**
   Computes the velocity induced by a source panel of unit strength.  If the velocity
   has been computed before, the cached velocity is returned.  If not, the quantity
   is computed and cached.
   
   @param[in]   other       Mesh on which the velocity is evaluated.
   @param[in]   other_panel Panel on which the velocity is evaluated.
   @param[in]   this_panel  Panel on which the vortex ring is located.
   
   @returns Velocity induced by the source panel.
*/
Vector3d
Mesh::source_unit_velocity(Mesh &other, int other_panel, int this_panel)
{
    // Retrieve inner matrix:
    map<int, std::vector<std::vector<Vector3d, Eigen::aligned_allocator<Vector3d> > > >::iterator it = source_unit_velocity_cache.find(other.id);
    if (it == source_unit_velocity_cache.end()) {
        cout << "Mesh " << id << ": Generating source unit velocity cache for mesh " << other.id << "." << endl;
        
        std::vector<std::vector<Vector3d, Eigen::aligned_allocator<Vector3d> > > inner;
        inner.reserve(other.n_panels());
        for (int i = 0; i < other.n_panels(); i++) {
            Vector3d x = other.panel_collocation_point(i, true);
            
            std::vector<Vector3d, Eigen::aligned_allocator<Vector3d> > inner_inner;
            inner_inner.reserve(n_panels());
            for (int j = 0; j < n_panels(); j++)        
                inner_inner.push_back(source_unit_velocity(x, j));
            inner.push_back(inner_inner);
        }
        
        source_unit_velocity_cache[other.id] = inner;
        
        return inner[other_panel][this_panel];
    } else {
        std::vector<std::vector<Vector3d, Eigen::aligned_allocator<Vector3d> > > &inner = it->second;
 
        return inner[other_panel][this_panel];
    }
}

/**
   Computes the velocity induced by a vortex ring of unit strength.  If the velocity
   has been computed before, the cached velocity is returned.  If not, the quantity
   is computed and cached.
   
   @param[in]   other       Mesh on which the velocity is evaluated.
   @param[in]   other_panel Panel on which the velocity is evaluated.
   @param[in]   this_panel  Panel on which the vortex ring is located.
   
   @returns Velocity induced by the vortex ring.
*/
Vector3d
Mesh::vortex_ring_unit_velocity(Mesh &other, int other_panel, int this_panel)
{
    // Retrieve inner matrix:
    map<int, std::vector<std::vector<Vector3d, Eigen::aligned_allocator<Vector3d> > > >::iterator it = vortex_ring_unit_velocity_cache.find(other.id);
    if (it == vortex_ring_unit_velocity_cache.end()) {
        cout << "Mesh " << id << ": Generating vortex ring unit velocity cache for mesh " << other.id << "." << endl;
        
        std::vector<std::vector<Vector3d, Eigen::aligned_allocator<Vector3d> > > inner;
        inner.reserve(other.n_panels());
        for (int i = 0; i < other.n_panels(); i++) {
            Vector3d x = other.panel_collocation_point(i, true);
            
            std::vector<Vector3d, Eigen::aligned_allocator<Vector3d> > inner_inner;
            inner_inner.reserve(n_panels());
            for (int j = 0; j < n_panels(); j++)        
                inner_inner.push_back(vortex_ring_unit_velocity(x, j));
            inner.push_back(inner_inner);
        }
        
        vortex_ring_unit_velocity_cache[other.id] = inner;
        
        return inner[other_panel][this_panel];
    } else {
        std::vector<std::vector<Vector3d, Eigen::aligned_allocator<Vector3d> > > &inner = it->second;
 
        return inner[other_panel][this_panel];
    }
}

/**
   Computes the velocity induced by a Ramasamy-Leishman vortex ring.  If the velocity
   has been computed before, the cached velocity is returned.  If not, the quantity
   is computed and cached.
   
   @param[in]   other       Mesh on which the velocity is evaluated.
   @param[in]   other_panel Panel on which the velocity is evaluated.
   @param[in]   this_panel  Panel on which the vortex ring is located.
   @param[in]   core_radii  Radii of the filaments forming the vortex ring.
   @param[in]   vorticity   Strength of vortex ring.
   
   @returns Velocity induced by the Ramasamy-Leishman vortex ring.
   
   @note See M. Ramasamy and J. G. Leishman, Reynolds Number Based Blade Tip Vortex Model, University of Maryland, 2005.
*/
Vector3d
Mesh::vortex_ring_ramasamy_leishman_velocity(Mesh &other, int other_panel, int this_panel, std::vector<double> core_radii, double vorticity)
{
    // Retrieve inner matrix:
    map<int, std::vector<std::vector<Vector3d, Eigen::aligned_allocator<Vector3d> > > >::iterator it = vortex_ring_ramasamy_leishman_velocity_cache.find(other.id);
    if (it == vortex_ring_ramasamy_leishman_velocity_cache.end()) {
        cout << "Mesh " << id << ": Generating vortex ring unit velocity cache for mesh " << other.id << "." << endl;
        
        std::vector<std::vector<Vector3d, Eigen::aligned_allocator<Vector3d> > > inner;
        inner.reserve(other.n_panels());
        for (int i = 0; i < other.n_panels(); i++) {
            Vector3d x = other.panel_collocation_point(i, true);
            
            std::vector<Vector3d, Eigen::aligned_allocator<Vector3d> > inner_inner;
            inner_inner.reserve(n_panels());
            for (int j = 0; j < n_panels(); j++)        
                inner_inner.push_back(vortex_ring_ramasamy_leishman_velocity(x, j, core_radii, vorticity));
            inner.push_back(inner_inner);
        }
        
        vortex_ring_ramasamy_leishman_velocity_cache[other.id] = inner;
        
        return inner[other_panel][this_panel];
    } else {
        std::vector<std::vector<Vector3d, Eigen::aligned_allocator<Vector3d> > > &inner = it->second;
 
        return inner[other_panel][this_panel];
    }
}

/**
   Invalidates all caches of mesh data including collocation points, normals,
   surface areas, diameters, influence coefficients, and unit velocities.
*/
void
Mesh::invalidate_cache()
{
    panel_collocation_point_cache[0].clear();
    panel_collocation_point_cache[1].clear();
    panel_normal_cache.clear();
    panel_surface_area_cache.clear();
    panel_diameter_cache.clear();
    doublet_influence_cache.clear();
    source_influence_cache.clear();
    source_unit_velocity_cache.clear();
    vortex_ring_unit_velocity_cache.clear();
    vortex_ring_ramasamy_leishman_velocity_cache.clear();
}
