//
// Vortexje -- Mesh.
//
// Copyright (C) 2012 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __MESH_HPP__
#define __MESH_HPP__

#include <vector>
#include <map>
#include <string>

#include <Eigen/Core>

#include <vortexje/parameters.hpp>

namespace Vortexje
{

class Mesh
{
public:
    // Constructor:
    int id;
    
    Mesh();
    Mesh(std::string file);
    
    // Destructor:
    void clear_node_panel_neighbors();
    
    ~Mesh();
    
    // Loading and saving of gmsh MSH files:
    bool load(std::string file);
    void save(std::string file,
              int node_offset = 0, int panel_offset = 0);
    void save(std::string file,
              std::vector<std::string> &view_names, std::vector<Eigen::VectorXd> &view_data,
              int node_offset = 0, int panel_offset = 0);
              
    // Mesh construction:
    int add_triangle(int node_a, int node_b, int node_c);
    int add_quadrangle(int node_a, int node_b, int node_c, int node_d);
    
    void compute_panel_neighbors();
    
    // Geometry:
    int n_nodes();
    int n_panels();
    
    std::vector<Eigen::Vector3d> nodes;
    std::vector<Eigen::Vector3d> node_deformation_velocities;
    std::vector<std::vector<int> *> node_panel_neighbors;
    std::vector<std::vector<int> > panel_nodes;
    std::vector<std::vector<int> > panel_neighbors;
    
    void rotate(Eigen::Vector3d axis, double angle);
    virtual void transform(Eigen::Matrix3d transformation);
    virtual void translate(Eigen::Vector3d translation);
    
    void rotate(Eigen::Vector3d axis, double angle, std::vector<Mesh*> &corotating_meshes);
    virtual void transform(Eigen::Matrix3d transformation, std::vector<Mesh*> &cotransforming_meshes);
    virtual void translate(Eigen::Vector3d translation, std::vector<Mesh*> &cotranslating_meshes);
    
    double distance_to_panel(Eigen::Vector3d x, int panel);
    virtual bool closest_panel(Eigen::Vector3d x, int &panel, double &distance);
    
    Eigen::Vector3d panel_collocation_point(int panel, bool below_surface);
    
    Eigen::Vector3d panel_normal(int panel);
    
    double panel_surface_area(int panel);
    
    double panel_diameter(int panel);
    
    Eigen::Vector3d panel_deformation_velocity(int panel);
    
    virtual Eigen::Vector3d close_to_body_point(int node);
    
    Eigen::Vector3d scalar_field_gradient(Eigen::VectorXd &scalar_field, int panel);
    
    // Panel influence expressions:
    double doublet_influence(Eigen::Vector3d x, int this_panel);
    double source_influence(Eigen::Vector3d x, int this_panel);
    
    Eigen::Vector3d source_unit_velocity(Eigen::Vector3d x, int this_panel);
    Eigen::Vector3d vortex_ring_unit_velocity(Eigen::Vector3d x, int this_panel);
    Eigen::Vector3d vortex_ring_ramasamy_leishman_velocity(Eigen::Vector3d x, int this_panel, std::vector<double> core_radii, double vorticity);
    
    // Panel influence expressions (cached):
    double doublet_influence(Mesh &other, int other_panel, int this_panel);
    double source_influence(Mesh &other, int other_panel, int this_panel);
    
    Eigen::Vector3d source_unit_velocity(Mesh &other, int other_panel, int this_panel);
    Eigen::Vector3d vortex_ring_unit_velocity(Mesh &other, int other_panel, int this_panel);
    Eigen::Vector3d vortex_ring_ramasamy_leishman_velocity(Mesh &other, int other_panel, int this_panel, std::vector<double> core_radii, double vorticity);
    
    // Caches:
    void invalidate_cache();
    
    std::vector<Eigen::Vector3d> panel_collocation_point_cache[2];
    std::vector<Eigen::Vector3d> panel_normal_cache;
    std::vector<double> panel_surface_area_cache;
    std::vector<double> panel_diameter_cache;
    std::map<int, std::vector<std::vector<double> > > doublet_influence_cache;
    std::map<int, std::vector<std::vector<double> > > source_influence_cache;
    std::map<int, std::vector<std::vector<Eigen::Vector3d> > > source_unit_velocity_cache;
    std::map<int, std::vector<std::vector<Eigen::Vector3d> > > vortex_ring_unit_velocity_cache;
    std::map<int, std::vector<std::vector<Eigen::Vector3d> > > vortex_ring_ramasamy_leishman_velocity_cache;
};

};

#endif // __MESH_HPP__
