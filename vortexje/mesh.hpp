//
// Vortexje -- Mesh.
//
// Copyright (C) 2012 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __MESH_HPP__
#define __MESH_HPP__

#include <map>
#include <string>

#include <Eigen/Core>
#include <Eigen/StdVector>

#include <vortexje/parameters.hpp>

namespace Vortexje
{

/**
   Mesh representation using node-panel, panel-node, and panel-panel data structures,
   and implements geometrical and singularity panel influence operations.
   
   @brief Mesh representation.
*/
class Mesh
{
public:
    /**
       Automatically generated mesh identification number.
    */
    int id;
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    Mesh();
    Mesh(std::string file);
    
    void clear_node_panel_neighbors();
    
    virtual ~Mesh();
    
    bool load(std::string file);
    void save(std::string file,
              int node_offset = 0, int panel_offset = 0);
    void save(std::string file,
              std::vector<std::string> &view_names, std::vector<Eigen::VectorXd> &view_data,
              int node_offset = 0, int panel_offset = 0);
              
    int add_triangle(int node_a, int node_b, int node_c);
    int add_quadrangle(int node_a, int node_b, int node_c, int node_d);
    
    void compute_panel_neighbors();
    
    int n_nodes();
    int n_panels();
    
    /**
       Node number to point map.
    */
    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > nodes;
    
    /**
       Node number to node deformation velocity map.
    */
    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > node_deformation_velocities;
    
    /**
       Node number to neigboring panel numbers map.
    */
    std::vector<std::vector<int> *> node_panel_neighbors;
    
    /**
       Panel number to comprising vertex numbers map.
    */
    std::vector<std::vector<int> > panel_nodes;
    
    /**
       Panel number to neighboring panel numbers map.
    */
    std::vector<std::vector<int> > panel_neighbors;
    
    void rotate(const Eigen::Vector3d &axis, double angle);
    virtual void transform(const Eigen::Matrix3d &transformation);
    virtual void translate(const Eigen::Vector3d &translation);
    
    void rotate(const Eigen::Vector3d &axis, double angle, std::vector<Mesh*> &corotating_meshes);
    virtual void transform(const Eigen::Matrix3d &transformation, std::vector<Mesh*> &cotransforming_meshes);
    virtual void translate(const Eigen::Vector3d &translation, std::vector<Mesh*> &cotranslating_meshes);
    
    double distance_to_panel(const Eigen::Vector3d &x, int panel);
    virtual bool closest_panel(const Eigen::Vector3d &x, int &panel, double &distance);
    
    Eigen::Vector3d panel_collocation_point(int panel, bool below_surface);
    
    Eigen::Vector3d panel_normal(int panel);
    
    double panel_surface_area(int panel);
    
    double panel_diameter(int panel);
    
    Eigen::Vector3d panel_deformation_velocity(int panel) const;
    
    virtual Eigen::Vector3d close_to_body_point(int node);
    
    Eigen::Vector3d scalar_field_gradient(const Eigen::VectorXd &scalar_field, int panel);
    
    double doublet_influence(const Eigen::Vector3d &x, int this_panel);
    double source_influence(const Eigen::Vector3d &x, int this_panel);
    
    Eigen::Vector3d source_unit_velocity(const Eigen::Vector3d &x, int this_panel);
    Eigen::Vector3d vortex_ring_unit_velocity(const Eigen::Vector3d &x, int this_panel);
    Eigen::Vector3d vortex_ring_ramasamy_leishman_velocity(const Eigen::Vector3d &x, int this_panel, std::vector<double> core_radii, double vorticity);
    
    double doublet_influence(Mesh &other, int other_panel, int this_panel);
    double source_influence(Mesh &other, int other_panel, int this_panel);
    
    Eigen::Vector3d source_unit_velocity(Mesh &other, int other_panel, int this_panel);
    Eigen::Vector3d vortex_ring_unit_velocity(Mesh &other, int other_panel, int this_panel);
    Eigen::Vector3d vortex_ring_ramasamy_leishman_velocity(Mesh &other, int other_panel, int this_panel, std::vector<double> core_radii, double vorticity);
    
    void invalidate_cache();
    
    /**
       Panel collocation point caches.
    */
    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > panel_collocation_point_cache[2];
    
    /**
       Panel normal vector cache.
    */
    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > panel_normal_cache;
    
    /**
       Panel surface area cache.
    */
    std::vector<double> panel_surface_area_cache;
    
    /**
       Panel diameter cache.
    */
    std::vector<double> panel_diameter_cache;
    
    /**
       Doublet influence coefficient cache.
    */
    std::map<int, std::vector<std::vector<double> > > doublet_influence_cache;
    
    /**
       Source influence coefficient cache.
    */
    std::map<int, std::vector<std::vector<double> > > source_influence_cache;
    
    /**
       Unit source panel induced velocity cache.
    */
    std::map<int, std::vector<std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > > > source_unit_velocity_cache;
    
    /**
       Unit vortex ring induced velocity cache.
    */
    std::map<int, std::vector<std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > > > vortex_ring_unit_velocity_cache;
    
    /**
       Ramasamy-Leishman vortex ring induced velocity cache.
    */
    std::map<int, std::vector<std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > > > vortex_ring_ramasamy_leishman_velocity_cache;
};

};

#endif // __MESH_HPP__
