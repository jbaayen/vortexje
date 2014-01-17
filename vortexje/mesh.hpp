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
    
    bool load(const std::string file);
    void save(const std::string file,
              int node_offset = 0, int panel_offset = 0) const;
    void save(const std::string file,
              const std::vector<std::string> &view_names, const std::vector<Eigen::VectorXd> &view_data,
              int node_offset = 0, int panel_offset = 0) const;
              
    int add_triangle(int node_a, int node_b, int node_c);
    int add_quadrangle(int node_a, int node_b, int node_c, int node_d);
    
    void compute_topology();
    void compute_geometry();
    
    int n_nodes() const;
    int n_panels() const;
    
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
    
    /**
       Panel number to collocation point map.
    */
    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > panel_collocation_points[2];
    
    /**
       Panel number to normal map.
    */
    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > panel_normals;
    
    /**
       Panel number to surface area map.
    */
    std::vector<double> panel_surface_areas;
    
    /**
       Panel number to diameter map.
    */
    std::vector<double> panel_diameters;
    
    void rotate(const Eigen::Vector3d &axis, double angle);
    virtual void transform(const Eigen::Matrix3d &transformation);
    virtual void translate(const Eigen::Vector3d &translation);
    
    void rotate(const Eigen::Vector3d &axis, double angle, std::vector<Mesh*> &corotating_meshes);
    virtual void transform(const Eigen::Matrix3d &transformation, std::vector<Mesh*> &cotransforming_meshes);
    virtual void translate(const Eigen::Vector3d &translation, std::vector<Mesh*> &cotranslating_meshes);
    
    double distance_to_panel(const Eigen::Vector3d &x, int panel) const;
    virtual bool closest_panel(const Eigen::Vector3d &x, int &panel, double &distance) const;
    
    Eigen::Vector3d panel_collocation_point(int panel, bool below_surface) const;
    
    Eigen::Vector3d panel_normal(int panel) const;
    
    double panel_surface_area(int panel) const;
    
    double panel_diameter(int panel) const;
    
    Eigen::Vector3d panel_deformation_velocity(int panel) const;
    
    virtual Eigen::Vector3d close_to_body_point(int node) const;
    
    Eigen::Vector3d scalar_field_gradient(const Eigen::VectorXd &scalar_field, int this_panel) const;
    
    double doublet_influence(const Eigen::Vector3d &x, int this_panel) const;
    double source_influence(const Eigen::Vector3d &x, int this_panel) const;
    
    Eigen::Vector3d source_unit_velocity(const Eigen::Vector3d &x, int this_panel) const;
    Eigen::Vector3d vortex_ring_unit_velocity(const Eigen::Vector3d &x, int this_panel) const;
    Eigen::Vector3d vortex_ring_ramasamy_leishman_velocity(const Eigen::Vector3d &x, int this_panel, const std::vector<double> core_radii, double vorticity) const;
    
    double doublet_influence(const Mesh &other, int other_panel, int this_panel) const;
    double source_influence(const Mesh &other, int other_panel, int this_panel) const;
    
    Eigen::Vector3d source_unit_velocity(const Mesh &other, int other_panel, int this_panel) const;
    Eigen::Vector3d vortex_ring_unit_velocity(const Mesh &other, int other_panel, int this_panel) const;
    Eigen::Vector3d vortex_ring_ramasamy_leishman_velocity(const Mesh &other, int other_panel, int this_panel, const std::vector<double> core_radii, double vorticity) const;
};

};

#endif // __MESH_HPP__
