//
// Vortexje -- Surface.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __SURFACE_HPP__
#define __SURFACE_HPP__

#include <memory>
#include <utility>
#include <map>
#include <string>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>

#include <vortexje/parameters.hpp>

namespace Vortexje
{

/**
   Surface representation using node-panel, panel-node, and panel-panel data structures.
   Implements geometrical and singularity panel influence operations.
   
   @brief Surface representation.
*/
class Surface
{
public:            
    /**
       Automatically generated surface identification number.
    */
    int id;
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    Surface();
    
    virtual ~Surface();
              
    int add_triangle(int node_a, int node_b, int node_c);
    int add_quadrangle(int node_a, int node_b, int node_c, int node_d);
    
    void compute_topology();
    void compute_geometry();
    
    void cut_panels(int panel_a, int panel_b);
    
    int n_nodes() const;
    int n_panels() const;

    /**
       Node number to point map.
    */
    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > nodes;
    
    /**
       Node number to neigboring panel numbers map.
    */
    std::vector<std::shared_ptr<std::vector<int> > > node_panel_neighbors;
    
    /**
       Panel number to comprising vertex numbers map.
    */
    std::vector<std::vector<int> > panel_nodes;
    
    /**
       Panel number to (edge number to (neighboring panel number, edge number)) map.
    */
    std::vector<std::map<int, std::pair<int, int> > > panel_neighbors;
    
    /**
       Panel number to comprising vertex points (in the panel coordinate system) map.
    */
    std::vector<std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > > panel_transformed_points;
    
    void rotate(const Eigen::Vector3d &axis, double angle);
    virtual void transform(const Eigen::Matrix3d &transformation);
    virtual void transform(const Eigen::Transform<double, 3, Eigen::Affine> &transformation);
    virtual void translate(const Eigen::Vector3d &translation);
    
    const Eigen::Vector3d &panel_collocation_point(int panel, bool below_surface) const;
    
    const Eigen::Vector3d &panel_normal(int panel) const;
    
    const Eigen::Transform<double, 3, Eigen::Affine> &panel_coordinate_transformation(int panel) const;
    
    double panel_surface_area(int panel) const;
    
    double panel_diameter(int panel) const;
    
    virtual void source_and_doublet_influence(const Eigen::Vector3d &x, int this_panel, double &source_influence, double &doublet_influence) const;
    
    double source_influence(const Eigen::Vector3d &x, int this_panel) const;
    double doublet_influence(const Eigen::Vector3d &x, int this_panel) const;
    
    virtual Eigen::Vector3d source_unit_velocity(const Eigen::Vector3d &x, int this_panel) const;
    virtual Eigen::Vector3d vortex_ring_unit_velocity(const Eigen::Vector3d &x, int this_panel) const;
    
    double doublet_influence(const std::shared_ptr<Surface> &other, int other_panel, int this_panel) const;
    double source_influence(const std::shared_ptr<Surface> &other, int other_panel, int this_panel) const;
    
    void source_and_doublet_influence(const std::shared_ptr<Surface> &other, int other_panel, int this_panel, double &source_influence, double &doublet_influence) const;
    
    Eigen::Vector3d source_unit_velocity(const std::shared_ptr<Surface> &other, int other_panel, int this_panel) const;
    Eigen::Vector3d vortex_ring_unit_velocity(const std::shared_ptr<Surface> &other, int other_panel, int this_panel) const;
    
protected:
    /**
       Panel number to collocation point map.
    */
    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > panel_collocation_points[2];
    
    /**
       Panel number to normal map.
    */
    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > panel_normals;
    
    /**
       Panel number to panel coordinate transformation map.
    */
    std::vector<Eigen::Transform<double, 3, Eigen::Affine>, Eigen::aligned_allocator<Eigen::Transform<double, 3, Eigen::Affine> > > panel_coordinate_transformations;
    
    /**
       Panel number to surface area map.
    */
    std::vector<double> panel_surface_areas;
    
    /**
       Panel number to diameter map.
    */
    std::vector<double> panel_diameters;
};

};

#endif // __SURFACE_HPP__
