//
// Vortexje -- Wing builder.
//
// Copyright (C) 2012 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __WING_BUILDER_HPP__
#define __WING_BUILDER_HPP__

#include <Eigen/Core>
#include <Eigen/StdVector>

#include <vortexje/wing.hpp>

namespace Vortexje
{

/**
   Class for constructions of wings.
   
   @brief Wing construction tools.
*/
class WingBuilder
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    WingBuilder(Wing &wing);
    
    /**
       Associated Wing.
    */
    Wing &wing;

    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > generate_naca_airfoil(double max_camber, double max_camber_dist, double max_thickness, bool finite_te_thickness, double chord, int n_points, int &trailing_edge_point_id) const;
    
    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > generate_clarky_airfoil(double chord, int n_points, int &trailing_edge_point_id) const;
    
    std::vector<int> add_points(const std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > &points, int trailing_edge_point_id);
    
    void connect_nodes(const std::vector<int> &first_nodes, const std::vector<int> &second_nodes,
                       int trailing_edge_point_id, int &trailing_edge_top_panel_id, int &trailing_edge_bottom_panel_id,
                       bool cyclic);
    
    std::vector<int> fill_airfoil(const std::vector<int> &airfoil_nodes, int trailing_edge_point_id, int z_sign);
};

};

#endif // __WING_BUILDER_HPP__
