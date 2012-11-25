//
// Vortexje -- Wing builder.
//
// Copyright (C) 2012 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __WING_BUILDER_HPP__
#define __WING_BUILDER_HPP__

#include <vector>

#include <Eigen/Core>

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
    WingBuilder(Wing &wing);
    
    /**
       Associated Wing.
    */
    Wing &wing;

    std::vector<Eigen::Vector3d> generate_naca_airfoil(double max_camber, double max_camber_dist, double max_thickness, double chord, int n_points, int &trailing_edge_point_id);
    
    std::vector<Eigen::Vector3d> generate_clarky_airfoil(double chord, int n_points, int &trailing_edge_point_id);
    
    std::vector<int> add_points(std::vector<Eigen::Vector3d> &points, int trailing_edge_point_id);
    
    /**
       Node connection mode. 
    */
    typedef enum {
        TRIANGLES_A,  /**< Two triangles per quadrangle, diagonal direction A. */
        TRIANGLES_B,  /**< Two triangles per quadrangle, diagonal direction B. */
        TRIANGLES_X,  /**< Four triangles per quadrangle. */
        QUADRANGLES   /**< Quadrangles only. */
    } ConnectNodesMode;
    
    void connect_nodes(std::vector<int> &first_nodes, std::vector<int> &second_nodes,
                       int trailing_edge_point_id, int &trailing_edge_top_panel_id, int &trailing_edge_bottom_panel_id,
                       bool cyclic, ConnectNodesMode mode);
    
    std::vector<int> fill_airfoil(std::vector<int> airfoil_nodes, int trailing_edge_point_id, int z_sign);
};

};

#endif // __WING_BUILDER_HPP__
