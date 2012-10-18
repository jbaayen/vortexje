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

class WingBuilder
{
public:
    // Constructor:
    WingBuilder(Wing &wing);

    // Wing construction:
    std::vector<Eigen::Vector3d> generate_naca_airfoil(double max_camber, double max_camber_dist, double max_thickness, double chord, int n_points, int &trailing_edge_point_id);
    
    std::vector<Eigen::Vector3d> generate_clarky_airfoil(double chord, int n_points, int &trailing_edge_point_id);
    
    std::vector<int> add_points(std::vector<Eigen::Vector3d> &points, int trailing_edge_point_id);
    
    typedef enum {
        TRIANGLES_A,
        TRIANGLES_B,
        TRIANGLES_X,
        QUADRANGLES
    } ConnectNodesMode;
    
    void connect_nodes(std::vector<int> &first_nodes, std::vector<int> &second_nodes,
                       int trailing_edge_point_id, bool cyclic, ConnectNodesMode mode);
    
    std::vector<int> fill_airfoil(std::vector<int> airfoil_nodes, int trailing_edge_point_id, int z_sign);
    
    // Properties:
    Wing &wing;
};

};

#endif // __WING_BUILDER_HPP__
