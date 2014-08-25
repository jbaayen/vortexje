//
// Vortexje -- PLY surface loader.
//
// Copyright (C) 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __PLY_SURFACE_LOADER_HPP__
#define __PLY_SURFACE_LOADER_HPP__

#include <string>

#include <vortexje/surface-loader.hpp>

namespace Vortexje
{

/**
   Polygon File Format/Stanford Triangle Format loader.
   
   @brief PLY surface loader.
*/
class PLYSurfaceLoader : public SurfaceLoader
{
public:
    const char *file_extension() const;
    
    bool load(std::shared_ptr<Surface> surface, const std::string &filename);
    
    // rply callbacks.
    void read_vertex_coordinate(int index, double value);
    void read_panel_node(int index, int length, int node);
    
private:
    // This class is not re-entrant.
    std::shared_ptr<Surface> surface;
    
    int current_panel;
    
    Eigen::Vector3d current_point;
    std::vector<int> current_panel_nodes;  
};

};

#endif // __PLY_SURFACE_LOADER_HPP__
