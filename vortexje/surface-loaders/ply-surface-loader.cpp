//
// Vortexje -- PLY surface loader.
//
// Copyright (C) 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <iostream>

#include <vortexje/surface-loaders/ply-surface-loader.hpp>

#include "rply/rply.h"

using namespace std;
using namespace Eigen;
using namespace Vortexje;

/**
   Returns the PLY file extension (".ply").
   
   @returns The PLY file extension (".ply").
*/
const char *
PLYSurfaceLoader::file_extension() const
{
    return ".ply";
}

// rply parser callbacks.
static int
vertex_cb(p_ply_argument argument)
{
    void *ptr;
    long index;
    ply_get_argument_user_data(argument, &ptr, &index);
    
    PLYSurfaceLoader *loader = static_cast<PLYSurfaceLoader *>(ptr);
    
    loader->read_vertex_coordinate(index, ply_get_argument_value(argument));
    
    return 1;
}

static int
face_cb(p_ply_argument argument)
{  
    void *ptr;
    ply_get_argument_user_data(argument, &ptr, NULL);
    
    PLYSurfaceLoader *loader = static_cast<PLYSurfaceLoader *>(ptr);
    
    long length, index;
    ply_get_argument_property(argument, NULL, &length, &index);
    
    if (index >= 0)
        loader->read_panel_node(index, length, (int) ply_get_argument_value(argument));
    
    return 1;
}

/**
   Loads and parses the contents of a PLY file into a surface->

   @param[in]   surface    Surface to load to.
   @param[in]   filename   Filename pointing to the PLY file to load.
   
   @returns true on success.
*/
bool
PLYSurfaceLoader::load(shared_ptr<Surface> surface, const string &filename)
{
    cout << "Surface " << surface->id << ": Loading from " << filename << "." << endl;
    
    this->surface = surface;
    
    current_panel = 0;
    
    // Load surface from PLY file.
    // We use "rply" from http://w3.impa.br/~diego/software/rply/, available under the MIT license.
    p_ply ply = ply_open(filename.c_str(), NULL, 0, NULL); 
    if (!ply)
        return false;
        
    if (!ply_read_header(ply))
        return false;
        
    ply_set_read_cb(ply, "vertex", "x", vertex_cb, this, 0);
    ply_set_read_cb(ply, "vertex", "y", vertex_cb, this, 1);
    ply_set_read_cb(ply, "vertex", "z", vertex_cb, this, 2);
    
    ply_set_read_cb(ply, "face", "vertex_indices", face_cb, this, 0);
    
    if (!ply_read(ply))
        return false;
        
    ply_close(ply);
    
    // Compute surface topology:
    surface->compute_topology();
    
    // Compute panel geometry:
    surface->compute_geometry();
    
    // Done:
    return true;
}

/**
   Internal interface method to rply library.
*/
void
PLYSurfaceLoader::read_vertex_coordinate(int index, double value)
{
    current_point(index) = value;
    
    if (index == 2) {
        // This is the last coordinate.  Process vertex.
        surface->nodes.push_back(current_point);      
            
        shared_ptr<vector<int> > neighbor_list = make_shared<vector<int> >();
        surface->node_panel_neighbors.push_back(neighbor_list);
    }
}

/**
   Internal interface method to rply library.
*/
void
PLYSurfaceLoader::read_panel_node(int index, int length, int node)
{
    current_panel_nodes.push_back(node);
    
    surface->node_panel_neighbors[node]->push_back(current_panel);
    
    if (index == length - 1) {
        // This is the last node.  Process panel.
        surface->panel_nodes.push_back(current_panel_nodes);
        
        current_panel_nodes.clear();
        
        current_panel++;
    }
}
