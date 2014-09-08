//
// Vortexje -- Gmsh surface writer.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <cmath>
#include <iostream>
#include <fstream>

#include <vortexje/surface-writers/gmsh-surface-writer.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

/**
   Returns the Gmsh MSH surface file extension (".msh").
   
   @returns The Gmsh MSH surface file extension (".msh").
*/
const char *
GmshSurfaceWriter::file_extension() const
{
    return ".msh";
}

/**
   Saves the given surface to a Gmsh MSH file, including data vectors associating numerical values to each panel.
  
   @param[in]   surface        Surface to write.
   @param[in]   filename       Destination filename.
   @param[in]   node_offset    Node numbering offset in output file.
   @param[in]   panel_offset   Panel numbering offset in output file.
   @param[in]   view_names     List of names of data vectors to be stored.
   @param[in]   view_data      List of data vectors to be stored.
   
   @returns true on success.
*/
bool
GmshSurfaceWriter::write(const shared_ptr<Surface> &surface, const string &filename, 
                         int node_offset, int panel_offset,
                         const std::vector<std::string> &view_names, const vector<MatrixXd, Eigen::aligned_allocator<MatrixXd> > &view_data)
{
    cout << "Surface " << surface->id << ": Saving to " << filename << "." << endl;
    
    // Save surface to gmsh file:
    ofstream f;
    f.open(filename.c_str());
    
    f << "$MeshFormat" << endl;
    f << "2.2 0 8" << endl; 
    f << "$EndMeshFormat" << endl;
    
    f << "$Nodes" << endl;
    f << surface->n_nodes() << endl;
    
    for (int i = 0; i < surface->n_nodes(); i++) {
        f << i + node_offset + 1;
        
        for (int j = 0; j < 3; j++) {
            f << ' ';
            f << surface->nodes[i](j);
        }
        
        f << endl;
    }
    
    f << "$EndNodes" << endl;
    f << "$Elements" << endl;
    f << surface->n_panels() << endl;
    
    for (int i = 0; i < surface->n_panels(); i++) {
        int element_type;
        switch (surface->panel_nodes[i].size()) {
        case 3:
            element_type = 2;
            break;
        case 4:
            element_type = 3;
            break;
        default:
            cerr << "Surface " << surface->id << ": Unknown polygon at panel " << i << "." << endl;
            continue;
        }
        
        f << i + panel_offset + 1;
        
        f << ' ';
        
        f << element_type;
        
        f << ' ';
        
        f << 0;
        
        for (int j = 0; j < (int) surface->panel_nodes[i].size(); j++) {
            f << ' ';
            f << surface->panel_nodes[i][j] + node_offset + 1;
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
        f << view_data[k].cols() << endl;
        f << surface->n_panels() << endl;
        
        for (int i = 0; i < surface->n_panels(); i++) {            
            f << i + panel_offset + 1;
            
            for (int j = 0; j < view_data[k].cols(); j++) {
                f << ' ';
            
                f << view_data[k](i, j);
            }
            
            f << endl;       
        }
        
        f << "$EndElementData" << endl;
    }
    
    f.close();
    
    // Done:
    return true;
}
