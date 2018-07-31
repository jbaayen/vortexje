//
// Vortexje -- VTK surface writer.
//
// Copyright (C) 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <cmath>
#include <iostream>
#include <fstream>

#include <vortexje/surface-writers/vtk-surface-writer.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

/**
   Returns the VTK surface file extension (".vtk").
   
   @returns The VTK surface file extension (".vtk").
*/
const char *
VTKSurfaceWriter::file_extension() const
{
    return ".vtk";
}

/**
   Saves the given surface to a VTK file, including data vectors associating numerical values to each panel.
  
   @param[in]   surface        Surface to write.
   @param[in]   filename       Destination filename.
   @param[in]   node_offset    Node numbering offset in output file.
   @param[in]   panel_offset   Panel numbering offset in output file.
   @param[in]   view_names     List of names of data vectors to be stored.
   @param[in]   view_data      List of data vectors to be stored.
   
   @returns true on success.
*/
bool
VTKSurfaceWriter::write(const std::shared_ptr<Surface> &surface, const string &filename, 
                        int node_offset, int panel_offset,
                        const std::vector<std::string> &view_names, const vector_aligned<MatrixXd> &view_data)
{
    cout << "Surface " << surface->id << ": Saving to " << filename << "." << endl;
    
    // Save surface to VTK file:
    ofstream f;
    f.open(filename.c_str());
    
    f << "# vtk DataFile Version 2.0" << endl;
    f << "FieldData" << endl;
    f << "ASCII" << endl;
    f << "DATASET UNSTRUCTURED_GRID" << endl;
    f << "POINTS " << surface->n_nodes() << " double" << endl;
    
    for (int i = 0; i < surface->n_nodes(); i++) {
        for (int j = 0; j < 3; j++) {
            if (j > 0)
                f << ' ';
            f << surface->nodes[i](j);
        }
        f << endl;
    }
    
    f << endl;
    
    int size = 0;
    for (int i = 0; i < surface->n_panels(); i++)
        size += surface->panel_nodes[i].size() + 1;
    
    f << "CELLS " << surface->n_panels() << " " << size << endl;
    
    for (int i = 0; i < surface->n_panels(); i++) {
        f << surface->panel_nodes[i].size();

        for (int j = 0; j < (int) surface->panel_nodes[i].size(); j++) {
            f << ' ';
            f << surface->panel_nodes[i][j];
        }
        
        f << endl;
    }
    
    f << endl;
    
    f << "CELL_TYPES " << surface->n_panels() << endl;
    
    for (int i = 0; i < surface->n_panels(); i++) {
        int cell_type;
        switch (surface->panel_nodes[i].size()) {
        case 3:
            cell_type = 5;
            break;
        case 4:
            cell_type = 9;
            break;
        default:
            cerr << "Surface " << surface->id << ": Unknown polygon at panel " << i << "." << endl;
            continue;
        }
        
        f << cell_type << endl;
    }
    
    f << endl;

    f << "CELL_DATA " << surface->n_panels() << endl;
    
    for (int k = 0; k < (int) view_names.size(); k++) {
        if (view_data[k].cols() == 1) {
            f << "SCALARS " << view_names[k] << " double 1" << endl;
            f << "LOOKUP_TABLE default" << endl;
        } else
            f << "VECTORS " << view_names[k] << " double" << endl;       
    
        for (int i = 0; i < surface->n_panels(); i++) {  
            for (int j = 0; j < view_data[k].cols(); j++) {
                if (j > 0)
                    f << ' ';
                    
                f << view_data[k](i, j);
            }
            
            f << endl;
        }
        
        f << endl;
    }
    
    f.close();
    
    // Done:
    return true;
}
