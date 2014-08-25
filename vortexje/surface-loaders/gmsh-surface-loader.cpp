//
// Vortexje -- Gmsh surface loader.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <cmath>
#include <iostream>
#include <fstream>

#include <vortexje/surface-loaders/gmsh-surface-loader.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

/**
   Returns the Gmsh MSH surface file extension (".msh").
   
   @returns The Gmsh MSH surface file extension (".msh").
*/
const char *
GmshSurfaceLoader::file_extension() const
{
    return ".msh";
}

/**
   Loads and parses the contents of a Gmsh MSH file into a surface->

   @param[in]   surface    Surface to load to.
   @param[in]   filename   Filename pointing to the Gmsh MSH file to load.
   
   @returns true on success.
*/
bool
GmshSurfaceLoader::load(shared_ptr<Surface> surface, const string &filename)
{
    cout << "Surface " << surface->id << ": Loading from " << filename << "." << endl;
    
    // Load surface from gmsh MSH file:
    ifstream f;
    f.open(filename.c_str());
    
    bool in_nodes = false, in_elements = false;
    int current_panel = 0;
    while (f.good()) {
        string line;
        getline(f, line);
        
        if (line[0] == '$') {
            if        (line == "$MeshFormat"   ) {
                getline(f, line);
                
                istringstream tokens(line);
                
                string version;
                int file_type, data_size;
                tokens >> version >> file_type >> data_size;
                
                if (version != "2.2" || file_type != 0 || data_size != 8) {
                    cerr << "Surface " << surface->id << ": Unknown data format in " << filename << "." << endl;
                    
                    f.close();
                    
                    return false;
                }
                
                getline(f, line);
                
            } else if (line == "$Nodes"        ) {
                getline(f, line);
                
                in_nodes = true;
              
            } else if (line == "$EndNodes"     ) {
                in_nodes = false;
                
            } else if (line == "$Elements"     ) {
                getline(f, line);
                
                in_elements = true;
                
            } else if (line == "$EndElements"  ) {
                in_elements = false;
                
            }
            
        } else if (in_nodes) {
            istringstream tokens(line);
            
            int node_number;
            double x, y, z;
            tokens >> node_number >> x >> y >> z;
            
            surface->nodes.push_back(Vector3d(x, y, z));            
            
            shared_ptr<vector<int> > neighbor_list = make_shared<vector<int> >();
            surface->node_panel_neighbors.push_back(neighbor_list);
            
        } else if (in_elements) {
            istringstream tokens(line);
            
            int single_number, single_type, number_of_tags, number_of_nodes;
            tokens >> single_number >> single_type >> number_of_tags;
            
            switch (single_type) {
            case 2:
                number_of_nodes = 3;
                break;
            case 3:
                number_of_nodes = 4;
                break;
            default:
                continue; // Only read triangles and quatrangles.
            }
            
            for (int i = 0; i < number_of_tags; i++) {
                int tag;
                tokens >> tag;
            }
            
            vector<int> single_panel_nodes;
            for (int i = 0; i < number_of_nodes; i++) {
                int node;
                tokens >> node;
                
                single_panel_nodes.push_back(node - 1);
                
                surface->node_panel_neighbors[node - 1]->push_back(current_panel);
            }

            surface->panel_nodes.push_back(single_panel_nodes);
            
            current_panel++;
        }  
    }
    
    f.close();
    
    // Compute surface topology:
    surface->compute_topology();
    
    // Compute panel geometry:
    surface->compute_geometry();
    
    // Done:
    return true;
}
