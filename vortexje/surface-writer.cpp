//
// Vortexje -- Surface writer base class.
//
// Copyright (C) 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <vortexje/surface-writer.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

/**
   Saves the given surface.
  
   @param[in]   surface    Surface to write.
   @param[in]   filename   Destination filename.
   
   @returns true on success.
*/
bool
SurfaceWriter::write(const std::shared_ptr<Surface> &surface, const std::string &filename)
{
    return write(surface, filename, 0, 0);
}

/**
   Saves the given surface.
  
   @param[in]   surface        Surface to write.
   @param[in]   filename       Destination filename.
   @param[in]   node_offset    Node numbering offset in output file.
   @param[in]   panel_offset   Panel numbering offset in output file.
   
   @returns true on success.
*/
bool
SurfaceWriter::write(const std::shared_ptr<Surface> &surface, const std::string &filename, 
                     int node_offset, int panel_offset)
{
    vector<string> empty_names;
    vector_aligned<MatrixXd> empty_data;
    
    return write(surface, filename, 0, 0, empty_names, empty_data);
}
