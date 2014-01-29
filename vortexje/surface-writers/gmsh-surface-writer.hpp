//
// Vortexje -- Gmsh surface writer.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __GMSH_SURFACE_WRITER_HPP__
#define __GMSH_SURFACE_WRITER_HPP__

#include <string>

#include <vortexje/surface-writer.hpp>

namespace Vortexje
{

namespace SurfaceWriters
{

/**
   Gmsh MSH file writer.
   
   @brief Gmsh surface writer.
*/
class GmshSurfaceWriter : public SurfaceWriter
{
public:
    const char *file_extension() const;
                    
    bool write(const Surface &surface, const std::string &filename,
               int node_offset, int panel_offset,
               const std::vector<std::string> &view_names, const std::vector<Eigen::VectorXd> &view_data);
};

};

};

#endif // __GMSH_SURFACE_WRITER_HPP__
