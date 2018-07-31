//
// Vortexje -- VTK surface writer.
//
// Copyright (C) 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __VTK_SURFACE_WRITER_HPP__
#define __VTK_SURFACE_WRITER_HPP__

#include <string>

#include <vortexje/surface-writer.hpp>

namespace Vortexje
{

/**
   VTK surface file writer.
   
   @brief VTK surface writer.
*/
class VTKSurfaceWriter : public SurfaceWriter
{
public:
    const char *file_extension() const;
       
    bool write(const std::shared_ptr<Surface> &surface, const std::string &filename,
               int node_offset, int panel_offset,
               const std::vector<std::string> &view_names, const vector_aligned<Eigen::MatrixXd> &view_data);
};

};

#endif // __VTK_SURFACE_WRITER_HPP__
