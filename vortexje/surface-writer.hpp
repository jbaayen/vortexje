//
// Vortexje -- Surface writer base class.
//
// Copyright (C) 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __SURFACE_WRITER_HPP__
#define __SURFACE_WRITER_HPP__

#include <string>

#include <vortexje/surface.hpp>

namespace Vortexje
{

/**
   Surface writer base class.
   
   @brief Surface writer base class.
*/
class SurfaceWriter
{
public:
    virtual const char *file_extension() const = 0;
    
    bool write(const Surface &surface, const std::string &filename);
    
    bool write(const Surface &surface, const std::string &filename,
               int node_offset, int panel_offset);
                    
    virtual bool write(const Surface &surface, const std::string &filename,
                       int node_offset, int panel_offset,
                       const std::vector<std::string> &view_names, const std::vector<Eigen::VectorXd> &view_data) = 0;
};

};

#endif // __SURFACE_WRITER_HPP__
