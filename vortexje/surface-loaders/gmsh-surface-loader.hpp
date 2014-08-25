//
// Vortexje -- Gmsh surface loader.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __GMSH_SURFACE_LOADER_HPP__
#define __GMSH_SURFACE_LOADER_HPP__

#include <string>

#include <vortexje/surface-loader.hpp>

namespace Vortexje
{

/**
   Gmsh MSH file loader.
   
   @brief Gmsh surface loader.
*/
class GmshSurfaceLoader : public SurfaceLoader
{
public:
    const char *file_extension() const;
    
    bool load(std::shared_ptr<Surface> surface, const std::string &filename);
};

};

#endif // __GMSH_SURFACE_LOADER_HPP__
