//
// Vortexje -- Surface loader base class.
//
// Copyright (C) 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __SURFACE_LOADER_HPP__
#define __SURFACE_LOADER_HPP__

#include <string>

#include <vortexje/surface.hpp>

namespace Vortexje
{

/**
   Surface loader base class.
   
   @brief Surface loader base class.
*/
class SurfaceLoader
{
public:
    virtual const char *file_extension() const = 0;
    
    virtual bool load(Surface &surface, const std::string &filename) = 0;
};

};

#endif // __SURFACE_LOADER_HPP__
