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
    /**
       Destructor.
    */
    virtual ~SurfaceLoader() {};
    
    /**
       Returns the appropriate file extension for this FieldWriter.
   
       @returns The file extension.
    */
    virtual const char *file_extension() const = 0;
    
    /**
       Loads and parses the contents of a file into a Surface.

       @param[in]   surface    Surface to load to.
       @param[in]   filename   Filename pointing to the file to load.
       
       @returns true on success.
    */
    virtual bool load(std::shared_ptr<Surface> surface, const std::string &filename) = 0;
};

};

#endif // __SURFACE_LOADER_HPP__
