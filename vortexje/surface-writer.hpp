//
// Vortexje -- Surface writer base class.
//
// Copyright (C) 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __SURFACE_WRITER_HPP__
#define __SURFACE_WRITER_HPP__

#include <memory>
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
    /**
       Destructor.
    */
    virtual ~SurfaceWriter() {};
    
    /**
       Returns the appropriate file extension for this SurfaceWriter.
   
       @returns The file extension.
    */
    virtual const char *file_extension() const = 0;
    
    bool write(const std::shared_ptr<Surface> &surface, const std::string &filename);
    
    bool write(const std::shared_ptr<Surface> &surface, const std::string &filename,
               int node_offset, int panel_offset);
      
    /**
       Saves the given surface to a file, including data vectors associating numerical values to each panel.
      
       @param[in]   surface        Surface to write.
       @param[in]   filename       Destination filename.
       @param[in]   node_offset    Node numbering offset in output file.
       @param[in]   panel_offset   Panel numbering offset in output file.
       @param[in]   view_names     List of names of data vectors to be stored.
       @param[in]   view_data      List of data vectors to be stored.
       
       @returns true on success.
    */               
    virtual bool write(const std::shared_ptr<Surface> &surface, const std::string &filename,
                       int node_offset, int panel_offset,
                       const std::vector<std::string> &view_names, const vector_aligned<Eigen::MatrixXd> &view_data) = 0;
};

};

#endif // __SURFACE_WRITER_HPP__
