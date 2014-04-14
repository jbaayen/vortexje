//
// Vortexje -- Surface writer base class.
//
// Copyright (C) 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __FIELD_WRITER_HPP__
#define __FIELD_WRITER_HPP__

#include <string>

#include <vortexje/solver.hpp>

namespace Vortexje
{

/**
   Field writer base class.
   
   @brief Field writer base class.
*/
class FieldWriter
{
public:
    /**
       Destructor.
    */
    virtual ~FieldWriter() {};
    
    /**
       Returns the appropriate file extension for this FieldWriter.
   
       @returns The file extension.
    */
    virtual const char *file_extension() const = 0;
    
    /**
       Logs the velocity vector field for a specified grid. 
      
       @param[in]   solver     Solver whose state to output.
       @param[in]   filename   Destination filename.
       @param[in]   x_min      Minimum X coordinate of grid.
       @param[in]   x_max      Maximum X coordinate of grid.
       @param[in]   y_min      Minimum Y coordinate of grid.
       @param[in]   y_max      Maximum Y coordinate of grid.
       @param[in]   z_min      Minimum Z coordinate of grid.
       @param[in]   z_max      Maximum Z coordinate of grid.
       @param[in]   dx         Grid step size in X-direction.
       @param[in]   dy         Grid step size in Y-direction.
       @param[in]   dz         Grid step size in Z-direction.
       
       @returns true on success.
    */
    virtual bool write_velocity_field(const Solver &solver,
                                      const std::string &filename,
                                      double x_min, double x_max,
                                      double y_min, double y_max,
                                      double z_min, double z_max,
                                      double dx, double dy, double dz) = 0;
             
    /**
       Logs the velocity potential scalar field for a specified grid.
      
       @param[in]   solver     Solver whose state to output.
       @param[in]   filename   Destination filename.
       @param[in]   x_min      Minimum X coordinate of grid.
       @param[in]   x_max      Maximum X coordinate of grid.
       @param[in]   y_min      Minimum Y coordinate of grid.
       @param[in]   y_max      Maximum Y coordinate of grid.
       @param[in]   z_min      Minimum Z coordinate of grid.
       @param[in]   z_max      Maximum Z coordinate of grid.
       @param[in]   dx         Grid step size in X-direction.
       @param[in]   dy         Grid step size in Y-direction.
       @param[in]   dz         Grid step size in Z-direction.
       
       @returns true on success.
    */
    virtual bool write_velocity_potential_field(const Solver &solver,
                                                const std::string &filename,
                                                double x_min, double x_max,
                                                double y_min, double y_max,
                                                double z_min, double z_max,
                                                double dx, double dy, double dz) = 0;
};

};

#endif // __FIELD_WRITER_HPP__
