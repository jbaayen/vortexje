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
       Returns the appropriate file extension for this FieldWriter.
   
       @returns The file extension.
    */
    virtual const char *file_extension() const = 0;
    
    /**
       Logs the velocity vector field.  The grid is the smallest box
       encompassing all nodes of all surfaces, expanded in the X, Y, and Z
       directions by the given margins.
      
       @param[in]   solver     Solver whose state to output.
       @param[in]   filename   Destination filename.
       @param[in]   dx         Grid step size in X-direction.
       @param[in]   dy         Grid step size in Y-direction.
       @param[in]   dz         Grid step size in Z-direction.
       @param[in]   x_margin   Grid expansion margin in X-direction.
       @param[in]   y_margin   Grid expansion margin in Y-direction.
       @param[in]   z_margin   Grid expansion margin in Z-direction.
       
       @returns true on success.
    */
    virtual bool write_velocity_field(const Solver &solver, const std::string &filename,
                                      double dx, double dy, double dz,
                                      double x_margin = 0.0, double y_margin = 0.0, double z_margin = 0.0) = 0;
             
    /**
       Logs the velocity potential scalar field.  The grid is the smallest box
       encompassing all nodes of all surfaces, expanded in the X, Y, and Z
       directions by the given margins.
      
       @param[in]   solver     Solver whose state to output.
       @param[in]   filename   Destination filename.
       @param[in]   dx         Grid step size in X-direction.
       @param[in]   dy         Grid step size in Y-direction.
       @param[in]   dz         Grid step size in Z-direction.
       @param[in]   x_margin   Grid expansion margin in X-direction.
       @param[in]   y_margin   Grid expansion margin in Y-direction.
       @param[in]   z_margin   Grid expansion margin in Z-direction.
       
       @returns true on success.
    */
    virtual bool write_velocity_potential_field(const Solver &solver, const std::string &filename,
                                                double dx, double dy, double dz,
                                                double x_margin = 0.0, double y_margin = 0.0, double z_margin = 0.0) = 0;
                                                
    void compute_field_envelope(const Solver &solver, double dx, double dy, double dz, double x_margin, double y_margin, double z_margin);
    
protected:
    // This class is not re-entrant, just like Solver is not.
    
    /**
       Minimum grid coordinate on X-axis.
    */
    double x_min;
    
    /**
       Maximum grid coordinate on X-axis.
    */
    double x_max;
    
    /**
       Minimum grid coordinate on Y-axis.
    */
    double y_min;
    
    /**
       Maximum grid coordinate on Y-axis.
    */
    double y_max;
    
    /**
       Minimum grid coordinate on Z-axis.
    */
    double z_min;
    
    /**
       Maximum grid coordinate on Z-axis.
    */
    double z_max;
    
    /**
       Grid step size in X-direction.
    */
    double dx;
    
    /**
       Grid step size in Y-direction.
    */
    double dy;
    
    /**
       Grid step size in Z-direction.
    */
    double dz;
    
    /**
       Number of grid points in X-direction.
    */
    int nx;
    
    /**
       Number of grid points in Y-direction.
    */
    int ny;
    
    /**
       Number of grid points in Z-direction.
    */
    int nz;
};

};

#endif // __FIELD_WRITER_HPP__
