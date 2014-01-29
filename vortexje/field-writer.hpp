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
    virtual const char *file_extension() const = 0;
    
    virtual bool write_velocity_field(const Solver &solver, const std::string &filename,
                                      double dx, double dy, double dz,
                                      double x_margin = 0.0, double y_margin = 0.0, double z_margin = 0.0) = 0;
                                     
    virtual bool write_velocity_potential_field(const Solver &solver, const std::string &filename,
                                                double dx, double dy, double dz,
                                                double x_margin = 0.0, double y_margin = 0.0, double z_margin = 0.0) = 0;
                                                
    void compute_field_envelope(const Solver &solver, double dx, double dy, double dz, double x_margin, double y_margin, double z_margin);
    
protected:
    // This class is not re-entrant, just like Solver is not.
    double x_min;
    double x_max;
    double y_min;
    double y_max;
    double z_min;
    double z_max;
    
    double dx;
    double dy;
    double dz;
    
    int nx;
    int ny;
    int nz;
};

};

#endif // __FIELD_WRITER_HPP__
