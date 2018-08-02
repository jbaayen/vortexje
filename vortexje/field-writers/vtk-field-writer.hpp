//
// Vortexje -- VTK field writer.
//
// Copyright (C) 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __VTK_FIELD_WRITER_HPP__
#define __VTK_FIELD_WRITER_HPP__

#include <string>
#include <fstream>

#include <vortexje/field-writer.hpp>
#include <vortexje/numeric-stream.hpp>

namespace Vortexje
{

/**
   VTK field file writer.
   
   @brief VTK field writer.
*/
class VTKFieldWriter : public FieldWriter
{
    nstream::serial_mode_e mode;
    nstream::float_type_e float_type;
public:
    VTKFieldWriter(nstream::serial_mode_e m = nstream::ASCII,
                   nstream::float_type_e ft = nstream::DOUBLE)
      : mode(m), float_type(ft) {}
    const char *file_extension() const;
    
    bool write_velocity_field(const Solver &solver,
                              const std::string &filename,
                              double x_min, double x_max,
                              double y_min, double y_max,
                              double z_min, double z_max,
                              double dx, double dy, double dz);
                              
    bool write_velocity_potential_field(const Solver &solver,
                                        const std::string &filename,
                                        double x_min, double x_max,
                                        double y_min, double y_max,
                                        double z_min, double z_max,
                                        double dx, double dy, double dz);
                                        
private:
    void write_preamble(std::ofstream &f, nstream::onstream &nf,
                        double x_min, double y_min, double z_min,
                        double dx, double dy, double dz,
                        int nx, int ny, int nz) const;

    const char *float_str() const;
};

};

#endif // __VTK_FIELD_WRITER_HPP__
