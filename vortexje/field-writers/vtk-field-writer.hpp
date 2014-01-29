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

namespace Vortexje
{

/**
   VTK field file writer.
   
   @brief VTK field writer.
*/
class VTKFieldWriter : public FieldWriter
{
public:
    const char *file_extension() const;
    
    bool write_velocity_field(const Solver &solver, const std::string &filename,
                              double dx, double dy, double dz,
                              double x_margin, double y_margin, double z_margin);
                              
    bool write_velocity_potential_field(const Solver &solver, const std::string &filename,
                                        double dx, double dy, double dz,
                                        double x_margin, double y_margin, double z_margin);
                                        
private:
    void write_preamble(std::ofstream &f) const;
};

};

#endif // __VTK_FIELD_WRITER_HPP__
