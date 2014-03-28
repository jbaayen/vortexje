//
// Vortexje -- VTK field writer.
//
// Copyright (C) 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <cmath>
#include <iostream>
#include <fstream>

#include <vortexje/field-writers/vtk-field-writer.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

/**
   Returns the VTK file extension (".vtk").
   
   @returns The VTK file extension (".vtk").
*/
const char *
VTKFieldWriter::file_extension() const
{
    return ".vtk";
}

/**
   Logs the velocity vector field into a VTK file.  The grid is the smallest
   box encompassing all nodes of all surfaces, expanded in the X, Y, and Z
   directions by the given margins.
  
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
bool
VTKFieldWriter::write_velocity_field(const Solver &solver, const std::string &filename,
                                     double x_min, double x_max,
                                     double y_min, double y_max,
                                     double z_min, double z_max,
                                     double dx, double dy, double dz)
{
    cout << "Solver: Computing and saving velocity vector field to " << filename << "." << endl;
    
    int nx = round((x_max - x_min) / dx) + 1;
    int ny = round((y_max - y_min) / dy) + 1;
    int nz = round((z_max - z_min) / dz) + 1;

    // Write output in VTK format:
    ofstream f;
    f.open(filename.c_str());
    
    write_preamble(f, x_min, y_min, z_min, dx, dy, dz, nx, ny, nz);
    
    double x, y, z;
    
    // Velocity vector field;    
    f << "VECTORS Velocity double" << endl;
    
    z = z_min;
    for (int i = 0; i < nz; i++) {
        y = y_min;
        for (int j = 0; j < ny; j++) {
            x = x_min;
            for (int k = 0; k < nx; k++) {
                Vector3d v = solver.velocity(Vector3d(x, y, z));
                
                f << v(0) << " " << v(1) << " " << v(2) << endl;
            
                x += dx;
            }
            
            y += dy;
        }
        
        z += dz;
    }
    
    // Close file:
    f.close();
    
    // Done:
    return true;
}

/**
   Logs the velocity potential scalar field into a VTK file.  The grid is the 
   smallest box encompassing all nodes of all surfaces, expanded in the X, Y,
   and Z directions by the given margins.
  
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
bool
VTKFieldWriter::write_velocity_potential_field(const Solver &solver, const std::string &filename,
                                               double x_min, double x_max,
                                               double y_min, double y_max,
                                               double z_min, double z_max,
                                               double dx, double dy, double dz)
{
    cout << "Solver: Computing and saving velocity potential field to " << filename << "." << endl;

    int nx = round((x_max - x_min) / dx) + 1;
    int ny = round((y_max - y_min) / dy) + 1;
    int nz = round((z_max - z_min) / dz) + 1;
    
    // Write output in VTK format:
    ofstream f;
    f.open(filename.c_str());
    
    write_preamble(f, x_min, y_min, z_min, dx, dy, dz, nx, ny, nz);
    
    double x, y, z;
    
    // Velocity potential field:    
    f << "SCALARS VelocityPotential double 1" << endl;
    f << "LOOKUP_TABLE default" << endl;
    
    z = z_min;
    for (int i = 0; i < nz; i++) {
        y = y_min;
        for (int j = 0; j < ny; j++) {
            x = z_min;
            for (int k = 0; k < nx; k++) {
                double p = solver.velocity_potential(Vector3d(x, y, z));
                
                f << p << endl;
            
                x += dx;
            }
            
            y += dy;
        }
        
        z += dz;
    }
    
    // Close file:
    f.close();
    
    // Done:
    return true;
}

/**
   Write preamble for VTK output file.
*/
void
VTKFieldWriter::write_preamble(ofstream &f,
                               double x_min, double y_min, double z_min,
                               double dx, double dy, double dz,
                               int nx, int ny, int nz) const
{
    f << "# vtk DataFile Version 2.0" << endl;
    f << "FieldData" << endl;
    f << "ASCII" << endl;
    f << "DATASET RECTILINEAR_GRID" << endl;
    f << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    f << "X_COORDINATES " << nx << " double" << endl;
    for (int i = 0; i < nx; i++) {
        if (i > 0)
            f << ' ';
        f << x_min + i * dx;
    }
    f << endl;
    f << "Y_COORDINATES " << ny << " double" << endl;
    for (int i = 0; i < ny; i++) {
        if (i > 0)
            f << ' ';
        f << y_min + i * dy;
    }
    f << endl;
    f << "Z_COORDINATES " << nz << " double" << endl;
    for (int i = 0; i < nz; i++) {
        if (i > 0)
            f << ' ';
        f << z_min + i * dz;
    }
    f << endl;
    
    f << endl;
    
    f << "POINT_DATA " << nx * ny * nz << endl;
}
