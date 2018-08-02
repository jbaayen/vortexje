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
#include <vortexje/numeric-stream.hpp>

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

const char *
VTKFieldWriter::float_str() const
{
    switch(float_type) {
    case nstream::FLOAT:
      return "FLOAT";
    default:
      return "DOUBLE";
    }
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
    int nx = round((x_max - x_min) / dx) + 1;
    int ny = round((y_max - y_min) / dy) + 1;
    int nz = round((z_max - z_min) / dz) + 1;

    // Compute velocity vector field:
    cout << "VTKFieldWriter: Computing velocity vector field." << endl;

    vector<Vector3d, Eigen::aligned_allocator<Vector3d> > velocities;
    velocities.resize(nx * ny * nz);
    
    int i;
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic, 1)
        for (i = 0; i < nx * ny * nz; i++) {
            double x, y, z;
            x = x_min + (i % (nx * ny)) % nx * dx;
            y = y_min + (i % (nx * ny)) / nx * dy;
            z = z_min + (i / (nx * ny)) * dz;
            velocities[i] = solver.velocity(Vector3d(x, y, z));
        }
    }

    // Write output in VTK format:
    cout << "VTKFieldWriter: Saving velocity vector field to " << filename << "." << endl;
    
    ofstream f;
    f.open(filename.c_str(), ios::binary);
    nstream::onstream nf(f, mode, float_type, nstream::BIGENDIAN);
    
    write_preamble(f, nf, x_min, y_min, z_min, dx, dy, dz, nx, ny, nz);
    
    // Velocity vector field;    
    f << "VECTORS Velocity " << float_str() << endl;
    
    vector<Vector3d, Eigen::aligned_allocator<Vector3d> >::const_iterator it;
    for (it = velocities.begin(); it != velocities.end(); it++) {
        Vector3d v = *it;
        nf << v(0) << " " << v(1) << " " << v(2) << endl;
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
    int nx = round((x_max - x_min) / dx) + 1;
    int ny = round((y_max - y_min) / dy) + 1;
    int nz = round((z_max - z_min) / dz) + 1;

    // Compute velocity vector field:
    cout << "VTKFieldWriter: Computing velocity potential field." << endl;

    vector<double> velocity_potentials;
    velocity_potentials.resize(nx * ny * nz);
    
    int i;
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic, 1)
        for (i = 0; i < nx * ny * nz; i++) {
            double x, y, z;
            x = x_min + (i % (nx * ny)) % nx * dx;
            y = y_min + (i % (nx * ny)) / nx * dy;
            z = z_min + (i / (nx * ny)) * dz;
            velocity_potentials[i] = solver.velocity_potential(Vector3d(x, y, z));
        }
    }

    // Write output in VTK format:
    cout << "VTKFieldWriter: Saving velocity potential field to " << filename << "." << endl;
    
    ofstream f;
    f.open(filename.c_str(), ios::binary);
    nstream::onstream nf(f, mode, float_type, nstream::BIGENDIAN);
    
    write_preamble(f, nf, x_min, y_min, z_min, dx, dy, dz, nx, ny, nz);
    
    // Velocity potential field:    
    f << "SCALARS VelocityPotential double 1" << endl;
    f << "LOOKUP_TABLE default" << endl;
    
    vector<double>::const_iterator it;
    for (it = velocity_potentials.begin(); it != velocity_potentials.end(); it++) {
        double p = *it;
        nf << p << endl;
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
VTKFieldWriter::write_preamble(ofstream &f, nstream::onstream &nf,
                               double x_min, double y_min, double z_min,
                               double dx, double dy, double dz,
                               int nx, int ny, int nz) const
{
    f << "# vtk DataFile Version 2.0" << endl;
    f << "FieldData" << endl;
    f << (mode==nstream::BINARY? "BINARY" : "ASCII") << endl;
    f << "DATASET RECTILINEAR_GRID" << endl;
    f << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    f << "X_COORDINATES " << nx << " " << float_str() << endl;
    for (int i = 0; i < nx; i++) {
        if (i > 0)
            nf << ' ';
        nf << x_min + i * dx;
    }
    f << endl;
    f << "Y_COORDINATES " << ny << " " << float_str() << endl;
    for (int i = 0; i < ny; i++) {
        if (i > 0)
            nf << ' ';
        nf << y_min + i * dy;
    }
    f << endl;
    f << "Z_COORDINATES " << nz << " " << float_str() << endl;
    for (int i = 0; i < nz; i++) {
        if (i > 0)
            nf << ' ';
        nf << z_min + i * dz;
    }
    f << endl;
    
    f << endl;
    
    f << "POINT_DATA " << nx * ny * nz << endl;
}
