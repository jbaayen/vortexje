//
// Vortexje -- Field writer base class.
//
// Copyright (C) 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <vortexje/field-writer.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

/**
   Computes the field envelope for a given Solver.
  
   @param[in]   solver     Solver whose state to output.
   @param[in]   x_margin   Grid expansion margin in X-direction.
   @param[in]   y_margin   Grid expansion margin in Y-direction.
   @param[in]   z_margin   Grid expansion margin in Z-direction.
*/
void
FieldWriter::compute_field_envelope(const Solver &solver, double dx, double dy, double dz, double x_margin, double y_margin, double z_margin)
{
    this->dx = dx;
    this->dy = dy;
    this->dz = dz;
    
    // Find extremal points:
    x_min = numeric_limits<double>::max();
    x_max = numeric_limits<double>::min();
    y_min = numeric_limits<double>::max();
    y_max = numeric_limits<double>::min();
    z_min = numeric_limits<double>::max();
    z_max = numeric_limits<double>::min();
    
    vector<Surface *> surfaces;
    for (int i = 0; i < (int) solver.bodies.size(); i++) {
        Body *body = solver.bodies[i];
  
        surfaces.insert(surfaces.end(), body->lifting_surfaces.begin(), body->lifting_surfaces.end());
        surfaces.insert(surfaces.end(), body->wakes.begin(), body->wakes.end());
        surfaces.insert(surfaces.end(), body->non_lifting_surfaces.begin(), body->non_lifting_surfaces.end());
    }
        
    for (int i = 0; i < (int) surfaces.size(); i++) {
        Surface *surface = surfaces[i];
        
        for (int j = 0; j < surface->n_nodes(); j++) {
            Vector3d &node = surface->nodes[j];
            
            if (node(0) < x_min)
                x_min = node(0);
            if (node(0) > x_max)
                x_max = node(0);
            if (node(1) < y_min)
                y_min = node(1);
            if (node(1) > y_max)
                y_max = node(1);
            if (node(2) < z_min)
                z_min = node(2);
            if (node(2) > z_max)
                z_max = node(2); 
        }
    }
    
    // Apply margin:
    x_min -= x_margin;
    x_max += x_margin;
    y_min -= y_margin;
    y_max += y_margin;
    z_min -= z_margin;
    z_max += z_margin;
    
    // Count points:
    nx = round((x_max - x_min) / dx) + 1;
    ny = round((y_max - y_min) / dy) + 1;
    nz = round((z_max - z_min) / dz) + 1;
}
