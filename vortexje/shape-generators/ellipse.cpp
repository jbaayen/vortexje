//
// Vortexje -- Ellipse generator.
//
// Copyright (C) 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <cmath>

#include <vortexje/shape-generators/ellipse.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje::ShapeGenerators;

/**
   Generates points tracing an ellipse.

   @param[in]   a          Semi-major axis (x-coordinate).
   @param[in]   b          Semi-minor axis (y-coordinate).
   @param[in]   n_points   Number of points to return.
   
   @returns List of points.
*/
vector<Vector3d, Eigen::aligned_allocator<Vector3d> >
Ellipse::generate(double a, double b, int n_points)
{
    vector<Vector3d, Eigen::aligned_allocator<Vector3d> > points;
    
    // Go in the clockwise (negative) direction, for consistency with the airfoil generators.
    double dt = -2 * M_PI / (double) n_points;
    for (int i = 0; i < n_points; i++) {
        double t = i * dt;
        
        Vector3d point(a * cos(t), b * sin(t), 0.0);

        points.push_back(point);
    }
    
    return points;
}
