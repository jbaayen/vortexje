//
// Vortexje -- NACA 4-digit series airfoil generator.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <cmath>
#include <iostream>

#include <vortexje/shape-generators/airfoils/naca4-airfoil-generator.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

static const double pi = 3.141592653589793238462643383279502884;

// Cosine rule:
static double
cosine_rule(int n_points, int i)
{
    return 0.5 * (1 - cos(pi * i / (double) n_points));
}

// NACA airfoil generation:
static double
compute_y_c(double x, double max_camber, double max_camber_dist, double max_thickness, double chord)
{
    if (x < max_camber_dist * chord)
        return max_camber * x / pow(max_camber_dist, 2) * (2 * max_camber_dist - x / chord);
    else
        return max_camber * (chord - x) / pow(1 - max_camber_dist, 2) * (1 + x / chord - 2 * max_camber_dist);
}

static double
compute_theta(double x, double max_camber, double max_camber_dist, double max_thickness, double chord)
{
    if (x < max_camber_dist * chord)
        return atan(max_camber / pow(max_camber_dist, 2) * (2 * max_camber_dist - x / chord) 
                    - max_camber * x / pow(max_camber_dist, 2) / chord);
    else
        return atan(-max_camber / pow(1 - max_camber_dist, 2) * (1 + x / chord - 2 * max_camber_dist) 
                    + max_camber * (chord - x) / pow(1 - max_camber_dist, 2) / chord);
}

static double
compute_y_t(double x, double max_camber, double max_camber_dist, double max_thickness, bool finite_te_thickness, double chord)
{ 
    double k5;
    if (finite_te_thickness) {
        if (x == chord)
            return 0.0; // Ensure a symmetric airfoil.
            
        k5 = -0.1015;
    } else
        k5 = -0.1036;
        
    return max_thickness / 0.2 * chord * (0.2969 * sqrt(x / chord) - 0.1260 * (x / chord)
                                          - 0.3516 * pow(x / chord, 2) + 0.2843 * pow(x / chord, 3) + k5 * pow(x / chord, 4));
}

/**
   Generates points tracing a 4-digit series NACA airfoil.

   @param[in]   max_camber              Maximum camber as a percentage of the chord (the first digit, divided by 100).
   @param[in]   max_camber_dist         Distance of maximum camber from the leading edge (the second digit, dividid by 10).
   @param[in]   max_thickness           Maximum thickness of the airfoil (the last two digits, divided by 100).
   @param[in]   finite_te_thickness     True to use a trailing edge with finite thickness.
   @param[in]   chord                   Chord length.
   @param[in]   n_points                Number of points to return.
   @param[out]  trailing_edge_point_id Index of the trailing edge node in the returned list.
   
   @returns List of points.
   
   @note See S. Yon, J. Katz, and A. Plotkin, Effect of Airfoil (Trailing-Edge) Thickness on the Numerical Solution of Panel Methods Based on the Dirichlet Boundary Condition, AIAA Journal, Vol. 30, No. 3, March 1992, for the issues that may arise when using an infinitely thin trailing edge.
*/
vector_aligned<Vector3d>
NACA4AirfoilGenerator::generate(double max_camber, double max_camber_dist, double max_thickness, bool finite_te_thickness, double chord, int n_points, int &trailing_edge_point_id)
{
    if (n_points % 2 == 1) {
        cerr << "NACA4::generate(): n_nodes must be even." << endl;
        exit(1);
    }
    
    vector_aligned<Vector3d> airfoil_points;
    
    // Add upper nodes:
    for (int i = 0; i < n_points / 2; i++) {
        double x     = chord * cosine_rule(n_points / 2, i);
        
        double y_c   = compute_y_c  (x, max_camber, max_camber_dist, max_thickness, chord);
        double theta = compute_theta(x, max_camber, max_camber_dist, max_thickness, chord);
        double y_t   = compute_y_t  (x, max_camber, max_camber_dist, max_thickness, finite_te_thickness, chord);
        
        Vector3d upper_point(x - y_t * sin(theta), y_c + y_t * cos(theta), 0.0);
        airfoil_points.push_back(upper_point);
    }
    
    // Add lower nodes:
    for (int i = 0; i < n_points / 2; i++) {
        double x     = chord * (1 - cosine_rule(n_points / 2, i));
        
        double y_c   = compute_y_c  (x, max_camber, max_camber_dist, max_thickness, chord);
        double theta = compute_theta(x, max_camber, max_camber_dist, max_thickness, chord);
        double y_t   = compute_y_t  (x, max_camber, max_camber_dist, max_thickness, finite_te_thickness, chord);
        
        Vector3d lower_point(x + y_t * sin(theta), y_c - y_t * cos(theta), 0.0);
        airfoil_points.push_back(lower_point);
    }
    
    // Done.
    trailing_edge_point_id = n_points / 2;
    
    return airfoil_points;
}
