//
// Vortexje -- Test NACA0012 lift- and drag coefficients against a reference.
//
// Copyright (C) 2013 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <iostream>

#include <vortexje/solver.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

#define TEST_TOLERANCE 1e-4

Vector2d
run_test(double alpha)
{
    // Set up parameters for simplest possible simulation:
    Parameters::unsteady_bernoulli = false;
    Parameters::convect_wake = false;
    
    // Load meshes:
    Mesh nolift_mesh;
    
    // Create wing:
    Vector3d location(0, 0, 0);
    Vector3d chord_direction(1, 0, 0);
    Vector3d top_direction(0, 1, 0);
    Vector3d span_direction(0, 0, 1);
    
    Mesh wing_mesh(string("naca0012.msh"));
    
    Wing wing(wing_mesh, location, chord_direction, top_direction, span_direction);
    
    // Rotate by angle of attack:
    wing.rotate(span_direction, -alpha);
    
    // Create mesh collection:
    Collection collection(string("test-wing"),
                          nolift_mesh);
    collection.add_wing(&wing);
    
    // Set up solver:
    Solver solver("test-naca0012-log");
    solver.add_collection(collection);
    
    Vector3d freestream_velocity(30, 0, 0);
    solver.set_freestream_velocity(freestream_velocity);
    
    double fluid_density = 1.2;
    solver.set_fluid_density(fluid_density);
    
    // Compute:
    solver.initialize_wakes(0);
    solver.update_coefficients(0);
    
    // Output lift and drag coefficients:
    Vector3d F_a = solver.aerodynamic_force(collection);
    
    double q = 0.5 * fluid_density * freestream_velocity.squaredNorm();
    
    double S = 4.5 * 0.75;
    
    double C_L = F_a(1) / (q * S);
    double C_D = F_a(0) / (q * S);
    
    return Vector2d(C_L, C_D);
}

int
main (int argc, char **argv)
{
    // Load reference values.
    std::vector<Vector3d, Eigen::aligned_allocator<Vector3d> > reference_results;
    
    ifstream f;
    f.open("test-naca0012-reference.txt");
    
    while (f.good()) {
        string line;
        getline(f, line);
        
        if (line.length() == 0)
            break;
        
        istringstream tokens(line);
        
        Vector3d reference_result;
        tokens >> reference_result(0) >> reference_result(1) >> reference_result(2);
        
        reference_results.push_back(reference_result);
    }
    
    f.close();
    
    // Compute new values, and compare.
    for (unsigned int i = 0; i < reference_results.size(); i++) {        
        Vector3d &reference_result = reference_results[i];
        
        Vector2d res = run_test(reference_result(0) / 180.0 * M_PI);
        
        if (fabs(res[0] - reference_result(1)) > TEST_TOLERANCE) {
            cerr << " *** TEST FAILED *** " << endl;
            cerr << " alpha = " << reference_result(0) << " deg" << endl;
            cerr << " C_L(ref) = " << reference_result(1) << endl;
            cerr << " C_L = " << res(0) << endl;
            cerr << " ******************* " << endl;
            
            exit(1);
        }
        
        if (fabs(res[1] - reference_result(2)) > TEST_TOLERANCE) {
            cerr << " *** TEST FAILED *** " << endl;
            cerr << " alpha = " << reference_result(0) << " deg" << endl;
            cerr << " C_D(ref) = " << reference_result(2) << endl;
            cerr << " C_D = " << res(1) << endl;
            cerr << " ******************* " << endl;
            
            exit(1);
        }
    }
    
    // Done.
    return 0;
}
