//
// Vortexje -- Ramasamy-Leishman wake model.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <cmath>

#include <vortexje/empirical-wakes/ramasamy-leishman-wake.hpp>
#include <vortexje/parameters.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

static const double pi = 3.141592653589793238462643383279502884;

// Default parameter values:
double RamasamyLeishmanWake::Parameters::fluid_kinematic_viscosity  = 15.68e-6;

double RamasamyLeishmanWake::Parameters::initial_vortex_core_radius = 0.07;

double RamasamyLeishmanWake::Parameters::min_vortex_core_radius     = 0.07;

double RamasamyLeishmanWake::Parameters::lambs_constant             = 1.25643;

double RamasamyLeishmanWake::Parameters::a_prime                    = 6.5e-5;

// Ramasamy-Leishman series data:
typedef struct {
    double vortex_reynolds_number;
    double a_1;
    double b_1;
    double a_2;
    double b_2;
    double b_3;
} ramasamy_leishman_data_row;

static ramasamy_leishman_data_row ramasamy_leishman_data[12] = { {    1, 1.0000, 1.2560, 0.0000, 0.00000, 0.0000},
                                                                 {  100, 1.0000, 1.2515, 0.0000, 0.00000, 0.0000},
                                                                 { 1000, 1.0000, 1.2328, 0.0000, 0.00000, 0.0000},
                                                                 {10000, 0.8247, 1.2073, 0.1753, 0.02630, 0.0000},
                                                                 {2.5e4, 0.5933, 1.3480, 0.2678, 0.01870, 0.2070},
                                                                 {4.8e4, 0.4602, 1.3660, 0.3800, 0.01380, 0.1674},
                                                                 {7.5e4, 0.3574, 1.3995, 0.4840, 0.01300, 0.1636},
                                                                 {  1e5, 0.3021, 1.4219, 0.5448, 0.01220, 0.1624},
                                                                 {2.5e5, 0.1838, 1.4563, 0.6854, 0.00830, 0.1412},
                                                                 {  5e5, 0.1386, 1.4285, 0.7432, 0.00580, 0.1144},
                                                                 {7.5e5, 0.1011, 1.4462, 0.7995, 0.00480, 0.1078},
                                                                 {  1e6, 0.0792, 1.4716, 0.8352, 0.00420, 0.1077} };
                                                                 
// Avoid having to divide by 4 pi all the time:
static const double one_over_4pi = 1.0 / (4 * pi);

/**
   Constructs an empty Ramasamy-Leishman wake.
   
   @param[in]   lifting_surface   Associated lifting surface.
*/
RamasamyLeishmanWake::RamasamyLeishmanWake(std::shared_ptr<LiftingSurface> lifting_surface): Wake(lifting_surface)
{
}

/**
   Adds new layer of wake panels.
*/
void
RamasamyLeishmanWake::add_layer()
{
    // Add layer:
    this->Wake::add_layer();
    
    // Add R-L data:
    if (n_panels() >= lifting_surface->n_spanwise_panels()) {
        for (int k = 0; k < lifting_surface->n_spanwise_panels(); k++) {
            int panel = n_panels() - lifting_surface->n_spanwise_panels() + k;
            
            // Add initial vortex core radii. 
            vector<double> panel_vortex_core_radii;
            for (int i = 0; i < 4; i++)
                panel_vortex_core_radii.push_back(RamasamyLeishmanWake::Parameters::initial_vortex_core_radius);
            vortex_core_radii.push_back(panel_vortex_core_radii);
            
            // Store base edge lengths.
            vector<double> edge_lengths;
            for (int i = 0; i < 4; i++) {
                int prev_idx;
                if (i == 0)
                    prev_idx = 3;
                else
                    prev_idx = i - 1;
                    
                const Vector3d &node_a = nodes[panel_nodes[panel][prev_idx]];
                const Vector3d &node_b = nodes[panel_nodes[panel][i]];
                
                Vector3d edge = node_b - node_a;
                edge_lengths.push_back(edge.norm());
            }
            base_edge_lengths.push_back(edge_lengths);
        }
    }
}

/**
   Computes the unit velocity induced by a Ramasamy-Leishman vortex ring.
   
   @param[in]   x            Point at which the velocity is evaluated.
   @param[in]   this_panel   Panel on which the vortex ring is located.
   
   @returns Unit velocity induced by the Ramasamy-Leishman vortex ring.
   
   @note See M. Ramasamy and J. G. Leishman, Reynolds Number Based Blade Tip Vortex Model, University of Maryland, 2005.
*/
Vector3d
RamasamyLeishmanWake::vortex_ring_unit_velocity(const Eigen::Vector3d &x, int this_panel) const
{
    if (this_panel >= n_panels() - lifting_surface->n_spanwise_panels()) {
        // This panel is contained in the latest row of wake panels.  To satisfy the Kutta condition
        // exactly, we use the unmodified vortex ring unit velocity here.
        return this->Surface::vortex_ring_unit_velocity(x, this_panel);
    }
     
    // Compute vortex Reynolds number:                          
    double vortex_reynolds_number = doublet_coefficients[this_panel] / RamasamyLeishmanWake::Parameters::fluid_kinematic_viscosity;
    
    // Interpolate Ramasamy-Leishman series values piecewise-linearly:
    int less_than_idx;
    for (less_than_idx = 0; less_than_idx < 12; less_than_idx++) {
        ramasamy_leishman_data_row &row = ramasamy_leishman_data[less_than_idx];
        
        if (vortex_reynolds_number < row.vortex_reynolds_number)
            break;
    }
    
    double a[3];
    double b[3];
    if (less_than_idx == 0) {
        a[0] = ramasamy_leishman_data[0].a_1;
        a[1] = ramasamy_leishman_data[0].a_2;
        
        b[0] = ramasamy_leishman_data[0].b_1;
        b[1] = ramasamy_leishman_data[0].b_2;
        b[2] = ramasamy_leishman_data[0].b_3;
    } else if (less_than_idx == 12) {
        a[0] = ramasamy_leishman_data[11].a_1;
        a[1] = ramasamy_leishman_data[11].a_2;
        
        b[0] = ramasamy_leishman_data[11].b_1;
        b[1] = ramasamy_leishman_data[11].b_2;
        b[2] = ramasamy_leishman_data[11].b_3;
    } else {
        double one_over_delta_vortex_reynolds_number =
            1.0 / (ramasamy_leishman_data[less_than_idx].vortex_reynolds_number - ramasamy_leishman_data[less_than_idx - 1].vortex_reynolds_number);
        double x = vortex_reynolds_number - ramasamy_leishman_data[less_than_idx - 1].vortex_reynolds_number;
        double slope;
        
        slope = (ramasamy_leishman_data[less_than_idx].a_1 - ramasamy_leishman_data[less_than_idx - 1].a_1) * one_over_delta_vortex_reynolds_number;
        a[0] = ramasamy_leishman_data[less_than_idx - 1].a_1 + slope * x;
        
        slope = (ramasamy_leishman_data[less_than_idx].a_2 - ramasamy_leishman_data[less_than_idx - 1].a_2) * one_over_delta_vortex_reynolds_number;
        a[1] = ramasamy_leishman_data[less_than_idx - 1].a_2 + slope * x;
        
        slope = (ramasamy_leishman_data[less_than_idx].b_1 - ramasamy_leishman_data[less_than_idx - 1].b_1) * one_over_delta_vortex_reynolds_number;
        b[0] = ramasamy_leishman_data[less_than_idx - 1].b_1 + slope * x;
        
        slope = (ramasamy_leishman_data[less_than_idx].b_2 - ramasamy_leishman_data[less_than_idx - 1].b_2) * one_over_delta_vortex_reynolds_number;
        b[1] = ramasamy_leishman_data[less_than_idx - 1].b_2 + slope * x;
        
        slope = (ramasamy_leishman_data[less_than_idx].b_3 - ramasamy_leishman_data[less_than_idx - 1].b_3) * one_over_delta_vortex_reynolds_number;
        b[2] = ramasamy_leishman_data[less_than_idx - 1].b_3 + slope * x;
    }
    
    a[2] = 1 - a[0] - a[1];
    
    // Compute velocity:
    Vector3d velocity(0, 0, 0);
    
    for (int i = 0; i < (int) panel_nodes[this_panel].size(); i++) {
        int previous_idx;
        if (i == 0)
            previous_idx = panel_nodes[this_panel].size() - 1;
        else
            previous_idx = i - 1;
            
        const Vector3d &node_a = nodes[panel_nodes[this_panel][previous_idx]];
        const Vector3d &node_b = nodes[panel_nodes[this_panel][i]];
        
        Vector3d r_0 = node_b - node_a;
        Vector3d r_1 = node_a - x;
        Vector3d r_2 = node_b - x;
        
        double r_0_norm = r_0.norm();
        double r_1_norm = r_1.norm();
        double r_2_norm = r_2.norm();
        
        Vector3d r_1xr_2 = r_1.cross(r_2);
        double r_1xr_2_sqnorm = r_1xr_2.squaredNorm();
        double r_1xr_2_norm = sqrt(r_1xr_2_sqnorm);
        
        if (r_0_norm < Vortexje::Parameters::inversion_tolerance ||
            r_1_norm < Vortexje::Parameters::inversion_tolerance ||
            r_2_norm < Vortexje::Parameters::inversion_tolerance ||
            r_1xr_2_sqnorm < Vortexje::Parameters::inversion_tolerance)
            continue;
            
        double d = r_1xr_2_norm / r_0_norm;
        
        double dr = pow(d / vortex_core_radii[this_panel][i], 2);
            
        double sum = 0;
        for (int j = 0; j < 3; j++)
            sum += a[j] * exp(-b[j] * dr);

        velocity += (1 - sum) * r_1xr_2 / r_1xr_2_sqnorm * r_0.dot(r_1 / r_1_norm - r_2 / r_2_norm);
    }

    return one_over_4pi * velocity;
}

/**
   Updates the Ramasamy-Leishman vortex ring core radii.
  
   @param[in]   panel   Panel number.
   @param[in]   dt      Time step size.
   
   @note See M. Ramasamy and J. G. Leishman, Reynolds Number Based Blade Tip Vortex Model, University of Maryland, 2005.
*/
void
RamasamyLeishmanWake::update_vortex_ring_radii(int panel, double dt)
{
    for (int i = 0; i < 4; i++) {
        int prev_idx;
        if (i == 0)
            prev_idx = 3;
        else
            prev_idx = i - 1;
            
        double vortex_reynolds_number = fabs(doublet_coefficients[panel]) / RamasamyLeishmanWake::Parameters::fluid_kinematic_viscosity;
        
        double t_multiplier = 4 * RamasamyLeishmanWake::Parameters::lambs_constant *
            (1 + vortex_reynolds_number * RamasamyLeishmanWake::Parameters::a_prime) *
             RamasamyLeishmanWake::Parameters::fluid_kinematic_viscosity;
        
        double t = (pow(vortex_core_radii[panel][i], 2) - pow(RamasamyLeishmanWake::Parameters::initial_vortex_core_radius, 2)) / t_multiplier;
        
        double vortex_core_size_0 = sqrt(pow(RamasamyLeishmanWake::Parameters::initial_vortex_core_radius, 2) + t_multiplier * (t + dt));
        
        Vector3d node_a = nodes[panel_nodes[panel][prev_idx]];
        Vector3d node_b = nodes[panel_nodes[panel][i]];
        
        Vector3d edge = node_b - node_a;
        
        double strain = (edge.norm() - base_edge_lengths[panel][i]) / base_edge_lengths[panel][i];
        
        vortex_core_radii[panel][i] = fmax(RamasamyLeishmanWake::Parameters::min_vortex_core_radius, vortex_core_size_0 / sqrt(1 + strain));
    }
}

/**
   Updates the Ramasamy-Leishman vortex ring core radii.
  
   @param[in]   dt   Time step size.
   
   @note See M. Ramasamy and J. G. Leishman, Reynolds Number Based Blade Tip Vortex Model, University of Maryland, 2005.
*/
void
RamasamyLeishmanWake::update_properties(double dt)
{
    int i;
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic, 1)
        for (i = 0; i < n_panels(); i++)
            update_vortex_ring_radii(i, dt);
    }
}
