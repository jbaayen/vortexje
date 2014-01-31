//
// Vortexje -- Lifting surface.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <vortexje/lifting-surface.hpp>
#include <vortexje/parameters.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

/**
   Constructs an empty LiftingSurface.
*/
LiftingSurface::LiftingSurface() : Surface()
{
    // Do nothing.
}

/**
   Returns the number of chordwise nodes.

   @returns The number of chordwise nodes.
*/
int 
LiftingSurface::n_chordwise_nodes() const
{
    return (int) upper_nodes.rows();
}

/**
   Returns the number of chordwise panels.
   
   @returns The number of chordwise panels.
*/
int 
LiftingSurface::n_chordwise_panels() const
{
    return (int) upper_panels.rows();
}

/**
   Returns the number of spanwise nodes.

   @returns The number of spanwise nodes.
*/
int 
LiftingSurface::n_spanwise_nodes() const
{
    return (int) upper_nodes.cols();
}

/**
   Returns the number of spanwise panels.
   
   @returns The number of spanwise panels.
*/
int 
LiftingSurface::n_spanwise_panels() const
{
    return (int) upper_panels.cols();
}

/**
   Returns the index'th trailing edge node.
  
   @param[in]   index   Trailing edge node index.
   
   @returns The node number of the index'th trailing edge node.
*/
int 
LiftingSurface::trailing_edge_node(int index) const
{
    return upper_nodes(upper_nodes.rows() - 1, index);
}

/**
   Returns the index'th upper trailing edge panel.
  
   @param[in]   index   Trailing edge panel index.
   
   @returns The panel number of the index'th upper trailing edge panel.
*/
int
LiftingSurface::trailing_edge_upper_panel(int index) const
{
    return upper_panels(upper_panels.rows() - 1, index);
}

/**
   Returns the index'th lower trailing edge panel.
  
   @param[in]   index   Trailing edge panel index.
   
   @returns The panel number of the index'th lower trailing edge panel.
*/
int
LiftingSurface::trailing_edge_lower_panel(int index) const
{
    return lower_panels(lower_panels.rows() - 1, index);
}
