//
// Vortexje -- Body.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __BODY_HPP__
#define __BODY_HPP__

#include <string>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>

#include <vortexje/surface.hpp>
#include <vortexje/lifting-surface.hpp>
#include <vortexje/wake.hpp>

namespace Vortexje
{

/**
   Container for surfaces.
   
   This class is designed for the application of kinematic operations to a body consisting of multiple surfaces.
   A typical application would be an airplane consisting of a fuselage and lifting surfaces.
   
   @brief Surface container.
*/
class Body
{
public:
    /**
       Body name.
    */
    std::string id;
    
    /**
       Data structure containing a non-lifting surface.
       
       @brief Surface data.
    */
    class SurfaceData {
    public:
        /**
           Constructs a SurfaceData object.
          
           @param[in]   surface          Surface object.
        */
        SurfaceData(std::shared_ptr<Surface> surface) :
            surface(surface) {}
        
        /**
           Associated surface object.
        */
        std::shared_ptr<Surface> surface;
    };
    
    /**
       Data structure grouping a lifting surface with its wake.
       
       @brief Lifting surface data.
    */
    class LiftingSurfaceData : public SurfaceData {
    public:
        /**
           Constructs a LiftingSurfaceData object.

           @param[in]   lifting_surface   Lifting surface object.
           @param[in]   wake              Wake object for this surface.
        */
        LiftingSurfaceData(std::shared_ptr<LiftingSurface> lifting_surface, std::shared_ptr<Wake> wake) :
            SurfaceData(lifting_surface), lifting_surface(lifting_surface), wake(wake) { }
        
        /**
           Associated lifting surface object.
        */
        std::shared_ptr<LiftingSurface> lifting_surface;
        
        /**
           Associated wake object.
        */
        std::shared_ptr<Wake> wake;
    };
    
    /**
       List of non-lifting surfaces.
    */
    std::vector<std::shared_ptr<SurfaceData> > non_lifting_surfaces;
    
    /**
       List of lifting surfaces.
    */
    std::vector<std::shared_ptr<LiftingSurfaceData> > lifting_surfaces;
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    Body(const std::string &id);
    
    ~Body();

    void add_non_lifting_surface(std::shared_ptr<Surface> surface);
    
    void add_lifting_surface(std::shared_ptr<LiftingSurface> lifting_surface);    
    void add_lifting_surface(std::shared_ptr<LiftingSurface> lifting_surface, std::shared_ptr<Wake> wake);
    
    /**
       Stitching of neighbour relationships between surfaces.
    */
    
    /**
       Data structure bundling a Surface, a panel ID, and an edge number.
       
       @brief Surface, panel ID, and edge bundle.
    */
    class SurfacePanelEdge {
    public:
        /**
           Constructor.
        */
        SurfacePanelEdge() {}
        
        /**
           Constructor.
           
           @param[in]   surface   Associated Surface object.
           @param[in]   panel     Panel ID.
        */
        SurfacePanelEdge(std::shared_ptr<Surface> surface, int panel, int edge) : surface(surface), panel(panel), edge(edge) { };
        
        /**
           Associated Surface object.
        */
        std::shared_ptr<Surface> surface;
        
        /**
           Panel ID.
        */
        int panel;
        
        /**
           Edge ID.
        */
        int edge;
        
        /**
           Equality operator.
        */
        bool operator==(const SurfacePanelEdge &other) const
        {
            return (surface.get() == other.surface.get()) && (panel == other.panel) && (edge == other.edge);
        }
    };
    
    void stitch_panels(std::shared_ptr<Surface> surface_a, int panel_a, int edge_a, std::shared_ptr<Surface> surface_b, int panel_b, int edge_b);
    
    std::vector<SurfacePanelEdge> panel_neighbors(const std::shared_ptr<Surface> &surface, int panel) const;
    
    std::vector<SurfacePanelEdge> panel_neighbors(const std::shared_ptr<Surface> &surface, int panel, int edge) const;
    
    /**
       Linear position of the entire body.
    */
    Eigen::Vector3d position;
    
    /**
       Linear velocity of the entire body.
    */
    Eigen::Vector3d velocity;
    
    /**
       Attitude (orientation) of the entire body.
    */
    Eigen::Quaterniond attitude;
    
    /**
       Rotational velocity of the entire body.
    */
    Eigen::Vector3d rotational_velocity;

    void set_position(const Eigen::Vector3d &position);
    void set_attitude(const Eigen::Quaterniond &attitude);
    
    void set_velocity(const Eigen::Vector3d &velocity);
    void set_rotational_velocity(const Eigen::Vector3d &rotational_velocity);
    
    Eigen::Vector3d panel_kinematic_velocity(const std::shared_ptr<Surface> &surface, int panel) const;
    
    Eigen::Vector3d node_kinematic_velocity(const std::shared_ptr<Surface> &surface, int node) const;
    
protected:
    /**
       Helper class to compare two SurfacePanelEdge objects:  first by surface, then by panel.  
       The edge numbers are not compared.
       
       @brief Helper class to compare two SurfacePanelEdge objects.
    */
    class CompareSurfacePanelEdge {
	public:
	    /**
	       Compare two SurfacePanelEdge objects: first by surface, then by panel.
	       
           @param[in]   a   First SurfacePanelEdge.
           @param[in]   b   Second SurfacePanelEdge.
	    */
		bool operator() (const SurfacePanelEdge a, const SurfacePanelEdge b) const {
		    if (a.surface->id == b.surface->id)
		        if (a.panel == b.panel)
		            return (a.edge < b.edge);
		        else
		            return (a.panel < b.panel);
		    else
		        return (a.surface->id < b.surface->id);
		}
    };
    
    /**
       List of stitches.
    */
    std::map<SurfacePanelEdge, SurfacePanelEdge, CompareSurfacePanelEdge> stitches;
};

};

#endif // __BODY_HPP__
