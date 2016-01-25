#ifndef SEGMENTCONTROLLER_H
#define SEGMENTCONTROLLER_H

#include "../../../../mepp/Polyhedron/polyhedron.h"

#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include "Correspondence_Polyhedron.h"

#include "libicp/icpPointToPoint.h"

#include <set>

#define CGAL_EIGEN3_ENABLED
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <boost/graph/graph_concepts.hpp>
//#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Surface_mesh_deformation.h>

#include "../components/Tools/Boolean_Operations/src/Boolean_Operations_Component.h"

//typedef CGAL::Exact_predicates_inexact_constructions_kernel simpleKernel;
typedef CGAL::Simple_cartesian<double> simpleKernel;
typedef CGAL::Polyhedron_3<simpleKernel,CGAL::Polyhedron_items_with_id_3> simplePolyhedron;

typedef boost::shared_ptr<Boolean_Operations_Component> Boolean_Operations_ComponentPtr;

typedef std::vector<Vertex_handle> border;

// For mesh deformation
typedef CGAL::Surface_mesh_deformation<simplePolyhedron> surface_mesh_deformation;
typedef boost::graph_traits<simplePolyhedron>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<simplePolyhedron>::halfedge_descriptor vertex_descriptor;
typedef std::map<vertex_descriptor, std::size_t>   Internal_vertex_map;
typedef std::map<halfedge_descriptor, std::size_t>     Internal_hedge_map;
typedef boost::associative_property_map<Internal_vertex_map>   Vertex_index_map;
typedef boost::associative_property_map<Internal_hedge_map>    Hedge_index_map;

struct snaxel
{
	
	snaxel * previous;
	Halfedge_handle he;
	snaxel * next;
};


// A model of SurfaceModelingWeights using a map of pre-computed weights
struct Weights_from_map
{
	typedef simplePolyhedron Halfedge_graph;
	Weights_from_map(std::map<halfedge_descriptor, double>* weight_map) : weight_map(weight_map)
	{}
	template<class VertexPointMap>
	double operator()(halfedge_descriptor e, Polyhedron& P,VertexPointMap vpm)
	{
		return (*weight_map)[e];
	}
	std::map<halfedge_descriptor,double>* weight_map;
};

class SegmentController
{

public:
	SegmentController(PolyhedronPtr p);
	~SegmentController();
	
	void cutSegments();
	
	double energyBorder(Halfedge_handle prev, Halfedge_handle curr, Halfedge_handle next,double alpha = 0.2, double beta = 0.8);
	
	void moveBorder(std::vector<Halfedge_handle> & snake);
	
	void getBorder();
	
	void optimizeBorders();
	
	snaxel * createSnake(std::vector<Halfedge_handle> & border);
	
	void moveSnake(snaxel * snake);
	
	void replaceSnaxel(snaxel * curr, Halfedge_handle candidate);
	
	void glueSegments(Viewer * v);
	
// 	//void sewSegments(Viewer * v);
	
	void fitSegments(Viewer * v, PolyhedronPtr target, PolyhedronPtr model, int & idBordTarget, int & idBordModel);
	
	void sewSegments(Viewer * v, PolyhedronPtr target, PolyhedronPtr model);
	
	void softICP(Viewer * v, PolyhedronPtr target, PolyhedronPtr model, double elasticity, double regionSize, int itermax);
	
	void fuseMeshes(Viewer *v, PolyhedronPtr target, PolyhedronPtr model);
	
	void fuseMeshes2(Viewer *v, PolyhedronPtr target, PolyhedronPtr model,  int idHETarget, int idHEModel);
	
	void joinSegments(Viewer * v);
	
	void unionSegments(Viewer * v);
	
	void addMainPart( PolyhedronPtr p);
	
	PolyhedronPtr fillHoles( PolyhedronPtr p);
	
	void addSegement( PolyhedronPtr p);
	
	void alignSegments( Viewer *v, PolyhedronPtr s, PolyhedronPtr t, int sourceFrameID, int targetFrameID);
	
	std::vector<PolyhedronPtr> m_mainPart;
	std::vector<PolyhedronPtr> m_parts;
	
	std::vector< std::vector<Halfedge_handle > > m_borders;
	std::vector< snaxel *> m_snakes;
	
	std::vector<std::vector<std::vector<Halfedge_handle> > > m_partsBorders;
	
private:
	
	PolyhedronPtr m_polyhedron;
};

template<class HDS>
class polyhedron_builder : public CGAL::Modifier_base<HDS> {
public:
	std::vector<double> &coords;
	std::vector<int>    &tris;
	polyhedron_builder( std::vector<double> &_coords, std::vector<int> &_tris ) : coords(_coords), tris(_tris) {}
	void operator()( HDS& hds) {
	typedef typename HDS::Vertex   Vertex;
	typedef typename Vertex::Point Point;
 
	// create a cgal incremental builder
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
	B.begin_surface( coords.size()/3, tris.size()/3 );
   
	// add the polyhedron vertices
	for( int i=0; i<(int)coords.size(); i+=3 )
	{
		B.add_vertex( Point( coords[i+0], coords[i+1], coords[i+2] ) );
	}
   
	// add the polyhedron triangles
	for( int i=0; i<(int)tris.size(); i+=3 )
	{
		B.begin_facet();
		B.add_vertex_to_facet( tris[i+0] );
		B.add_vertex_to_facet( tris[i+1] );
		B.add_vertex_to_facet( tris[i+2] );
		B.end_facet();
	}
   
	// finish up the surface
        B.end_surface();
    }
};

void visitVertexSelection(Halfedge_around_vertex_circulator h,
			  std::map<Vertex_iterator,bool> & isSelected, 
			  std::map<Vertex_iterator,int> & cc, int nbcc,
			  std::vector<std::set<Vertex_handle> > & ccVertices);

void visitVertexSelection(Halfedge_around_vertex_circulator h,
			  std::map<Vertex_iterator,bool> & isSelected, 
			  std::map<Vertex_iterator,int> & cc, int nbcc);
 
void visitBorder(Halfedge_handle h, std::map<Halfedge_handle,bool> & isVisited, std::vector<Halfedge_handle> & border);

bool isColoredBlack(Vertex_handle v);

void colorMesh(PolyhedronPtr p, float r, float g, float b);

void copyDescriptor(PolyhedronPtr s, PolyhedronPtr d);

void collectVertsAndFaces(PolyhedronPtr p, std::vector<double> & coords, std::vector<int> & faces);

Point3d toWorld(Viewer * v, int p, Point3d lcp);
Point3d toPartFrame(Viewer * v, int p, Point3d wcp);

double L2Dist(std::vector<double>& descr1, std::vector<double>& descr2);

simplePolyhedron * convertToSimplePolyhedron(PolyhedronPtr p);
PolyhedronPtr convertToEnrichedPolyhedron(simplePolyhedron * p);

void test_softICP_SVD();

#endif // SEGMENTCONTROLLER_H
