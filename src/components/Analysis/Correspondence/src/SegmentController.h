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
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

#include "../components/Tools/Boolean_Operations/src/Boolean_Operations_Component.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel simpleKernel;
typedef CGAL::Polyhedron_3<simpleKernel,CGAL::Polyhedron_items_with_id_3> simplePolyhedron;

typedef boost::shared_ptr<Boolean_Operations_Component> Boolean_Operations_ComponentPtr;

typedef std::vector<Vertex_handle> border;

class SegmentController
{

public:
	SegmentController(PolyhedronPtr p);
	~SegmentController();
	
	void cutSegments();
	
	void glueSegments(Viewer * v);
	
	//void sewSegments(Viewer * v);
	
	void joinSegments(Viewer * v);
	
	void unionSegments(Viewer * v);
	
	void addMainPart( PolyhedronPtr p);
	
	PolyhedronPtr fillHoles( PolyhedronPtr p);
	
	void addSegement( PolyhedronPtr p);
	
	void alignSegments( Viewer *v, PolyhedronPtr s, PolyhedronPtr t, int sourceFrameID, int targetFrameID);
	
	std::vector<PolyhedronPtr> m_mainPart;
	std::vector<PolyhedronPtr> m_parts;
	
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

Point3d toWorld(Viewer * v, int p, Point3d lcp);
Point3d toPartFrame(Viewer * v, int p, Point3d wcp);

simplePolyhedron * convertToSimplePolyhedron(PolyhedronPtr p);
PolyhedronPtr convertToEnrichedPolyhedron(simplePolyhedron * p);

#endif // SEGMENTCONTROLLER_H
