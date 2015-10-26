#ifndef SEGMENTCONTROLLER_H
#define SEGMENTCONTROLLER_H

#include "../../../../mepp/Polyhedron/polyhedron.h"

#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include "Correspondence_Polyhedron.h"

#define CGAL_EIGEN3_ENABLED
#include "CGAL/Surface_mesh_deformation.h"
#include <set>

typedef CGAL::Simple_cartesian<double>	simpleKernel;
typedef CGAL::Polyhedron_3<simpleKernel,CGAL::Polyhedron_items_with_id_3> simplePolyhedron;
typedef CGAL::Surface_mesh_deformation<simplePolyhedron> surface_mesh_deformation;


class SegmentController
{

public:
	SegmentController(PolyhedronPtr p);
	~SegmentController();
	
	void cutSegments();
	
	void glueSegments();
	
	void fillHoles(Viewer * v, PolyhedronPtr p);
	
	void addMainPart( PolyhedronPtr p);
	
	void addSegement( PolyhedronPtr p);
	
	std::vector<PolyhedronPtr> m_mainPart;
	std::vector<PolyhedronPtr> m_parts;
	
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
 
bool isColoredBlack(Vertex_handle v);

void colorMesh(PolyhedronPtr p, float r, float g, float b);

void copyDescriptor(PolyhedronPtr s, PolyhedronPtr d);

simplePolyhedron * convertToSimplePolyhedron(PolyhedronPtr p);
PolyhedronPtr convertToEnrichedPolyhedron(simplePolyhedron * p);

#endif // SEGMENTCONTROLLER_H
