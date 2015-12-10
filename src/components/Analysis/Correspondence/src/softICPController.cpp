#include "softICPController.h"
#include <CGAL/Kernel/interface_macros.h>
#include <CGAL/Sphere_3.h>


softICPController::softICPController(PolyhedronPtr m1, PolyhedronPtr m2) : m_polyhedron1(m1), m_polyhedron2(m2)
{

}

void softICPController::buildTreeStructure()
{
	/*	First, identify the snapping region 
	 *	on the two meshes.
	 */
	getSnappingRegion();
}

void softICPController::getSnappingRegion()
{
	double R = 0.05; //Region Size
	
	// First, normalize boundary edges order
	m_polyhedron1->normalize_border();
	m_polyhedron2->normalize_border();
	
	// Then we find the two border halfedges that are the closest to each other
	double distMin = std::numeric_limits<double>::max();
	Halfedge_handle border1,border2;
	
	for(auto b1 = m_polyhedron1->border_halfedges_begin(); b1!=m_polyhedron1->halfedges_end();++b1)
	{
		Point3d p1 = b1->vertex()->point();
		Halfedge_handle closestToP1;
		double distMinToP1 = std::numeric_limits<double>::max();
		for(auto b2 = m_polyhedron2->border_halfedges_begin(); b2!=m_polyhedron2->halfedges_end();++b2)
		{
			Point3d p2 = b2->vertex()->point();
			double dist = CGAL::squared_distance(p1,p2);
			if(dist<distMinToP1)
			{
				distMinToP1 = dist;
				closestToP1 = b2;
			}
		}
		if(distMinToP1<distMin)
		{
			border1 = b1;
			border2 = closestToP1;
		}
	}
	// Border 1 and border2 are on the two closest boundary loops
	
	// Define the snapping region
	Vertex_handle v1 = border1->vertex();
	semanticDescr & sd1 = v1->getSemantic();
	Vertex_handle v2 = border2->vertex();
	semanticDescr & sd2 = v2->getSemantic();
	
	// Add the two borders to the snapping region
	/*Halfedge_iterator it = border1;
	do{
		m_treeStructure1.m_halfedges.insert(it);
	}
	while(++it!=border1);
	it = border2;
	do{
		m_treeStructure2.m_halfedges.insert(it);
	}
	while(++it!=border2);*/
	// find the set of closest point to the other shape's boundary loop (euclidean distance)
	/*double euclidean_radius = 1.0;
	Tree tree(m_polyhedron1->facets_begin(),m_polyhedron1->facets_end());
	tree.accelerate_distance_queries();
	// Utiliser structure acceleratrice
	Point3d query = border2->vertex()->point();*/
	
	double squared_euclidean_radius = 0.01;

	
	Point3d center = border2->vertex()->point();
	for(auto pHalfedge = m_polyhedron1->halfedges_begin();pHalfedge!=m_polyhedron1->halfedges_end();++pHalfedge)
	{
 		Point3d p = pHalfedge->vertex()->point();
		if( CGAL::squared_distance(center,p) < squared_euclidean_radius)
		{	// point is candidate for semantic descriptor testing
			pHalfedge->vertex()->color(0,1,0);
			semanticDescr & heDescr = pHalfedge->vertex()->getSemantic();
			
			Halfedge_iterator it = border1;
			
			do{
				it->vertex()->color(1,1,0);
				semanticDescr & borderDescr = it->vertex()->getSemantic();
				std::cout << heDescr.size();
				std::cout << borderDescr.size() << std::endl;
				
				if(L2Dist(heDescr,borderDescr)<R)
				{
					m_treeStructure1.m_halfedges.insert(pHalfedge);
					pHalfedge->vertex()->color(1,0,0);
				}
				it = it->next();
			}
			while(it!=border1);
		}
	}
	////////////////////////////////////////////////////////////////////
	
	
	
}
