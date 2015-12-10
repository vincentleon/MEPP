#ifndef SOFTICPCONTROLLER_H
#define SOFTICPCONTROLLER_H

#include "Correspondence_Polyhedron.h"
#include <set>
#include "../components/Analysis/Correspondence/src/SegmentController.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////

typedef std::vector<double> semanticDescr;

class deformationNode
{
public:
	std::set<Halfedge_handle> m_halfedges;
	std::vector<deformationNode *> m_childrenNodes;
};

class softICPController
{
public:
	
	// Constructor
	softICPController(PolyhedronPtr poly1, PolyhedronPtr poly2);
	
	void buildTreeStructure();

private:
	
	PolyhedronPtr m_polyhedron1;
	PolyhedronPtr m_polyhedron2;
	
	deformationNode m_treeStructure1;
	deformationNode m_treeStructure2;
	
	void getSnappingRegion();

};

#endif // SOFTICPCONTROLLER_H
