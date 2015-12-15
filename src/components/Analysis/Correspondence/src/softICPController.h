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
	Halfedge_handle m_rep;
	std::map<Halfedge_handle,int > m_cluster;
};

class softICPController
{
public:
	
	// Constructor
	softICPController(PolyhedronPtr poly1, PolyhedronPtr poly2);
	
	void buildTreeStructure(const int sizeOfTree);
	


private:
	
	std::map<Halfedge_handle,Halfedge_handle> m_Phi; // the snapping Region correspondence
	
	PolyhedronPtr m_polyhedron1;
	PolyhedronPtr m_polyhedron2;
	
	deformationNode m_treeStructure1;
	deformationNode m_treeStructure2;
	
	void getSnappingRegion();
	
	void computeSnappingRegionCorrespondence();
	
	void kMedoidMain(deformationNode * root, int k);
	
	void colorCluster(deformationNode * root);
	
	void cluster(deformationNode * root);
	
	void updateRep(deformationNode * root);
	
	void hierarchicalBuild(deformationNode * root, const int sizeOfTree, int level);

};

double computePhiDistance(Vertex_handle v1, Vertex_handle v2, double w1 = 0.4, double w2=0.2, double w3=0.4);

#endif // SOFTICPCONTROLLER_H
