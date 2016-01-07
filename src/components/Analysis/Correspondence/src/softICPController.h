#ifndef SOFTICPCONTROLLER_H
#define SOFTICPCONTROLLER_H

#include "Correspondence_Polyhedron.h"
#include <set>
#include "../components/Analysis/Correspondence/src/SegmentController.h"

#include "CGAL/Linear_algebraCd.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////

typedef std::vector<double> semanticDescr;
typedef CGAL::Linear_algebraCd<double>::Matrix myMatrix;
typedef CGAL::Linear_algebraCd<double>::Vector myVector;

class deformationNode
{
public:
	//std::set<Halfedge_handle> m_halfedges;
	std::set<Vertex_handle> m_vertices;
	std::vector<deformationNode *> m_childrenNodes;
	deformationNode * m_parent;
 	int m_clusterID; // clusterID regarding Parents (-1 if root)
	
	Vertex_handle m_rep;
	std::map<Vertex_handle,int > m_cluster;
	
	myMatrix m_distanceMatrix;
};

class pointTransformation
{
public:
	qglviewer::Quaternion Q;
	qglviewer::Vec T;
};

class softICPController
{
public:
	
	// Constructor
	softICPController(PolyhedronPtr poly1, PolyhedronPtr poly2);
	
	void buildTreeStructure(const int sizeOfTree, double R=0.1, double squared_euclidean_radius=0.05);

	void snapRegions(double R, unsigned elasticity);
	
	pointTransformation computeTransformation(std::vector<Vertex_handle> & N, std::vector<Vertex_handle> & phiN);

	void applyTransformation(Vertex_handle p, pointTransformation & ti, int iter, const int itermax);
	
private:
	
	std::map<Vertex_handle,Vertex_handle> m_Phi; // the snapping Region correspondence
	std::map<Vertex_handle,double> m_distToLoop;
	
	PolyhedronPtr m_polyhedron1;
	PolyhedronPtr m_polyhedron2;
	
	deformationNode m_treeStructure1;
	deformationNode m_treeStructure2;
	
	void getSnappingRegion(double R, double squared_euclidean_radius);
	
	void computeSnappingRegionCorrespondence(bool order);
	
	void kMedoidMain(deformationNode * root, int k);
	
	void computeMatrixDistance(deformationNode * root);
	
	void colorCluster(deformationNode * root);
	
	void cluster(deformationNode * root);
	
	void updateRep(deformationNode * root);
	
	void hierarchicalBuild(deformationNode * root, int sizeOfTree, int level,int k);
	
	std::vector<Vertex_handle> getNeighborhood(Vertex_handle p, double R, unsigned iter, unsigned elasticity, bool order);
	
	std::vector<Vertex_handle> getCorrespondingNeighborhood( std::vector<Vertex_handle> & N);
	
	
	
	void displayNode( deformationNode * root);
};

double computePhiDistance(Vertex_handle v1, Vertex_handle v2, double w1 = 0.4, double w2=0.2, double w3=0.4);
//double computePhiDistance(Vertex_handle v1, Vertex_handle v2, double w1 = 0.7, double w2=0.3, double w3=0.0);


#endif // SOFTICPCONTROLLER_H
