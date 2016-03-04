#ifndef SOFTICPCONTROLLER_H
#define SOFTICPCONTROLLER_H

#include <set>
#include "../components/Analysis/Correspondence/src/Correspondence_Component.h"
#include "../components/Analysis/Correspondence/src/SegmentController.h"
#include "../components/Analysis/Correspondence/src/mepp_component_Correspondence_plugin.hxx"

#include "../components/Analysis/Correspondence/src/libicp/matrix.h"

#include "CGAL/Linear_algebraCd.h"

#include "../components/Analysis/Correspondence/src/geodesic/geodesic_algorithm_dijkstra_alternative.h"

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
	
	std::vector<deformationNode *> m_closest; // the closest nodes
	
	std::vector<double> m_distClosest; // the distance to the closest nodes
	
	int m_level; // level in the tree
};

typedef std::vector<double>::const_iterator myiter;

struct ordering {
	bool operator ()(std::pair<deformationNode*, double> const& a, std::pair<deformationNode*, double> const& b)
	{
		return (a.second) < (b.second);
	}
};

class pointTransformation
{
public:
	qglviewer::Quaternion Q;
	Matrix R;
	qglviewer::Vec T;
	double S;
	qglviewer::Vec mu_t;
};

class softICPController
{
public:
	
	// Constructor
	softICPController(PolyhedronPtr poly1, PolyhedronPtr poly2);
	
	void buildTreeStructure(const int sizeOfTree, double R=0.1, double squared_euclidean_radius=0.05);
	
	void buildFullTreeStructure(const int sizeOfTree, double factor);

	void snapRegions(double R, double elasticity, int itermax=5, int treeDepth=2);
	
	pointTransformation computeTransformation(std::vector<Vertex_handle> & N, std::vector<Vertex_handle> & phiN);

	void applyTransformation(Vertex_handle p, pointTransformation & ti, int iter, const int itermax);
	
	void remesh(Viewer * v);
	
private:
	
	geodesic::Mesh m_g1;
	geodesic::Mesh m_g2;
	
	
	double m_R;
	
	std::map< std::pair<Vertex_handle,Vertex_handle>, double> m_distanceMemo;
	
	std::map< Vertex_handle, bool > m_isSnappingRegion;
	
	std::map< Vertex_handle, bool > m_isFactoredRegion;
	
	bool m_stop_for_debug;
	
	std::vector< std::vector<deformationNode*> > m_levelNodes1;
	std::vector< std::vector<deformationNode*> > m_levelNodes2;
	
	std::vector < Matrix * > m_distanceMatrices1;
	std::vector < Matrix * > m_distanceMatrices2;	
	
	std::map<Vertex_handle,Vertex_handle> m_Phi; // the snapping Region correspondence
	std::map<Vertex_handle,double> m_distToLoop; // distance to the mesh's b-loop
	std::map<Vertex_handle,double> m_distToSnap; // distance to the snapping region (inner loop)
		
	std::vector<Vertex_handle> m_sr1;
	std::vector<Vertex_handle> m_sr2;
	
	PolyhedronPtr m_polyhedron1;
	PolyhedronPtr m_polyhedron2;
	
	deformationNode m_treeStructure1;
	deformationNode m_treeStructure2;
	
	std::vector<Halfedge_handle> m_loop1;
	std::vector<Halfedge_handle> m_loop2;
	
	double getDistance(Vertex_handle v1, Vertex_handle v2, double w1 = 0.4, double w2=0.2, double w3=0.4);
	
	void getSnappingRegionAABB();
	
	void getSnappingRegion(const double factor);
	
	void getSnappingRegionGeo();
	
	void getSnappingRegionGeo(const double factor);
	
	void computeSnapRegionDistances();
	
	void computeSnappingRegionCorrespondence(bool order);
	
	void kMedoidMain(deformationNode * root, int k=4, int clusterIter=3);
	
	void computeMatrixDistance(deformationNode * root);
	
	void colorCluster(deformationNode * root);
	
	void colorLastClusters(deformationNode * root);
	
	void cluster(deformationNode * root);
	
	void updateRep(deformationNode * root);
	
	void hierarchicalBuild(deformationNode * root, int sizeOfTree, int level,int k);
	
	std::vector<Vertex_handle> getNeighborhoodOld(Vertex_handle p, double R, unsigned iter, unsigned itermax, double elasticity, bool order);
	
	std::vector<Vertex_handle> getNeighborhoodNoTree(Vertex_handle p, double R, unsigned iter, unsigned itermax, double elasticity, bool order);
	
	std::vector<Vertex_handle> getNeighborhood(Vertex_handle p, double R, unsigned iter, unsigned itermax, double elasticity, bool order);
	
	std::vector<Vertex_handle> getCorrespondingNeighborhood( std::vector<Vertex_handle> & N);
	
	double centerDistance(Facet_handle f);
	
	double computeDelta(bool order);
	
	void displayNode( deformationNode * root);
	
	void subdivideLeaves( deformationNode * root);
	
	void computeClosest(  deformationNode * root, int m, std::vector< std::vector<deformationNode*> > & levelNodes, std::vector < Matrix * > & distanceMatrices  );
	
	void computeClosestGeo(  deformationNode * root, int m, std::vector< std::vector<deformationNode*> > & levelNodes, std::vector < Matrix * > & distanceMatrices, geodesic::Mesh & g);
	
	void storeNodeLevel( deformationNode * root, int level, std::vector< std::vector<deformationNode*> > & levelNodes);
	
	void fixBorder();
	
	void gaussianSmooth();
	
	void initClosest(std::vector< std::vector<deformationNode*> > & levelNodes, std::vector < Matrix * > & distanceMatrices, int m);
	
	PolyhedronPtr getMeshOutsideSR();
	
	PolyhedronPtr getMeshInsideSR(PolyhedronPtr p);
	
	PolyhedronPtr buildSRMesh(std::set<Point3d> & border);
	
	void remeshSR(std::vector<double> & coords, std::vector<int> & tris, int vertexOffset);
	
	void findCouples(PolyhedronPtr meshA, PolyhedronPtr meshB);
	
	void floodConstruct(PolyhedronPtr meshA, PolyhedronPtr meshB);
	
	void rmCouple(Facet_handle A, Facet_handle B);
	
	void computeIntersections();
	
	void interTriangleTriangle(Facet_handle A, Facet_handle B);
	
	void cutIntersectedFacets(PolyhedronPtr meshA, PolyhedronPtr meshB);
	
	PolyhedronPtr stitchAndSmooth(PolyhedronPtr outMesh, PolyhedronPtr inMesh, std::set<Point3d> & borders);
	
	void finalTransform(bool order);
	
	void moveToCorrespondence(bool order);
};

double computePhiDistance(Vertex_handle v1, Vertex_handle v2, double w1 = 0.4, double w2=0.2, double w3=0.4);

void initGeodesicMesh(PolyhedronPtr p, geodesic::Mesh * g);

bool sphere_clip_vector(Point3d &O, double r,const Point3d &P, Vector &V);

#endif // SOFTICPCONTROLLER_H
