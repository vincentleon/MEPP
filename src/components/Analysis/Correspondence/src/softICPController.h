#ifndef SOFTICPCONTROLLER_H
#define SOFTICPCONTROLLER_H

#include "Correspondence_Polyhedron.h"
#include <set>
#include "../components/Analysis/Correspondence/src/SegmentController.h"
#include "../components/Analysis/Correspondence/src/mepp_component_Correspondence_plugin.hxx"

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


   struct Triangle_Cut {
		/*! \brief true if the facet belongs to the first polyhedron*/
		bool								Facet_from_A;
		/*! \brief An exact vector giving the direction of the normal*/
		Vector_exact						norm_dir;
		/*! \brief A list of segments (the intersections with the facets of the other polyhedron)*/
                std::vector<std::vector<InterId> >	CutList;
		/*! \brief A list of points (when the intersection is a point)*/
		std::set<InterId>					PtList;
		/*! \brief The list of the intersections*/
		std::map<HalfedgeId, InterId>		RefInter;

                /*! \brief Default constructor*/
                Triangle_Cut() {}
                /*! \brief Constructor
                 \param V : The normal direction
                 \param ffA : Must be true if the facet belongs to the first polyhedron*/
                Triangle_Cut(Vector_exact V, bool ffA) { norm_dir=V; Facet_from_A=ffA; } // MT
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

	void snapRegions(double R, double elasticity, int itermax=5, int treeDepth=2);
	
	pointTransformation computeTransformation(std::vector<Vertex_handle> & N, std::vector<Vertex_handle> & phiN);

	void applyTransformation(Vertex_handle p, pointTransformation & ti, int iter, const int itermax);
	
	PolyhedronPtr remesh();
	
private:
	
	double m_R;
	
	std::map< std::pair<Vertex_handle,Vertex_handle>, double> m_distanceMemo;
	
	std::map< Vertex_handle, bool > m_isSnappingRegion;
	
	bool m_stop_for_debug;
	
	std::vector< std::vector<deformationNode*> > m_levelNodes1;
	std::vector< std::vector<deformationNode*> > m_levelNodes2;
	
	std::vector < Matrix * > m_distanceMatrices1;
	std::vector < Matrix * > m_distanceMatrices2;	
	
	std::map<Vertex_handle,Vertex_handle> m_Phi; // the snapping Region correspondence
	std::map<Vertex_handle,double> m_distToLoop;
	
	PolyhedronPtr m_polyhedron1;
	PolyhedronPtr m_polyhedron2;
	
	deformationNode m_treeStructure1;
	deformationNode m_treeStructure2;
	
	std::vector<Halfedge_handle> m_loop1;
	std::vector<Halfedge_handle> m_loop2;
	
	std::vector<Triangle_Cut> m_inter_tri;
	
	std::map<Facet_handle, std::set<Facet_handle> > m_Couples;
	
	double getDistance(Vertex_handle v1, Vertex_handle v2, double w1 = 0.4, double w2=0.2, double w3=0.4);
	
	void getSnappingRegionAABB();
	
	void getSnappingRegionOld(double R, double squared_euclidean_radius);
	
	void computeSnappingRegionCorrespondence(bool order);
	
	void kMedoidMain(deformationNode * root, int k=4, int clusterIter=3);
	
	void computeMatrixDistance(deformationNode * root);
	
	void colorCluster(deformationNode * root);
	
	void colorLastClusters(deformationNode * root);
	
	void cluster(deformationNode * root);
	
	void updateRep(deformationNode * root);
	
	void hierarchicalBuild(deformationNode * root, int sizeOfTree, int level,int k);
	
	std::vector<Vertex_handle> getNeighborhoodOld(Vertex_handle p, double R, unsigned iter, unsigned itermax, double elasticity, bool order);
	std::vector<Vertex_handle> getNeighborhood(Vertex_handle p, double R, unsigned iter, unsigned itermax, double elasticity, bool order);
	
	std::vector<Vertex_handle> getCorrespondingNeighborhood( std::vector<Vertex_handle> & N);
	
	double computeDelta(bool order);
	
	void displayNode( deformationNode * root);
	
	void subdivideLeaves( deformationNode * root);
	
	void computeClosest(  deformationNode * root, int m, std::vector< std::vector<deformationNode*> > & levelNodes, std::vector < Matrix * > & distanceMatrices  );
	
	void storeNodeLevel( deformationNode * root, int level, std::vector< std::vector<deformationNode*> > & levelNodes);
	
	void fixBorder();
	
	void gaussianSmooth();
	
	void initClosest(std::vector< std::vector<deformationNode*> > & levelNodes, std::vector < Matrix * > & distanceMatrices, int m);
	
	void getMeshOutsideSR(std::vector<double> & coords, std::vector<int> & tris, PolyhedronPtr p, int vertexOffset);
	
	void getMeshInsideSR(vector< double >& coords, vector< int >& tris, PolyhedronPtr p, int vertexOffset);
	
	PolyhedronPtr buildSRMesh(std::vector<double> & coords, std::vector<int> & tris, int vertexOffset, std::set<Point3d> & border);
	
	void remeshSR(std::vector<double> & coords, std::vector<int> & tris, int vertexOffset);
	
	void findCouples(PolyhedronPtr meshA, PolyhedronPtr meshB);
	
	void floodConstruct(PolyhedronPtr meshA, PolyhedronPtr meshB);
	
	void rmCouple(Facet_handle A, Facet_handle B);
	
	void computeIntersections();
	
	void interTriangleTriangle(Facet_handle A, Facet_handle B);
	
	void cutIntersectedFacets(PolyhedronPtr meshA, PolyhedronPtr meshB);
};

double computePhiDistance(Vertex_handle v1, Vertex_handle v2, double w1 = 0.4, double w2=0.2, double w3=0.4);
//double computePhiDistance(Vertex_handle v1, Vertex_handle v2, double w1 = 0.7, double w2=0.3, double w3=0.0);


#endif // SOFTICPCONTROLLER_H
