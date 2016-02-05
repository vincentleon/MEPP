#include "softICPController.h"
#include <CGAL/Kernel/interface_macros.h>
#include <CGAL/Sphere_3.h>
#include "time.h"
#include <list>
#include "omp.h"

// For quick distance computation (Snapping region)
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>

// For re-meshing ( snapping region )
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>

#include <CGAL/convex_hull_3_to_polyhedron_3.h>

#include <CGAL/Advancing_front_surface_reconstruction.h>

#include <CGAL/Subdivision_method_3.h>

typedef CGAL::Simple_cartesian<double> AABB_Kernel;
typedef CGAL::AABB_polyhedron_triangle_primitive<AABB_Kernel,Polyhedron> AABB_Primitive;
typedef CGAL::AABB_traits<AABB_Kernel, AABB_Primitive> AABB_Traits;
typedef CGAL::AABB_tree<AABB_Traits> AABB_Tree;
typedef AABB_Tree::Object_and_primitive_id Object_and_primitive_id;
typedef AABB_Tree::Point_and_primitive_id Point_and_primitive_id;

typedef CGAL::Exact_predicates_inexact_constructions_kernel delaunayKernel;
typedef CGAL::Triangulation_3<AABB_Kernel>      Triangulation;
typedef Triangulation::Cell_handle    Cell_handle;
typedef Triangulation::Vertex_handle  VoronoiVertex_handle;
typedef Triangulation::Locate_type    Locate_type;
typedef Triangulation::Point          Point;
typedef AABB_Kernel::Triangle_3 Triangle;

typedef CGAL::cpp11::array<std::size_t,3> RFacet;


struct Perimeter {
  double bound;
  Perimeter(double bound)
    : bound(bound)
  {}
  // The point type that will be injected here will be
  // CGAL::Exact_predicates_inexact_constructions_kernel::Point_3
  template <typename Point>
  bool operator()(const Point& p, const Point& q, const Point& r) const
  {
    // bound == 0 is better than bound < infinity
    // as it avoids the distance computations
    if(bound == 0){
      return false;
    }
    double d  = sqrt(squared_distance(p,q));
    if(d>bound) return true;
    d += sqrt(squared_distance(p,r)) ;
    if(d>bound) return true;
    d+= sqrt(squared_distance(q,r));
    return d>bound;
  }
};

struct DistanceSurface {
	AABB_Tree & m_tree;
	std::set<Point3d> & m_borders;
	DistanceSurface(AABB_Tree & tree, std::set<Point3d> & borders) : m_tree(tree), m_borders(borders){}
	
	
	template <typename Point>
	bool operator()(const Point& p, const Point& q, const Point& r) const
	{	
		
		Point3d pp(p.x(),p.y(),p.z());
		Point3d qq(q.x(),q.y(),q.z());
		Point3d rr(r.x(),r.y(),r.z());
		
		
		return m_borders.count(pp) && m_borders.count(qq) && m_borders.count(rr);
	}
	
};


std::stack<clock_t> timer_stack;

void timer_tic() {
  timer_stack.push(clock());
}

void timer_toc() {
    std::cout << "Time elapsed: "
              << ((double)(clock() - timer_stack.top())) / CLOCKS_PER_SEC
              << std::endl;
    timer_stack.pop();
}

softICPController::softICPController(PolyhedronPtr m1, PolyhedronPtr m2) : m_polyhedron1(m1), m_polyhedron2(m2), m_stop_for_debug(false)
{

}

// MEMOIZE DISTANCES
double softICPController::getDistance(Vertex_handle v1, Vertex_handle v2, double w1, double w2, double w3)
{
	std::pair<Vertex_handle,Vertex_handle> key(v1,v2);
	std::pair<Vertex_handle,Vertex_handle> skey(v2,v1);
	
	auto itr = m_distanceMemo.find(key);
	if(itr == m_distanceMemo.end())
	{
		double dist = computePhiDistance(v1,v2,w1,w2,w3);
		m_distanceMemo.insert(make_pair(key,dist));
		m_distanceMemo.insert(make_pair(skey,dist));
		return dist;
	}
	else
	{
		return itr->second;
	}
}

void softICPController::computeSnappingRegionCorrespondence(bool order)
{
	/*for every vertex in the snapping region, find its correspondence 
	 * that minimizes the weighted sum
	 */
	
	m_Phi.clear(); // recompute correspondence at each iteration
	
	//std::set<Halfedge_handle> * s1; 
	//std::set<Halfedge_handle> * s2; 
	std::set<Vertex_handle> * s1; 
	std::set<Vertex_handle> * s2; 
	
	if(order)
	{
		s1 = &m_treeStructure1.m_vertices;
		s2 = &m_treeStructure2.m_vertices;
	}
	else // reverse order between meshes
	{
		s2 = &m_treeStructure1.m_vertices;
		s1 = &m_treeStructure2.m_vertices;
	}
	
	
	for(auto it = s1->begin(); it!=s1->end();++it)
	{
		double distMin = std::numeric_limits<double>::max();
		Vertex_handle bestCorres;
		for(auto c = s2->begin(); c!=s2->end(); ++c)
		{
			//double dist = computePhiDistance((*it),(*c));
			//double dist = computePhiDistance((*it),(*c),0.7,0.3,0.0);
			double dist = computePhiDistance((*it),(*c),0.5,0.1,0.4);
			if(dist < distMin)
			{
				bestCorres = *c;
				distMin = dist;
			}
		}
		m_Phi[(*it)] = bestCorres;
	}	
}

void softICPController::computeClosest( deformationNode * root, int m, std::vector< std::vector<deformationNode*> > & levelNodes, std::vector < Matrix * > & distanceMatrices )
{
	// Compute full matrix distance at each level
	for(unsigned l=0; l<levelNodes.size();++l)
	{
		int sizeLevel = levelNodes[l].size();
		distanceMatrices.push_back(new Matrix(sizeLevel,sizeLevel));
		Matrix * mat = distanceMatrices[l];
		for(unsigned p=0;p<sizeLevel;++p)
		{
			deformationNode* nodeP = levelNodes[l][p];
			for(unsigned q=0;q<p;++q)
			{
				deformationNode* nodeQ = levelNodes[l][q];
					double dist = computePhiDistance(nodeP->m_rep,nodeQ->m_rep,1.0,0.0,0.0);
				//double dist = computePhiDistance(nodeP->m_rep,nodeQ->m_rep,0.5,0.1,0.4);
				mat->val[p][q] = dist;
				mat->val[q][p] = dist;
			}
		}
	}
}

void softICPController::initClosest(std::vector< std::vector<deformationNode*> > & levelNodes, std::vector < Matrix * > & distanceMatrices, int m)
{
	for(unsigned l=0;l<levelNodes.size();++l)
	{
		int sizeLevel = levelNodes[l].size();
		Matrix * mat = distanceMatrices[l];
		for(unsigned i=0;i<sizeLevel;++i)
		{
			//sort levelNodes according to m_distanceMatrices[m][l]
			std::vector<std::pair<deformationNode*, double> > order(sizeLevel);
			int j=0;
			for(auto it=levelNodes[l].begin(); it!=levelNodes[l].end();++it,++j)
			{
				order[j] = std::make_pair(*it,mat->val[i][j]);
			}
			std::sort(order.begin(), order.end(), ordering());
			for(unsigned c=0;(c<=m && c<order.size());++c)
			{
				levelNodes[l][i]->m_closest.push_back(order[c].first);
				levelNodes[l][i]->m_distClosest.push_back(order[c].second);
			}
			//std::cout << levelNodes[l][i]->m_closest.size() << std::endl;
		}
		//std::cout << std::endl;
	}
	
	//std::cout << "initClosest Ok" << std::endl;
}

void softICPController::storeNodeLevel(deformationNode* root, int level, std::vector< std::vector<deformationNode*> > & levelNodes)
{
	if(levelNodes.size()<=level){ levelNodes.resize(level+1);}
	levelNodes[level].push_back(root);
	for(unsigned c=0;c<root->m_childrenNodes.size();++c)
	{
		storeNodeLevel(root->m_childrenNodes[c],level+1,levelNodes);
	}
}


double softICPController::computeDelta(bool order)
{
	double delta = 0.0;
	std::set<Vertex_handle> * s1; 
	
	if(order)
	{
		s1 = &m_treeStructure1.m_vertices;
	}
	else // reverse order between meshes
	{
		s1 = &m_treeStructure2.m_vertices;
	}
	for(auto it = s1->begin(); it!=s1->end(); ++it)
	{
		delta += CGAL::squared_distance((*it)->point(),m_Phi[(*it)]->point());
	}
	return (sqrt(delta)/s1->size());
}


void softICPController::buildTreeStructure(const int sizeOfTree,double R, double squared_euclidean_radius)
{
	// Identify the snapping region on the two meshes
	//std::cout << "\n\t\tget Snapping region ";
	//timer_tic();
	//getSnappingRegionOld(R,squared_euclidean_radius);
	getSnappingRegionAABB();
	//timer_toc();
	R = m_R;
	
	//std::cout << "\t\thierarchicalBuild " << std::endl;
	//timer_tic();
	//#pragma omp parallel sections
	//{
	//	#pragma omp section
	//	{
			hierarchicalBuild(&m_treeStructure1,sizeOfTree,0,4);
	//	}
	//	#pragma omp section
	//	{
			hierarchicalBuild(&m_treeStructure2,sizeOfTree,0,4);
	//	}
	//}
	//timer_toc();
	
	// Init. root nodes properly
	m_treeStructure1.m_parent = NULL;
	m_treeStructure2.m_parent = NULL;
	m_treeStructure1.m_clusterID = -1;
	m_treeStructure2.m_clusterID = -1;
}

void softICPController::colorLastClusters(deformationNode* root)
{
	if(root->m_childrenNodes.size() == 0)
	{
		deformationNode * parent = root->m_parent;
		colorCluster(parent);
	}
	else
	{
		for(unsigned c = 0; c < root->m_childrenNodes.size(); ++c)
		{
			colorLastClusters(root->m_childrenNodes[c]);
		}
	}
}

void softICPController::hierarchicalBuild(deformationNode* root, const int sizeOfTree, int level, int k)
{
	root->m_level = level;
	
	if(level != sizeOfTree)
	{
		if(root->m_vertices.size() > k)
		{
			kMedoidMain(root,k);
			for(unsigned c=0;c<root->m_childrenNodes.size();++c)
			{	
				hierarchicalBuild(root->m_childrenNodes[c],sizeOfTree,level+1,k);
			}
		}
	}
}

void softICPController::displayNode(deformationNode* root)
{
	std::cout << root->m_childrenNodes.size() << " how many children nodes" << std::endl;
	std::cout << root->m_cluster.size() << " size of patch vertices IDs" <<  std::endl;
	std::cout << root->m_clusterID << " id wrt the parent node" <<  std::endl;
	std::cout << root->m_vertices.size()<< " how many vertices in this patch" << std::endl;	
	std::cout << std::endl;
	for(unsigned c=0;c <root->m_childrenNodes.size(); ++c)
	{
		displayNode(root->m_childrenNodes[c]);
	}
}




void softICPController::colorCluster(deformationNode* root)
{
	if(root->m_childrenNodes.size() != 4)
	{
		return;
	}
	for(auto h = root->m_childrenNodes[0]->m_vertices.begin(); 
	    h != root->m_childrenNodes[0]->m_vertices.end();
		++h)
	{
		(*h)->color(1,0,0);
	}
	for(auto h = root->m_childrenNodes[1]->m_vertices.begin(); 
	    h != root->m_childrenNodes[1]->m_vertices.end();
		++h)
	{
		(*h)->color(1,1,0);
	}
	for(auto h = root->m_childrenNodes[2]->m_vertices.begin(); 
	    h != root->m_childrenNodes[2]->m_vertices.end();
		++h)
	{
		(*h)->color(0,0,1);
	}
	for(auto h = root->m_childrenNodes[3]->m_vertices.begin(); 
	    h != root->m_childrenNodes[3]->m_vertices.end();
		++h)
	{
		(*h)->color(0,1,0);
	}
}

void softICPController::getSnappingRegionAABB()
{
	
	// First, normalize boundary edges order
	m_polyhedron1->normalize_border();
	m_polyhedron2->normalize_border();
	
	// Then we find the two border halfedges that are the closest to each other
	double distMin = std::numeric_limits<double>::max();
	Halfedge_handle border1,border2;
	
	//#pragma omp parallel for private(b1)
	for(auto b1 = m_polyhedron1->border_halfedges_begin(); b1!=m_polyhedron1->halfedges_end();++b1)
	{
		//b1->vertex()->color(1,0,0);
		Point3d p1 = b1->vertex()->point();
		Halfedge_iterator closestToP1;
		double distMinToP1 = std::numeric_limits<double>::max();
		for(auto b2 = m_polyhedron2->border_halfedges_begin(); b2!=m_polyhedron2->halfedges_end();++b2)
		{
			//b2->vertex()->color(0,0,1);
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
			distMin = distMinToP1;
		}
	}
	// Border 1 and border2 are on the two closest boundary loops
	
	std::vector<Vertex_handle> closestToB1;
	std::vector<Vertex_handle> closestToB2;
	
	std::vector<Vertex_handle> loop1;
	std::vector<Vertex_handle> loop2;
	
	AABB_Tree tree1(m_polyhedron1->facets_begin(),m_polyhedron1->facets_end());
	AABB_Tree tree2(m_polyhedron2->facets_begin(),m_polyhedron2->facets_end());
	
	tree1.accelerate_distance_queries();
	tree2.accelerate_distance_queries();
	
	// Find the set of closest points to the other shapes B-loop
	Halfedge_handle it = border1;
	do{
		loop1.push_back(it->vertex());
		m_loop1.push_back(it);
		// find the vertex closest to it->vertex() in m_polyhedron2
		Point_and_primitive_id pp = tree2.closest_point_and_primitive(it->vertex()->point());
		Facet_iterator f_nearest = pp.second;
		closestToB1.push_back(f_nearest->facet_begin()->vertex());
		
		// progress alongside border
		Halfedge_handle next;
		Halfedge_around_vertex_circulator hC = it->vertex()->vertex_begin();
		do
		{
			if(hC->opposite()->is_border_edge() && hC->opposite()!=it->opposite())
			{
				next = hC->opposite();
				break;
			}
			
		}
		while(++hC!=it->vertex()->vertex_begin());
		
		it = next;
	}
	while(it!=border1);

	it = border2;
	do{
		loop2.push_back(it->vertex());
		m_loop2.push_back(it);
		// find the vertex closest to it->vertex() in m_polyhedron2
		Point_and_primitive_id pp = tree1.closest_point_and_primitive(it->vertex()->point());
		Facet_iterator f_nearest = pp.second;
		closestToB2.push_back(f_nearest->facet_begin()->vertex());
		
		// progress alongside border
		Halfedge_handle next;
		Halfedge_around_vertex_circulator hC = it->vertex()->vertex_begin();
		do
		{
			if(hC->opposite()->is_border_edge() && hC->opposite()!=it->opposite())
			{
				next = hC->opposite();
				break;
			}
			
		}
		while(++hC!=it->vertex()->vertex_begin());
		it = next;
	}
	while(it!=border2);
	
	// compute for each point in the two sets, compute the distance to the other shapes b-loop
	double R = 0.0;
	for(int i=0;i<closestToB2.size();++i)
	{	
		Point3d p = closestToB2[i]->point();
		double closest = std::numeric_limits<double>::max();
		for(int b=0;b<loop1.size();++b)
		{
			Point3d q = loop1[b]->point();
			double dist = CGAL::squared_distance(p,q);
			if(dist<closest)
			{
				closest = dist;
			}
		}
		if(closest>R){R = closest;}
	}
	for(int i=0;i<closestToB1.size();++i)
	{
		Point3d p = closestToB1[i]->point();
		double closest = std::numeric_limits<double>::max();
		for(int b=0;b<loop2.size();++b)
		{
			Point3d q = loop2[b]->point();
			double dist = CGAL::squared_distance(p,q);
			if(dist<closest)
			{
				closest = dist;
			}
		}
		if(closest>R){R = closest;}
	}
	m_R = sqrt(R);
	
	// R is the maximum of these distances
	// The SR is the submesh within R distance to the boundary loop
	for(auto pVertex = m_polyhedron1->vertices_begin();
	    pVertex != m_polyhedron1->vertices_end();
		++pVertex)
	{
		Point3d p = pVertex->point();
		double distToLoop = std::numeric_limits<double>::max();
		for(auto b=loop1.begin();b!=loop1.end();++b)
		{
			double dist = CGAL::squared_distance(p,(*b)->point());
			if(dist<distToLoop)
			{
				distToLoop = dist;
			}
		}
		distToLoop = sqrt(distToLoop);
		if(distToLoop < m_R)
		{
			//pVertex->color(1.0,1.0,0.0);
			m_distToLoop[pVertex] = distToLoop;
			m_treeStructure1.m_vertices.insert(pVertex);
			m_isSnappingRegion[pVertex] = true;
		}
	}
	
	for(auto pVertex = m_polyhedron2->vertices_begin();
	    pVertex != m_polyhedron2->vertices_end();
		++pVertex)
	{
		Point3d p = pVertex->point();
		double distToLoop = std::numeric_limits<double>::max();
		for(auto b=loop2.begin();b!=loop2.end();++b)
		{
			double dist = CGAL::squared_distance(p,(*b)->point());
			if(dist<distToLoop)
			{
				distToLoop = dist;
			}
		}
		distToLoop = sqrt(distToLoop);
		if(distToLoop < m_R)
		{
			//pVertex->color(1.0,1.0,0.0);
			m_distToLoop[pVertex] = distToLoop;
			m_treeStructure2.m_vertices.insert(pVertex);
			m_isSnappingRegion[pVertex] = true;
		}
	}
}


void softICPController::getSnappingRegionOld(double R, double squared_euclidean_radius)
{
	// First, normalize boundary edges order
	m_polyhedron1->normalize_border();
	m_polyhedron2->normalize_border();
	
	// Then we find the two border halfedges that are the closest to each other
	double distMin = std::numeric_limits<double>::max();
	Halfedge_handle border1,border2;
	
	//#pragma omp parallel for private(b1)
	for(auto b1 = m_polyhedron1->border_halfedges_begin(); b1!=m_polyhedron1->halfedges_end();++b1)
	{
		b1->vertex()->color(1,0,0);
		Point3d p1 = b1->vertex()->point();
		Halfedge_iterator closestToP1;
		double distMinToP1 = std::numeric_limits<double>::max();
		for(auto b2 = m_polyhedron2->border_halfedges_begin(); b2!=m_polyhedron2->halfedges_end();++b2)
		{
			b2->vertex()->color(0,0,1);
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
			distMin = distMinToP1;
		}
	}
	// Border 1 and border2 are on the two closest boundary loops
	
	// Define the snapping region
	Vertex_handle v1 = border1->vertex();
	semanticDescr & sd1 = v1->getSemantic();
	Vertex_handle v2 = border2->vertex();
	semanticDescr & sd2 = v2->getSemantic();
	
	Point3d center = border2->vertex()->point();
	
		/*Halfedge_handle it = border1;
		do{
			Halfedge_handle next;
			Halfedge_around_vertex_circulator hC = it->vertex()->vertex_begin();
			do
			{
				if(hC->opposite()->is_border_edge() && hC->opposite()!=it->opposite())
				{
					next = hC->opposite();
					break;
				}
				
			}
			while(++hC!=it->vertex()->vertex_begin());
			it->vertex()->color(1,1,0);
			it = next;
		}
		while(it!=border1);*/
	//#pragma omp parallel for private(pHalfedge)
	for(auto pHalfedge = m_polyhedron1->halfedges_begin();pHalfedge!=m_polyhedron1->halfedges_end();++pHalfedge)
	{
 		Point3d p = pHalfedge->vertex()->point();
		//if( CGAL::squared_distance(center,p) < squared_euclidean_radius)
		//{	// point is candidate for semantic descriptor testing
			//pHalfedge->vertex()->color(0,1,0);
			semanticDescr & heDescr = pHalfedge->vertex()->getSemantic();
			
			Halfedge_handle prec = border1;
			Halfedge_handle it = border1;
			
			int bc = 0;
			
			do{
				semanticDescr & borderDescr = it->vertex()->getSemantic();
				//double distToLoop = L2Dist(heDescr,borderDescr);
				double distToLoop = CGAL::squared_distance(pHalfedge->vertex()->point(),it->vertex()->point());
				distToLoop = sqrt(distToLoop);
				
				if(distToLoop<R)
				{
					m_treeStructure1.m_vertices.insert(pHalfedge->vertex());
					//pHalfedge->vertex()->color(1,1,0);
					if(!m_distToLoop.count(pHalfedge->vertex()))
					{
						m_distToLoop[pHalfedge->vertex()] = distToLoop;
					}
					else if (m_distToLoop[pHalfedge->vertex()] > distToLoop)
					{
						m_distToLoop[pHalfedge->vertex()] = distToLoop;
					}
				}
				
				//progress alongside border
				Halfedge_handle next;
				Halfedge_around_vertex_circulator hC = it->vertex()->vertex_begin();
				do
				{
					if(hC->opposite()->is_border_edge() && hC->opposite()!=it->opposite())
					{
						next = hC->opposite();
						break;
					}
					
				}
				while(++hC!=it->vertex()->vertex_begin());
				
				it = next;
				}
			while(it!=border1);
		//}
	}
	////////////////////////////////////////////////////////////////////
	center = border1->vertex()->point();
	for(auto pHalfedge = m_polyhedron2->halfedges_begin();pHalfedge!=m_polyhedron2->halfedges_end();++pHalfedge)
	{
 		Point3d p = pHalfedge->vertex()->point();
		//if( CGAL::squared_distance(center,p) < squared_euclidean_radius)
		//{	// point is candidate for semantic descriptor testing
			//pHalfedge->vertex()->color(0,1,0);
			semanticDescr & heDescr = pHalfedge->vertex()->getSemantic();
			
			Halfedge_handle it = border2;
			Halfedge_handle prec = border2;
			
			do{
				//it->vertex()->color(1,1,0);
				semanticDescr & borderDescr = it->vertex()->getSemantic();
				
				//double distToLoop = L2Dist(heDescr,borderDescr);
				double distToLoop = CGAL::squared_distance(pHalfedge->vertex()->point(),it->vertex()->point());
				distToLoop = sqrt(distToLoop);
				
				if(distToLoop<R)
				{
					m_treeStructure2.m_vertices.insert(pHalfedge->vertex());
					//pHalfedge->vertex()->color(1,1,0);
					if(!m_distToLoop.count(pHalfedge->vertex()))
					{
						m_distToLoop[pHalfedge->vertex()] = distToLoop;
					}
					else if (m_distToLoop[pHalfedge->vertex()] > distToLoop)
					{
						m_distToLoop[pHalfedge->vertex()] = distToLoop;
					}
				}
				//progress alongside border
				
				Halfedge_handle next;
				Halfedge_around_vertex_circulator hC = it->vertex()->vertex_begin();
				do
				{
					if(hC->opposite()->is_border_edge() && hC->opposite()!=it->opposite())
					{
						next = hC->opposite();
						break;
					}
					
				}
				while(++hC!=it->vertex()->vertex_begin());
				//it->vertex()->color(1,1,0);
				it = next;
			}
			while(it!=border2);
		//}
	}	
}

void softICPController::cluster(deformationNode* root)
{
	std::map<Vertex_handle,bool> isRep;
	for(unsigned r=0;r<root->m_childrenNodes.size();++r)
	{
		isRep[root->m_childrenNodes[r]->m_rep] = true;
	}
	
	for(auto v = root->m_vertices.begin(); v!=root->m_vertices.end();++v)
	{
		if(isRep[*v]){continue;} // dont cluster the rep
		// For each point in the set, find the closest rep
		double distMin = std::numeric_limits<double>::max();
		for( int clust = 0; clust < root->m_childrenNodes.size(); ++clust)
		{
 			Vertex_handle & rep = root->m_childrenNodes[clust]->m_rep;
			
			//double dist = computePhiDistance(rep,*v);
			//double dist = computePhiDistance(rep,*v,0.7,0.3,0.0);
			double dist = computePhiDistance(rep,*v,0.5,0.1,0.4);
			if(dist<distMin)
			{
				distMin = dist;
				root->m_cluster[*v] = clust;
			}
		}
		
	}
	// Just check that rep cluster has not been changed
}

void softICPController::subdivideLeaves(deformationNode* root)
{
	unsigned sizeChildren = root->m_childrenNodes.size();
	if(sizeChildren == 0) //This is a leaf noe
	{
		//std::cout<< "size of leaf : " << root->m_vertices.size() << std::endl;
	}
	else
	{
		for(unsigned l=0;l<sizeChildren;++l)
		{
			subdivideLeaves(root->m_childrenNodes[l]);
		}
	}
}

void softICPController::updateRep(deformationNode* root)
{
	std::vector<Vertex_handle> candidateRep;
	std::vector<double> distMin(root->m_childrenNodes.size(),std::numeric_limits<double>::max());
	
	candidateRep.resize(root->m_childrenNodes.size());
	
	for(unsigned c = 0; c < root->m_childrenNodes.size();++c)
	{
		double eMin = std::numeric_limits<double>::max();
		for(auto v = root->m_vertices.begin(); v!=root->m_vertices.end();++v)
		{
			if(root->m_cluster[*v] != c){continue;}
			// compute sum of squared distance to the rest of the cluster
			double e = 0.0;
			for(auto g = root->m_vertices.begin(); g!=root->m_vertices.end();++g)
			{
				if(root->m_cluster[*g] !=c){continue;}
				{
					//double dist = computePhiDistance((*v),(*g));
					//double dist = computePhiDistance(*v,*g,0.7,0.3,0.0);
					double dist = computePhiDistance(*v,*g,0.5,0.1,0.4);
					e += dist*dist;
				}
			}
			if(e < eMin)
			{
				root->m_childrenNodes[c]->m_rep = *v;
				eMin = e;
			}
		}
	}
}

void softICPController::kMedoidMain(deformationNode * root, int k, int clusterIter)
{
	srand(time(NULL));
	
	// Compute the first medoid 
	std::set<int> seq; // make sure we have unique random rep
	//std::cout << "\n\t\t\t Init kmedoid rep";
	//timer_tic();
	for(unsigned c=0; c<k; ++c)
	{
		root->m_childrenNodes.push_back(new deformationNode);
		root->m_childrenNodes.back()->m_parent = root;
		root->m_childrenNodes.back()->m_clusterID = c;
		
		int r = -1;
		bool cond = true;
		while(cond)
		{
			r = (rand() / (double)RAND_MAX ) * (root->m_vertices.size() -1); 
			if(seq.count(r)==0)
			{
				seq.insert(r);
				cond = false;
			}
		}
		auto it = root->m_vertices.begin();
		
		for(unsigned v=0;v<=r;++v){++it;}
		root->m_childrenNodes.back()->m_rep = *it;
		root->m_cluster[root->m_childrenNodes.back()->m_rep]=c;
	}
	//timer_toc();
	
	//std::cout << "\n\t\t\t The kmedoid loop";
	//timer_tic();
	// The kMedoidloop
	cluster(root);
	for(unsigned i=0;i<clusterIter;++i)
	{
		updateRep(root);
		cluster(root);
	}
	//timer_toc();
	
	//std::cout << "\n\t\t\t collect verts";
	// Finally put halfedges in the right patch (cluster)
	//timer_tic();
	for(auto it = root->m_cluster.begin(); it!=root->m_cluster.end(); ++it)
	{
		root->m_childrenNodes[it->second]->m_vertices.insert(it->first);
	}
	
	
}

void softICPController::computeMatrixDistance(deformationNode* root)
{
	int nbcluster = root->m_childrenNodes.size();

	root->m_distanceMatrix = myMatrix(nbcluster,nbcluster);
	myMatrix & m = root->m_distanceMatrix;
	
	for(unsigned i=0; i<nbcluster;++i)
	{
		Vertex_handle rep1 = root->m_childrenNodes[i]->m_rep;
		for(unsigned j=0;j<i;++j)
		{
			Vertex_handle rep2 = root->m_childrenNodes[j]->m_rep;
			//double dist = computePhiDistance(rep1,rep2,0.7,0.3,0.0);
			double dist = computePhiDistance(rep1,rep2,0.5,0.1,0.4);
			m[i][j] = dist;
			m[j][i] = dist;
		}
	}
}




void softICPController::snapRegions(double R, double elasticity, int itermax, int treeDepth)
{
	std::cout << "Total time : ";
	timer_tic();
	bool convergenceCriterion = false;
	bool order = false;
	int iter = 1;

	//std::cout << "omp_get_max_threads(): " << omp_get_max_threads() << std::endl;
	
	// Select snapping region and create patches on surface
	//std::cout << "buildTreeStructure : "; 
	//timer_tic();
	buildTreeStructure(treeDepth,R);
	//timer_toc();
	
	R = m_R;
	//std::cout << "R : " << m_R;
	
	//std::cout << "init closest patches";
	//timer_tic();
	storeNodeLevel(&m_treeStructure1,0,m_levelNodes1);
	storeNodeLevel(&m_treeStructure2,0,m_levelNodes2);
	computeClosest(&m_treeStructure1,4,m_levelNodes1,m_distanceMatrices1);
	computeClosest(&m_treeStructure1,4,m_levelNodes2,m_distanceMatrices2);
	initClosest(m_levelNodes1,m_distanceMatrices1,4);
	initClosest(m_levelNodes2,m_distanceMatrices2,4);
	//timer_toc();
	
	
	
	//std::cout << "convergence loop";
	//timer_tic();
	
	double delta = 0.0;
	while(!convergenceCriterion && !m_stop_for_debug)
	{
		// swap meshes at each iteration
		order = !order; 
		
		PolyhedronPtr ma;
		PolyhedronPtr mb;
		
		convergenceCriterion = (iter == itermax);
		
		
		
		// find correspondence between the 2 snapping regions
		computeSnappingRegionCorrespondence(order);
		//double oldDelta = delta;
		/*delta = computeDelta(order);
		if( delta > oldDelta && iter!=1)
		{
			m_stop_for_debug = true;
		}*/
		if(!m_stop_for_debug)
		{
			deformationNode * treeStructure;
			if(order)
			{
				treeStructure = &m_treeStructure1;
			}
			else
			{
				treeStructure = &m_treeStructure2;
				iter++;
			}
			
			// Compute a transformation for each point
			unsigned int nbP = 0;
			std::vector <pointTransformation> transf;
			//#pragma omp parallel for private(it,nbP,m_Phi) default(none) shared(transf)
			for(auto it=treeStructure->m_vertices.begin(); it!=treeStructure->m_vertices.end() && !m_stop_for_debug;++it)
			{
				Vertex_handle pVertex = *it;
				// Get neighborhood N around pVertex
				std::vector<Vertex_handle> N = getNeighborhood(pVertex,R,iter,itermax,elasticity,order);
				// Get corresponding neighborhood
				std::vector<Vertex_handle> phiN = getCorrespondingNeighborhood(N);
				// Compute and store transformation
				pointTransformation ti = computeTransformation(N,phiN);
				
				transf.push_back(ti);
				nbP++;
			}
			
			// Apply the transformation for each point, but scale it (iter/itermax)
			int i=0;
			//#pragma omp parallel for private(i,it) default(none) shared(transf)
			int sizeS = treeStructure->m_vertices.size();
			for(auto it=treeStructure->m_vertices.begin(); it!=treeStructure->m_vertices.end() && (i<nbP) ;++it)
			//for(i=0;i<sizeS;++i)
			{
				Vertex_handle pVertex = *it;
				applyTransformation(pVertex,transf[i],iter,itermax);
				i++;
			}
		}
		//if(m_stop_for_debug)
		//{
		//	break;
		//}
	}
	//timer_toc();
	//std::cout << "fixBorder";
	//timer_tic();
	fixBorder();
	for(unsigned s=0;s<4;++s){gaussianSmooth();}
 	timer_toc();	
	//std::cout << "#of iterations : " << iter << "/" << itermax << std::endl;
}
																																
void softICPController::applyTransformation(Vertex_handle p, pointTransformation & ti, int iter, int itermax)
{
	//double li = iter/(double)itermax;
	double li = 1.0;
	ti.T = li*ti.T;
	ti.Q.setAxisAngle(ti.Q.axis(),li*ti.Q.angle());
	
	Point3d pos = p->point();
	qglviewer::Vec V = ti.Q * qglviewer::Vec(pos.x(),pos.y(),pos.z()) + ti.T;
	p->point() = CGAL::ORIGIN + Vector(V[0],V[1],V[2]);
	
	double dist = sqrt(CGAL::squared_distance(p->point(),m_Phi[p]->point()));
}

vector< Vertex_handle > softICPController::getNeighborhood(Vertex_handle p, double R, unsigned int iter,unsigned itermax, double elasticity, bool order)
{
	double distToLoop = m_distToLoop[p];
	//distToLoop = std::min(distToLoop,R-distToLoop);
	// compute the size of the local neighborhood
	double radius = 0.0;
	if(distToLoop != 0.0)
	{
		double expo = (elasticity)/(distToLoop);
		//double expo = ((iter)*elasticity)/(distToLoop);
		radius = R * exp(-expo*expo);
		
	}
	
	//std::cout << radius << std::endl;
	
	deformationNode * node;
	deformationNode * containsP;
	if(order){node = &m_treeStructure1;}
	else{node = &m_treeStructure2;}
	
	int nbClust = node->m_childrenNodes.size();
	
	while(node->m_childrenNodes.size()!=0)
	{
		int clust = node->m_cluster[p];
		node = node->m_childrenNodes[clust];
	}
	
	//std::cout << iter << " " << radius << std::endl;
	
	std::vector<Vertex_handle> N;
 	containsP = node; // leaf node
	
	bool getout = false;
	
	while(!getout)
	{
		if(node->m_distClosest.size() == 0)
		{
			for(auto it = node->m_vertices.begin();it!=node->m_vertices.end();++it)
				{
					N.push_back(*it);
				}
			getout = true;
			break;
		}
		double maxDist = node->m_distClosest.back(); // the m furthest distance
		
		if( radius < maxDist ) // collect all the patches whose distance is smaller than radius
		{
			for(unsigned c=0;c<node->m_distClosest.size();++c)
			{
				if(node->m_distClosest[c]<radius) 
				{
					deformationNode * patch = node->m_closest[c];
					//std::cout << patch->m_vertices.size() << " " << patch->m_clusterID << std::endl;
					for(auto it = patch->m_vertices.begin();it!=patch->m_vertices.end();++it)
					{
						N.push_back(*it);
					}
				}
			}
			
			getout = true;
		}
		//else // go up one level in the hierarchy and repeat the process
		if ((radius >=maxDist || N.size() < 3) )
		{
			N.clear();
			if(node->m_clusterID == -1 ) // can't go up one level
			{
				for(auto it = node->m_vertices.begin();it!=node->m_vertices.end();++it)
				{
					N.push_back(*it);
				}
				getout = true;
			}
			else{
				node = node->m_parent;
				getout = false;
			}
		}
	}
	return N;
}



vector< Vertex_handle > softICPController::getNeighborhoodOld(Vertex_handle p, double R, unsigned int iter,unsigned itermax, double elasticity, bool order)
{
	double distToLoop = m_distToLoop[p];
	
	// compute the size of the local neighborhood
	double radius = 0.0;
	//distToLoop = std::min(distToLoop,R-distToLoop);
	if(distToLoop != 0.0)
	{
		double expo = iter*elasticity/(distToLoop);
		radius = R * exp(-expo*expo);
		//std::cout << distToLoop << " " << expo << " " << expo*expo << " " << radius << std::endl;
	}
	
	deformationNode * node;
	deformationNode * containsP;
	if(order) {node = &m_treeStructure1;}
	else{node = &m_treeStructure2;}
	
	int nbClust = node->m_childrenNodes.size();
	
	while(node->m_childrenNodes.size()!=0)
	{
		int clust = node->m_cluster[p];
		node = node->m_childrenNodes[clust];
	}
	
	std::vector<Vertex_handle> N;
 	containsP = node; // leaf node
	node = node->m_parent; // parent node, to access matrix
	
	std::cout << "radius : " << radius << std::endl;
	
	bool getout = false;
	
	while(!getout)
	{
		// compute the maximum distance to containsP
		double maxDist = 0.0;
		for(unsigned k=0;k<nbClust;++k)
		{
// 			//double dist = node->m_distanceMatrix[k][containsP->m_clusterID];
			Vertex_handle rep = node->m_childrenNodes[k]->m_rep;
			double dist = computePhiDistance(p,rep,0.7,0.3,0.0);
			
			if(dist>maxDist)
			{
				maxDist = dist;
			}
		}
		
		if( radius < maxDist ) // collect all the patches whose distance is smaller than radius
		{
			for(unsigned k=0;k<nbClust;++k)
			{
				//double dist = node->m_distanceMatrix[k][containsP->m_clusterID];
				Vertex_handle rep = node->m_childrenNodes[k]->m_rep;
				double dist = computePhiDistance(p,rep,0.7,0.3,0.0);
				
				if(dist < radius)
				{
					for(auto it = node->m_childrenNodes[k]->m_vertices.begin();it!=node->m_childrenNodes[k]->m_vertices.end();++it)
					{
						N.push_back(*it);
					}
				}
			}
			
			getout = true;
		}
		//else // go up one level in the hierarchy and repeat the process
		if ((radius >=maxDist || N.size() < 3) )
		{
			N.clear();
			if(node->m_clusterID == -1 )
			{
				for(auto it = node->m_vertices.begin();it!=node->m_vertices.end();++it)
				{
					N.push_back(*it);
				}
				getout = true;
			}
			else{
				node = node->m_parent;
				getout = false;
			}
		}
		
		/*if(getout)
		{
			// Check for intermediate node : 
			// NOT root AND NOT leaf
			//if(node->m_clusterID !=-1 && node->m_childrenNodes.size()!=0)
			if(radius !=0 && radius < 0.001)
			{
				m_stop_for_debug = true;
				for(unsigned v = 0; v < N.size(); ++v)
				{
					N[v]->color(0,1,1);
					m_Phi[N[v]]->color(0,0,0);
					
				}
				p->color(1,0,0);
				m_Phi[p]->color(0,1,0);
				std::cout << "stop at iteration #"<< iter << std::endl;
			}
		}*/
	}
	
	/*bool isIn = false;
	for(unsigned v = 0; v < N.size(); ++v)
	{
		if(N[v] == p)
		{
			isIn = true;
		}
	}*/
	
	/*if(N.size() > 100)
	{
		for(unsigned v = 0; v < N.size(); ++v)
		{
			N[v]->color(0,1,1);
			m_Phi[N[v]]->color(0,0,0);
		}
	}*/
	
	
	/*for(unsigned v = 0; v < N.size(); ++v)
	{
		N[v]->color(0,1,1);
		m_Phi[N[v]]->color(0,0,0);
		p->color(0,0,0);
		m_Phi[p]->color(1,0,0);
	}*/
	
	return N;
}

std::vector<Vertex_handle> softICPController::getCorrespondingNeighborhood( std::vector<Vertex_handle> & N)
{
	std::vector<Vertex_handle> cN;
	for(auto it = N.begin(); it!= N.end(); ++it)
	{
		cN.push_back(m_Phi[*it]);
	}
	return cN;
}

void softICPController::fixBorder()
{
	this->computeSnappingRegionCorrespondence(true);
	for(auto it = m_treeStructure1.m_vertices.begin();
	    it!= m_treeStructure1.m_vertices.end();
		++it)
	    {
		if(m_distToLoop[*it]==0.0)
		{
			Point3d p = m_Phi[*it]->point();
			(*it)->point() = p;
		}
	    }
	this->computeSnappingRegionCorrespondence(false);
	for(auto it = m_treeStructure2.m_vertices.begin();
	    it!= m_treeStructure2.m_vertices.end();
		++it)
	    {
		if(m_distToLoop[*it]==0.0)
		{
			Point3d p = m_Phi[*it]->point();
			(*it)->point() = p;
		}
	    }
}

void softICPController::gaussianSmooth()
{
	std::vector<Point3d> newP;
	for(auto it = m_treeStructure1.m_vertices.begin();
	    it!= m_treeStructure1.m_vertices.end();
		++it)
	{
		Vertex_handle pVertex = *it;
		Halfedge_around_vertex_circulator hC = pVertex->vertex_begin();
		double sx = 0.0;
		double sy = 0.0;
		double sz = 0.0;
		do
		{
			Point3d n = hC->opposite()->vertex()->point();
			sx = sx + n.x();
			sy = sy + n.y();
			sz = sz + n.z();
		}while(++hC!=pVertex->vertex_begin());
		sx/=pVertex->vertex_degree();
		sy/=pVertex->vertex_degree();
		sz/=pVertex->vertex_degree();
		Point3d smoothed(sx,sy,sz);
		newP.push_back(smoothed);
	}
	int i=0;
	for(auto it = m_treeStructure1.m_vertices.begin();
	    it!= m_treeStructure1.m_vertices.end();
		++it,++i)
	{
		Vertex_handle pVertex = *it;
		pVertex->point() = newP[i];
		//pVertex->color(1.0,1.0,0.0);
	}
	newP.clear();
	for(auto it = m_treeStructure2.m_vertices.begin();
	    it!= m_treeStructure2.m_vertices.end();
		++it)
	{
		Vertex_handle pVertex = *it;
		Halfedge_around_vertex_circulator hC = pVertex->vertex_begin();
		double sx = 0.0;
		double sy = 0.0;
		double sz = 0.0;
		do
		{
			Point3d n = hC->opposite()->vertex()->point();
			sx = sx + n.x();
			sy = sy + n.y();
			sz = sz + n.z();
		}while(++hC!=pVertex->vertex_begin());
		sx/=pVertex->vertex_degree();
		sy/=pVertex->vertex_degree();
		sz/=pVertex->vertex_degree();
		Point3d smoothed(sx,sy,sz);
		newP.push_back(smoothed);
	}
	i=0;
	for(auto it = m_treeStructure2.m_vertices.begin();
	    it!= m_treeStructure2.m_vertices.end();
		++it,++i)
	{
		Vertex_handle pVertex = *it;
		pVertex->point() = newP[i];
		//pVertex->color(1.0,1.0,0.0);
	}
}

pointTransformation softICPController::computeTransformation(std::vector<Vertex_handle> & N, std::vector<Vertex_handle> & phiN)
{
	pointTransformation pi;
	
	// Compute Rotation with svd minimization
	Matrix p_m(N.size(),3);
	Matrix p_t(N.size(),3);
	
	// Compute Translation with center-of-mass alignment
	Matrix mu_m(1,3);
	Matrix mu_t(1,3);
	double mum0 = 0.0, mum1 = 0.0, mum2 = 0.0;
	double mut0 = 0.0, mut1 = 0.0, mut2 = 0.0;
	
	auto itt=phiN.begin();
	int i=0;
	int sizeN = N.size();
	
	//#pragma omp parallel for private(i,itt) default(none) shared(p_m,p_t,N,phiN)
	//for(auto it=N.begin();it!=N.end();++it)
	
	//std::cout << omp_get_num_threads() << std::endl;
	
	//#pragma omp parallel for private(i) default(none) shared(p_m,p_t,N,phiN,sizeN) reduction(+:mum0,mum1,mum2, mut0,mut1,mut2)
	for(i=0;i<sizeN;++i)
	{	
		Point3d np = N[i]->point();
		p_t.val[i][0] = np.x(); mut0 +=p_t.val[i][0];
		p_t.val[i][1] = np.y(); mut1 +=p_t.val[i][1];
		p_t.val[i][2] = np.z(); mut2 +=p_t.val[i][2];
		
		// get nearest point
		Point3d npp = phiN[i]->point();
		p_m.val[i][0] = npp.x(); mum0 += p_m.val[i][0];
		p_m.val[i][1] = npp.y(); mum1 += p_m.val[i][1];
		p_m.val[i][2] = npp.z(); mum2 += p_m.val[i][2];
	}
	
	mu_m.val[0][0] = mum0;
	mu_m.val[0][1] = mum1;
	mu_m.val[0][2] = mum2;
	
	mu_t.val[0][0] = mut0;
	mu_t.val[0][1] = mut1;
	mu_t.val[0][2] = mut2;
	
	mu_m = mu_m/(double)N.size();
	mu_t = mu_t/(double)N.size();
	
	Matrix q_m = p_m - Matrix::ones(N.size(),1)*mu_m;
	Matrix q_t = p_t - Matrix::ones(N.size(),1)*mu_t;
	
	// compute rotation matrix R and translation vector t
	Matrix H = ~q_t*q_m;
	Matrix U,W,V;
	H.svd(U,W,V);
	Matrix R = V*~U;
	
	if(R.det()<0)
	{
		Matrix B = Matrix::eye(3);
		B.val[2][2] = R.det();
		R = V*B*~U;
	}
	
	Matrix t = ~mu_m - R*~mu_t;
	
	double rData[3][3];
	for(unsigned i=0;i<3;++i)
	{
		for(unsigned j=0;j<3;++j)
		{
			rData[i][j] = R.val[i][j];
		}
	}
	pi.Q.setFromRotationMatrix(rData);
	pi.T.setValue(t.val[0][0],t.val[1][0],t.val[2][0]);
	
	return pi;
}

PolyhedronPtr softICPController::remesh()
{
	PolyhedronPtr sp(new Polyhedron);
	
	std::vector<double> coords;
	std::vector<int> tris;
	std::cout << "remesh :"<<std::endl;
	//getMeshOutsideSR(coords,tris,m_polyhedron1,0);
	//getMeshOutsideSR(coords,tris,m_polyhedron2,coords.size()/3);
	//std::cout << "getMeshOutsideSR ok :"<<std::endl;
	
	getMeshOutsideSR(coords,tris,m_polyhedron1,0);
	getMeshOutsideSR(coords,tris,m_polyhedron2,coords.size()/3);
	
	
	//buildSRMesh(coords,tris,coords.size()/3);
	std::cout << "buildSRMesh ok :"<<std::endl;
	polyhedron_builder<HalfedgeDS> builder(coords,tris);
	
	std::cout << "call builder :"<<std::endl;
	sp->delegate(builder);
	sp->compute_normals();
	
	std::set<Point3d> border;
	sp->normalize_border();
	for(auto h = sp->border_halfedges_begin(); h!=sp->halfedges_end();++h)
	{
		border.insert(h->vertex()->point());
	}
	
	coords.clear(); tris.clear();
	
	buildSRMesh(coords,tris,0,border);
	PolyhedronPtr poly(new Polyhedron);
	polyhedron_builder<HalfedgeDS> finalbuilder(coords,tris);
	poly->delegate(finalbuilder);
	//std::cout << sp->size_of_vertices() << std::endl;
//	isolatedVertices_remover<HalfedgeDS> verticeRemover;
//	sp->delegate(verticeRemover);
//	std::cout << sp->size_of_vertices() << std::endl;
	
	poly->compute_normals();
	//sp->compute_normals();

	return poly;
	//return sp;
}

void softICPController::getMeshOutsideSR(vector< double >& coords, vector< int >& tris, PolyhedronPtr p, int vertexOffset)
{
	p->set_index_vertices();
	for(auto pVertex = p->vertices_begin(); pVertex!=p->vertices_end(); ++pVertex)
	{
			Point3d p = pVertex->point();
			coords.push_back(p.x());coords.push_back(p.y());coords.push_back(p.z());
	}
	
	for(auto pFacet = p->facets_begin(); pFacet!=p->facets_end(); ++pFacet)
	{	
		int inside = 0;
		Halfedge_around_facet_circulator hC = pFacet->facet_begin();
		int vertsIds[3];
		int vId=0;
		do
		{
			Vertex_handle v = hC->vertex();
			if(m_isSnappingRegion[v]){inside++;}
			vertsIds[vId] = v->tag()+vertexOffset;
			vId++;
		}
		while(++hC!=pFacet->facet_begin());
		if(inside<3)
		{
			tris.push_back(vertsIds[0]);tris.push_back(vertsIds[1]);tris.push_back(vertsIds[2]);
		}
	}
}

void softICPController::getMeshInsideSR(vector< double >& coords, vector< int >& tris, PolyhedronPtr p, int vertexOffset)
{
	p->set_index_vertices();
	for(auto pVertex = p->vertices_begin(); pVertex!=p->vertices_end(); ++pVertex)
	{
			Point3d p = pVertex->point();
			coords.push_back(p.x());coords.push_back(p.y());coords.push_back(p.z());
	}
	
	for(auto pFacet = p->facets_begin(); pFacet!=p->facets_end(); ++pFacet)
	{	
		int inside = 0;
		Halfedge_around_facet_circulator hC = pFacet->facet_begin();
		int vertsIds[3];
		int vId=0;
		do
		{
			Vertex_handle v = hC->vertex();
			if(m_isSnappingRegion[v]){inside++;}
			vertsIds[vId] = v->tag()+vertexOffset;
			vId++;
		}
		while(++hC!=pFacet->facet_begin());
		if(inside!=0)
		{
			tris.push_back(vertsIds[0]);tris.push_back(vertsIds[1]);tris.push_back(vertsIds[2]);
		}
	}
}


PolyhedronPtr softICPController::buildSRMesh(vector< double >& coords, vector< int >& tris, int vertexOffset, std::set<Point3d> & border)
{
	// construction from a list of points
	std::list<Point> L;
	for(auto p = m_treeStructure1.m_vertices.begin(); p!= m_treeStructure1.m_vertices.end(); ++p)
	{
		Point3d pp = (*p)->point();
		//L.push_front(std::make_pair<Point,unsigned>(Point(pp.x(),pp.y(),pp.z()),(*p)->tag()));
		L.push_front(Point(pp.x(),pp.y(),pp.z()));
	}
	
	for(auto p = m_treeStructure2.m_vertices.begin(); p!= m_treeStructure2.m_vertices.end(); ++p)
	{
		Point3d pp = (*p)->point();
		//L.push_front(std::make_pair<Point,unsigned>(Point(pp.x(),pp.y(),pp.z()),(*p)->tag()));
		L.push_front(Point(pp.x(),pp.y(),pp.z()));
	}
	
	Triangulation T(L.begin(),L.end());
	if(T.is_valid()){ std::cout << "triangulation is valid" << std::endl;}
	else{ std::cout << "triangulation is not valid" << std::endl;}
	
	std::vector<RFacet> facets;
	
	AABB_Tree tree(m_polyhedron1->facets_begin(),m_polyhedron1->facets_end());
	tree.insert(m_polyhedron2->facets_begin(),m_polyhedron2->facets_end());
	tree.build();
	
	Perimeter perimeter(0.0);
	DistanceSurface distSurf(tree,border);
	
	CGAL::advancing_front_surface_reconstruction(L.begin(),L.end(),std::back_inserter(facets),distSurf);
	
	for(auto v = L.begin(); v!=L.end();++v)
	{
		auto p = (*v);
		coords.push_back(p.x());coords.push_back(p.y());coords.push_back(p.z());
	}
	
	for(auto f = facets.begin(); f!=facets.end(); ++f)
	{
		//tris.push_back(*f[0]);tris.push_back(*f[1]);tris.push_back(*f[2]);
		tris.push_back((*f)[0]);tris.push_back((*f)[1]);tris.push_back((*f)[2]);
	}
}

/*void softICPController::findCouples(PolyhedronPtr meshA, PolyhedronPtr meshB)
{

	Facet_iterator pFacet =	NULL;
	std::list<AABB_Tree::Primitive_id> primitives;
	std::list<Triangle> triangles;

	unsigned hEID = 0;
	unsigned fID = 0;
	
	AABB_Tree tree;
	
	//The AABB-tree is built on the polyhedron with the less number of facets
	if(meshA->size_of_facets() < meshB->size_of_facets())
	{
		for(pFacet = meshA->facets_begin(); pFacet != meshA->facets_end(); pFacet++) 
		{triangles.push_back(Triangle(pFacet));}
		tree.rebuild(triangles.begin(),triangles.end());
		
		//collision test with each facet of the second polyhedron
		for (pFacet = meshB->facets_begin(); pFacet != meshB->facets_end(); pFacet++)
		{
			tree.all_intersected_primitives(Triangle(pFacet), std::back_inserter(primitives));
			if(primitives.size() !=0)
			{
				
			}
		}
	}
	
	else
	{
		for(pFacet = meshB->facets_begin(); pFacet != meshB->facets_end(); pFacet++)
		{triangles.push_back(Triangle(pFacet));}
		tree.rebuild(triangles.begin(),triangles.end());
		
		//collision test with each facet of the first polyhedron
		for(pFacet = meshA->facets_begin(); pFacet != meshA->facets_end(); pFacet++)
		{
			tree.all_intersected_primitives(Triangle(pFacet), std::back_inserter(primitives));
			if(primitives.size() !=0)
			{
				
			}
		}
	}
}

void softICPController::remeshSR(std::vector<double> & coords, std::vector<int> & tris, int vertexOffset)
{
	PolyhedronPtr subMeshA(new Polyhedron);
	PolyhedronPtr subMeshB(new Polyhedron);
	
	std::vector<double> coords;
	std::vector<int> tris;
	
	getMeshInsideSR(coords,tris,m_polyhedron1,0);
	polyhedron_builder<HalfedgeDS> builder(coords,tris);
	subMeshA->delegate(builder);
	subMeshA->compute_normals();
	coords.clear();tris.clear;
	getMeshInsideSR(coords,tris,m_polyhedron2,0);
	subMeshB->delegate(builder);
	subMeshB->compute_normals();
	
}

void softICPController::floodConstruct(PolyhedronPtr mesh1, PolyhedronPtr mesh2)
{
	Facet_iterator pFacet =	NULL;
	std::list<AABB_Tree::Primitive_id> primitives;
	std::list<Triangle> triangles;

	std::list<Facet_handle> solution;
	
	std::map<Facet_handle, bool> isVisited;
	
	AABB_Tree tree;
	
	unsigned heID = 0;
	unsigned faceID = 0;

	if(mesh1->size_of_facets() > mesh2->size_of_facets())
	{
		//Build the AABB-tree on the other surface
		//for(pFacet = mesh2->facets_begin(); pFacet != mesh2->facets_end(); pFacet++) {triangles.push_back(Triangle(pFacet));}
		//tree.rebuild(triangles.begin(),triangles.end());
		
		//Get all the facets adjacent to the border
		std::list<Facet_handle> front1;
		std::list<Facet_handle> front2;
		
		
		for(auto it = m_loop1.begin(); it!=m_loop1.end();++it)
		{	
			solution.push_back((*it)->opposite()->face());
			front1.push_back((*it)->opposite()->face());
			isVisited[front1.back()] = true;
		}
		
		for(auto it = m_loop2.begin(); it!=m_loop2.end();++it)
		{
			solution.push_back((*it)->opposite()->face());
			front2.push_back((*it)->opposite()->face());
			isVisited[front2.back()] = true;
		}
		
		for(auto pFacet = front1.begin(); pFacet != front1.end(); ++pFacet)
		{
			triangles.push_back(Triangle(pFacet));
		}
		tree.rebuild(triangles.begin(),triangles.end());
		
		bool stopCondition = false;
		
		while(!stopCondition)
		{
			int front1Size = front1.size();
			int front2Size = front2.size();
			for(unsigned i = 0; i< front1Size();++i)
			{
				Facet_handle f = front1.front();
				Halfedge_around_facet_circulator h = f->facet_begin();
				do
				{
					Facet_handle nF = h->opposite()->face();
					if(!isVisited[nF])
					{
						front1.push_back(nF);
					}
				}
				while(++h!=f->facet_begin());
				front1.pop_front();
			}// front1 has been updated
			
			for(unsigned i = 0; i< front2Size();++i)
			{
				Facet_handle f = front2.front();
				Halfedge_around_facet_circulator h = f->facet_begin();
				do
				{
					Facet_handle nF = h->opposite()->face();
					if(!isVisited[nF])
					{
						front2.push_back(nF);
					}
				}
				while(++h!=f->facet_begin());
				front2.pop_front();
			}// front2 has been updated
			
			// Check distance // intersection betwen front1 and 2
			for(unsigned i=0;i<front1.size();++i)
			{
				triangles.push_back(Triangle(front1[i]));
			}
			tree.rebuild();
			for(auto it = front2.begin(); it!= front2.end(); ++it)
			{
				if(tree.do_intersect(*it))
				{
					front2.erase(it);
				}
			}
		}
		
		
		//Flood while there is no intersection
		
		while(!front1.empty())
		{
			Facet_handle f = front1.front();
			tree.all_intersected_primitives(Triangle(f),std::back_inserter(primitives));
			if(primitives.size()!=0)
			{
				m_inter_tri.push_back(Triangle_Cut(Compute_Normal_direction(pFacet->facet_begin()), false));
				do
				{
					
					front2.push_back(primitives.back()->facet());
					isVisited[primitives.back()];
					m_inter_tri.push_back(Triangle_Cut(Compute_Normal_direction(primitives.back()->facet()->facet_begin()), true));
					m_Couples[primitives.back()->facet()].insert(pFacet);
					primitives.pop_back();
				}while(primitives.size()!=0);
				
				for(unsigned i=0;i<front1.size();++i)
				{
					front2.push_back(tree.closest_point_and_primitive(front1[i]));
				}
				
			}
			else
			{
				solution.push_back(f);
				Halfedge_around_facet_circulator h = f->facet_begin();
				do
				{
					Facet_handle nF = h->opposite()->face();
					if(!isVisited[nF])
					{
						front1.push_back(nF);
						isVisited[nF] = true;
					}
				}
				while(++h!=f->facet_begin());
			}
			front1.pop_front();
		}
		
		
		
		
	
	}
}

void softICPController::computeIntersections()
{
	while(!m_Couples.empty())
	{
		Facet_handle fA, fB;
		fA = m_Couples.begin()->first;
		fB = *m_Couples[fA].begin();
		interTriangleTriangle(fA, fB);
		rmCouple(fA, fB);
	}
}*/

/*void softICPController::cutIntersectedFacets(PolyhedronPtr meshA, PolyhedronPtr meshB)
{
	Triangle_Cut TriCut;
	Halfedge_handle he;

	//every intersected facet is triangulated if at least one of the intersections is a segment
	for(FacetId Facet = 0 ; Facet != m_inter_tri.size() ; ++Facet)
	{
		if(!m_inter_tri[Facet].CutList.empty())
		{
			TriCut = m_inter_tri[Facet];
			he = Facet_Handle[Facet]->facet_begin();
			bool IsExt[3];
			//creation of a triangulation
			Triangulation<Exact_Kernel> T(he, TriCut.norm_dir);
			//add the list of intersection points (only happens in case of intersection of two edges)
			for(std::set<InterId>::iterator i = TriCut.PtList.begin();i != TriCut.PtList.end();++i)
			{
				T.add_new_pt(InterPts[*i], (unsigned long &)*i);    // MT: ajout cast
			}
			//add the intersection segments
			for(int i = 0;i!=(int)TriCut.CutList.size();++i)
			{
				T.add_segment(InterPts[TriCut.CutList[i][0]], InterPts[TriCut.CutList[i][1]], TriCut.CutList[i][0], TriCut.CutList[i][1]);
			}
			//get the triangles of the triangulation thay belong to the result
			//and determine if the three neighboring facets belongs to the result (using IsExt[3])
			vector<vector<unsigned long> > Tri_set = T.get_triangles((m_BOOP == MINUS && !TriCut.Facet_from_A)?true:false, IsExt);
			//add these triangles to the result
			ppbuilder.add_triangle(Tri_set, he);

			//update the tags
			Facet_Handle[Facet]->IsOK = true;
			if(IsExt[0]) he->opposite()->facet()->IsExt = true;
			if(IsExt[1]) he->next()->opposite()->facet()->IsExt = true;
			if(IsExt[2]) he->next()->next()->opposite()->facet()->IsExt = true;
		}
	}
}*/

/*void softICPController::rmCouple(Facet_handle A, Facet_handle B)
{
	if(m_Couples[A].count(B) != 0) m_Couples[A].erase(B);
	if(m_Couples[A].empty()) m_Couples.erase(A);
}*/

double computePhiDistance(Vertex_handle v1, Vertex_handle v2, double w1, double w2, double w3)
{
	double dist = w1 * sqrt(CGAL::squared_distance(v1->point(),v2->point()));
	dist += w2*acos(v1->normal()*v2->normal());
	dist += w3*L2Dist(v1->getSemantic(),v2->getSemantic());
	return dist;
}

	