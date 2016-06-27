#include "softICPController.h"
#include <Tools/Tools_Polyhedron.h>
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

#include <CGAL/IO/output_surface_facets_to_polyhedron.h>


#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#include <CGAL/Subdivision_method_3.h>

#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

#include <CGAL/bounding_box.h>

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits.h>


#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double> AABB_Kernel;
typedef CGAL::AABB_polyhedron_triangle_primitive<AABB_Kernel,Polyhedron> AABB_Primitive;
typedef CGAL::AABB_traits<AABB_Kernel, AABB_Primitive> AABB_Traits;
typedef CGAL::AABB_tree<AABB_Traits> AABB_Tree;
typedef AABB_Tree::Object_and_primitive_id Object_and_primitive_id;
typedef AABB_Tree::Point_and_primitive_id Point_and_primitive_id;
typedef AABB_Kernel::Triangle_3 Triangle;

typedef CGAL::Exact_predicates_inexact_constructions_kernel delaunayKernel;
typedef CGAL::Triangulation_3<AABB_Kernel>      Triangulation;
typedef Triangulation::Cell_handle    Cell_handle;
typedef Triangulation::Vertex_handle  VoronoiVertex_handle;
typedef Triangulation::Locate_type    Locate_type;
typedef Triangulation::Point          Point;


typedef CGAL::cpp11::array<std::size_t,3> RFacet;

typedef CGAL::Simple_cartesian<double> cartesianK;
typedef CGAL::Search_traits_3<cartesianK> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree sTree;

namespace PMP = CGAL::Polygon_mesh_processing;

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
	double m_r;
	DistanceSurface(AABB_Tree & tree, std::set<Point3d> & borders, double r) : m_tree(tree), m_borders(borders), m_r(r)
	{
		std::cout << "m_tree info : " << m_tree.size() << std::endl;
	}
	
	template <typename Point>
	bool operator()(const Point& p, const Point& q, const Point& r) const
	{	
		Point3d pp(p.x(),p.y(),p.z());
		Point3d qq(q.x(),q.y(),q.z());
		Point3d rr(r.x(),r.y(),r.z());
		Triangle t(pp,qq,rr);
		
		bool border = m_borders.count(pp) && m_borders.count(qq) && m_borders.count(rr); 
		int borderC = m_borders.count(pp) +  m_borders.count(qq) + m_borders.count(rr); 
		
		double pq = sqrt(CGAL::squared_distance(p,q));
		double qr = sqrt(CGAL::squared_distance(q,r));
		double rp = sqrt(CGAL::squared_distance(r,p));
		
		/*double ratio = false;
		if( pq/qr > 3.0 || pq/qr < 0.33){ratio = true;}
		if( pq/rp > 3.0 || pq/rp < 0.33){ratio = true;}
		if( qr/rp > 3.0 || qr/rp < 0.33){ratio = true;}
		
		if(ratio && borderC > 0) {return true;}*/
		
		//bool tree = m_tree.do_intersect(t);
		
		/*if(pq > m_r/5.0){return true;}
		if(qr > m_r/5.0){return true;}
		if(rp > m_r/5.0){return true;}*/
		
		//std::cout <<  "do intersect ? " << tree << "\n";
		
		//return border || ratio || tree;
		//return (borderC>2);
		return border;
		//return false;
	}
	
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel delaunayKernel;
typedef CGAL::Advancing_front_surface_reconstruction<CGAL::Default, DistanceSurface> Reconstruction;
typedef Reconstruction::Triangulation_3 RTriangulation_3;
typedef Reconstruction::Outlier_range Outlier_range;
typedef Reconstruction::Boundary_range Boundary_range;
typedef Reconstruction::Vertex_on_boundary_range Vertex_on_boundary_range;
typedef Reconstruction::Triangulation_data_structure_2 TDS_2;

Halfedge_handle create_center_vertex( PolyhedronPtr p, Facet_iterator f) {
    Vector vec( 0.0, 0.0, 0.0);
    std::size_t order = 0;
    auto h = f->facet_begin();
    
    do {
        vec = vec + ( h->vertex()->point() - CGAL::ORIGIN);
        ++ order;
    } while ( ++h != f->facet_begin());
    CGAL_assertion( order >= 3); // guaranteed by definition of polyhedron
    Point3d center =  CGAL::ORIGIN + (vec / static_cast<double>(order));
    Halfedge_handle new_center = p->create_center_vertex( f->halfedge());
    new_center->vertex()->point() = center;
    return new_center;
}




softICPController::softICPController(PolyhedronPtr m1, PolyhedronPtr m2, Viewer * v)	 : m_polyhedron1(m1), m_polyhedron2(m2), m_v(v), m_stop_for_debug(false)
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
		if(!m_isSnappingRegion[*it]){continue;}
		double distMin = std::numeric_limits<double>::max();
		Vertex_handle bestCorres;
		for(auto c = s2->begin(); c!=s2->end(); ++c)
		{
			if(!m_isSnappingRegion[*c]){continue;}
			//double dist = computePhiDistance((*it),(*c));
			//double dist = computePhiDistance((*it),(*c),1.0,0.0,0.0);
			double dist = computePhiDistance((*it),(*c),0.7,0.3,0.0);
			//double dist = computePhiDistance((*it),(*c),0.5,0.1,0.4);
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
				//double dist = computePhiDistance(nodeP->m_rep,nodeQ->m_rep,0.5,0.1,0.4);
				double dist = computePhiDistance(nodeP->m_rep,nodeQ->m_rep,1.0,0.0,0.0);
				mat->val[p][q] = dist;
				mat->val[q][p] = dist;
			}
		}
	}
}

void softICPController::computeClosestGeo( deformationNode * root, int m, std::vector< std::vector<deformationNode*> > & levelNodes, std::vector < Matrix * > & distanceMatrices, geodesic::Mesh & g)
{	
	geodesic::GeodesicAlgorithmDijkstraAlternative * geoAlg = new geodesic::GeodesicAlgorithmDijkstraAlternative(&g);
	
	// for l == 0 
	distanceMatrices.push_back(new Matrix(1,1));
	
	for(unsigned l=1; l<levelNodes.size();++l)
	{
		int sizeLevel = levelNodes[l].size();
		distanceMatrices.push_back(new Matrix(sizeLevel,sizeLevel));
		Matrix * mat = distanceMatrices.back();
		
		for(unsigned p=0;p<sizeLevel;++p)
		{
			int index = levelNodes[l][p]->m_rep->tag();
			
			std::vector<geodesic::SurfacePoint> sources;
			geodesic::SurfacePoint gp(&g.vertices()[index]);
			sources.push_back(gp);
		
			geoAlg->propagate(sources);
			
			for(unsigned q=0;q<p;++q)
			{
				int idx = levelNodes[l][q]->m_rep->tag();
				double dist;
				geodesic::SurfacePoint gq(&g.vertices()[idx]);
				geoAlg->best_source(gq,dist);
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
	getSnappingRegionAABB();
	m_polyhedron1->compute_normals();
	m_polyhedron2->compute_normals();
	
	hierarchicalBuild(&m_treeStructure1,sizeOfTree,0,4);
	hierarchicalBuild(&m_treeStructure2,sizeOfTree,0,4);
	
	// Init. root nodes properly
	m_treeStructure1.m_parent = NULL;
	m_treeStructure2.m_parent = NULL;
	m_treeStructure1.m_clusterID = -1;
	m_treeStructure2.m_clusterID = -1;
}

void softICPController::buildFullTreeStructure(const int sizeOfTree, const double factor)
{
	
	//getSnappingRegion(factor);
	std::cout << "getSnappingRegion geodesic: " << std::endl;
	//timer_tic();
	getSnappingRegionGeo(factor);
	//timer_toc();
	
	m_polyhedron1->compute_normals();
	m_polyhedron2->compute_normals();
	
	//std::cout << "structure 1 : " << m_treeStructure1.m_vertices.size() << std::endl;
	//std::cout << "structure 2 : " << m_treeStructure2.m_vertices.size() << std::endl;
	
	std::cout << "Hierarchical build : " << std::endl;
	//timer_tic();
	//hierarchicalBuild(&m_treeStructure1,sizeOfTree,0,4);
	//hierarchicalBuild(&m_treeStructure2,sizeOfTree,0,4);
	//timer_toc();
	
	// Init. root nodes properly
	m_treeStructure1.m_parent = NULL;
	m_treeStructure2.m_parent = NULL;
	m_treeStructure1.m_clusterID = -1;
	m_treeStructure2.m_clusterID = -1;
	
	/*colorLastClusters(&m_treeStructure1);
 	colorLastClusters(&m_treeStructure2);*/
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
	
	//std::cout << "level : " << level << " nbVertices : " << root->m_vertices.size() << std::endl;
	
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
		std::cout << "size of leaf : " << root->m_childrenNodes.size() << std::endl;
		return;
	}
	double r = rand()/(double)RAND_MAX;
	double g = rand()/(double)RAND_MAX;
	double b = rand()/(double)RAND_MAX;
	for(auto h = root->m_childrenNodes[0]->m_vertices.begin(); 
	    h != root->m_childrenNodes[0]->m_vertices.end();
		++h)
	{
		
		(*h)->color(r,g,b);
	}
	r = rand()/(double)RAND_MAX;
	g = rand()/(double)RAND_MAX;
	b = rand()/(double)RAND_MAX;
	for(auto h = root->m_childrenNodes[1]->m_vertices.begin(); 
	    h != root->m_childrenNodes[1]->m_vertices.end();
		++h)
	{
		(*h)->color(r,g,b);
	}
	r = rand()/(double)RAND_MAX;
	g = rand()/(double)RAND_MAX;
	b = rand()/(double)RAND_MAX;
	for(auto h = root->m_childrenNodes[2]->m_vertices.begin(); 
	    h != root->m_childrenNodes[2]->m_vertices.end();
		++h)
	{
		(*h)->color(r,g,b);
	}
	r = rand()/(double)RAND_MAX;
	g = rand()/(double)RAND_MAX;
	b = rand()/(double)RAND_MAX;
	for(auto h = root->m_childrenNodes[3]->m_vertices.begin(); 
	    h != root->m_childrenNodes[3]->m_vertices.end();
		++h)
	{
		(*h)->color(r,g,b);
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
	
	const double double_max = numeric_limits<double>::max();
	
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
// 			//pVertex->color(1.0,1.0,0.0);
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
	
	Facet_iterator B1 = m_polyhedron1->facets_begin();
	Facet_iterator E1 = m_polyhedron1->facets_end(); --E1;
	
	Facet_iterator B2 = m_polyhedron2->facets_begin();
	Facet_iterator E2 = m_polyhedron2->facets_end(); --E2;

	Facet_iterator f = B1;
	do {
		int countSR = 0;
		auto hC = f->facet_begin();
		do
		{
			if(m_isSnappingRegion[hC->vertex()]){countSR++;}
		}
		while(++hC!=f->facet_begin());
		if(countSR == 3)
		{
			double distToLoop = double_max;
			Halfedge_handle h = create_center_vertex(m_polyhedron1,f);
			Vertex_handle p = h->vertex();
			m_isSnappingRegion[p] = true;
			m_treeStructure1.m_vertices.insert(p);
			for(auto b=loop1.begin();b!=loop1.end();++b)
			{
				double dist = CGAL::squared_distance(p->point(),(*b)->point());
				if(dist<distToLoop)
				{
					distToLoop = dist;
				}
			}
			distToLoop = sqrt(distToLoop);
			m_distToLoop[h->vertex()] = distToLoop;
		}
	} while ( f++ != E1);
	
	Facet_iterator g = B2;
	do {
		int countSR = 0;
		auto hC = g->facet_begin();
		do
		{
			if(m_isSnappingRegion[hC->vertex()]){countSR++;}
		}
		while(++hC!=g->facet_begin());
		if(countSR == 3)
		{
			double distToLoop = double_max;
			Halfedge_handle h = create_center_vertex(m_polyhedron2,g);
			Vertex_handle p = h->vertex();
			m_isSnappingRegion[p] = true;
			m_treeStructure2.m_vertices.insert(p);
			
			for(auto b=loop2.begin();b!=loop2.end();++b)
			{
				double dist = CGAL::squared_distance(p->point(),(*b)->point());
				if(dist<distToLoop)
				{
					distToLoop = dist;
				}
			}
			distToLoop = sqrt(distToLoop);
			m_distToLoop[p] = distToLoop;
		}
	} while ( g++ != E2);
}

void softICPController::getSnappingRegion(const double factor)
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


	const double double_max = numeric_limits<double>::max();

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

	double factoredR = m_R;

	if(factor>1.0)
	{
		factoredR*=factor;
	}

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
			pVertex->color(1.0,1.0,0.0);
			m_distToLoop[pVertex] = distToLoop;
			m_treeStructure1.m_vertices.insert(pVertex);
			m_isSnappingRegion[pVertex] = true;
			m_isFactoredRegion[pVertex] = true;
		}
		//else if(distToLoop < factoredR)
		else
		{
			pVertex->color(1.0,0.0,0.0);
			m_distToLoop[pVertex] = distToLoop;
			m_treeStructure1.m_vertices.insert(pVertex);
			m_isFactoredRegion[pVertex] = true;
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
			pVertex->color(1.0,1.0,0.0);
			m_distToLoop[pVertex] = distToLoop;
			m_treeStructure2.m_vertices.insert(pVertex);
			m_isSnappingRegion[pVertex] = true;
			m_isFactoredRegion[pVertex] = true;
		}
		//else if(distToLoop < factoredR)
		else
		{
			pVertex->color(1.0,0.0,0.0);
			m_distToLoop[pVertex] = distToLoop;
			m_treeStructure2.m_vertices.insert(pVertex);
			m_isFactoredRegion[pVertex] = true;
		}

	}

	Facet_iterator B1 = m_polyhedron1->facets_begin();
	Facet_iterator E1 = m_polyhedron1->facets_end(); --E1;

	Facet_iterator B2 = m_polyhedron2->facets_begin();
	Facet_iterator E2 = m_polyhedron2->facets_end(); --E2;

	Facet_iterator f = B1;
	do {
		int countFR = 0;
		int countSR = 0;
		auto hC = f->facet_begin();
		do
		{
			//if(m_isFactoredRegion[hC->vertex()]){countFR++;}
			if(m_isSnappingRegion[hC->vertex()]){countSR++;}
		}
		while(++hC!=f->facet_begin());
		if(countSR == 3)
		{
			double distToLoop = double_max;
			Halfedge_handle h = create_center_vertex(m_polyhedron1,f);
			Vertex_handle p = h->vertex();
			m_isFactoredRegion[p] = true;
			m_treeStructure1.m_vertices.insert(p);
			for(auto b=loop1.begin();b!=loop1.end();++b)
			{
				double dist = CGAL::squared_distance(p->point(),(*b)->point());
				if(dist<distToLoop)
				{
					distToLoop = dist;
				}
			}
			distToLoop = sqrt(distToLoop);
			m_distToLoop[h->vertex()] = distToLoop;
			if(countSR == 3){m_isSnappingRegion[p]=true;}
		}
	} while ( f++ != E1);

	Facet_iterator g = B2;
	do {
		int countFR = 0;
		int countSR = 0;
		auto hC = g->facet_begin();
		do
		{
			//if(m_isFactoredRegion[hC->vertex()]){countFR++;}
			if(m_isSnappingRegion[hC->vertex()]){countSR++;}
		}
		while(++hC!=g->facet_begin());
		if(countSR == 3)
		{
			double distToLoop = double_max;
			Halfedge_handle h = create_center_vertex(m_polyhedron2,g);
			Vertex_handle p = h->vertex();
			m_isSnappingRegion[p] = true;
			m_treeStructure2.m_vertices.insert(p);

			for(auto b=loop2.begin();b!=loop2.end();++b)
			{
				double dist = CGAL::squared_distance(p->point(),(*b)->point());
				if(dist<distToLoop)
				{
					distToLoop = dist;
				}
			}
			distToLoop = sqrt(distToLoop);
			m_distToLoop[p] = distToLoop;
			if(countSR == 3){m_isSnappingRegion[p]=true;}
		}
	} while ( g++ != E2);
}

void softICPController::getSnappingRegionGeo()
{
	//std::cout << "Before initGeodesicMesh" << std::endl;
	// initialize geodesic graph and algorithm
	initGeodesicMesh(m_polyhedron1,&m_g1);
	
	//std::cout << "After initGeodesicMesh" << std::endl;
	
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


	const double double_max = numeric_limits<double>::max();

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
	
	geodesic::GeodesicAlgorithmDijkstraAlternative * alg1 = new geodesic::GeodesicAlgorithmDijkstraAlternative(&m_g1);
	// compute for each point in the two sets the geodesic distance to the shape's b-loop
	std::vector<geodesic::SurfacePoint> sources1;
	std::vector<geodesic::SurfacePoint> sources2;
	
	for(unsigned i=0;i<loop1.size();++i)
	{
		int index = loop1[i]->tag();
		geodesic::SurfacePoint gp(&m_g1.vertices()[index]);
		sources1.push_back(gp);
	}
	alg1->propagate(sources1);
	// For every vertex, geodesic distance to the boundary loop
	for(auto pVertex = m_polyhedron1->vertices_begin(); pVertex!=m_polyhedron1->vertices_end(); ++pVertex)
	{
		geodesic::SurfacePoint gp(&m_g1.vertices()[pVertex->tag()]);
		double distance;
		alg1->best_source(gp,distance);
		m_distToLoop[pVertex] = distance;
		m_treeStructure1.m_vertices.insert(pVertex);
	}
	//delete alg1;
	
	initGeodesicMesh(m_polyhedron2,&m_g2);
	geodesic::GeodesicAlgorithmDijkstraAlternative * alg2 = new geodesic::GeodesicAlgorithmDijkstraAlternative(&m_g2);
	for(unsigned i=0;i<loop2.size();++i)
	{
		int index = loop2[i]->tag();
		geodesic::SurfacePoint gp(&m_g2.vertices()[index]);
		sources2.push_back(gp);
	}
	alg2->propagate(sources2);
	for(auto pVertex = m_polyhedron2->vertices_begin(); pVertex!=m_polyhedron2->vertices_end(); ++pVertex)
	{
		geodesic::SurfacePoint gp(&m_g2.vertices()[pVertex->tag()]);
		double distance;
		alg2->best_source(gp,distance);
		m_distToLoop[pVertex] = distance;
		m_treeStructure2.m_vertices.insert(pVertex);
	}
	//delete alg2;

	
	// m_R is the largest distance among the vertices in the snapping region
	double sR = 0.0;
	for(int i=0; i<closestToB2.size();i++)
	{
		double dist = m_distToLoop[closestToB2[i]];
		if(dist > sR) {sR = dist;}
	}
	for(int i=0; i<closestToB1.size();i++)
	{
		double dist = m_distToLoop[closestToB1[i]];
		if(dist < sR) {sR = dist;}
	}
	m_R = sqrt(sR);
	
	//The snapping region is the submesh within R geodesic distance to the b-loop
	for(auto pVertex = m_polyhedron1->vertices_begin();
		pVertex != m_polyhedron1->vertices_end();
		++pVertex)
	{
		double dist = m_distToLoop[pVertex];
		if(dist < m_R)
		{
			m_isSnappingRegion[pVertex]=true;
			m_sr1.push_back(pVertex);
			pVertex->color(1,0,0);
		}
	}
	for(auto pVertex = m_polyhedron2->vertices_begin();
		pVertex != m_polyhedron2->vertices_end();
		++pVertex)
	{
		double dist = m_distToLoop[pVertex];
		if(dist < m_R)
		{
			m_isSnappingRegion[pVertex]=true;
			m_sr2.push_back(pVertex);
			pVertex->color(1,0,0);
		}
	}
	computeSnapRegionDistances();
	
	delete alg1;
	delete alg2;
}

void softICPController::getSnappingRegionGeo(const double factor)
{
	initGeodesicMesh(m_polyhedron1,&m_g1);
	
	//std::cout << "After initGeodesicMesh" << std::endl;
	
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


	const double double_max = numeric_limits<double>::max();

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
	
	geodesic::GeodesicAlgorithmDijkstraAlternative * alg1 = new geodesic::GeodesicAlgorithmDijkstraAlternative(&m_g1);
	// compute for each point in the two sets the geodesic distance to the shape's b-loop
	std::vector<geodesic::SurfacePoint> sources1;
	std::vector<geodesic::SurfacePoint> sources2;
	
	for(unsigned i=0;i<loop1.size();++i)
	{
		int index = loop1[i]->tag();
		geodesic::SurfacePoint gp(&m_g1.vertices()[index]);
		sources1.push_back(gp);
	}
	alg1->propagate(sources1);
	// For every vertex, geodesic distance to the boundary loop
	for(auto pVertex = m_polyhedron1->vertices_begin(); pVertex!=m_polyhedron1->vertices_end(); ++pVertex)
	{
		geodesic::SurfacePoint gp(&m_g1.vertices()[pVertex->tag()]);
		double distance;
		alg1->best_source(gp,distance);
		m_distToLoop[pVertex] = distance;
// 		//m_treeStructure1.m_vertices.insert(pVertex);
	}
	//delete alg1;
	
	initGeodesicMesh(m_polyhedron2,&m_g2);
	geodesic::GeodesicAlgorithmDijkstraAlternative * alg2 = new geodesic::GeodesicAlgorithmDijkstraAlternative(&m_g2);
	for(unsigned i=0;i<loop2.size();++i)
	{
		int index = loop2[i]->tag();
		geodesic::SurfacePoint gp(&m_g2.vertices()[index]);
		sources2.push_back(gp);
	}
	alg2->propagate(sources2);
	for(auto pVertex = m_polyhedron2->vertices_begin(); pVertex!=m_polyhedron2->vertices_end(); ++pVertex)
	{
		geodesic::SurfacePoint gp(&m_g2.vertices()[pVertex->tag()]);
		double distance;
		alg2->best_source(gp,distance);
		m_distToLoop[pVertex] = distance;
		//m_treeStructure2.m_vertices.insert(pVertex);
	}
	//delete alg2;

	
	// m_R is the largest distance among the vertices in the snapping region
	double sR = 0.0;
	for(int i=0; i<closestToB2.size();i++)
	{
		double dist = m_distToLoop[closestToB2[i]];
		if(dist > sR) {sR = dist;}
	}
	for(int i=0; i<closestToB1.size();i++)
	{
		double dist = m_distToLoop[closestToB1[i]];
		if(dist > sR) {sR = dist;}
	}
	m_R = sR;//*1.10;
	
	//m_R = sR*1.1;
	
	double factoredR = m_R;
	if(factor>1)
	{
		factoredR*=factor;
	}
		
	
	//The snapping region is the submesh within R geodesic distance to the b-loop
	for(auto pVertex = m_polyhedron1->vertices_begin();
		pVertex != m_polyhedron1->vertices_end();
		++pVertex)
	{
		double dist = m_distToLoop[pVertex];
		if(dist <= m_R)
		{
			m_isSnappingRegion[pVertex]=true;
			m_sr1.push_back(pVertex);
			pVertex->color(1,0,0);
			m_treeStructure1.m_vertices.insert(pVertex);
		}
		else if( dist < factoredR)
		{
			m_treeStructure1.m_vertices.insert(pVertex);
			m_isFactoredRegion[pVertex] = true;
			pVertex->color(1,1,0);
		}
	}
	for(auto pVertex = m_polyhedron2->vertices_begin();
		pVertex != m_polyhedron2->vertices_end();
		++pVertex)
	{
		double dist = m_distToLoop[pVertex];
		if(dist <= m_R)
		{
			m_treeStructure2.m_vertices.insert(pVertex);
			m_isSnappingRegion[pVertex] = true;
			m_sr2.push_back(pVertex);
			pVertex->color(1,0,0);
		}
		else if( dist < factoredR)
		{
			m_treeStructure2.m_vertices.insert(pVertex);
			m_isFactoredRegion[pVertex] = true;
			pVertex->color(1,1,0);
		}
	}	
	
	//Divide the factored region
	/*int nbSubdiv = 0;
	while(nbSubdiv<2){
	Facet_iterator B1 = m_polyhedron1->facets_begin();
	Facet_iterator E1 = m_polyhedron1->facets_end(); --E1;
	
	Facet_iterator B2 = m_polyhedron2->facets_begin();
	Facet_iterator E2 = m_polyhedron2->facets_end(); --E2;

	Facet_iterator f = B1;
	do {
		int countER = 0;
		auto hC = f->facet_begin();
		do
		{
			if(m_isSnappingRegion[hC->vertex()] || m_isFactoredRegion[hC->vertex()]){countER++;}
		}
		while(++hC!=f->facet_begin());
		if(countER == 3)
		{
			double distToLoop = centerDistance(f);
			Halfedge_handle h = create_center_vertex(m_polyhedron1,f);
			Vertex_handle p = h->vertex();
			
			m_treeStructure1.m_vertices.insert(p);
			m_distToLoop[p] = distToLoop;
			
			if(distToLoop <= m_R)
			{
				m_isSnappingRegion[p] = true;
				m_sr1.push_back(p);
				p->color(1,0,0);
			}
			else if( distToLoop < factoredR )
			{
				m_isFactoredRegion[p] = true;
				p->color(1,1,0);
			}
		}
	} while ( f++ != E1);
	
	Facet_iterator g = B2;
	do {
		int countER = 0;
		auto hC = g->facet_begin();
		do
		{
			if(m_isSnappingRegion[hC->vertex()] || m_isFactoredRegion[hC->vertex()]){countER++;}
		}
		while(++hC!=g->facet_begin());
		if(countER == 3)
		{
			double distToLoop = centerDistance(g);
			Halfedge_handle h = create_center_vertex(m_polyhedron2,g);
			Vertex_handle p = h->vertex();
			
			m_treeStructure2.m_vertices.insert(p);
			m_distToLoop[p] = distToLoop;
			
			if(distToLoop <= m_R)
			{
				m_isSnappingRegion[p] = true;
				m_sr2.push_back(p);
				p->color(1,0,0);
			}
			else if( distToLoop < factoredR )
			{
				m_isFactoredRegion[p] = true;
				p->color(1,1,0);
			}
		}
	} while ( g++ != E2); 
	nbSubdiv++;
	}*/
	
	computeSnapRegionDistances();
	m_polyhedron1->compute_normals();
	m_polyhedron2->compute_normals();
	
	delete alg1;
	delete alg2;
}

void softICPController::computeSnapRegionDistances()
{
	std::vector<geodesic::SurfacePoint> sources1;
	std::vector<geodesic::SurfacePoint> sources2;
	
	geodesic::GeodesicAlgorithmExact * alg1 = new geodesic::GeodesicAlgorithmExact(&m_g1);
	for(unsigned i=0;i<m_sr1.size();++i)
	{
		int index = m_sr1[i]->tag();
		geodesic::SurfacePoint gp(&m_g1.vertices()[index]);
		sources1.push_back(gp);
	}
	alg1->propagate(sources1);
	for(auto pVertex = m_polyhedron1->vertices_begin();
		pVertex!=m_polyhedron1->vertices_end();
		++pVertex)
	{
		if(!m_isSnappingRegion[pVertex] && m_isFactoredRegion[pVertex])
		{
			double dist;
			geodesic::SurfacePoint gp(&m_g1.vertices()[pVertex->tag()]);
			alg1->best_source(gp,dist);
			m_distToSnap[pVertex] = dist;
		}
	}
	
	
	geodesic::GeodesicAlgorithmExact * alg2 = new geodesic::GeodesicAlgorithmExact(&m_g2);
	for(unsigned i=0;i<m_sr2.size();++i)
	{
		int index = m_sr2[i]->tag();
		geodesic::SurfacePoint gp(&m_g2.vertices()[index]);
		sources2.push_back(gp);
	}
	alg2->propagate(sources2);
	for(auto pVertex = m_polyhedron2->vertices_begin();
		pVertex!=m_polyhedron2->vertices_end();
		++pVertex)
	{
		if(!m_isSnappingRegion[pVertex] && m_isFactoredRegion[pVertex])
		{
			double dist;
			geodesic::SurfacePoint gp(&m_g2.vertices()[pVertex->tag()]);
			alg2->best_source(gp,dist);
			m_distToSnap[pVertex] = dist;
		}
	}
	
	
	sources1.clear();
	for(auto pVertex = m_polyhedron1->vertices_begin(); pVertex!=m_polyhedron1->vertices_end();++pVertex)
	{
		if(!m_isSnappingRegion[pVertex] && m_isFactoredRegion[pVertex])
		{
			int index = pVertex->tag();
			geodesic::SurfacePoint gp(&m_g1.vertices()[index]);
			sources1.push_back(gp);
		}
	}
	alg1->propagate(sources1);
	for(unsigned v = 0; v < m_sr1.size(); ++v)
	{
		Vertex_handle pVertex = m_sr1[v];
		if(m_isSnappingRegion[pVertex] || m_isFactoredRegion[pVertex])
		{
			double dist;
			geodesic::SurfacePoint gp(&m_g1.vertices()[pVertex->tag()]);
			alg1->best_source(gp,dist);
			m_distToSnap[pVertex] = dist;
		}
	}
	
	sources2.clear();
	for(auto pVertex = m_polyhedron2->vertices_begin(); pVertex!=m_polyhedron2->vertices_end();++pVertex)
	{
		if(!m_isSnappingRegion[pVertex] && m_isFactoredRegion[pVertex])
		{
			int index = pVertex->tag();
			geodesic::SurfacePoint gp(&m_g2.vertices()[index]);
			sources2.push_back(gp);
		}
	}
	alg2->propagate(sources2);
	for(unsigned v = 0; v < m_sr2.size(); ++v)
	{
		Vertex_handle pVertex = m_sr2[v];
		if(m_isSnappingRegion[pVertex] || m_isFactoredRegion[pVertex])
		{
			double dist;
			geodesic::SurfacePoint gp(&m_g2.vertices()[pVertex->tag()]);
			alg2->best_source(gp,dist);
			m_distToSnap[pVertex] = dist;
		}
	}
	
	
	delete alg1;
	delete alg2;
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
			//double dist = computePhiDistance(rep,*v,0.5,0.1,0.4);
			double dist = computePhiDistance(rep,*v,1.0,0.0,0.0);
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
					//double dist = computePhiDistance(*v,*g,0.5,0.1,0.4);
					double dist = computePhiDistance(*v,*g,1.0,0.0,0.0);
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
			//double dist = computePhiDistance(rep1,rep2,0.5,0.1,0.4);
			double dist = computePhiDistance(rep1,rep2,1.0,0.0,0.0);
			m[i][j] = dist;
			m[j][i] = dist;
		}
	}
}


void softICPController::snapRegions(double R, double elasticity, int itermax, int treeDepth)
{
	std::cout << "itermax :" << itermax << std::endl;
	std::cout << "Total time : ";
	//_tic();
	bool convergenceCriterion = false;
	bool order = false;
	int iter = 0;
	
	// Select snapping region and create patches on surface
	//std::cout << "buildTreeStructure : "; 
	//timer_tic();
	
	//buildTreeStructure(treeDepth,R);
	
	buildFullTreeStructure(treeDepth,R);
	//timer_toc();
	
	//R = m_R*1.20;
	R = m_R;	
	
	//storeNodeLevel(&m_treeStructure1,0,m_levelNodes1);
	//storeNodeLevel(&m_treeStructure2,0,m_levelNodes2);	
	//std::cout << "computeClosestGeo" << std::endl;
	//timer_tic();
	//computeClosestGeo(&m_treeStructure1,7,m_levelNodes1,m_distanceMatrices1,m_g1);
	//computeClosestGeo(&m_treeStructure2,7,m_levelNodes2,m_distanceMatrices2,m_g2);
	//timer_toc();
	//std::cout << "end computeClosestGeo" << std::endl;
	
	//initClosest(m_levelNodes1,m_distanceMatrices1,7);
	//initClosest(m_levelNodes2,m_distanceMatrices2,7);

	//colorCluster(&m_treeStructure1);
	//colorCluster(&m_treeStructure2);
	
	//colorLastClusters(&m_treeStructure1);
	//colorLastClusters(&m_treeStructure2);
	
	int id =0; 
	for(auto pVertex =m_polyhedron1->vertices_begin();pVertex!=m_polyhedron1->vertices_end();++pVertex)
	{
		pVertex->tag() = id;
		id++;
	}
	id =0; 
	for(auto pVertex =m_polyhedron2->vertices_begin();pVertex!=m_polyhedron2->vertices_end();++pVertex)
	{
		pVertex->tag() = id;
		id++;
	}
	
	
	/*ofstream source, corresp,transform;
	source.open("/home/leon/source.txt");
	corresp.open("/home/leon/corresp.txt");			
	transform.open("/home/leon/transform.txt");*/
	
	while(!convergenceCriterion && !m_stop_for_debug)
	{
		// swap meshes at each iteration
		order = !order; 
		
		PolyhedronPtr ma;
		PolyhedronPtr mb;
		
		convergenceCriterion = (iter == itermax-1);
		
		// find correspondence between the 2 snapping regions
		computeSnappingRegionCorrespondence(order);
	
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
			}
			iter++;
			std::cout << "iteration: " << iter << std::endl;
			
			// Compute a transformation for each point
			unsigned int nbP = 0;
			std::vector <pointTransformation> transf;
			
			int csr = 0;
			for(auto it=treeStructure->m_vertices.begin(); it!=treeStructure->m_vertices.end() && !m_stop_for_debug;++it)
			{
				
				Vertex_handle pVertex = *it;
				
				/*if(iter == 10 && m_distToLoop[pVertex] == 0.0)
				{
					m_stop_for_debug = true;
				}*/
				
				// Get neighborhood N around pVertex
				std::vector<Vertex_handle> N = getNeighborhoodNoTree(pVertex,R,iter,itermax,elasticity,order);
				//std::vector<Vertex_handle> N = getNeighborhood(pVertex,R,iter,itermax,elasticity,order);
		
				// Get corresponding neighborhoods
				std::vector<Vertex_handle> phiN = getCorrespondingNeighborhood(N);
				
				// Compute and store transformation
				pointTransformation ti;
				
				/*if(m_stop_for_debug)
				{
					ti.S = 1.0;
					ti.T = Vec(0.0,0.0,0.0);
					ti.R = Matrix::eye(3);
					ti.Q = Quaternion();
					
				}*/
				//else 
				//{
				ti = computeTransformationNew(N,phiN);
				//}
					
				transf.push_back(ti);
				nbP++;
			}
			
			// Apply the transformation for each point, but scale it (iter/itermax)
			int i=0;
			for(auto it=treeStructure->m_vertices.begin(); it!=treeStructure->m_vertices.end() && (i<nbP) ;++it)
			{	
				Vertex_handle pVertex = *it;
				applyTransformation(pVertex,transf[i],iter,itermax);
				i++;
			}
			//if(!m_stop_for_debug)
			//{
			finalTransform(order,iter,itermax);
			//}
			//m_v->recreateListsAndUpdateGL();
			//SleeperThread::msleep(300);
		}
	}
	//moveBorder(order);
	//moveBorder(!order);
	//finalTransform(order,1,1);
	//finalTransform(!order,1,1);
	//timer_toc();
}
																																
void softICPController::applyTransformation(Vertex_handle p, pointTransformation & ti, int iter, int itermax)
{
	double li = iter/(double)itermax;
	
	//double scale = li*ti.S + (1-li);
	
	//Vec scaledTranslation = li*ti.T;
	
	Matrix sR = ti.R;
	double rData[3][3];
	for(unsigned i=0;i<3;++i)
	{
		for(unsigned j=0;j<3;++j)
		{
 			rData[i][j] = sR.val[i][j];
		}
	}
	
	//Quaternion scaledQ;
	Quaternion q;
	q.setFromRotationMatrix(rData);
	//ti.Q.setFromRotationMatrix(rData);
	//std::cout << "qu : " << ti.Q.axis() << " " << ti.Q.angle() << std::endl;	
	//Quaternion unitQ;
	//scaledQ = Quaternion::slerp(unitQ,ti.Q,li);
	
	//std::cout << "translation : " << scaledTranslation << std::endl;
	//std::cout << "rotation : " << scaledQ.axis() << " " << scaledQ.angle() << std::endl;
	
	Point3d pos = p->point();
	 
	qglviewer::Vec V = q*qglviewer::Vec(pos.x(),pos.y(),pos.z())*ti.S + ti.T;
	Point3d newPos = CGAL::ORIGIN + Vector(V[0],V[1],V[2]);
	double intX,intY,intZ;
	intX = newPos.x()*li + pos.x()*(1-li);	
	intY = newPos.y()*li + pos.y()*(1-li);
	intZ = newPos.z()*li + pos.z()*(1-li);
	Point3d interP(intX,intY,intZ);
	p->point() = interP;
}

std::vector< Vertex_handle > softICPController::getNeighborhood(Vertex_handle p, double R, unsigned int iter,unsigned itermax, double elasticity, bool order)
{
	double distToLoop = m_distToLoop[p];
	
	double radius = R/20.0;
	
	double elasticityR = elasticity*R;
	
	// compute the size of the local neighborhood
	if(distToLoop != 0.0)
	{
		double expo = ((iter/(double)itermax)*elasticityR)/(distToLoop);
		radius += R * exp(-expo*expo);
	}
	
	
	double ratio = radius / R;
	
	if(!m_isSnappingRegion[p]) // if the vertex is outside the snapping region, add distance to the snapping region
	{
		//radius += 2*m_distToSnap[p];
		radius += m_distToLoop[p];
	}
	
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
	
	std::vector<Vertex_handle> N;
 	containsP = node; // leaf node
	
	bool getout = false;
	
	while(!getout)
	{
		if(node->m_distClosest.size() == 0)
		{
			for(auto it = node->m_vertices.begin();it!=node->m_vertices.end();++it)
				{
					if(m_isSnappingRegion[*it])
					{N.push_back(*it);}
				}
			getout = true;
			break;
		}
		double maxDist = node->m_distClosest.back(); // the m furthest distance
		
		if( radius <= maxDist ) // collect all the patches whose distance is smaller than radius
		{
			for(unsigned c=0;c<node->m_distClosest.size();++c)
			{
				if(node->m_distClosest[c]<=radius) 
				{
					deformationNode * patch = node->m_closest[c];
					//std::cout << patch->m_vertices.size() << " " << patch->m_clusterID << std::endl;
					for(auto it = patch->m_vertices.begin();it!=patch->m_vertices.end();++it)
					{
						if(m_isSnappingRegion[*it])
						{
							N.push_back(*it);
						}
					}
				}
			}
			
			getout = true;
		}
		
		// go up one level in the hierarchy and repeat the process
		if ((radius >maxDist || N.size() < 4) )
		{
			N.clear();
			if(node->m_clusterID == -1 ) // can't go up one level
			{
				for(auto it = node->m_vertices.begin();it!=node->m_vertices.end();++it)
				{
					if(m_isSnappingRegion[*it])
					{
						N.push_back(*it);
					}
				}
				getout = true;
			}
			else{
				//if( radius == 0){std::cout << "radius = 0 " << N.size() << std::endl;}
				node = node->m_parent;
				getout = false;
			}
		}
	}
	
	//std::cout << iter << " " << order << " " << distToLoop << " " <<  N.size () << " " << radius << std::endl;
	
	
	
	return N;
}

std::vector< Vertex_handle > softICPController::getNeighborhoodNoTree(Vertex_handle p, double R, unsigned int iter,unsigned itermax, double elasticity, bool order)
{
	
	double distToLoop = m_distToLoop[p];
	
	double radius = R/10.0;
	//double radius = R/40.0;
	
	double elasticityR = elasticity*R;
	
	// compute the size of the local neighborhood
	if(distToLoop != 0.0)
	{
		double expo = ((iter/(double)itermax)*elasticityR)/(distToLoop);
		radius += R * exp(-expo*expo);
	}
	
	if(!m_isSnappingRegion[p]) // if the vertex is outside the snapping region, add distance to the snapping region
	{
		radius += m_distToSnap[p];
 		//radius += 2*m_distToSnap[p];
		//radius += m_distToLoop[p];
	}
	else if(m_stop_for_debug)
	{
		std::cout << "radius : " << radius << std::endl;
		std::cout << "R : " << R << std::endl;
	}
	
	// find the geodesic neighborhood of size 'radius'
	std::vector<Vertex_handle> N;
	std::set<Vertex_handle> vertices;
	Point3d O = p->point();
	std::stack<Vertex_handle> S;
	int nbLoop = 1;
	while(nbLoop <= 20){
		S.push(p);
		vertices.insert(p);
		while(!S.empty())
		{
			Vertex_handle v = S.top();
			S.pop();
			Point3d P =v->point();
			
			Halfedge_around_vertex_circulator h = v->vertex_begin();
			Halfedge_around_vertex_circulator pHalfedgeStart = h;
			CGAL_For_all(h,pHalfedgeStart)
			{
				Point3d p1 = h->vertex()->point();
				Point3d p2 = h->opposite()->vertex()->point();
				Vector V = (p2-p1);
				//if( v == p || V * (P - O) > 0.0)
				//{
					bool isect = sphere_clip_vectorG(O,radius, P, V);
					if(!isect)
					{
						Vertex_handle w = h->opposite()->vertex();
						if(vertices.find(w) == vertices.end())
						{
							
							vertices.insert(w);
							S.push(w);
						}
					}
				//}
			}
		}
		
		for(auto v = vertices.begin(); v!=vertices.end(); ++v)
		{
			
			if(m_isSnappingRegion[*v])
			{
				N.push_back(*v);
			}
		}
		
		/*if(N.size() == 0)
		{
			for(auto v = vertices.begin(); v!=vertices.end(); ++v)
			{
				//(*v)->color(0.0,1.0,1.0);
			}
		}*/
		
		std::set<Point3d> corresp;
		for(unsigned c=0;c<N.size();++c)
		{
			corresp.insert(m_Phi[N[c]]->point());
		}
		
		if( N.size() == 0)
		{
			radius+= R/10.0;
			//std::cout << "size is zero" << std::endl;
			//std::cout << "radius : " << radius << std::endl;
		}
		else if( corresp.size() < 3 && nbLoop!=20 ) 
		{
			//std::cout << "nbLoop : " << nbLoop << " radius : " << radius << std::endl;
			N.clear();
			corresp.clear();
			radius+=R/10.0;
			nbLoop++;
		}
		else 
		{
			break;
		}
	}
	if(N.size() == 0 )
	{
		//std::cout << "N.size() = 0  _ iter=" << nbLoop << std::endl;
		p->color(1.0,0.0,1.0);
		//m_stop_for_debug = true;
	}
	
	if(m_stop_for_debug)
	{
		for(auto it = N.begin(); it != N.end(); ++it)
		{
			Vertex_handle v = (*it);
			v->color(1.0,1.0,1.0);
			m_Phi[v]->color(0.0,1.0,1.0);
		}
		std::cout << " N.size : " << N.size() << std::endl;
		std::cout << " radius : " << radius << std::endl;	
	}
	
	return N;
}

void softICPController::moveToCorrespondence(bool order)
{
	std::vector<Vertex_handle> * s;
	PolyhedronPtr m;
	if(order)
	{
		s = &m_sr1;
		m = m_polyhedron1;
	}
	else // reverse order between meshes
	{
		s = &m_sr2;
		m = m_polyhedron2;
	}
	
	computeSnappingRegionCorrespondence(order);
	
	for(auto v = s->begin(); v!= s->end(); ++v)
	{
		Point3d p = m_Phi[*v]->point();
		(*v)->point() = p;
	}
}

void softICPController::moveBorder(bool order)
{
	std::vector<Vertex_handle> * s;
	PolyhedronPtr m;
	if(order)
	{
		s = &m_sr1;
		m = m_polyhedron1;
	}
	else // reverse order between meshes
	{
		s = &m_sr2;
		m = m_polyhedron2;
	}
	
	computeSnappingRegionCorrespondence(order);
	
	for(auto v = s->begin(); v!= s->end(); ++v)
	{
		if(m_distToLoop[*v]==0.0)
		{
			Point3d p = m_Phi[*v]->point();
			(*v)->point() = p;
		}
	}
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

void softICPController::finalTransform(bool order, int iter, int itermax)
{
	std::vector<Vertex_handle> * s;
	PolyhedronPtr m;
	if(order)
	{
		s = &m_sr1;
		m = m_polyhedron1;
	}
	else // reverse order between meshes
	{
		s = &m_sr2;
		m = m_polyhedron2;
	}
	
	//computeSnappingRegionCorrespondence(order);
	
	std::vector<Vertex_handle> N;
	
	// Get all possible vertices (biggest support region)
	for(auto v = s->begin(); v!=s->end();++v)
	{
		N.push_back(*v);
	}
	
	std::vector<Vertex_handle> phiN = getCorrespondingNeighborhood(N);
	
	pointTransformation pT = computeTransformationNew(N,phiN);
	
	// Apply the final rigid transformation to all points that are not in the extended region
	for(auto pVertex=m->vertices_begin();pVertex!=m->vertices_end();++pVertex)
	{
		
		if(!m_isSnappingRegion[pVertex] && !m_isFactoredRegion[pVertex])
		{
			pVertex->color(0,0,1);
			applyTransformation(pVertex,pT,iter,itermax);
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

pointTransformation softICPController::computeTransformationNew(std::vector<Vertex_handle> & N, std::vector<Vertex_handle> & phiN)
{
	pointTransformation pi;
	
	//Center data
	Matrix p_m(N.size(),3);
	Matrix p_t(N.size(),3);

	Matrix mu_m(1,3);
	Matrix mu_t(1,3);
	double mum0 = 0.0, mum1 = 0.0, mum2 = 0.0;
	double mut0 = 0.0, mut1 = 0.0, mut2 = 0.0;
	
	double distM = 0.0;
	double distT = 0.0;
	
	int i=0;
	int sizeN = N.size();
	
	std::vector<Point3d> nPts;
	std::vector<Point3d> pnPts;

	std::set<Point3d> corres;
	
	for(i=0;i<sizeN;++i)
	{	
		Point3d npp =  N[i]->point();
		nPts.push_back(npp);
		p_m.val[i][0] = npp.x(); mum0 += p_m.val[i][0];
		p_m.val[i][1] = npp.y(); mum1 += p_m.val[i][1];
		p_m.val[i][2] = npp.z(); mum2 += p_m.val[i][2];
		
		Point3d np = phiN[i]->point();
		pnPts.push_back(np);
		corres.insert(np);
		p_t.val[i][0] = np.x(); mut0 +=p_t.val[i][0];
		p_t.val[i][1] = np.y(); mut1 +=p_t.val[i][1];
		p_t.val[i][2] = np.z(); mut2 +=p_t.val[i][2];
	}
	
	//std::cout << sizeN << " " << corres.size() << std::endl;
	
	if(N.size() == 1 && corres.size() == 1)
	{
		pi.Q = Quaternion();
		pi.R = Matrix::eye(3);
		pi.S = 1.0;
		
		double tx,ty,tz;
		tx = phiN[0]->point().x() - N[0]->point().x();
		ty = phiN[0]->point().y() - N[0]->point().y(); 
		tz = phiN[0]->point().z() - N[0]->point().z(); 
		pi.T.setValue(tx,ty,tz);
		return pi;
	}
	else if (N.size() == 0)
	{
		pi.Q = Quaternion();
		pi.R = Matrix::eye(3);
		pi.S = 1.0;
		pi.T.setValue(0.0,0.0,0.0);
		return pi;
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
	
	//Scale
	double scale = 1.0;
	double scale_m = q_m.l2norm()/q_m.m;
	double scale_t = q_t.l2norm()/q_t.m;
	scale = scale_t/scale_m;
	q_m = q_m * scale;
	
	//Rotation
	Matrix H = ~q_m * q_t;
	Matrix T;
	
	Matrix U,W,V;
	H.svd(U,W,V);
	
	Matrix R = V*~U;
	double detR = R.det();
	if(detR<0)
	{
		Matrix B = Matrix::eye(3);
		B.val[2][2] = detR;
		R = V*B*~U;
	}
	double rData[3][3];
	for(unsigned i=0;i<3;++i)
	{
		for(unsigned j=0;j<3;++j)
		{
			rData[i][j] = R.val[i][j];
		}
	}	
	
	//Translation 
	Matrix t = ~mu_t - (R*scale)*~mu_m;
	
	pi.Q.setFromRotationMatrix(rData);
	pi.R = R;
	pi.S = scale;
	pi.T.setValue(t.val[0][0],t.val[1][0],t.val[2][0]);
	
	return pi;
}


pointTransformation softICPController::computeTransformation(std::vector<Vertex_handle> & N, std::vector<Vertex_handle> & phiN, Vertex_handle p, int iter, const int itermax)
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
	
	double distM = 0.0;
	double distT = 0.0;
	
	int i=0;
	int sizeN = N.size();
	
	std::vector<Point3d> nPts;
	std::vector<Point3d> pnPts;
	
	
	for(i=0;i<sizeN;++i)
	{	
		Point3d npp =  phiN[i]->point();
		nPts.push_back(npp);
		p_m.val[i][0] = npp.x(); mum0 += p_m.val[i][0];
		p_m.val[i][1] = npp.y(); mum1 += p_m.val[i][1];
		p_m.val[i][2] = npp.z(); mum2 += p_m.val[i][2];
		
		Point3d np = N[i]->point();
		pnPts.push_back(np);
		p_t.val[i][0] = np.x(); mut0 +=p_t.val[i][0];
		p_t.val[i][1] = np.y(); mut1 +=p_t.val[i][1];
		p_t.val[i][2] = np.z(); mut2 +=p_t.val[i][2];
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
	
	for(unsigned i=0;i<sizeN;++i)
	{
		distM += q_m.getMat(i,0,i,2).l2norm();
		distT += q_t.getMat(i,0,i,2).l2norm();
	}
	
	distM/=N.size();
	distT/=N.size();
	
	// compute rotation matrix R and translation vector 
	Matrix H = ~q_t*q_m;
	//Matrix H = ~q_m*q_t;
	Matrix U,W,V;
	H.svd(U,W,V);
	Matrix R = V*~U;
	
	double detR = R.det();
	
	//std::cout << "detR : " << detR << std::endl;
	
	if(detR<0)
	{
		Matrix B = Matrix::eye(3);
		B.val[2][2] = detR;
		R = V*B*~U;
	}	
	
	double singMean = 0.0;
	for(unsigned i=0; i<W.m;++i)
	{
		singMean+= abs(W.val[i][i]);
	}
	singMean/=W.m;

	Matrix t = ~mu_m - R*~mu_t;

	double rData[3][3];
	for(unsigned i=0;i<3;++i)
	{
		for(unsigned j=0;j<3;++j)
		{
			rData[i][j] = R.val[i][j];
		}
	}
	
	double scaling = distT/distM;
	
	pi.R = R;
	pi.Q.setFromRotationMatrix(rData);
	pi.T.setValue(t.val[0][0],t.val[1][0],t.val[2][0]);
	pi.S = singMean ;
	
	if(m_isSnappingRegion[p])
	{
		double li = iter/(double)itermax;
		
		double scale = li*pi.S + (1-li);
	
		Vec scaledTranslation = li*pi.T;
	
		Matrix sR = pi.R;
		double rData[3][3];
		for(unsigned i=0;i<3;++i)
		{
			for(unsigned j=0;j<3;++j)
			{
				rData[i][j] = sR.val[i][j];
			}
		}
		
		Quaternion scaledQ,res;
		scaledQ.setFromRotationMatrix(rData);	
		Quaternion unitQ;
		res = Quaternion::slerp(unitQ,scaledQ,li);
	
		Point3d pos = m_Phi[p]->point();
		qglviewer::Vec VV = res*qglviewer::Vec(pos.x(),pos.y(),pos.z())*scaling + scaledTranslation;
		Point3d q  = CGAL::ORIGIN + Vector(VV[0],VV[1],VV[2]);
		
		double distTransform = CGAL::squared_distance(pos,q);
	}
	
	return pi;
}

void softICPController::remeshTri(Viewer *v)
{
	// for both meshes, collect all points in the snapping region
	std::vector<Point3d> points;
	for(auto pVertex = m_polyhedron1->vertices_begin();
	    pVertex != m_polyhedron1->vertices_end();
		++pVertex)
	{
		if(m_isSnappingRegion[pVertex]||m_isFactoredRegion[pVertex])
		{
			points.push_back(pVertex->point());
		}
	}
	for(auto pVertex = m_polyhedron2->vertices_begin();
	    pVertex != m_polyhedron2->vertices_end();
		++pVertex)
	{
		if(m_isSnappingRegion[pVertex]||m_isFactoredRegion[pVertex])
		{
			points.push_back(pVertex->point());
		}
	}
	// Smooth the point set of the snapping region
	this->smoothPointSet(points,2);
	
	//Affect the new values to the vertices
	unsigned u = 0;
	for(auto pVertex = m_polyhedron1->vertices_begin();
	    pVertex != m_polyhedron1->vertices_end();
		++pVertex)
	{
		if(m_isSnappingRegion[pVertex]||m_isFactoredRegion[pVertex])
		{
			pVertex->point() = points[u];
			u++;
		}
	}
	for(auto pVertex = m_polyhedron2->vertices_begin();
	    pVertex != m_polyhedron2->vertices_end();
		++pVertex)
	{
		if(m_isSnappingRegion[pVertex]||m_isFactoredRegion[pVertex])
		{
			pVertex->point() = points[u];
			u++;
		}
	}
	
	v->recreateListsAndUpdateGL();
	
	// Now fill the borders of the meshes for boolean union
	
	PolyhedronPtr p1 = fillHoles(m_polyhedron1);
	PolyhedronPtr p2 = fillHoles(m_polyhedron2);
	
	v->getScenePtr()->add_polyhedron(p1);
	v->getScenePtr()->add_polyhedron(p2);
	v->getScenePtr()->todoIfModeSpace(v,0.0);
	
}

PolyhedronPtr softICPController::fillHolesSimple(PolyhedronPtr p)
{
	while(!p->is_closed())
	{
		p->fill_hole(p->border_halfedges_begin()->opposite());
		p->normalize_border();
	}
	p->triangulate();
	return p;
}


PolyhedronPtr softICPController::fillHoles(PolyhedronPtr p)
{
	p->normalize_border();
	
	std::vector<double> coords;
	std::vector<int> tris;
	
	p->set_index_vertices();
	for(auto po = p->points_begin(); po != p->points_end(); ++po)
	{
		coords.push_back(po->x());
		coords.push_back(po->y());
		coords.push_back(po->z());
	}
	for(auto f = p->facets_begin(); f!=p->facets_end(); ++f)
	{
		
		auto h = f->facet_begin();
		do
		{
			tris.push_back(h->vertex()->tag());
		}
		while(++h!=f->facet_begin());
		
	}
	
	constructPolyhedron * fullMesh(new constructPolyhedron);
	polyhedron_builder<constructPolyhedron::HalfedgeDS> builder(coords,tris);
	isolatedVertices_remover<constructPolyhedron::HalfedgeDS> isolatedVerts;
	fullMesh->delegate(builder);
	fullMesh->delegate(isolatedVerts);
	fullMesh->normalize_border();
	
	std::cout << " reconstructed fullMesh ok " << std::endl;
	
	
	  bool intersecting = CGAL::Polygon_mesh_processing::does_self_intersect(*fullMesh,
     CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point,* fullMesh)));

	  std::cout
    << (intersecting ? "There are self-intersections." : "There is no self-intersection.")
    << std::endl;
			
      std::vector<std::pair<constructPolyhedron::Facet_handle,constructPolyhedron::Facet_handle > > intersected_tris;
     CGAL::Polygon_mesh_processing::self_intersections(*fullMesh,
    std::back_inserter(intersected_tris),
       CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, *fullMesh)));
	
	for(unsigned u = 0 ; u < intersected_tris.size(); ++u)
	{
		std::cout << u << std::endl;
		constructPolyhedron::Facet_handle f = intersected_tris[u].first;
		std::cout << "ok" << std::endl;
		
		if (f!= constructPolyhedron::Facet_handle() && f->halfedge() != constructPolyhedron::Halfedge_handle())
		{
			if( ! f->halfedge()->is_border())
			{
				std::cout << "before erase facet" << std::endl;
				fullMesh->erase_facet(f->halfedge());
			}
		}
	}
	
	std::cout << "remove self_intersections " << std::endl;
     
	std::vector<constructPolyhedron::Facet_handle>  patch_facets;
	std::vector<constructPolyhedron::Vertex_handle> patch_vertices;
	for(auto h = fullMesh->halfedges_begin(); h!=fullMesh->halfedges_end(); ++h)
	{
		if(h->is_border())
		{
		CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(
			*fullMesh,
			h,
			std::back_inserter(patch_facets),
			std::back_inserter(patch_vertices),
			CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, *fullMesh)).
			geom_traits(constructionK())); 
		}
	}
	std::cout << " triangulate hole ok " << std::endl;
	
	PolyhedronPtr result = convertToEnrichedPolyhedron(fullMesh);
	isolatedVertices_remover<HalfedgeDS> ivr;
	result->delegate(ivr);
	result->compute_normals();
	
	return result;
}

void softICPController::remeshNew(Viewer *v)
{
	std::cout << "start remesh" << std::endl;
	PolyhedronPtr outMesh = getMeshOutsideSR();
	for(auto p=outMesh->vertices_begin();p!=outMesh->vertices_end();++p)
	{
		p->color(1.0,0.0,0.0);
	}
	outMesh->normalize_border();
	outMesh->compute_normals();
	
	std::cout << "get outside mesh done " << std::endl;
	
	v->getScenePtr()->add_polyhedron(outMesh);
		
	// collect border points for surface Reconstruction filtering
	std::set<Point3d> border;
	for(auto h = outMesh->border_halfedges_begin(); h!=outMesh->halfedges_end();++h)
	{
		border.insert(h->vertex()->point());
	}
	
	PolyhedronPtr srMesh = buildSRMesh(border,outMesh);
	for(auto p=srMesh->vertices_begin();p!=srMesh->vertices_end();++p)
	{
		p->color(0.0,1.0,0.0);
	}
	srMesh->compute_normals();
	v->getScenePtr()->add_polyhedron(srMesh);
	
	
	//v->getScenePtr()->setVisible(0,false);
	//v->getScenePtr()->setVisible(1,false);
	
	srMesh->normalize_border();
	
	std::vector<double> coords;
	std::vector<int> tris;
	
	outMesh->set_index_vertices();
	srMesh->set_index_vertices();
	
	for(auto po = outMesh->points_begin(); po != outMesh->points_end(); ++po)
	{
		coords.push_back(po->x());
		coords.push_back(po->y());
		coords.push_back(po->z());
	
	}
	
	
	int vertOffset = outMesh->size_of_vertices();
	for(auto f = outMesh->facets_begin(); f!=outMesh->facets_end(); ++f)
	{
		
		auto h = f->facet_begin();
		do
		{
			tris.push_back(h->vertex()->tag());
		}
		while(++h!=f->facet_begin());
		
	}
	
	constructPolyhedron * fullOutMesh(new constructPolyhedron);
	polyhedron_builder<constructPolyhedron::HalfedgeDS> builder(coords,tris);
	fullOutMesh->delegate(builder);


	// fill holes that are not part of the main mesh borders
	std::vector<constructPolyhedron::Facet_handle>  patch_facets;
	std::vector<constructPolyhedron::Vertex_handle> patch_vertices;
	for(auto h = fullOutMesh->halfedges_begin(); h!=fullOutMesh->halfedges_end(); ++h)
	{
		auto p = h->vertex()->point();
		Point3d pp(p.x(),p.y(),p.z());
		if(h->is_border() )
		{
		CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(
			*fullOutMesh,
			h,
			std::back_inserter(patch_facets),
			std::back_inserter(patch_vertices),
			CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, *fullOutMesh)).
			geom_traits(constructionK())); 
		}
	}
	
	coords.clear();
	tris.clear();
	
	
	for(auto po = srMesh->points_begin(); po != srMesh->points_end(); ++po)
	{
		coords.push_back(po->x());
		coords.push_back(po->y());
		coords.push_back(po->z());
	
	}
	for(auto f = srMesh->facets_begin(); f!=srMesh->facets_end(); ++f)
	{
		
		auto h = f->facet_begin();
		do
		{
			tris.push_back(h->vertex()->tag());
		}
		while(++h!=f->facet_begin());
		
	}
	
	constructPolyhedron * fullSRMesh(new constructPolyhedron);
	polyhedron_builder<constructPolyhedron::HalfedgeDS> builderSR(coords,tris);
	fullSRMesh->delegate(builderSR);


	// fill holes that are not part of the main mesh borders
	std::vector<constructPolyhedron::Facet_handle>  patch_facetsSR;
	std::vector<constructPolyhedron::Vertex_handle> patch_verticesSR;
	for(auto h = fullSRMesh->halfedges_begin(); h!=fullSRMesh->halfedges_end(); ++h)
	{
		auto p = h->vertex()->point();
		Point3d pp(p.x(),p.y(),p.z());
		if(h->is_border() )
		{
		CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(
			*fullSRMesh,
			h,
			std::back_inserter(patch_facetsSR),
			std::back_inserter(patch_verticesSR),
			CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, *fullSRMesh)).
			geom_traits(constructionK())); 
		}
	}
	
	PolyhedronPtr resultOut = convertToEnrichedPolyhedron(fullOutMesh);
	PolyhedronPtr resultSR = convertToEnrichedPolyhedron(fullSRMesh);
	
	isolatedVertices_remover<HalfedgeDS> ivr;
	resultOut->delegate(ivr);
	resultOut->compute_normals();
	
	resultSR->delegate(ivr);
	resultSR->compute_normals();
	
	resultOut->set_index_vertices();
	resultSR->set_index_vertices();
	
	
	for(unsigned u=0; u < v->getScenePtr()->get_nb_polyhedrons(); ++u )
	{
		v->getScenePtr()->setVisible(u,false);
	}
	
	v->getScenePtr()->add_polyhedron(resultOut);
	v->getScenePtr()->add_polyhedron(resultSR);
	v->getScenePtr()->todoIfModeSpace(v,0.0);
	
	tris.clear();
	coords.clear();
	for(auto p = resultOut->points_begin(); p!= resultOut->points_end(); ++p)
	{
		coords.push_back(p->x());
		coords.push_back(p->y());
		coords.push_back(p->z());
	}
	
	
	//extract CC from result Out	
	const int nb = resultOut->calc_nb_components(); // now each face is tagged with #CC
	resultOut->tag_facets(0);
	std::cout << "nb components : " << nb << std::endl;
	int c = 1;
	for(auto pFacet = resultOut->facets_begin();
		    pFacet != resultOut->facets_end();
			++pFacet)
	{
		
		if(pFacet->tag() == 0)
		{
			resultOut->tag_component(pFacet,0,c);
			c++;
		}
	}
	std::cout << std::endl;
	
	//resultOut->tag_facets(0);
	m_toMerge.push_back(resultSR);
	for(unsigned cc = 1; cc <= nb; ++cc)
	{
		//std::cout << "cc : " << cc << std::endl;
		for(auto pFacet = resultOut->facets_begin();
		    pFacet != resultOut->facets_end();
			++pFacet)
		{
			//std::cout << "facet tag : " << pFacet->tag() << std::endl;
			if(pFacet->tag() == cc)
			{
				Halfedge_around_facet_circulator h = pFacet->facet_begin();
				do
				{
					tris.push_back(h->vertex()->tag());
				}
				while(++h!=pFacet->facet_begin());
			}
		}
		//std::cout << "size of tris : " << tris.size() << std::endl;
		polyhedron_builder<HalfedgeDS> partBuilder(coords,tris);
		PolyhedronPtr mpart(new Polyhedron);
		mpart->delegate(partBuilder);
		isolatedVertices_remover<HalfedgeDS> partIsoRem;
		mpart->delegate(partIsoRem);
		mpart->compute_bounding_box();
		m_toMerge.push_back(mpart);
		std::cout << "size of m_toMerge : " << m_toMerge.size() << std::endl;
		tris.clear();
		std::cout << "creating outMesh cc " << cc << std::endl;
	}
	std::cout << "fin remesh New" << std::endl;
}


void softICPController::remesh(Viewer * v)
{	
	PolyhedronPtr outMesh = getMeshOutsideSR();
	for(auto p=outMesh->vertices_begin();p!=outMesh->vertices_end();++p)
	{
		p->color(1.0,0.0,0.0);
	}
	outMesh->normalize_border();
	outMesh->compute_normals();
	v->getScenePtr()->add_polyhedron(outMesh);
		
	// collect border points for surface Reconstruction filtering
	std::set<Point3d> border;
	for(auto h = outMesh->border_halfedges_begin(); h!=outMesh->halfedges_end();++h)
	{
		border.insert(h->vertex()->point());
	}
	
	PolyhedronPtr srMesh = buildSRMesh(border,outMesh);
	for(auto p=srMesh->vertices_begin();p!=srMesh->vertices_end();++p)
	{
		p->color(0.0,1.0,0.0);
	}
	outMesh->compute_normals();
	v->getScenePtr()->add_polyhedron(srMesh);
	
	
	v->getScenePtr()->setVisible(0,false);
	v->getScenePtr()->setVisible(1,false);
	
	srMesh->normalize_border();
	/*for(auto h = srMesh->border_halfedges_begin(); h!=srMesh->halfedges_end();++h)
	{
		border.erase(h->vertex()->point());
	}*/
	
	PolyhedronPtr completeMesh = stitchAndSmooth(outMesh,srMesh,border);
	
	completeMesh->compute_normals();
	v->getScenePtr()->add_polyhedron(completeMesh);
	v->getScenePtr()->todoIfModeSpace(v,0.0);
}

PolyhedronPtr softICPController::getMeshOutsideSR()
{
	PolyhedronPtr OutMesh(new Polyhedron);
	
	std::vector<double> coords;
	std::vector<int> tris;
	PolyhedronPtr p = m_polyhedron1;
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
			if(m_isSnappingRegion[v]||m_isFactoredRegion[v]){inside++;}
			vertsIds[vId] = v->tag();
			vId++;
		}
		while(++hC!=pFacet->facet_begin());
		if(inside<3)
		{
			tris.push_back(vertsIds[0]);tris.push_back(vertsIds[1]);tris.push_back(vertsIds[2]);
		}
	}
	
	int vertOffset = coords.size()/3;
	p = m_polyhedron2;
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
			if(m_isSnappingRegion[v]||m_isFactoredRegion[v]){inside++;}
			vertsIds[vId] = v->tag() + vertOffset ;
			vId++;
		}
		while(++hC!=pFacet->facet_begin());
		if(inside<3)
		{
			tris.push_back(vertsIds[0]);tris.push_back(vertsIds[1]);tris.push_back(vertsIds[2]);
		}
	}
	
	polyhedron_builder<HalfedgeDS> builder(coords,tris);
	OutMesh->delegate(builder);
	OutMesh->normalize_border();
	return OutMesh;
}

PolyhedronPtr softICPController::getMeshInsideSR(PolyhedronPtr p)
{
	PolyhedronPtr InMesh(new Polyhedron);
	
	std::vector<double> coords;
	std::vector<int> tris;
	
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
			if(m_isSnappingRegion[v]||m_isFactoredRegion[v]){inside++;}
			vertsIds[vId] = v->tag();
			vId++;
		}
		while(++hC!=pFacet->facet_begin());
		if(inside>=2)
		{
			tris.push_back(vertsIds[0]);tris.push_back(vertsIds[1]);tris.push_back(vertsIds[2]);
		}
	}
	
	polyhedron_builder<HalfedgeDS> builder(coords,tris);
	isolatedVertices_remover<HalfedgeDS> isolatedRemover;
	InMesh->delegate(builder);
	InMesh->delegate(isolatedRemover);
	
	return InMesh;
}

double softICPController::centerDistance(Facet_handle f)
{
	double center_dist = 0.0;
	auto h = f->facet_begin();
	do {
		center_dist+=m_distToLoop[h->vertex()];
	} while ( ++h != f->facet_begin());
	center_dist/=3.0;
	return center_dist;
}

PolyhedronPtr softICPController::buildSRMesh(std::set<Point3d> & border, PolyhedronPtr outsideMesh)
{
	std::vector<double> coords;
	std::vector<int> tris;
	PolyhedronPtr m1 = getMeshInsideSR(m_polyhedron1);
	PolyhedronPtr m2 = getMeshInsideSR(m_polyhedron2);
	std::list<delaunayKernel::Point_3> L;
	
	std::vector<Point3d> pts;
	std::vector<Point3d> temp;
	
	for(auto p = m1->vertices_begin(); p!=m1->vertices_end();++p)
	{
		if(border.count(p->point())==0)
		{ 
			pts.push_back(p->point());
		}
		else
		{
			temp.push_back(p->point());
		}
	}
	
	for(auto p = m2->vertices_begin(); p!=m2->vertices_end();++p)
	{
		if(border.count(p->point())==0)
		{
			pts.push_back(p->point());
		}
		else
		{
			temp.push_back(p->point());
		}
	}
	
	
	for(int s=0;s<2;++s)
	{
		smoothPointSet(pts,8);
	}
	
	for(auto p = temp.begin(); p!=temp.end(); ++p)
	{
		pts.push_back(*p);
	}
	
	/*for(auto p = border.begin(); p!= border.end(); ++p)
	{
		pts.push_back(*p);
	}*/
	
	for(auto p = pts.begin(); p!=pts.end(); ++p)
	{
		L.push_front(delaunayKernel::Point_3(p->x(),p->y(),p->z()));
	}
		
	std::vector<RFacet> facets;
	AABB_Tree tree(outsideMesh->facets_begin(),outsideMesh->facets_end());
	tree.build();

	DistanceSurface pdistSurf(tree,border,m_R);
	CGAL::advancing_front_surface_reconstruction(L.begin(),L.end(),std::back_inserter(facets),pdistSurf,1.0,0.0005);

	//CGAL::advancing_front_surface_reconstruction(L.begin(),L.end(),std::back_inserter(facets),pdistSurf);
	for(auto v = L.begin(); v!=L.end();++v)
	{
		auto p = (*v);
		coords.push_back(p.x());coords.push_back(p.y());coords.push_back(p.z());
	}
	
	for(auto f = facets.begin(); f!=facets.end(); ++f)
	{
		tris.push_back((*f)[0]);tris.push_back((*f)[1]);tris.push_back((*f)[2]);
	}
	
	PolyhedronPtr reconstructedMesh(new Polyhedron);
	polyhedron_builder<HalfedgeDS> builder(coords,tris);
	reconstructedMesh->delegate(builder);
	//reconstructedMesh->keep_largest_connected_components(1);
	reconstructedMesh->compute_normals();
	
	return reconstructedMesh;
}



PolyhedronPtr softICPController::stitchAndSmooth(PolyhedronPtr outMesh, PolyhedronPtr inMesh, std::set<Point3d> & borders)
{
	// merge w/polyhedronbuilder
	std::vector<double> coords;
	std::vector<int> tris;
	
	outMesh->set_index_vertices();
	inMesh->set_index_vertices();
	
	std::vector<Point3d> newP;
	std::set<Point3d> sr;
	sr.insert(inMesh->points_begin(),inMesh->points_end());
	
	std::set<Point3d> inMeshBorder;
	
	
	
	for(auto po = outMesh->points_begin(); po != outMesh->points_end(); ++po)
	{
		coords.push_back(po->x());
		coords.push_back(po->y());
		coords.push_back(po->z());
	
	}
	for(auto po = inMesh->points_begin(); po != inMesh->points_end(); ++po)
	{
		coords.push_back(po->x());
		coords.push_back(po->y());
		coords.push_back(po->z());
	
	}
	for(auto it = newP.begin(); it!= newP.end(); ++it)
	{
		coords.push_back(it->x());
		coords.push_back(it->y());
		coords.push_back(it->z());
	}
	
	int vertOffset = outMesh->size_of_vertices();
	for(auto f = outMesh->facets_begin(); f!=outMesh->facets_end(); ++f)
	{
		
		auto h = f->facet_begin();
		do
		{
			tris.push_back(h->vertex()->tag());
		}
		while(++h!=f->facet_begin());
		
	}
	
	for(auto f = inMesh->facets_begin(); f!=inMesh->facets_end(); ++f)
	{
		
		auto h = f->facet_begin();
		do
		{
			tris.push_back(h->vertex()->tag()+vertOffset);
		}
		while(++h!=f->facet_begin());
	}
	
	constructPolyhedron * stitchedMesh(new constructPolyhedron);
	polyhedron_builder<constructPolyhedron::HalfedgeDS> builder(coords,tris);
	stitchedMesh->delegate(builder);

	// stitch non-connected vertices
	//CGAL::Polygon_mesh_processing::stitch_borders(*stitchedMesh);
	
	// fill holes that are not part of the main mesh borders
	std::vector<constructPolyhedron::Facet_handle>  patch_facets;
	std::vector<constructPolyhedron::Vertex_handle> patch_vertices;
	for(auto h = stitchedMesh->halfedges_begin(); h!=stitchedMesh->halfedges_end(); ++h)
	{
		auto p = h->vertex()->point();
		Point3d pp(p.x(),p.y(),p.z());
		if(h->is_border() && borders.count(pp)==0){
		CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(
			*stitchedMesh,
			h,
			std::back_inserter(patch_facets),
			std::back_inserter(patch_vertices),
			CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, *stitchedMesh)).
			geom_traits(constructionK())); 
		}
	}
	
	CGAL::Polygon_mesh_processing::stitch_borders(*stitchedMesh);
	
	// Return result as Enriched Polyhedron
	PolyhedronPtr result = convertToEnrichedPolyhedron(stitchedMesh);
	
	isolatedVertices_remover<HalfedgeDS> ivr;
	result->delegate(ivr);
	
	// Laplacian smoothing
	int nbSmooth = 0;
	/*while(nbSmooth < 10)
	{
		for(auto p = inMesh->vertices_begin(); p!=inMesh->vertices_end();++p)
		{
			if(sr.count(p->point()) && !borders.count(p->point())){
			double sx = 0.0;
			double sy = 0.0;
			double sz = 0.0;
			Halfedge_around_vertex_circulator hC = p->vertex_begin();
			do
			{
				Point3d n = hC->opposite()->vertex()->point();
				sx = sx + n.x();
				sy = sy + n.y();
				sz = sz + n.z();
				
				
			}while(++hC != p->vertex_begin());
			sx/=p->vertex_degree();
			sy/=p->vertex_degree();
			sz/=p->vertex_degree();
			Point3d smoothed(sx,sy,sz);
			newP.push_back(smoothed);
			
			}
		}
		int i =0;
		for(auto p = inMesh->vertices_begin(); p!=inMesh->vertices_end();++p)
		{
			if(sr.count(p->point()) && !borders.count(p->point()))
			{
				sr.erase(p->point());
				p->point() = newP[i];
				sr.insert(newP[i]);
				i++;
			}
		}
		nbSmooth++;
	}*/
	
	//if (!CGAL::Polygon_mesh_processing::is_outward_oriented(*result)	)
	//CGAL::Polygon_mesh_processing::reverse_face_orientations(*result);
	
	result->compute_normals();
	return result;
}

void softICPController::smoothPointSet(std::vector<Point3d> & points,const unsigned k)
{
	sTree tree(points.begin(),points.end());
	
	for(unsigned i=0;i<points.size();++i)
	{
		Neighbor_search tsearch(tree, points[i], k);
		
		double x=0.0;
		double y=0.0;
		double z=0.0;
		for(Neighbor_search::iterator it = tsearch.begin(); it != tsearch.end(); ++it)
		{
			x+=it->first.x();
			y+=it->first.y();
			z+=it->first.z();
		}
		x/=k; y/=k; z/=k;
		points[i] = Point3d(x,y,z);
	}
}

std::list<PolyhedronPtr> & softICPController::getPolyhedronToMerge()
{
	return m_toMerge;
}


double computePhiDistance(Vertex_handle v1, Vertex_handle v2, double w1, double w2, double w3)
{
	double dist = w1 * sqrt(CGAL::squared_distance(v1->point(),v2->point()));
	dist += w2*acos(v1->normal()*v2->normal());
	dist += w3*L2Dist(v1->getSemantic(),v2->getSemantic());
	return dist;
}

void initGeodesicMesh(PolyhedronPtr p, geodesic::Mesh * g)
{
// 	//g = new geodesic::Mesh;
	std::vector<double> points;
	std::vector<int> faces;
	
	p->set_index_vertices();

	for(Vertex_iterator pVertex = p->vertices_begin();
		pVertex != p->vertices_end();
		++pVertex)
	{
		points.push_back(pVertex->point().x()); // add the vertex to the list of vertices
		points.push_back(pVertex->point().y());
		points.push_back(pVertex->point().z());
	}

	for( Facet_iterator pFacet = p->facets_begin();
		pFacet!=p->facets_end();
		++pFacet)
	{
		// Use halfedge to visit vertices around the face
		Halfedge_around_facet_circulator  pHalfedge = pFacet->facet_begin();
		do
		{
			faces.push_back(pHalfedge->vertex()->tag());
		}
		while(++pHalfedge != pFacet->facet_begin());
	}
	
	// initialize geodesic graph and algorithm
        g->initialize_mesh_data(points,faces);
	std::cout << "initialized geodesic mesh" << std::endl;
}



bool sphere_clip_vectorG(Point3d &O, double r,const Point3d &P, Vector &V)
    {
        Vector W = P - O ;
        double a = (V*V);
        double b = 2.0 * V * W ;
        double c = (W*W) - r*r ;
        double delta = b*b - 4*a*c ;
        if (delta < 0) {
            // Should not happen, but happens sometimes (numerical precision)

            return true ;
        }

        double t = (- b + ::sqrt(delta)) / (2.0 * a) ;
        if (t < 0.0) {
            // Should not happen, but happens sometimes (numerical precision)

            return true ;
        }

        if (t >= 1.0) {
            // Inside the sphere

            return false ;
        }

        V=V*t;

        return true ;
    }