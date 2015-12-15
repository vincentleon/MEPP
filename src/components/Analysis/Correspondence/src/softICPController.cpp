#include "softICPController.h"
#include <CGAL/Kernel/interface_macros.h>
#include <CGAL/Sphere_3.h>
#include "time.h"


softICPController::softICPController(PolyhedronPtr m1, PolyhedronPtr m2) : m_polyhedron1(m1), m_polyhedron2(m2)
{

}

void softICPController::computeSnappingRegionCorrespondence()
{
	/*for every vertex in the snapping region, find its correspondence 
	 * that minimizes the weighted sum
	 */
	std::set<Halfedge_handle> & s1 = m_treeStructure1.m_halfedges;
	std::set<Halfedge_handle> & s2 = m_treeStructure2.m_halfedges;
	
	for(auto it = s1.begin(); it!=s1.end();++it)
	{
		double distMin = std::numeric_limits<double>::max();
		Halfedge_handle bestCorres;
		for(auto c = s2.begin(); c!=s2.end(); ++c)
		{
			double dist = computePhiDistance((*it)->vertex(),(*c)->vertex());
			if(dist < distMin)
			{
				bestCorres = *c;
				distMin = dist;
			}
		}
	}
	
	
}


void softICPController::buildTreeStructure(const int sizeOfTree)
{
	/*	First, identify the snapping region 
	 *	on the two meshes.
	 */
	getSnappingRegion();
	
//#pragma omp parallel sections num_threads(2)
//{	
	kMedoidMain(&m_treeStructure1,4);
	colorCluster(&m_treeStructure1);
//#pragma omp section
	kMedoidMain(&m_treeStructure2,4);
	colorCluster(&m_treeStructure2);
//}
	
	/*deformationNode * base = &m_treeStructure1;
	std::cout << "base size : " << base->m_childrenNodes.size() << std::endl;
	for(unsigned level=0;level<sizeOfTree-1;++level)
	{
		std::cout << " level : " << level << std::endl;
		for(unsigned c=0; c<base->m_childrenNodes.size() ;++c)
		{
			std::cout << "\t children : " << c << std::endl;
			// cluster nodes in lower level
			kMedoidMain(base->m_childrenNodes[c],4);
			//colorCluster(base->m_childrenNodes[c]);
		}
	}
	std::cout << std::endl;
	base = &m_treeStructure2;
	for(unsigned level=0;level<sizeOfTree-1;++level)
	{
		for(unsigned c=0; c<base->m_childrenNodes.size() ;++c)
		{
			// cluster nodes in lower level
			kMedoidMain(base->m_childrenNodes[c],4);
			//colorCluster(base->m_childrenNodes[c]);
		}
	}
	std::cout << std::endl;*/
	
	
	/*	Compute the correspondence 
	 * 	Phi between the two snapping regions
	 */
	computeSnappingRegionCorrespondence();
}

/*void softICPController::hierarchicalBuild(deformationNode* root, const int sizeOfTree, int level)
{
	if(level != sizeOfTree)
	{
		for(unsigned c=0;c<root->m_childrenNodes.size();++c)
		{
			hierarchicalBuild(root->m_childrenNodes[c],sizeOfTree,level+1);
		}
	}
	kMedoidMain(root,4);
}*/



void softICPController::colorCluster(deformationNode* root)
{
	if(root->m_childrenNodes.size() != 4)
	{
		std::cout << " cluster number != 4 " << std::endl;
		return;
	}
	for(auto h = root->m_childrenNodes[0]->m_halfedges.begin(); 
	    h != root->m_childrenNodes[0]->m_halfedges.end();
		++h)
	{
		(*h)->vertex()->color(1,0,0);
	}
	for(auto h = root->m_childrenNodes[1]->m_halfedges.begin(); 
	    h != root->m_childrenNodes[1]->m_halfedges.end();
		++h)
	{
		(*h)->vertex()->color(1,1,0);
	}
	for(auto h = root->m_childrenNodes[2]->m_halfedges.begin(); 
	    h != root->m_childrenNodes[2]->m_halfedges.end();
		++h)
	{
		(*h)->vertex()->color(0,0,1);
	}
	for(auto h = root->m_childrenNodes[3]->m_halfedges.begin(); 
	    h != root->m_childrenNodes[3]->m_halfedges.end();
		++h)
	{
		(*h)->vertex()->color(0,1,0);
	}
}


void softICPController::getSnappingRegion()
{
	double R = 0.1; //Region Size
	
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
			distMin = distMinToP1;
		}
	}
	// Border 1 and border2 are on the two closest boundary loops
	
	// Define the snapping region
	Vertex_handle v1 = border1->vertex();
	semanticDescr & sd1 = v1->getSemantic();
	Vertex_handle v2 = border2->vertex();
	semanticDescr & sd2 = v2->getSemantic();
	
	double squared_euclidean_radius = 0.05;

	
	Point3d center = border2->vertex()->point();
	for(auto pHalfedge = m_polyhedron1->halfedges_begin();pHalfedge!=m_polyhedron1->halfedges_end();++pHalfedge)
	{
 		Point3d p = pHalfedge->vertex()->point();
		if( CGAL::squared_distance(center,p) < squared_euclidean_radius)
		{	// point is candidate for semantic descriptor testing
			//pHalfedge->vertex()->color(0,1,0);
			semanticDescr & heDescr = pHalfedge->vertex()->getSemantic();
			
			Halfedge_iterator it = border1;
			
			do{
				//it->vertex()->color(1,1,0);
				semanticDescr & borderDescr = it->vertex()->getSemantic();
				//std::cout << heDescr.size();
				//std::cout << borderDescr.size() << std::endl;
				
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
	center = border1->vertex()->point();
	for(auto pHalfedge = m_polyhedron2->halfedges_begin();pHalfedge!=m_polyhedron2->halfedges_end();++pHalfedge)
	{
 		Point3d p = pHalfedge->vertex()->point();
		if( CGAL::squared_distance(center,p) < squared_euclidean_radius)
		{	// point is candidate for semantic descriptor testing
			//pHalfedge->vertex()->color(0,1,0);
			semanticDescr & heDescr = pHalfedge->vertex()->getSemantic();
			
			Halfedge_iterator it = border2;
			
			do{
				//it->vertex()->color(1,1,0);
				semanticDescr & borderDescr = it->vertex()->getSemantic();
				//std::cout << heDescr.size();
				//std::cout << borderDescr.size() << std::endl;
				
				if(L2Dist(heDescr,borderDescr)<R)
				{
					m_treeStructure2.m_halfedges.insert(pHalfedge);
					pHalfedge->vertex()->color(1,0,0);
				}
				it = it->next();
			}
			while(it!=border2);
		}
	}
}

void softICPController::cluster(deformationNode* root)
{
	for(auto h = root->m_halfedges.begin(); h!=root->m_halfedges.end();++h)
	{
		// For each point in the set, find the closest rep
		double distMin = std::numeric_limits<double>::max();
		for( int clust = 0; clust < root->m_childrenNodes.size(); ++clust)
		{
 			Halfedge_handle rep = root->m_childrenNodes[clust]->m_rep;
			double dist = computePhiDistance(rep->vertex(),(*h)->vertex(),0.0,0.0,1.0);
			if(dist<distMin)
			{
				distMin = dist;
				root->m_cluster[*h] = clust;
			}
		}
	}
}

void softICPController::updateRep(deformationNode* root)
{
	std::vector<Halfedge_handle> candidateRep;
	std::vector<double> distMin(root->m_childrenNodes.size(),std::numeric_limits<double>::max());
	
	candidateRep.resize(root->m_childrenNodes.size());
	
	for(auto h = root->m_halfedges.begin(); h!=root->m_halfedges.end();++h)
	{
		int clust = root->m_cluster[*h];
		/* consider that h is a candidate medoid
		*	Compute the sum of squared distances
		* 	to every point in the set
		* */
		double e = 0.0;
		for(auto g = root->m_halfedges.begin(); g!=root->m_halfedges.end();++g)
		{
			if(root->m_cluster[*g]!=clust){continue;}
			double dist = computePhiDistance((*h)->vertex(),(*g)->vertex(),0.0,0.0,1.0);
			e += dist*dist;
		}
		if(e<distMin[clust])
		{
			distMin[clust] = e;
			root->m_childrenNodes[clust]->m_rep = *h;
		}
	}
}

void softICPController::kMedoidMain(deformationNode * root, int k)
{
	srand(time(NULL));
	// Compute the first medoid 
	//updateRep(root);
	// Create the k clusters and randomly initialize medoids
	for(unsigned c=0; c<k; ++c)
	{
		root->m_childrenNodes.push_back(new deformationNode);
		int r = (rand() / RAND_MAX) * root->m_halfedges.size();
		auto it = root->m_halfedges.begin();
		for(unsigned v=0;v<r;++v){++it;}
		root->m_childrenNodes.back()->m_rep = *it;
	}
	// The kMedoidloop
	for(unsigned i=0;i<10;++i)
	{
		cluster(root);
		for(unsigned c=0;c<k;++c)
		{
			updateRep(root);
		}
	}
	// Finally put halfedges in the right patch (cluster)
	for(auto it = root->m_cluster.begin(); it!=root->m_cluster.end(); ++it)
	{
		root->m_childrenNodes[it->second]->m_halfedges.insert(it->first);
	}
	std::cout << "Display cluster info : " << std::endl;
	for(unsigned l=0;l<k;++l)
	{
		std::cout << "\t cluster #"<<l<<"  size : " << root->m_childrenNodes[l]->m_halfedges.size() << std::endl;
	}
}

double computePhiDistance(Vertex_handle v1, Vertex_handle v2, double w1, double w2, double w3)
{
	double dist = w1 * sqrt(CGAL::squared_distance(v1->point(),v2->point()));
	dist += w2*acos(v1->normal()*v2->normal());
	dist += w3*L2Dist(v1->getSemantic(),v2->getSemantic());
	return dist;
}
