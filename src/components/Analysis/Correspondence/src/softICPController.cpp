#include "softICPController.h"
#include <CGAL/Kernel/interface_macros.h>
#include <CGAL/Sphere_3.h>
#include "time.h"


softICPController::softICPController(PolyhedronPtr m1, PolyhedronPtr m2) : m_polyhedron1(m1), m_polyhedron2(m2)
{

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
			double dist = computePhiDistance((*it),(*c));
			if(dist < distMin)
			{
				bestCorres = *c;
				distMin = dist;
			}
		}
		m_Phi[(*it)] = bestCorres;
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
	//computeSnappingRegionCorrespondence();
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
					m_treeStructure1.m_vertices.insert(pHalfedge->vertex());
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
					m_treeStructure2.m_vertices.insert(pHalfedge->vertex());
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
	for(auto v = root->m_vertices.begin(); v!=root->m_vertices.end();++v)
	{
		// For each point in the set, find the closest rep
		double distMin = std::numeric_limits<double>::max();
		for( int clust = 0; clust < root->m_childrenNodes.size(); ++clust)
		{
 			Vertex_handle rep = root->m_childrenNodes[clust]->m_rep;
			double dist = computePhiDistance(rep,*v,0.0,0.0,1.0);
			if(dist<distMin)
			{
				distMin = dist;
				root->m_cluster[*v] = clust;
			}
		}
	}
}

void softICPController::updateRep(deformationNode* root)
{
	std::vector<Vertex_handle> candidateRep;
	std::vector<double> distMin(root->m_childrenNodes.size(),std::numeric_limits<double>::max());
	
	candidateRep.resize(root->m_childrenNodes.size());
	
	for(auto v = root->m_vertices.begin(); v!=root->m_vertices.end();++v)
	{
		int clust = root->m_cluster[*v];
		/* consider that h is a candidate medoid
		*	Compute the sum of squared distances
		* 	to every point in the set
		* */
		double e = 0.0;
		for(auto g = root->m_vertices.begin(); g!=root->m_vertices.end();++g)
		{
			if(root->m_cluster[*g]!=clust){continue;}
			double dist = computePhiDistance((*v),(*g),0.0,0.0,1.0);
			e += dist*dist;
		}
		if(e<distMin[clust])
		{
			distMin[clust] = e;
			root->m_childrenNodes[clust]->m_rep = *v;
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
		root->m_childrenNodes.back()->m_parent = root;
		root->m_childrenNodes.back()->m_clusterID = c;
		int r = (rand() / RAND_MAX) * root->m_vertices.size();
		auto it = root->m_vertices.begin();
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
		root->m_childrenNodes[it->second]->m_vertices.insert(it->first);
	}
	std::cout << "Display cluster info : " << std::endl;
	for(unsigned l=0;l<k;++l)
	{
		std::cout << "\t cluster #"<<l<<"  size : " << root->m_childrenNodes[l]->m_vertices.size() << std::endl;
	}
	// compute distances between the rep of all clusters
	computeMatrixDistance(root);
}

void softICPController::computeMatrixDistance(deformationNode* root)
{
	int nbcluster = root->m_childrenNodes.size();
	
	
	myMatrix & m  = root->m_distanceMatrix;
	m(nbcluster,nbcluster);
	
	for(unsigned i=0; i<nbcluster;++i)
	{
		Vertex_handle rep1 = root->m_childrenNodes[i]->m_rep;
		for(unsigned j=0;j<i;++j)
		{
			Vertex_handle rep2 = root->m_childrenNodes[j]->m_rep;
			double dist = computePhiDistance(rep1,rep2,0.0,0.0,1.0);
			m[i][j] = dist;
			m[j][i] = dist;
		}
	}
}


void softICPController::snapRegions(double R, unsigned elasticity)
{
	bool convergenceCriterion = false;
	bool order = false;
	int iter = 0;
	const int itermax = 5;
	
	while(!convergenceCriterion)
	{
		// swap meshes at each iteration
		order = !order; 
		iter++;
		PolyhedronPtr ma;
		PolyhedronPtr mb;
		
		convergenceCriterion = (iter == itermax);
		
		if(order)
		{
			ma = m_polyhedron1;
			mb = m_polyhedron2;
		}
		else
		{
			mb = m_polyhedron1;
			ma = m_polyhedron2;
		}
		
		// find correspondence between the 2 snapping regions
		computeSnappingRegionCorrespondence(order);
		// for each point in ma
		for(auto pVertex = ma->vertices_begin();
		    pVertex!=ma->vertices_end();
			++pVertex)
		{
 			std::set<Vertex_handle> N = getNeighborhood(pVertex,R,iter,elasticity,order);
			std::set<Vertex_handle> phiN = getCorrespondingNeighborhood(N);
			pointTransformation ti = computeTransformation(N,phiN);
			// apply only part of the transformation
			applyTransformation(pVertex,ti,iter,itermax);
		}
	}
}

void softICPController::applyTransformation(Vertex_handle p, pointTransformation & ti, int iter, int itermax)
{
	double li = iter/itermax;
	ti.T = li*ti.T;
	ti.Q.setAxisAngle(ti.Q.axis(),li*ti.Q.angle());
	Point3d pos = p->point();
	qglviewer::Vec V = ti.Q * qglviewer::Vec(pos.x(),pos.y(),pos.z()) + ti.T;
	p->point() = CGAL::ORIGIN + Vector(V[0],V[1],V[2]);
}

set< Vertex_handle > softICPController::getNeighborhood(Vertex_handle p, double R, unsigned int iter, unsigned int elasticity, bool order)
{
	double distToLoop = 0;
	
	// compute the size of the local neighborhood
	double expo = (iter*elasticity)/distToLoop;
	double radius = exp(-expo*expo);
	
	// find the leaf node that contains p
	deformationNode * node;
	deformationNode * containsP;
	if(order) {node = &m_treeStructure1;}
	else{node = &m_treeStructure2;}
	while(node->m_childrenNodes.size()!=0)
	{
		int clust = node->m_cluster[p];
		node = node->m_childrenNodes[clust];
	}
	containsP = node->m_parent; // this is the leaf node that contains p
	
	std::set<Vertex_handle> N;
	
	while(true)
	{
		
		int nbClust = node->m_childrenNodes.size();
		node = containsP->m_parent; // parent node of the one that contains p
		
		if(node == NULL)
		{
			N.insert(node->m_vertices.begin(),node->m_vertices.end());
			break;
		} // use the full snapping region as neigborhood if root == m_treeStructureX
		
		// compute the maximum distance to containsP
		double maxDist = 0.0;
		for(unsigned k=0;k<nbClust;++k)
		{
			double dist = node->m_distanceMatrix[k][containsP->m_clusterID];
		}
		if( radius < maxDist ) // collect all the patches whose distance is smaller than radius
		{
			for(unsigned k=0;k<nbClust;++k)
			{
				double dist = node->m_distanceMatrix[k][containsP->m_clusterID];
				if(dist < radius)
				{
					N.insert(node->m_childrenNodes[k]->m_vertices.begin(),node->m_childrenNodes[k]->m_vertices.end());
				}
			}	
			break;
		}
		else // go up one level in the hierarchy and repeat the process
		{
			node = node->m_parent;
			containsP = containsP->m_parent;
		}
	}
	
	return N;
}

std::set<Vertex_handle> softICPController::getCorrespondingNeighborhood( std::set<Vertex_handle> & N)
{
	std::set<Vertex_handle> cN;
	for(auto it = N.begin(); it!= N.end(); ++it)
	{
		cN.insert(m_Phi[*it]);
	}
	return cN;
}

pointTransformation softICPController::computeTransformation(std::set<Vertex_handle> & N, std::set<Vertex_handle> & phiN)
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
#pragma omp parallel for private(i,itt) default(none) shared(p_m,p_t,N,phiN)
	for(auto it=N.begin();it!=N.end();++it,++itt,i++)
	{
		Point3d np = (*it)->point();
		p_t.val[i][0] = np.x(); mut0 +=p_t.val[i][0];
		p_t.val[i][1] = np.y(); mut1 +=p_t.val[i][1];
		p_t.val[i][2] = np.z(); mut2 +=p_t.val[i][2];
		
		// get nearest point
		Point3d npp = (*itt)->point();
		p_m.val[i][0] = npp.x(); mum0 += p_m.val[i][0];
		p_m.val[i][1] = npp.y(); mum1 += p_m.val[i][0];
		p_m.val[i][2] = npp.z(); mum2 += p_m.val[i][0];
	}
	
	mu_m = mu_m/(double)N.size();
	mu_t = mu_t/(double)N.size();
	
	// substract mean (translation)
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


double computePhiDistance(Vertex_handle v1, Vertex_handle v2, double w1, double w2, double w3)
{
	double dist = w1 * sqrt(CGAL::squared_distance(v1->point(),v2->point()));
	dist += w2*acos(v1->normal()*v2->normal());
	dist += w3*L2Dist(v1->getSemantic(),v2->getSemantic());
	return dist;
}
