///////////////////////////////////////////////////////////////////////////
// Author: Vincent Leon
// Year: 2015
// Universite Lille 1, CRIStAL
///////////////////////////////////////////////////////////////////////////
#include <mepp_config.h>
#ifdef BUILD_component_Correspondence

#include <iostream>
#include <sstream>
#include <algorithm>
#include <stdlib.h>
#include <stack>
	

#include "../../../../mepp/mepp_component.h"


#include "Correspondence_Component.h"
#include "Correspondence_Polyhedron.h"
#include "geodesic/geodesic_algorithm_exact.h"


#define Malloc(type,n) (type *)malloc((n)*sizeof(type))


Correspondence_Component::Correspondence_Component(Viewer* v, PolyhedronPtr p) : mepp_component(v,p), m_segCtr(p)
{
	componentName = "Correspondence_Component";
	init = 1;
	p->set_index_vertices();
}

void Correspondence_Component::initParameters(int nbLabel, int meshId,std::string meshDir)
{
	m_Shape.m_meshID = meshId;
	m_nbLabel = nbLabel;
	m_Shape.initFaceLabelsAndSegments(meshDir);
}


void Correspondence_Component::initGeodesicMesh(PolyhedronPtr p)
{
	//structure for geodesic computation
	std::vector<double> points;
	std::vector<unsigned> faces;
	std::vector<Vertex_handle> alreadyvisited;
	
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
        m_gmesh.initialize_mesh_data(points,faces);
        m_geoAlg = new geodesic::GeodesicAlgorithmExact(&m_gmesh);
}

void Correspondence_Component::learnDescriptor(PolyhedronPtr p,std::string meshDir)
{
	m_Shape.initFaceLabelsAndSegments(meshDir);
	std::cout << "initFaceLabelsAndSegments" << std::endl;
	computeDescriptorAllVertices(p);
	std::cout << "computeDescriptorAllVertices" << std::endl;
	saveDescriptor(p,meshDir);
}


void Correspondence_Component::saveDescriptor(PolyhedronPtr p,std::string meshDir)
{
	std::ofstream file;
	std::stringstream ss;
	
	ss<<meshDir<<"/"<<m_Shape.m_meshID<<".semantic";
	std::cout << "saving file : " << ss.str() << std::endl;
	file.open(ss.str().c_str());
	for(Vertex_iterator pVertex = p->vertices_begin();
	    pVertex!=p->vertices_end(); ++pVertex)
	{
		const std::vector<double> & localDescr = pVertex->getSemantic();
		for(unsigned l=0;l<localDescr.size();++l)
		{
			file<<localDescr[l];
			if(l!=localDescr.size()-1)
			{
				file<<" ";
			}
		}
		file<<"\n";
	}
	file.close();
}

bool Correspondence_Component::readDescriptor(PolyhedronPtr p,std::string meshDir, bool normalize)
{
	Point3d bb = Point3d(p->xmin(),p->ymin(),p->zmin());
	std::ifstream file;
	std::stringstream ss;
	ss<<meshDir<<"/"<<m_Shape.m_meshID<<".semantic";
	file.open(ss.str().c_str());
	Vertex_iterator pVertex = p->vertices_begin();
	if(file)
	{
		std::string line;
		//int i = 0;
		while(getline(file,line))
		{
			std::vector<double> localDescr;
			std::istringstream iss(line);
			do
			{
				std::string fs;
				iss >> fs;
				if(fs == "" || fs == " "){continue;}
				localDescr.push_back(atof(fs.c_str()));
			} while(iss);


			Point3d pt = pVertex->point();

			//localDescr.push_back(std::sqrt(CGAL::squared_distance(p,bb)));
			/*localDescr.push_back( std::abs(pt.x() - p->xmin()) );
			localDescr.push_back( std::abs(pt.y() - p->ymin()) );
			localDescr.push_back( std::abs(pt.z() - p->zmin()) );*/
			
			pVertex->setSemantic(localDescr);
			++pVertex;
		}
	}
	else
	{
		std::cout << "Impossible d'ouvrir le fichier "<<ss.str()<<" !\n";
		return false;
	}
	if(normalize)
	{
		initMaxVector(p);
		for(Vertex_iterator pVertex = p->vertices_begin();
			pVertex!=p->vertices_end();++pVertex)
		{
			std::vector<double> descr = pVertex->getSemantic();
			this->normalize(descr);
			pVertex->setSemantic(descr);
		}
	}
	return true;
}

void Correspondence_Component::showDescriptor(PolyhedronPtr p, int dim)
{
	//initMaxVector(p);
	for(Vertex_iterator pVertex = p->vertices_begin();
	    pVertex!=p->vertices_end();++pVertex)
	{
		if(pVertex->getSemantic().size() == 0)
		{
			pVertex->color(0,0,0);
		}
		else{
			std::vector<double> localDescr = pVertex->getSemantic();
			//this->normalize(localDescr);
			//localDescr[dim] /= m_maxVector[dim];
			pVertex->color(localDescr[dim],0.5,1.0-localDescr[dim]);
		}
	}
}


void Correspondence_Component::initMaxVector(PolyhedronPtr p)
{
	if(m_maxVector.size()!=0){m_maxVector.clear();}
	for(unsigned l=0;l<m_nbLabel;++l)
	{
		double maxV = 0.0;
		for(Vertex_iterator pVertex = p->vertices_begin();
		pVertex!=p->vertices_end();
		++pVertex)
		{ 
			if(pVertex->getSemantic().size() == 0)
			{
				continue;
			}
			double val = pVertex->getSemantic()[l];
			if(val>maxV && val!=std::numeric_limits<double>::infinity()){maxV = val;}
		}
		m_maxVector.push_back(maxV);
	}
}

void Correspondence_Component::normalize(vector< double >& descr)
{
	for(unsigned i=0;i<descr.size();++i)
	{
		if(m_maxVector[i]!=std::numeric_limits<double>::infinity())
		{
			descr[i]/=m_maxVector[i];
		}
	}
}

void Correspondence_Component::computeGeodesicDistancesAllVertices(PolyhedronPtr p)
{
	for(unsigned label=0; label<m_nbLabel;++label)
	{
		std::vector<geodesic::SurfacePoint> sources;
		int facetC = 0;
		for(Facet_iterator pFacet = p->facets_begin();
			pFacet!=p->facets_end();++pFacet)
		{
			int faceLabel = m_Shape.m_faceLabels[facetC];
			if(faceLabel == label)
			{
				geodesic::Face sourceFace = (m_gmesh.faces()[facetC]);
				sources.push_back(sourceFace.adjacent_vertices()[0]);
			}
			facetC++;
		}

		if(sources.size() != 0)
		{
			m_geoAlg->propagate(sources);
		}
		else
		{
			for(Vertex_iterator pVertex = p->vertices_begin();
			    pVertex!=p->vertices_end();++pVertex)
			{
				pVertex->pushSemantic(std::numeric_limits<double>::infinity());
			}

		}

		for(Vertex_iterator pVertex = p->vertices_begin();
		    pVertex!=p->vertices_end();++pVertex)
		{
			
			int pIndex = pVertex->tag();
			geodesic::SurfacePoint pt(&m_gmesh.vertices()[pIndex]);
			double distance;
			unsigned best_source = m_geoAlg->best_source(pt,distance);
			if(distance == geodesic::GEODESIC_INF)
			{
				pVertex->pushSemantic(std::numeric_limits<double>::infinity());
			}
			else
			{
				pVertex->pushSemantic(distance);
			}
		}
	}
}


void Correspondence_Component::computeDescriptorAllVertices(PolyhedronPtr p)
{
	this->initGeodesicMesh(p);
	std::cout << "initgeodesicMesh" << std::endl;
	this->computeGeodesicDistancesAllVertices(p);
	std::cout << "computeGeodesicDistancesAllVertices" << std::endl;
	
	this->initMaxVector(p);
	for(Vertex_iterator pVertex = p->vertices_begin();
		    pVertex!=p->vertices_end();++pVertex)
	{
		std::vector<double> descr = pVertex->getSemantic();
		this->normalize(descr);
		pVertex->setSemantic(descr);
	}
}

vector<double> & Correspondence_Component::getClosetVertexDescriptor(PolyhedronPtr p, Point3d pickedPoint)
{
	double distMin = std::numeric_limits<double>::max();
	Vertex_iterator cVertex;
	Facet_iterator cFacet;
	int i = 0;
	int min = 0;

		for(Vertex_iterator pVertex = p->vertices_begin();
			pVertex!=p->vertices_end();
			++pVertex)
		{
			Point3d pt = pVertex->point();
			double dist = pickedPoint.x() - pt.x();
			dist = dist*dist;
			dist += (pickedPoint.y() - pt.y())*(pickedPoint.y()-pt.y());
			dist += (pickedPoint.z() - pt.z())*(pickedPoint.z()-pt.z());
			dist = sqrt(dist);
			if(dist < distMin)
			{
				min = i;
				cVertex = pVertex;
				distMin = dist;
			}
			i++;
		}
	std::vector<double> & localDescr = cVertex->getSemantic();
	//std::cout << "closest has been found" << distMin << std::endl;
	return localDescr;
}

void Correspondence_Component::compareDescriptorWithSVM(PolyhedronPtr p, unsigned SVM_mode)
{
	myVector mu = myVector(m_nbLabel);
	myVector sig = myVector(m_nbLabel);
	
	unsigned resSize = 0;
	
	double radius = m_patchRadius;
	
	std::cout << "radius : " << radius << "nbCandidates : "<< m_nbCandidates << std::endl;
	
	std::vector<std::pair<Vertex_handle,double> > candidateVertices;
	
	for(Vertex_iterator pVertex = p->vertices_begin();
	    pVertex!=p->vertices_end();++pVertex)
	{
		std::vector<double> & localdescr = pVertex->getSemantic();
		double dist = L2Dist(localdescr,m_centreDescriptor);
		candidateVertices.push_back(std::make_pair(pVertex,dist));
	}
	std::sort(candidateVertices.begin(), candidateVertices.end(), VertDistOrdering());

	// Pour chacun des n candidats, construire un patch et classifier les sommets de ce patch
	//std::map<Vertex_handle, double> labelMap;
	std::vector<Vertex_handle> selVerts;
	
	for(unsigned i = 0; i < m_nbCandidates; ++i)
	{
		std::set<Vertex_handle> patch;
		getInRadius(p,candidateVertices[i].first,radius,patch);
		
		
		std::vector<Vertex_handle> patchV(patch.begin(), patch.end()); 
		Enriched_kernel::Vector_3 U,V,N;
		
		if(SVM_mode != 2)
		{
			computePatchBasis(patchV,U,V,N);
		}
		
		unsigned nbFeatures = m_nbLabel + 3;
		if(SVM_mode == 1)
		{
			nbFeatures = 3;
		}
		else if (SVM_mode == 2)
		{
			nbFeatures = m_nbLabel;
		}
		
		for(unsigned f = 0; f < nbFeatures ; ++f)
		{
			std::cout << featureMeans[f] << std::endl;
		}
		
		for(auto it = patch.begin(); it!= patch.end(); ++it)
		{
			
			Vertex_handle pVertex = *it;
			//pVertex->color(0.0,1.0,0.0);
			
			std::vector<double> & descr = pVertex->getSemantic();
			struct svm_node * x = Malloc(struct svm_node,nbFeatures+1);
	
			if(SVM_mode!=1)
			{
				unsigned l = 0;
				for(l=0;l<m_nbLabel;++l)
				{
					x[l].index = l;
					x[l].value = descr[l];
				}
			}
			
			
			if (SVM_mode != 2)
			{
				Enriched_kernel::Vector_3 cv = candidateVertices[i].first->point()-pVertex->point();
				Enriched_kernel::Vector_3 projcv_n = (cv*N)*N;
				Enriched_kernel::Vector_3 projcv_uv = cv - projcv_n;
				
				double d = sqrt(cv.squared_length());
				double h = sqrt(projcv_n.squared_length());
				double normU = sqrt(U.squared_length());
				double normProjcv_uv = sqrt(projcv_uv.squared_length());
				double a = acos((U*projcv_uv)/(normU*normProjcv_uv + 0.0001)); 
			
				if(SVM_mode == 0)
				{
					x[m_nbLabel].index = m_nbLabel;
					x[m_nbLabel].value = d;
					x[m_nbLabel+1].index = m_nbLabel+1;
					x[m_nbLabel+1].value = h;
					x[m_nbLabel+2].index = m_nbLabel+2;
					x[m_nbLabel+2].value = a;
				}
				else if(SVM_mode == 1)
				{
					x[0].index = 0;
					x[0].value = d;
					x[1].index = 1;
					x[1].value = h;
					x[2].index = 2;
					x[2].value = a;
				}
			}

			
			x[nbFeatures].index = -1;
			
			for(unsigned f = 0; f < nbFeatures; ++f)
			{
				x[f].value = (x[f].value - featureMeans[f])/(2*featureSTD[f]+0.001);
			}
			
			double predictedLabel = svm_predict(m_svmModel,x);
			
			//free(x);
			
			if(predictedLabel == 1)
			{
				resSize++;
				pVertex->color(0.29,0.71,0.85);
				
				for(unsigned l=0;l<descr.size();++l)
				{
					mu[l]+=descr[l];
				}
				selVerts.push_back(pVertex);
			}
			else
			{
				pVertex->color(0.5,0.5,0.5);
			}
		}
		candidateVertices[i].first->color(0.0,0.0,0.0);
	}
	std::cout << "end test SVM " << std::endl;
	std::cout << " restSize : " << resSize << std::endl;
	mu = (1.0/resSize)*mu;
	
	for(unsigned vid = 0; vid < selVerts.size(); ++vid)
	{
		Vertex_handle vs = selVerts[vid];
		std::vector<double > & descr = vs->getSemantic();
		for(unsigned l=0;l<descr.size();++l)
		{	
			double diff = descr[l]-mu[l];
			sig[l]+= diff*diff;
		}
	}
	sig = (1.0/resSize)*sig;
	for(unsigned l=0;l<sig.dimension();++l)
	{
		sig[l] = sqrt(sig[l]);
	}
	
	std::cout << "mean : " << mu << std::endl;
	std::cout << "stdv " << sig << std::endl;
}

void Correspondence_Component::getInRadius(PolyhedronPtr p, Vertex_handle c, double radius, set< Vertex_handle >& vertices)
{
	Point3d O = c->point();
	std::stack<Vertex_handle> S;
	
	S.push(c);
	vertices.insert(c);
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
				bool isect = sphere_clip_vector(O,radius, P, V);
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
}


void Correspondence_Component::compareDescriptorToEllipse(PolyhedronPtr p)
{
	std::cout << "Compare Descriptor to Ellipse : " << p->pName << std::endl;
	std::vector<double> centreDescr = m_centreDescriptor;
	//normalize(centreDescr);
	myVector mu = myVector(centreDescr.size());
	myVector sig = myVector(centreDescr.size());
	
	unsigned resSize = 0;
	
	for(Vertex_iterator pVertex = p->vertices_begin();
		    pVertex!=p->vertices_end();++pVertex)
	{
		std::vector<double> localDescr = pVertex->getSemantic();
		//normalize(localDescr);
		double eqEll = 0.0;
		for(unsigned l=0;l<localDescr.size();++l)
		{
			double diff = localDescr[l]- centreDescr[l];
			eqEll += (diff*diff/(m_ellipse[l]*m_ellipse[l]));
		}
		if(eqEll<=1)
		{
			++resSize;
			pVertex->color(0.29,0.71,0.85);
			
			std::vector<double > & descr = pVertex->getSemantic();
			for(unsigned l=0;l<descr.size();++l)
			{
				mu[l]+=descr[l];
			}
		}
		else
		{
			pVertex->color(0.5,0.5,0.5);
		}
	}
	mu = (1.0/resSize)*mu;
	
	for(Vertex_iterator pVertex = p->vertices_begin();
		    pVertex!=p->vertices_end();++pVertex)
	{
		std::vector<double> localDescr = pVertex->getSemantic();
		//normalize(localDescr);
		double eqEll = 0.0;
		for(unsigned l=0;l<localDescr.size();++l)
		{
			double diff = localDescr[l]- centreDescr[l];
			eqEll += (diff*diff/(m_ellipse[l]*m_ellipse[l]));
		}
		if(eqEll<=1)
		{
			//++resSize;
			
			std::vector<double > & descr = pVertex->getSemantic();
			for(unsigned l=0;l<descr.size();++l)
			{
				double diff = descr[l]-mu[l];
				sig[l]+= diff*diff;
			}
		}
	}
	sig = (1.0/resSize)*sig;
	for(unsigned l=0;l<sig.dimension();++l)
	{
		sig[l] = sqrt(sig[l]);
	}
	
	
	std::cout << "mean : " << mu << std::endl;
	std::cout << "stdv " << sig << std::endl;
	
}

void Correspondence_Component::compareDescriptorToEllipseRotation(PolyhedronPtr p)
{
	myMatrix EllipseMat(m_nbLabel,m_nbLabel);
	for(unsigned i=0;i<m_nbLabel;++i)
	{
		for(unsigned j=0;j<m_nbLabel;++j)
		{
			EllipseMat[i][j] = m_ellipse[i*m_nbLabel+j];
		}
	}
	
	myVector mu(m_centreDescriptor.size());
	for(unsigned l=0;l<m_centreDescriptor.size();++l)
	{
		mu[l] = m_centreDescriptor[l];
	}
	
	for(Vertex_iterator pVertex = p->vertices_begin();
		    pVertex!=p->vertices_end();++pVertex)
	{
		std::vector<double> localDescr = pVertex->getSemantic();
		double eqEll = 0.0;
		myVector x(localDescr.size());
		for(unsigned l=0;l<localDescr.size();++l)
		{
 			x[l] = localDescr[l];
		}
			
		myVector distToMean = x -mu;
		myMatrix distToMeanT = Linear_algebraCd<double>::transpose(distToMean);
		
		myVector res = distToMeanT*(EllipseMat*distToMean);
		eqEll = res[0];
		
		if(eqEll<=1)
		{
			pVertex->color(0,0,0);
		}
		else
		{
			pVertex->color(0.5,0.5,0.5);
		}
	}
}


void Correspondence_Component::compareDescriptorToGaussian(PolyhedronPtr p)
{	
	myVector & mu = m_mu;	
	
	myMatrix I,D;
	myVector C;
	std::vector<int> vv;
	double dd;

	double & det = m_gaussianDeterminant;
	myMatrix & inverse = m_inverseGaussianMatrix;
	
	double power = pow(2*3.14159,m_nbLabel);
	double norm = 1.0/(sqrt(power)*det);
	
	std::vector<double> dist;
		
	for(Vertex_iterator pVertex = p->vertices_begin();
	    pVertex!=p->vertices_end();++pVertex)
	{
		std::vector<double> localDescr = pVertex->getSemantic();
		myVector x(localDescr.size());
		for(unsigned l=0;l<localDescr.size();++l)
		{
 			x[l] = localDescr[l];
		}
		
		myVector distToMean = x -mu;
		myMatrix distToMeanT = Linear_algebraCd<double>::transpose(distToMean);
		
		myVector mahalanobis = (distToMeanT * inverse * distToMean);
		
		dist.push_back(sqrt(mahalanobis[0]));
		//	std::cout << sqrt(mahalanobis[0]) << " " <<std::flush;
	}
	int v = 0;
	for(Vertex_iterator pVertex = p->vertices_begin();
	    pVertex!=p->vertices_end();++pVertex)
	{
 		float m = dist[v];
		if( m < m_threshold )
		{
			//pVertex->color(m/m_threshold,m/m_threshold,m/m_threshold);
			pVertex->color(0,0,0);
		}
		else
		{
			pVertex->color(1.0,1.0,1.0);
		}
		v++;
	}
}

void Correspondence_Component::compareDescriptorToSVM(PolyhedronPtr p)
{
	for(Vertex_iterator pVertex=p->vertices_begin();
	    pVertex!=p->vertices_end();++pVertex)
	{
		std::vector<double> & descr = pVertex->getSemantic();
		struct svm_node * x = Malloc(struct svm_node,m_nbLabel+1);
		for(unsigned l=0;l<m_nbLabel;++l)
		{
			x[l].index = l;
			x[l].value = descr[l];
		}
		x[m_nbLabel].index = -1;
		double predictedLabel = svm_predict(m_svmModel,x);
		free(x);
		if(predictedLabel == 1)
		{
			pVertex->color(0,0,0);
		}
		else
		{
			pVertex->color(1,1,1);
		}
	}
}

Vertex_handle Correspondence_Component::getSelectionCenter()
{
	double scoreMin = std::numeric_limits<double>::max();
	for(int i = 0; i<m_selection.size(); ++i)
	{
		double score = 0;
		Point3d p1 = m_selection[i]->point();

		for(int j=0;j<m_selection.size();++j)
		{
			double dist = L2Dist(m_selection[i]->getSemantic(),m_selection[j]->getSemantic());
			score+=dist*dist;
		}
		if(score<scoreMin)
		{
			scoreMin=score;
			m_centreSelection = m_selection[i];
			m_centreDescriptor = m_centreSelection->getSemantic();
		}
	}
	//m_centreSelection->color(0,0,1);
	return m_centreSelection;
}

Vertex_handle Correspondence_Component::getFurtherFromSelectionCenter()
{
	double scoreMax = 0;
	Vertex_handle furtherFromCenter;
	std::vector<double> & descrCenter = m_centreSelection->getSemantic();

	for(int i=0; i<m_selection.size();++i)
	{
		double score = 0.0;
		std::vector<double> & lDescr = m_selection[i]->getSemantic();
		double dist = L2Dist(descrCenter,lDescr);
		score = dist*dist;
		if( score > scoreMax )
		{
			scoreMax = score;
			furtherFromCenter = m_selection[i];
		}
	}
	//furtherFromCenter->color(0,0,1);
	return furtherFromCenter;
}

void Correspondence_Component::tagSelectionVertices(PolyhedronPtr p)
{
	std::list<Vertex_handle> q;
	for(unsigned i=0; i<m_selection.size();++i)
	{
		Vertex_handle v = m_selection[i];
		m_tag[v] = 1;
		q.push_back(v);
		v->color(1,0,0);
	}
	while(!q.empty())
	{
		Vertex_handle s = q.front(); q.pop_front();

		if(m_tag[s] == 4){break;}
		// For each neighbor of S
		Halfedge_around_vertex_circulator he = s->vertex_begin();
		do
		{
			Vertex_handle n = he->opposite()->vertex();
			if(m_tag[n]==0)
			{
				q.push_back(n);
				m_tag[n] = m_tag[s] + 1;
				m_selection.push_back(n);
			}

		}
		while(++he!=s->vertex_begin());
	}

	for(unsigned i=0; i<m_selection.size();++i)
	{
		Vertex_handle v = m_selection[i];
		if( m_tag[v] == 2 )
				{
			v->color(1,1,0);
				}
		else if (m_tag[v] == 3 )
		{
			v->color(0,1,1);
		}
	}
}

void Correspondence_Component::readSelectionBasedOnColor(PolyhedronPtr p)
{
	for(Vertex_iterator pVertex = p->vertices_begin();
	    pVertex!=p->vertices_end();
	    ++pVertex)
	    {
		int r = pVertex->color(0);
		int g = pVertex->color(1);
		int b = pVertex->color(2);
		
		if( r==1 && g==0 && b==0 )
		{
			m_selection.push_back(pVertex);
		}
	    }
	this->tagSelectionVertices(p);
}


void Correspondence_Component::initializeEllipsoid(PolyhedronPtr p)
{
	m_centreSelection = getSelectionCenter();
	Vertex_handle extremumSelection = getFurtherFromSelectionCenter();
	m_isolineValue = L2Dist(m_centreSelection->getSemantic(),extremumSelection->getSemantic());
	m_ellipse = std::vector<double>(m_nbLabel,m_isolineValue);
}

void Correspondence_Component::computeEllipseParameters(PolyhedronPtr p)
{
	std::vector<double> & descr = m_selection[0]->getSemantic();
	
	myVector mu = myVector(descr.size());
	myVector sig = myVector(descr.size());
	
	for(unsigned l=0;l<descr.size();++l)
	{
		mu[l] = descr[l];
 	}
	
	for(unsigned s=1;s<m_selection.size();++s)
	{
		descr = m_selection[s]->getSemantic();
		for(unsigned l=0;l<descr.size();++l)
		{
			mu[l]+=descr[l];
		}
	}
	mu = (1.0/m_selection.size())*mu;
	
	for(unsigned s=0;s<m_selection.size();++s)
	{
		descr = m_selection[s]->getSemantic();
		for(unsigned l=0;l<descr.size();++l)
		{
			double diff = descr[l] - mu[l];
			sig[l] += diff*diff;
		}
	}
	sig = (1.0/m_selection.size())*sig;
	for(unsigned l=0;l<sig.dimension();++l)
	{
		sig[l] = sqrt(sig[l]);
	}
	std::cout << mu << std::endl;
	std::cout << sig << std::endl;
	
	m_ellipse = std::vector<double>(m_nbLabel,m_isolineValue);
	
	std::vector<double> ell = m_ellipse;
	
	int dim = m_ellipse.size();
	
	nlopt::opt opt(nlopt::LN_COBYLA,dim);
	
	void * data = &(*this);
	
	opt.set_min_objective(objectiveFun,data);
	
	std::vector<double> lBounds(dim,0.0);
	opt.set_lower_bounds(lBounds);
	opt.set_xtol_rel(1e-3);
	double minf;
	
	nlopt::result res = opt.optimize(ell,minf);
	
	m_ellipse = ell;
}

void Correspondence_Component::computeEllipseParametersRotation(PolyhedronPtr p)
{
	m_ellipse = std::vector<double>(m_nbLabel*m_nbLabel,0.0);
	
	for(int i=0;i<m_nbLabel;++i)
	{
		m_ellipse[i*m_nbLabel+i] = 1/(m_isolineValue*m_isolineValue);
	}
		
	std::vector<double> ell = m_ellipse;
	
	int dim = m_ellipse.size();
	
	nlopt::opt opt(nlopt::LN_COBYLA,dim);
	
	void * data = &(*this);
	
	opt.set_min_objective(objectiveFun,data);
	
	std::vector<double> lBounds(dim,0.0);
	opt.set_lower_bounds(lBounds);
	opt.set_xtol_rel(1e-2);
	double minf;
	
	nlopt::result res = opt.optimize(ell,minf);
	
	m_ellipse = ell;
}

void Correspondence_Component::computeGaussianParameters(PolyhedronPtr p)
{	
	m_centreSelection = getSelectionCenter();
	
	Vertex_handle extremumSelection = getFurtherFromSelectionCenter();
	std::vector<double> extrDescr = extremumSelection->getSemantic();
	
	myMatrix sig(m_nbLabel,m_nbLabel);
	
	for(unsigned s=0;s<m_selection.size();++s)
	{
		Vertex_handle v = m_selection[s];
		std::vector<double> localDescr = v->getSemantic();
		for(unsigned i=0;i<m_nbLabel;++i)
		{
			double t1 = localDescr[i] - m_centreDescriptor[i];
			for(unsigned j=0;j<m_nbLabel;++j)
			{
				sig[i][j] += t1 * (localDescr[j] - m_centreDescriptor[j]);
			}
		}
	}
	m_gaussianMatrix = (1.0/(m_selection.size()-1)) * sig;
	
	/*for(Vertex_iterator pVertex = p->vertices_begin();pVertex!=p->vertices_end();++pVertex)
	{
		std::vector<double> localDescr = pVertex->getSemantic();
		for(unsigned i=0;i<m_nbLabel;++i)
		{
			double t1 = localDescr[i] - m_centreDescriptor[i];
			for(unsigned j=0;j<m_nbLabel;++j)
			{
				sig[i][j] += t1 * (localDescr[j] - m_centreDescriptor[j]);
			}
		}
	}
	m_gaussianMatrix = (1.0/(p->size_of_vertices()-1)) * sig;*/
	myMatrix I,D;	
	myVector C;
	std::vector<int> vv;
	double dd;
	
	
	m_inverseGaussianMatrix = Linear_algebraCd<double>::inverse(m_gaussianMatrix,dd);
	m_gaussianDeterminant = Linear_algebraCd<double>::determinant(m_gaussianMatrix,I,D,vv,C);

	m_mu = myVector(m_centreDescriptor.size());
	
	for(unsigned l=0;l<m_centreDescriptor.size();++l)
	{
		m_mu[l] = m_centreDescriptor[l];
	}
	
	// threshold optimization
	
	/*std::vector<double> threshold(1,0.0);
	nlopt::opt opt(nlopt::LN_COBYLA,1);
	void * data = &(*this);
	opt.set_min_objective(objectiveFunGaussian,data);
	std::vector<double> lBounds(1,0.0);
	opt.set_lower_bounds(lBounds);
	opt.set_xtol_rel(1e-2);
	double minf;
	nlopt::result res = opt.optimize(threshold,minf);
	m_threshold = threshold[0];*/
	
	Vertex_handle v = getFurtherFromSelectionCenter();
	std::vector<double> localDescr = v->getSemantic();
	myVector x(localDescr.size());
	for(unsigned l=0;l<localDescr.size();++l)
	{
		x[l] = localDescr[l];
	}
		
	myVector distToMean = x - m_mu;
	for(unsigned s=0;s<localDescr.size();++s)
	{
		std::cout<< x[s] << " ";
	}
	std::cout << std::endl;
	
	for(unsigned s=0;s<localDescr.size();++s)
	{
		std::cout<< m_mu[s] << " ";
	}
	std::cout << std::endl;
	for(unsigned s=0;s<localDescr.size();++s)
	{
		std::cout<< distToMean[s] << " ";
	}
	std::cout << std::endl;
	
	
	myMatrix distToMeanT = Linear_algebraCd<double>::transpose(distToMean);
	myVector mahalanobis = ((distToMeanT * m_inverseGaussianMatrix )* distToMean);
	
	m_threshold = sqrt(mahalanobis[0]);
	std::cout << "\nthreshold : " << m_threshold << std::endl;
}

void Correspondence_Component::learnSVMPatch(PolyhedronPtr p, unsigned SVM_mode)
{
	p->compute_normals();
	m_centreSelection = getSelectionCenter();
	Vertex_handle furthest = getFurtherFromSelectionCenter();
	m_patchRadius = sqrt(CGAL::squared_distance(m_centreSelection->point(),furthest->point()));
	m_nbCandidates = 3;
	
	std::map<Vertex_handle,bool> isSelected;
	for(unsigned s = 0; s<m_selection.size(); ++s)
	{
		isSelected[m_selection[s]] = true;
		m_selection[s]->color(1.0,1.0,0.0);
	}
	
	/* Increase number of samples for SVM */
	// -> increase resolution
	/*Facet_iterator bF = p->facets_begin();
	Facet_iterator eF = p->facets_end(); --eF;

	Facet_iterator f = bF;
	do {
		int countSel = 0;
		auto hC = f->facet_begin();
		do
		{
			if(m_tag[hC->vertex()]==1)
			{countSel++;}
			
		}
		while(++hC!=f->facet_begin());
		if(countSel == 3)
		{
			//std::cout << "Create new vertex" << std::endl;
			Halfedge_handle h = create_center_vertex_with_descriptor(p,f,m_nbLabel);
			m_selection.push_back(h->vertex());
			//isSelected[h->vertex()] = true;
			h->vertex()->color(0.0,1.0,0.0);
			m_tag[h->vertex()] = 1;
		}
	} while ( ++f != eF);*/
	
	p->compute_normals();
	//std::cout << "getSelectionCenter()" << std::endl;
	// Define the reference frame
	//	1- Compute the average normal on the patch
	Enriched_kernel::Vector_3 N;
	for(unsigned s=0;s<m_selection.size();++s)
	{
		Enriched_kernel::Vector_3 n = m_selection[s]->normal();
		N = N + n;
	}
	N = (1.0/m_selection.size())*N;
	
	//std::cout << "averagePatchNormal" << std::endl;
	//	2- Compute the two principal directions of curvature on the patch
	//Enriched_kernel::Vector_3 U,V;
	
	long double ppMatrix_sum[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
	
	for(unsigned s=0; s<m_selection.size();++s)
	{
		principal_curv(m_selection[s],ppMatrix_sum,m_selection.size());	
	}
	//Eigen Values / Vector of ppMatrix_sum
	//std::cout << "principal curv. per vert" << std::endl;
	
	std::cout << "m_selection size : " << m_selection.size() << std::endl;
	
	Eigen::Matrix<long double,3,3> mat;
	mat.setZero();
 
	mat << ppMatrix_sum[0][0], ppMatrix_sum[0][1], ppMatrix_sum[0][2],
	ppMatrix_sum[1][0], ppMatrix_sum[1][1], ppMatrix_sum[1][2],
	ppMatrix_sum[2][0], ppMatrix_sum[2][1], ppMatrix_sum[2][2];
	
	Eigen::EigenSolver<Eigen::Matrix<long double, 3, 3> > es(mat);
	
	auto eigVal = es.eigenvectors();
	
	double Ux,Uy,Uz,Vx,Vy,Vz;
	Ux = eigVal(0,1).real(); Uy = eigVal(1,1).real(); Uz = eigVal(2,1).real();
	Vx = eigVal(0,2).real(); Vy = eigVal(1,2).real(); Vz = eigVal(2,2).real();
	Enriched_kernel::Vector_3 U(Ux,Uy,Uz);
	Enriched_kernel::Vector_3 V(Vx,Vy,Vz);
	std::cout << "Learning patch : U, V" << std::endl;
	std::cout << Ux << " " << Uy << " " << Uz << std::endl;
	std::cout << Vx << " " << Vy << " " << Vz << std::endl;
	
	// learn the svm classifier
	svm_problem selectionProblem;
	selectionProblem.l = m_selection.size();
	std::vector<double> labels;
	std::vector<svm_node*> inputs;
	selectionProblem.x = Malloc(struct svm_node *, selectionProblem.l);
	
	unsigned nbFeatures = m_nbLabel+3;
	if(SVM_mode == 1)
	{
		nbFeatures = 3;
	}
	else if(SVM_mode == 2)
	{
		nbFeatures = m_nbLabel;
	}
	
	for(unsigned s=0; s<m_selection.size();++s)
	{
		Vertex_handle v = m_selection[s];
		Enriched_kernel::Vector_3 cv = m_centreSelection->point()-v->point();
		Enriched_kernel::Vector_3 projcv_n = (cv*N)*N;
		Enriched_kernel::Vector_3 projcv_uv = cv - projcv_n;
		
 		double d = sqrt(cv.squared_length());
		double h = sqrt(projcv_n.squared_length());
		double normU = sqrt(U.squared_length());
		double normProjcv_uv = sqrt(projcv_uv.squared_length());
		double a = acos((U*projcv_uv)/(normU*normProjcv_uv + 0.0001)); 
		
		/*std::cout << " a : " << a << std::endl;
		std::cout << " normProjcv_uv : " << normProjcv_uv << std::endl;
		std::cout << " arg : " << (U*projcv_uv)/(normU*normProjcv_uv + 0.0001) << std::endl;*/
		//std::cout << " num : " << U << projcv_uv << std::endl;
		//std::cout << " denom : " << normU << normProjcv_uv << std::endl;
		
		std::vector<double> & descr = m_selection[s]->getSemantic();
		if(m_tag[m_selection[s]] == 1)
		{
 			labels.push_back(1);
			m_selection[s]->color(1.0,1.0,0.0);
		}
		else
		{
			labels.push_back(0);
			m_selection[s]->color(0.0,1.0,1.0);
		}
		
		selectionProblem.x[s] = Malloc(struct svm_node, nbFeatures+1);
		if(SVM_mode != 1)
		{
			for(unsigned l=0;l<m_nbLabel;++l)
			{
				selectionProblem.x[s][l].index = l;
				selectionProblem.x[s][l].value = descr[l];
			}
		}
		if(SVM_mode == 0)
		{
			selectionProblem.x[s][m_nbLabel].index = m_nbLabel;
			selectionProblem.x[s][m_nbLabel].value = d;
			selectionProblem.x[s][m_nbLabel+1].index = m_nbLabel+1;
			selectionProblem.x[s][m_nbLabel+1].value = h;
			selectionProblem.x[s][m_nbLabel+2].index = m_nbLabel+2;
			selectionProblem.x[s][m_nbLabel+2].value = a;
		}
		else if(SVM_mode == 1)
		{
			selectionProblem.x[s][0].index = 0;
			selectionProblem.x[s][0].value = d;
			selectionProblem.x[s][1].index = 1;
			selectionProblem.x[s][1].value = h;
			selectionProblem.x[s][2].index = 2;
			selectionProblem.x[s][2].value = a;
		}
		
		selectionProblem.x[s][nbFeatures].index = -1;
	}
	selectionProblem.y = labels.data();
	
	/*unsigned uu;
	std::cin >> uu;*/
	
	/*for(unsigned s=0; s < m_selection.size();++s)
	{
		for(unsigned f = 0; f < nbFeatures; ++f)
		{
			std::cout << selectionProblem.x[s][f].value <<",";
		}
		if(labels[s]==1.0)
		{
			std::cout<<"selected"<<std::endl;
		}
		else
		{
			std::cout<<"not"<<std::endl;
		}
	}*/
	
	featureMeans.resize(nbFeatures,0.0);
	featureSTD.resize(nbFeatures,0.0);
	
	for(unsigned f = 0; f < nbFeatures; ++f)
	{
		double mean = 0.0;
		for(unsigned s = 0; s < m_selection.size(); ++s)
		{
			mean += selectionProblem.x[s][f].value;
		}
		featureMeans[f] = mean / m_selection.size();
		
	}
	
	for(unsigned f = 0; f < nbFeatures; ++f)
	{
		double std = 0.0;
		for(unsigned s = 0; s < m_selection.size(); ++s)
		{
			double dev = selectionProblem.x[s][f].value - featureMeans[f];
			selectionProblem.x[s][f].value = dev;
			std += dev*dev;
		}
		featureSTD[f] = sqrt(std/m_selection.size());
		/*if(f == (nbFeatures-1))
		{
			std::cout << std << std::endl;
		}*/
	}
	
	std::map<Vertex_handle,svm_node*> normalizedFeature;
	
	// normalize all features
	for(unsigned f = 0; f < nbFeatures; ++f)
	{
		for(unsigned s = 0; s < m_selection.size(); ++s)
		{
			selectionProblem.x[s][f].value /= (2*featureSTD[f]+0.001);
			
		}
		
		
	}
	
	for(unsigned s = 0; s < m_selection.size(); ++s)
	{
		normalizedFeature[m_selection[s]] = selectionProblem.x[s];
	}
	
	//std::cout << "all features have been normalized" << std::endl;
	
	// Find best gamma parameter as median(median(dist(p,neigh(p))))
	std::vector<double> medV;
	for(unsigned s = 0; s < m_selection.size(); ++s)
	{
		Vertex_handle v = m_selection[s];
		
		std::vector<double> descriptor;
		for(unsigned f = 0; f < nbFeatures; ++f)
		{
			descriptor.push_back(selectionProblem.x[s][f].value);
		}
		
		std::vector<double> distNeigh;
		
		auto hC = v->vertex_begin();
		do
		{
			Vertex_handle nV = hC->opposite()->vertex();
			if(normalizedFeature.count(nV) == 0){continue;}
			svm_node * nD = normalizedFeature[nV];
			std::vector<double> nDescriptor;
			for(unsigned f = 0; f < nbFeatures; ++f)
			{
				nDescriptor.push_back(nD[f].value);
			}
			double dist = L2Dist(descriptor,nDescriptor);
			distNeigh.push_back(dist);
			
		}while(++hC!=v->vertex_begin());
		
		double medN = 0;
		if(distNeigh.size() != 0)
		{
			std::sort(distNeigh.begin(),distNeigh.end());
			double medN = distNeigh[(int)(distNeigh.size()/2.0)];
			medV.push_back(medN);
		}
		
	}
	std::sort(medV.begin(),medV.end());
	double median = medV[(int)(medV.size()/2.0)];
	
	// Set parameters for SVM classifier
	struct svm_parameter param;
	param.svm_type = C_SVC;
	param.kernel_type = GAUSSIAN;
	param.gamma = median;
	param.C = pow(10,-5);
	
	param.eps = 0.5;
	param.cache_size = 100;
	param.shrinking = 1;
	param.probability = 0;
	param.nr_weight = 0;
	param.weight_label = NULL;
	param.weight = NULL;
	
	// grid search for svm parameters
	/*int bestC = 0;
	int bestG = 0;
	int maxScore = 0;
	for(int cExp = -30; cExp <=-5; cExp+=2)
	{
		for(int gExp = -30; gExp <= -15; gExp+=2)
		{
			param.gamma = pow(2,gExp);
			param.C = pow(2,cExp);
			m_svmModel = svm_train(&selectionProblem,&param);
			int score = 0;
			for(unsigned i = 0;i<labels.size();++i)
			{
				double predictedLabel = svm_predict(m_svmModel,selectionProblem.x[i]);
				if(predictedLabel == labels[i])
				{	
					score++;
				}
			}
			if(score >= maxScore)
			{
				bestC = cExp;
				bestG = gExp;
				maxScore = score;
			}
		}
	}
	std::cout << bestC<< " " << bestG << std::endl;
	param.gamma = pow(2,bestG);
	param.C = pow(2,bestC)
	;*/
	
	m_svmModel = svm_train(&selectionProblem,&param);
	
	int score = 0;
	for(unsigned i = 0; i < labels.size();++i)
	{
		double predictedLabel = svm_predict(m_svmModel,selectionProblem.x[i]);
		//std::cout << "predicted label : " << predictedLabel << std::endl;
		if(predictedLabel == labels[i])
		{
			score++;
			//m_selection[i]->color(0.0,0.0,1.0);
		}
		else
		{
			//m_selection[i]->color(1.0,0.0,1.0);
		}
		
		if(predictedLabel == 1)
			{
				m_selection[i]->color(1.0,1.0,0.0);
			}
			else
			{
				m_selection[i]->color(0.0,1.0,1.0);
			}
		//m_selection[i]->color(predictedLabel,0.0,1.0);
		
	}
	std::cout << "Learning set score : " << score/(double)labels.size() << std::endl;
	for(unsigned s=0;s<m_selection.size();++s)
	{
		free(selectionProblem.x[s]);
	}
	free(selectionProblem.x);
	
	std::cout << "SVM MODEL HAS BEEN TRAINED" << std::endl;
	std::cout << "Size of training set : " << m_selection.size() << std::endl;	
}

void Correspondence_Component::learnSVMClassifier(PolyhedronPtr p)
{	
	m_centreSelection = getSelectionCenter();
	Vertex_handle furthest = getFurtherFromSelectionCenter();
	
	m_patchRadius = 1.5*sqrt(CGAL::squared_distance(m_centreSelection->point(),furthest->point()));
	m_nbCandidates = 3;
	
	svm_problem selectionProblem;
	selectionProblem.l = m_selection.size();
	std::vector<double> labels;
	std::vector<svm_node*> inputs;
	selectionProblem.x = Malloc(struct svm_node *, selectionProblem.l);
	
	for(unsigned s=0;s<m_selection.size();++s)
	{
		std::vector<double> & descr = m_selection[s]->getSemantic();
		if(m_tag[m_selection[s]] == 1)
		{
 			labels.push_back(1);
			m_selection[s]->color(1.0,1.0,0.0);
		}
		else
		{
			labels.push_back(0);
			m_selection[s]->color(0.0,1.0,1.0);
		}
		
		selectionProblem.x[s] = Malloc(struct svm_node, m_nbLabel+1);
		for(unsigned l=0;l<m_nbLabel;++l)
		{
			selectionProblem.x[s][l].index = l;
			selectionProblem.x[s][l].value = descr[l];
		}
		selectionProblem.x[s][m_nbLabel].index = -1;
	}
	selectionProblem.y= labels.data();
	
	struct svm_parameter param;
	param.svm_type = C_SVC;
	param.kernel_type = GAUSSIAN;
	param.degree = 3;
	param.gamma = 1.0/m_nbLabel;	// 1/num_features
	param.coef0 = 0;
	param.nu = 2e-5;
	param.cache_size = 100;
	param.C = pow(2,5);
	param.eps = 0.5;
	param.p = 0.1;
	param.shrinking = 1;
	param.probability = 0;
	param.nr_weight = 0;
	param.weight_label = NULL;
	param.weight = NULL;
	
	param.C = pow(2,17);
	
	// grid search for svm parameters
	int bestC = 0;
	int bestG = 0;
	int maxScore = 0;
	/*for(int cExp = -5; cExp <=30; cExp+=2)
	{
		for(int gExp = -3; gExp <= 30; gExp+=2)
		{
			param.gamma = pow(2,gExp);
			param.C = pow(2,cExp);
			m_svmModel = svm_train(&selectionProblem,&param);
			int score = 0;
			for(unsigned i = 0;i<labels.size();++i)
			{
				double predictedLabel = svm_predict(m_svmModel,selectionProblem.x[i]);
				if(predictedLabel == labels[i])
				{	
					score++;
				}
			}
			if(score > maxScore)
			{
				bestC = cExp;
				bestG = gExp;
				maxScore = score;
			}
		}
	}
	std::cout << bestC<< " " << bestG << std::endl;
	param.gamma = pow(2,bestG);
	param.C = pow(2,bestC);*/
	param.gamma = pow(2,1);
	param.C = pow(2,17);
	m_svmModel = svm_train(&selectionProblem,&param);
	for(unsigned s=0;s<m_selection.size();++s)
	{
		free(selectionProblem.x[s]);
	}
	free(selectionProblem.x);
}

// Energy Functions 

double Correspondence_Component::computeEnergy(const std::vector<double> & ellipse)
{
	double energy = 0.0;
	
	std::vector<double> cDescr = m_centreDescriptor;
	//normalize(cDescr);
	
	for(unsigned i=0;i<m_selection.size();++i)
	{
		Vertex_handle v = m_selection[i];
		std::vector<double> localDescr = v->getSemantic();
//  		//normalize(localDescr);
		double eqEll = 0.0;
		for(unsigned l=0;l<localDescr.size();++l)
		{
			double diff = localDescr[l] - cDescr[l];
			eqEll += (diff*diff/(ellipse[l]*ellipse[l]));
		}
		bool inEllipse = (eqEll<=1);

		bool inSelection = (m_tag[v] <=2);
 
		if(inEllipse != inSelection)
		{
			double sqrtM = std::sqrt(eqEll) - 1;
			energy+=sqrtM*sqrtM;
		}
	}
	return energy;
}

double Correspondence_Component::computeEnergyRotation(const vector< double >& ellipse)
{
	double energy = 0.0;
	
	myMatrix EllipseMat(m_nbLabel,m_nbLabel);
	for(unsigned i=0;i<m_nbLabel;++i)
	{
		for(unsigned j=0;j<m_nbLabel;++j)
		{
			EllipseMat[i][j] = ellipse[i*m_nbLabel+j];
		}
	}
	
	std::vector<double> cDescr = m_centreSelection->getSemantic();
	
	myVector mu(m_centreDescriptor.size());
	for(unsigned l=0;l<m_centreDescriptor.size();++l)
	{
		mu[l] = m_centreDescriptor[l];
	}
	
	for(unsigned i=0;i<m_selection.size();++i)
	{
		Vertex_handle v = m_selection[i];
		std::vector<double> localDescr = v->getSemantic();
		double eqEll = 0.0;
		
		myVector x(localDescr.size());
		for(unsigned l=0;l<localDescr.size();++l)
		{
 			x[l] = localDescr[l];
		}
			
		myVector distToMean = x -mu;
		myMatrix distToMeanT = Linear_algebraCd<double>::transpose(distToMean);
		
		myVector res = distToMeanT*(EllipseMat*distToMean);
		eqEll = res[0];
		
		bool inEllipse = (eqEll<=1);

		bool inSelection = (m_tag[v] <=2);

		if(inEllipse != inSelection)
		{
			double sqrtM = std::sqrt(eqEll) - 1;
			energy+=sqrtM*sqrtM;
		}
	}
	return energy;
}

double Correspondence_Component::computeEnergyGaussian(const vector< double >& threshold)
{
	double energy = 0.0;
	
	myVector & mu = m_mu;

	double power = pow(2*3.14159,m_nbLabel);
	double norm = 1.0/(sqrt(power)*m_gaussianDeterminant);
	
	std:vector<double> dist;
	
	for(unsigned i=0;i<m_selection.size();++i)
	{
		Vertex_handle v = m_selection[i];
		std::vector<double> localDescr = v->getSemantic();
		myVector x(localDescr.size());
		for(unsigned l=0;l<localDescr.size();++l)
		{
 			x[l] = localDescr[l];
		}
			
		myVector distToMean = x - mu;
		myMatrix distToMeanT = Linear_algebraCd<double>::transpose(distToMean);
		myVector mahalanobis = (distToMeanT * m_inverseGaussianMatrix) * (distToMean);
		
		dist.push_back(sqrt(mahalanobis[0]));
	}
	double distMax = *std::max_element(dist.begin(),dist.end());
	//double distMax = 1.0;
	
	int f = 0;
	for(unsigned i=0;i<m_selection.size();++i)
	{
		Vertex_handle v = m_selection[i];
		float m = dist[f]/distMax;
		
		bool inGaussian = (m<=threshold[0]);

		bool inSelection = (m_tag[v] <=2);

		if(inGaussian && !inSelection)
		{
			double incr= m;
			energy+=incr*incr;
		}
		f++;
	}
	return energy;
}


// Getters and setters

std::vector< double > Correspondence_Component::getCentreDescr() const
{
	return m_centreDescriptor;
}

std::vector< double > Correspondence_Component::getEllipse() const
{
	return m_ellipse;
}

myMatrix Correspondence_Component::getMatrix() const
{
	return m_gaussianMatrix;
}

myVector Correspondence_Component::getVector() const
{
	return m_mu;
}


myMatrix Correspondence_Component::getInverseMatrix() const
{
	return m_inverseGaussianMatrix;
}

double Correspondence_Component::getDeterminant() const
{
	return m_gaussianDeterminant;
}

double Correspondence_Component::getThreshold() const
{
	return m_threshold;
}

svm_model* Correspondence_Component::getSVM() const
{
	return m_svmModel;
}


void Correspondence_Component::setCentreDescriptor(const std::vector<double> & centreDescr)
{
	m_centreDescriptor = centreDescr;
}


void Correspondence_Component::setEllipse(const vector<double> & ellipse)
{
	m_ellipse = ellipse;
}

void Correspondence_Component::setMatrix(const myMatrix& m)
{
	m_gaussianMatrix = m;
}

void Correspondence_Component::setVector(const myVector& v)
{
	m_mu = v;
}

void Correspondence_Component::setSVM(svm_model* svm)
{
	m_svmModel = svm;
}


void Correspondence_Component::setInverseMatrix(const myMatrix& m)
{
	m_inverseGaussianMatrix = m;
}

void Correspondence_Component::setDeterminant(double det)
{
	m_gaussianDeterminant = det;
}

void Correspondence_Component::setThreshold(double thresh)
{
	m_threshold = thresh;
}

void Correspondence_Component::scaleMesh(Polyhedron::Iso_cuboid bbox, PolyhedronPtr p)
{
	double dx = bbox.xmax() - bbox.xmin();
	double dy = bbox.ymax() - bbox.ymin();
	double dz = bbox.zmax() - bbox.zmin();
	
	Polyhedron::Iso_cuboid nbox = p->bbox();
	double nx = nbox.xmax() - nbox.xmin();
	double ny = nbox.ymax() - nbox.ymin();
	double nz = nbox.zmax() - nbox.zmin();

	double scale = 1.0;
	
	if(dx > std::max(dy,dz))
	{
		scale = dx/nx;
	}
	if(dy > std::max(dx,dz))
	{
		scale = dy/ny;
	}
	else
	{
		scale = dz/nz;
	}
	for(Vertex_iterator pVertex = p->vertices_begin();
		pVertex!=p->vertices_end();++pVertex)
	{
		Point3d vp = pVertex->point();
		Point3d nvp(vp.x()*scale,vp.y()*scale,vp.z()*scale);
		pVertex->point() = nvp;
	}
	p->compute_bounding_box();
	p->compute_normals();
}

double Correspondence_Component::getRadius() const
{
	return m_patchRadius;
}
int Correspondence_Component::getNbCandidates() const
{
	return m_nbCandidates;
}

void Correspondence_Component::setRadius(double radius)
{
	m_patchRadius = radius;
}

void Correspondence_Component::setNbCandidates(int nbCandidates)
{
	m_nbCandidates = nbCandidates;
}

//// NON-MEMBER FUNCTIONS
/*double L2Dist(std::vector<double>& descr1, std::vector<double>& descr2)
{
	double dist = 0.0;
	for(unsigned el = 0; el<descr1.size(); ++el)
	{
		double diff = descr1[el] - descr2[el];
		if(descr1[el]==std::numeric_limits<double>::infinity() || descr2[el]==std::numeric_limits<double>::infinity())
		{continue;}
		dist += diff*diff;
	}
	return sqrt(dist);
}*/

double objectiveFun(const std::vector<double> & ellipse, std::vector<double> & grad, void *data)
{
	Correspondence_Component * correspondenceComp = (Correspondence_Component *) data;
	double energy = correspondenceComp->computeEnergy(ellipse);
	return energy;
}

double objectiveFunRotation(const std::vector<double> & ellipse, std::vector<double> & grad, void *data)
{
	Correspondence_Component * correspondenceComp = (Correspondence_Component *) data;
	double energy = correspondenceComp->computeEnergyRotation(ellipse);
	return energy;
}

double objectiveFunGaussian(const std::vector<double> & threshold, std::vector<double> & grad, void *data)
{
	Correspondence_Component * correspondenceComp = (Correspondence_Component *) data;
	double energy = correspondenceComp->computeEnergyGaussian(threshold);
	return energy;
}

void vector_times_transpose_multi(long double pVector[3],
                                               long double ppMatrix[3][3],
                                               double coeff)
{
  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      ppMatrix[i][j] = coeff * pVector[i] * pVector[j];
}

//**********************************************
// add two matrices
//**********************************************
void addM(long double pMatrix[3][3],
                       long double pMatrixSum[3][3])
{
  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      pMatrixSum[i][j] += pMatrix[i][j];
}

//**********************************************
// fix sine
//**********************************************
double fix_sin(double sine)
{
  if (sine >= 1)
    return 3.14159/2;
  else
    if (sine <= -1)
      return -3.14159/2;
    else
      return std::asin(sine);
}

 double areaFace(Facet_handle &f)
	{
		Halfedge_around_facet_circulator pHalfedge = f->facet_begin();
		Point3d P = pHalfedge->vertex()->point();
		Point3d Q = pHalfedge->next()->vertex()->point();
		Point3d R = pHalfedge->next()->next()->vertex()->point();

		Vector PQ=Q-P;
                //Vector PR=R-P; // MT
		Vector QR=R-Q;


		Vector normal = CGAL::cross_product(PQ,QR);
		double area=0.5*sqrt(normal*normal);

		return area;

	}

void principal_curv(Vertex_handle pVertex, long double ppMatrix_sum[3][3], int size)
{

double area=0;

  // iterate over all edges
  Halfedge_around_vertex_circulator pHalfedge = pVertex->vertex_begin();
  Halfedge_around_vertex_circulator pHalfedgeStart = pHalfedge;
  CGAL_For_all(pHalfedge,pHalfedgeStart)
  {

    // build edge vector and comput its norm
	Point3d p1 = pHalfedge->vertex()->point();
	Point3d p2 = pHalfedge->opposite()->vertex()->point();
	Vector edge = (p1-p2);
	double len_edge = std::sqrt(edge*edge);
	if (len_edge == 0) // avoid divide by zero
	continue;

	// compute (signed) angle between two incident faces, if exists
	Facet_handle pFacet1 = pHalfedge->facet();
	Facet_handle pFacet2 = pHalfedge->opposite()->facet();
	CGAL_assertion(pFacet1 != pFacet2);
	if (pFacet1 == NULL || pFacet2 == NULL)
	continue; // border edge

	area+=areaFace(pFacet1);

	Vector normal1 = pFacet1->normal();
	Vector normal2 = pFacet2->normal();

	double sine = (CGAL::cross_product(normal1,normal2)*edge)/len_edge;
	double beta = fix_sin(sine);

	// compute edge * edge^t * coeff, and add it to current matrix
	long double pVector_edge[3] = {edge.x(),edge.y(),edge.z()};
	long double ppMatrix[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
	
	
	double factor = pow(10,size);
	
	
	vector_times_transpose_multi(pVector_edge,ppMatrix,beta/(len_edge*factor));

	addM(ppMatrix,ppMatrix_sum);
  }
  
	for (int i=0;i<3;i++)
		for (int j=0;j<3;j++)
			ppMatrix_sum[i][j]/=area;
}

void computePatchBasis(std::vector<Vertex_handle> & sel, Enriched_kernel::Vector_3 &  U, Enriched_kernel::Vector_3 & V, Enriched_kernel::Vector_3 & N)
{
	for(unsigned s=0;s<sel.size();++s)
	{
		Enriched_kernel::Vector_3 n = sel[s]->normal();
		N = N + n;
	}
	N = (1.0/sel.size())*N;
	
	//std::cout << "averagePatchNormal" << std::endl;
	//	2- Compute the two principal directions of curvature on the patch
	//Enriched_kernel::Vector_3 U,V;
	
	long double ppMatrix_sum[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
	
	for(unsigned s=0; s<sel.size();++s)
	{
		principal_curv(sel[s],ppMatrix_sum,sel.size());
	}
	//Eigen Values / Vector of ppMatrix_sum
	//std::cout << "principal curv. per vert" << std::endl;
	
	Eigen::Matrix<long double,3,3> mat;
	mat.setZero();
	mat << ppMatrix_sum[0][0], ppMatrix_sum[0][1], ppMatrix_sum[0][2],
	ppMatrix_sum[1][0], ppMatrix_sum[1][1], ppMatrix_sum[1][2],
	ppMatrix_sum[2][0], ppMatrix_sum[2][1], ppMatrix_sum[2][2];
	
	Eigen::EigenSolver<Eigen::Matrix<long double, 3, 3> > es(mat);
	
	auto eigVal = es.eigenvectors();
	
	double Ux,Uy,Uz,Vx,Vy,Vz;
	Ux = eigVal(0,1).real(); Uy = eigVal(1,1).real(); Uz = eigVal(2,1).real();
	Vx = eigVal(0,2).real(); Vy = eigVal(1,2).real(); Vz = eigVal(2,2).real();
	
	std::cout << "test patch : U, V" << std::endl;
	std::cout << Ux << " " << Uy << " " << Uz << std::endl;
	std::cout << Vx << " " << Vy << " " << Vz << std::endl;
	
	U = Enriched_kernel::Vector_3(Ux,Uy,Uz);
	V = Enriched_kernel::Vector_3(Vx,Vy,Vz);
}

Halfedge_handle create_center_vertex_with_descriptor( PolyhedronPtr p, Facet_iterator f, int nbLabel) {
    Vector vec( 0.0, 0.0, 0.0);
    std::vector<double> descr(nbLabel,0.0);
    std::size_t order = 0;
    auto h = f->facet_begin();
    
    do {
	std::vector<double> nDescr = h->vertex()->getSemantic();
        vec = vec + ( h->vertex()->point() - CGAL::ORIGIN);
        ++ order;
	for(unsigned l = 0; l <nbLabel; ++l)
	{
		descr[l] = nDescr[l];
	}
    } while ( ++h != f->facet_begin());
    CGAL_assertion( order >= 3); // guaranteed by definition of polyhedron
    Point3d center =  CGAL::ORIGIN + (vec / static_cast<double>(order));
    for(unsigned l = 0; l < nbLabel; ++l)
    {
	    descr[l]/=order;
    }
    Halfedge_handle new_center = p->create_center_vertex( f->halfedge());
    new_center->vertex()->setSemantic(descr);
    new_center->vertex()->point() = center;
    new_center->vertex()->color(1.0,0.0,0.0);
    return new_center;
}


#endif
