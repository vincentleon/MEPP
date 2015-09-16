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
	

#include "../../../../mepp/mepp_component.h"

#include "Correspondence_Component.h"
#include "Correspondence_Polyhedron.h"
#include "geodesic/geodesic_algorithm_exact.h"

Correspondence_Component::Correspondence_Component(Viewer* v, PolyhedronPtr p) : mepp_component(v,p)
{
	componentName = "Correspondence_Component";
	init = 1;
	p->set_index_vertices();
}

void Correspondence_Component::initParameters(int nbLabel, int meshId)
{
	m_Shape.m_meshID = meshId;
	m_nbLabel = nbLabel;
	m_Shape.initFaceLabelsAndSegments();
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

void Correspondence_Component::saveDescriptor(PolyhedronPtr p)
{
	std::ofstream file;
	std::stringstream ss;
	ss<<m_Shape.m_meshID<<".semantic";
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

void Correspondence_Component::readDescriptor(PolyhedronPtr p)
{
	Point3d bb = Point3d(p->xmin(),p->ymin(),p->zmin());
	std::ifstream file;
	std::stringstream ss;
	ss<<"/home/leon/datasetHuman/"<<m_Shape.m_meshID<<".semantic";
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
	}
	initMaxVector(p);
}

void Correspondence_Component::showDescriptor(PolyhedronPtr p, int dim)
{
	for(Vertex_iterator pVertex = p->vertices_begin();
	    pVertex!=p->vertices_end();++pVertex)
	{
		std::vector<double> & localDescr = pVertex->getSemantic();
		pVertex->color(localDescr[dim],0.5,1.0-localDescr[dim]);
	}
}


void Correspondence_Component::initMaxVector(PolyhedronPtr p)
{
	for(unsigned l=0;l<m_nbLabel;++l)
	{
		double maxV = 0.0;
		for(Vertex_iterator pVertex = p->vertices_begin();
		pVertex!=p->vertices_end();
		++pVertex)
		{
			double val = pVertex->getSemantic()[l];
			if(val>maxV){maxV = val;}
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
			descr[i]/m_maxVector[i];
		}
	}
}

void Correspondence_Component::computeGeodesicDistancesAllVertices(PolyhedronPtr p)
{
	for(unsigned label=0; label>m_nbLabel;++label)
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
			pVertex->pushSemantic(distance);
		}
	}
}


void Correspondence_Component::computeDescriptorAllVertices(PolyhedronPtr p)
{
	this->initGeodesicMesh(p);
	this->computeGeodesicDistancesAllVertices(p);
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
	std::cout << "closest has been found" << distMin << std::endl;
	return localDescr;
}

void Correspondence_Component::compareDescriptorToEllipse(PolyhedronPtr p)
{
	std::vector<double> centreDescr = m_centreDescriptor;
	normalize(centreDescr);

	for(Vertex_iterator pVertex = p->vertices_begin();
		    pVertex!=p->vertices_end();++pVertex)
	{
		std::vector<double> localDescr = pVertex->getSemantic();
		normalize(localDescr);
		double eqEll = 0.0;
		for(unsigned l=0;l<localDescr.size();++l)
		{
			double diff = localDescr[l]- centreDescr[l];
			eqEll += (diff*diff/(m_ellipse[l]*m_ellipse[l]));
		}
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
	myVector mu(m_centreDescriptor.size());
	for(unsigned l=0;l<m_centreDescriptor.size();++l)
	{
		mu[l] = m_centreDescriptor[l];
	}
	
	myMatrix I,D;
	myVector C;
	std::vector<int> vv;
	double dd;
		
	std::cout << std::endl;
	std::cout << "Matrix : " << std::endl;
	for(unsigned i = 0;i<m_nbLabel;++i)
	{
		for(unsigned j = 0;j<m_nbLabel;++j)
		{
			std::cout << m_gaussianMatrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
	double det = Linear_algebraCd<double>::determinant(m_gaussianMatrix,I,D,vv,C);
 	std::cout << "determinant : " << det << std::endl;
	
	myMatrix inverse = Linear_algebraCd<double>::inverse(m_gaussianMatrix,dd);
	
	double power = pow(2*3.14159,m_nbLabel);
	double norm = 1.0/(sqrt(power)*det);
	
	std:vector<double> dist;
	
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
		myVector mahalanobis = (distToMeanT * inverse) * (distToMean);
		
		std::cout << sqrt(mahalanobis[0]) << std::endl;
		
		/*if(sqrt(mahalanobis[0])<10)
		{
			pVertex->color(0.5,0.5,0.5);
		}
		else
		{
			pVertex->color(0,0,0);
// 		}*/
		dist.push_back(log(sqrt(mahalanobis[0])));
	}
	//double distMax = *std::max_element(dist.begin(),dist.end());
	double distMax = 1.0;
	int v = 0;
	for(Vertex_iterator pVertex = p->vertices_begin();
	    pVertex!=p->vertices_end();++pVertex)
	{
		float m = dist[v]/distMax;
		pVertex->color(m,0,0);
		v++;
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
		v->color(0,1,0);
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
	//this->readSelectionBasedOnColor(p);
	m_centreSelection = getSelectionCenter();
	//m_centreSelection->color(1,0,0);
	Vertex_handle extremumSelection = getFurtherFromSelectionCenter();
	//extremumSelection->color(1,0,0);
	m_isolineValue = L2Dist(m_centreSelection->getSemantic(),extremumSelection->getSemantic());
	
	m_ellipse = std::vector<double>(m_nbLabel,m_isolineValue);
}

void Correspondence_Component::computeEllipseParameters(PolyhedronPtr p)
{
	m_ellipse = std::vector<double>(m_nbLabel,m_isolineValue);
	
	std::vector<double> ell = m_ellipse;
	
	int dim = m_ellipse.size();
	
	nlopt::opt opt(nlopt::LN_COBYLA,dim);
	
	void * data = &(*this);
	
	opt.set_min_objective(objectiveFun,data);
	
	std::vector<double> lBounds(dim,0.0);
	opt.set_lower_bounds(lBounds);
	opt.set_xtol_rel(1e-4);
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
}


double Correspondence_Component::computeEnergy(const std::vector<double> & ellipse)
{
	double energy = 0.0;
	
	std::vector<double> cDescr = m_centreDescriptor;
	normalize(cDescr);
	
	for(unsigned i=0;i<m_selection.size();++i)
	{
		Vertex_handle v = m_selection[i];
		std::vector<double> localDescr = v->getSemantic();
		normalize(localDescr);
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



//// NON-MEMBER FUNCTIONS
double L2Dist(std::vector<double>& descr1, std::vector<double>& descr2)
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
}

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


#endif
