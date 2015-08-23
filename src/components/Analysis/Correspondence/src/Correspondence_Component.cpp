///////////////////////////////////////////////////////////////////////////
// Author: Vincent Leon
// Year: 2015
// Universite Lille 1, CRIStAL
//
///////////////////////////////////////////////////////////////////////////
#include <mepp_config.h>
#ifdef BUILD_component_Correspondence

#include <iostream>
#include <sstream>

#include "../../../../mepp/mepp_component.h"

#include "Correspondence_Component.h"
#include "Correspondence_Polyhedron.h"
#include "geodesic/geodesic_algorithm_exact.h"

Correspondence_Component::Correspondence_Component(Viewer* v, PolyhedronPtr p) : mepp_component(v,p)
{
	
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
		points.push_back(pVertex->point().x());// add the vertex to the list of vertices
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
	ss<<m_Shape.m_meshID<<".semantic";
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


			Point3d p = pVertex->point();

			//localDescr.push_back(std::sqrt(CGAL::squared_distance(p,bb)));
			localDescr.push_back( std::abs(p.x() - p->xmin()) );
			localDescr.push_back( std::abs(p.y() - p->ymin()) );
			localDescr.push_back( std::abs(p.z() - p->zmin()) );
			pVertex->setSemantic(localDescr);
			++pVertex;
		}
	}
	else
	{
		std::cout << "Impossible d'ouvrir le fichier "<<ss.str()<<" !\n";
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
		if(max[i]!=std::numeric_limits<double>::infinity())
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
			    pVertex!=p->vertices_end();+pVertex)
			{
				pVertex->pushSemantic(std::numeric_limits<double>::infinity());
			}
			
		}
		
		for(Vertex_iterator pVertex = p->vertices_begin();
		    pVertex!=p->vertices_end();++pVertex)
		{
			int pIndex = pVertex->getIndex();
			geodesic::SurfacePoint p(&m_gmesh->vertices()[pIndex];
			double distance;
			unsigned best_source = m_geoAlg->best_source(p,distance);
			pVertex->pushSemantic(distance);
		
		}
	}
}


void Correspondence_Component::computeDescriptorAllVertices(PolyhedronPtr p)
{
	this->initGeodesicMesh();
	this->computeGeodesicDistancesAllVertices();
	this->initMaxVector();
	for(Vertex_iterator pVertex = p->vertices_begin();
		    pVertex!=p->vertices_end();++pVertex)
	{
		std::vector<double> descr = pVertex->getSemantic();
		this->normalize(descr);
		pVertex->setSemantic(descr);
	}
}

vector<double> Correspondence_Component::getClosetVertexDescriptor(PolyhedronPtr p, std::vector<double> & pt)
{
	//TODO: improve picking
}

void Correspondence_Component::compareToDescrEllipse(PolyhedronPtr p, vector< double >& ellipse)
{
	std::vector<double> centreDescr = m_centreSelection->getSemantic();
	
	for(Vertex_iterator pVertex = p->vertices_begin();
		    pVertex!=p->vertices_end();++pVertex)
	{
		std::vector<double> localDescr = pVertex->getSemantic();
		double eqEll = 0.0;
		for(unsigned l=0;l<localDescr.size();++l)
		{
			double diff = localDescr[l]- centreDescr[l];
			eqEll += (diff*diff/(ellipse[l]*ellipse[l]));
		}
		if(eqEll<=1)
		{
			pVertex->color(0,255,0);
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
			
			//double dist = L2Dist(m_Semantics[selection[i]],m_Semantics[selection[j]]);
			double dist = L2Dist(m_selection[i]->getSemantic(),m_selection[j]->getSemantic());
			score+=dist*dist;
			
		}
		if(score<scoreMin)
		{
			scoreMin=score;
			m_centreSelection = m_selection[i];
		}
	}
	m_centreSelection->color(0,0,255);
	return m_centreSelection;
	}
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
	furtherFromCenter->color(0,255,0);
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
			v->color(255,255,0);
				}
		else if (m_tag[v] == 3 )
		{
			v->color(0,255,255);
		}
	}
}

void Correspondence_Component::readSelectionBasedOnColor(PolyyhedronPtr p)
{
	for(Vertex_iterator pVertex = p->vertices_begin();
	    pVertex!=p->vertices_end();
	    ++pVertex)
	    {
		int r = pVertex->color(0);    
		int g = pVertex->color(1);
		int b = pVertex->color(2);
		
		if( r==255 && g==0 && b==0 )
		{
			m_selection.push_back(pVertex);
		}
	    }
	this->tagSelectionVertices(p);
}



void Correspondence_Component::initializeEllipsoid(PolyhedronPtr p)
{
	this->readSelectionBasedOnColor(p);
	m_centreSelection = getSelectionCenter();
	Vertex_handle extremumSelection = getFurtherFromSelectionCenter();
	m_isolineValue = L2Dist(m_centreSelection->getSemantic(),extremumSelection->getSemantic());
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


#endif