#include "SegmentController.h"
#include <set>


#include <CGAL/Polyhedron_items_with_id_3.h>
#include <../../src/components/Tools/Boolean_Operations/src/Boolean_Operations_Component.h>

SegmentController::SegmentController(PolyhedronPtr p) : m_polyhedron(p)
{}

SegmentController::~SegmentController()
{}

void SegmentController::addMainPart(PolyhedronPtr p)
{
	m_mainPart.push_back(p);
}
	
void SegmentController::addSegement(PolyhedronPtr p)
{
	m_parts.push_back(p);
}

void SegmentController::cutSegments()
{
	std::vector<Vertex_handle> selection;
	std::map<Vertex_iterator,int>  ccMap;
	std::map<Vertex_iterator,bool> isSelected;
	
	int nbcc = 1;
	int mainCC = -1;
	
	for(Facet_iterator pFacet = m_polyhedron->facets_begin();
	    pFacet!=m_polyhedron->facets_end();++pFacet)
	{
		bool facetBelongs = false;
		int nbBelongs = 0;
		Halfedge_around_facet_circulator h = pFacet->facet_begin();
		do
		{
			Vertex_handle v = h->vertex();
			facetBelongs = facetBelongs || isColoredBlack(v);
			if(isColoredBlack(v))
			{
				nbBelongs++;
			}
			
		}
		while(++h!=pFacet->facet_begin());
		//if(facetBelongs)
		if(nbBelongs>=2)
		{
			do
			{
				isSelected[h->vertex()] = true;
				selection.push_back(h->vertex());
			}
			while(++h!=pFacet->facet_begin());	
		}
	}
	
	for(Vertex_iterator pVertex = m_polyhedron->vertices_begin();
	    pVertex!=m_polyhedron->vertices_end();++pVertex)
	{
		if(ccMap[pVertex] == 0) // this vertex has not been visited
		{
			if(isSelected[pVertex]) // part of segment
			{
				nbcc = nbcc + 1; // create a new connected component
				ccMap[pVertex] = nbcc;
				//pVertex->color(0,1,0);
				Halfedge_around_vertex_circulator h = pVertex->vertex_begin();
				visitVertexSelection(h,isSelected,ccMap,nbcc);
			}
			else // part of main mesh
			{
				mainCC = mainCC - 1;
				ccMap[pVertex] = mainCC;
				//pVertex->color(0,0,1);
				Halfedge_around_vertex_circulator h = pVertex->vertex_begin();
				visitVertexSelection(h,isSelected,ccMap,mainCC);
			}
		}
	}
	
	std::set<int> labelSet;
	
	for(auto it = ccMap.begin(); it!= ccMap.end(); ++it)	
	{
		if(it->second >0 )
		{
			it->first->color(1,1,0);
		}
		else if (it->second<0)
		{
			it->first->color(0,1,1);
		}
		labelSet.insert(it->second);
	}
	// All vertices have been assigned with a cc number
	// Now build each cc accordingly
	
	std::vector<double> coords;
	for(Vertex_iterator pVertex = m_polyhedron->vertices_begin();
	    pVertex != m_polyhedron->vertices_end();++pVertex)
	{
		coords.push_back(pVertex->point().x());
		coords.push_back(pVertex->point().y());
		coords.push_back(pVertex->point().z());
	}// collect all vertices for builder
	
	m_polyhedron->set_index_vertices();
	for(auto it = labelSet.begin();it!=labelSet.end();++it)
	{
		if(*it == 0){continue;}
		std::vector<int> tris;
		for(Facet_iterator pFacet = m_polyhedron->facets_begin();
		pFacet!=m_polyhedron->facets_end();++pFacet)
		{
			bool facetBelongs = false;
			int nbBelongs = 0;
			Halfedge_around_facet_circulator h = pFacet->facet_begin();
			do
			{
				facetBelongs = facetBelongs || (ccMap[h->vertex()] == *it);
				if((ccMap[h->vertex()] == *it))
				{
					nbBelongs++;
				}
			}
			while(++h!=pFacet->facet_begin());
			h = pFacet->facet_begin();
			//if(facetBelongs)
			if(nbBelongs>=2)
			{
				do
				{
					tris.push_back(h->vertex()->tag());
					
				}
				while(++h!=pFacet->facet_begin());
			}
		}
	
		PolyhedronPtr p(new Polyhedron);
		polyhedron_builder<HalfedgeDS> builder(coords,tris);
		p->delegate(builder);
		copyDescriptor(m_polyhedron,p);
		
		if(p->size_of_facets() != 0)
		{
			if(*it>0)
			{	
				m_parts.push_back(p);
			}
			else 
			{	
				m_mainPart.push_back(p);
			}
		}
		p->compute_bounding_box();
		p->compute_normals();
		p->calc_nb_components();
	}
}

void SegmentController::glueSegments()
{
	// Init the indices of the halfedges and the vertices
	for(unsigned p=0; p<m_parts.size();++p)
	{
		PolyhedronPtr part = m_parts[p]; 
	
		simplePolyhedron * mesh = convertToSimplePolyhedron(part);
		mesh->keep_largest_connected_components(1); //remove isolated vertices
		mesh->normalize_border();
		// compute correspondence
		std::map<simplePolyhedron::Vertex_handle,simplePolyhedron::Vertex_handle> corresp;
		for(unsigned m=0;m<m_mainPart.size();++m)
		{
			simplePolyhedron * mainPart = convertToSimplePolyhedron(m_mainPart[m]);
			mainPart->keep_largest_connected_components(1); //remove isolated vertices
			auto h = mesh->halfedges_begin();
			for(;h!=mesh->halfedges_end();++h)
			{
				if(h->is_border())
				{
					double distMax = std::numeric_limits<double>::max();
					for(auto pVertex = mainPart->vertices_begin(); pVertex!= mainPart->vertices_end();++pVertex)
					{
						double dist = CGAL::squared_distance(pVertex->point(),h->vertex()->point());
						if(dist < distMax)
						{
							corresp[h->vertex()] = pVertex;
							distMax = dist;
						}
					}
				}
			}
		}
		
		CGAL::set_halfedgeds_items_id(*mesh);
		
		// Create a deformation object
		surface_mesh_deformation deform_mesh(*mesh);
		// Definition of the Region Of Interest
		/*for(simplePolyhedron::Vertex_iterator it =  mesh->vertices_begin(); it!= mesh->vertices_end(); ++it)
		{
			deform_mesh.insert_roi_vertex(it);
		}*/
		// Select control vertices
		
		auto h = mesh->halfedges_begin();
		for(;h!=mesh->halfedges_end();++h)
		{
			if(h->is_border())
			{
				deform_mesh.insert_roi_vertex(h->vertex());
				auto n = h->vertex()->vertex_begin();
				do
				{
					deform_mesh.insert_roi_vertex(n->opposite()->vertex());
				}
				while(++n!=h->vertex()->vertex_begin());
				deform_mesh.insert_control_vertex(h->vertex());
			}
		}
		
		bool is_matrix_factorization_OK = deform_mesh.preprocess();
		if(!is_matrix_factorization_OK)
		{
			std::cout << "preprocess fail" << std::endl;
		}
		
		h = mesh->halfedges_begin();
		for(;h!=mesh->halfedges_end();++h)
		{
			if(h->is_border())
			{
				deform_mesh.set_target_position(h->vertex(),corresp[h->vertex()]->point());
			} 
		}
		//deform_mesh.set_iterations(10);
		//deform_mesh.set_tolerance(0.0);
		deform_mesh.deform();
		m_parts[p] = convertToEnrichedPolyhedron(mesh);
	}
}

/* Fuse the meshes
	Look for border halfedge, find the closest vertex and fill the hole*/
void SegmentController::fillHoles(Viewer* v, PolyhedronPtr p)
{
	Boolean_Operations_Component booleanOp(v,p);
	for(unsigned m=0;m<m_parts.size();++m)
	{
		booleanOp.Boolean_Union(p,m_parts[m],p);
	}
	for(unsigned m=1;m<m_mainPart.size();++m)
	{
		booleanOp.Boolean_Union(p,m_mainPart[m],p);
	}
}

void visitVertexSelection(Halfedge_around_vertex_circulator h, std::map<Vertex_iterator,bool> & isSelected,
			  std::map<Vertex_iterator,int> & cc, int nbcc)
    {
	Vertex_iterator hVertex = h->opposite()->vertex();
	
	if(cc[hVertex] != 0)// neighbour has been visited
	{
		return;
	}
	else if (nbcc >0 && isSelected[hVertex])
	{
		cc[hVertex] = nbcc;
		//hVertex->color(0,1,0);
		Halfedge_around_vertex_circulator hNext = hVertex->vertex_begin();
		do
		{
			visitVertexSelection(hNext,isSelected,cc,nbcc);
		}
		while(++hNext!=hVertex->vertex_begin());
	}
	else if( nbcc <0 && !isSelected[hVertex] )
	{
		cc[hVertex] = nbcc;
		//hVertex->color(0,0,1);
		Halfedge_around_vertex_circulator hNext = hVertex->vertex_begin();
		do
		{
			visitVertexSelection(hNext,isSelected,cc,nbcc);
		}
		while(++hNext!=hVertex->vertex_begin());	
	}
}

bool isColoredBlack(Vertex_handle v)
{
	return ( v->color(0) == 0.0 && v->color(1) == 0.0 && v->color(2) == 0.0);
}

void colorMesh(PolyhedronPtr p, float r, float g, float b)
{
	for(Vertex_iterator pVertex = p->vertices_begin();
	    pVertex!= p->vertices_end();
		++pVertex)
	{
		pVertex->color(r,g,b);
	}
}

void copyDescriptor(PolyhedronPtr s, PolyhedronPtr d)
{
	Vertex_iterator dVertex = d->vertices_begin();
	for(Vertex_iterator pVertex = s->vertices_begin();
	    pVertex!= s->vertices_end();
		++pVertex)
	    {
		dVertex->setSemantic(pVertex->getSemantic());   
		++dVertex;  
	    }
}

// Make a copy in a Polyhedron_3 data structure, because EnrichedPolyhedron and boost::graph_traits crash
simplePolyhedron * convertToSimplePolyhedron(PolyhedronPtr p)
{
	std::vector<double> coords;
	for(Vertex_iterator pVertex = p->vertices_begin();
	    pVertex != p->vertices_end();++pVertex)
	{
		coords.push_back(pVertex->point().x());
		coords.push_back(pVertex->point().y());
		coords.push_back(pVertex->point().z());
	}// collect all vertices for builder
	p->set_index_vertices();
	std::vector<int> tris;
	for(Facet_iterator pFacet = p->facets_begin();
	pFacet!=p->facets_end();++pFacet)
	{
		Halfedge_around_facet_circulator h = pFacet->facet_begin();
		do
		{
			tris.push_back(h->vertex()->tag());
		}
		while(++h!=pFacet->facet_begin());
	}
	
	simplePolyhedron* sp(new simplePolyhedron);
	polyhedron_builder<simplePolyhedron::HalfedgeDS> builder(coords,tris);
	sp->delegate(builder);
	return sp;
}

PolyhedronPtr convertToEnrichedPolyhedron(simplePolyhedron * p)
{
	std::vector<double> coords;
	for(auto pVertex = p->vertices_begin();
	    pVertex != p->vertices_end();++pVertex)
	{
		coords.push_back(pVertex->point().x());
		coords.push_back(pVertex->point().y());
		coords.push_back(pVertex->point().z());
	}// collect all vertices for builder
	
	std::vector<int> tris;
	for(auto pFacet = p->facets_begin();
	pFacet!=p->facets_end();++pFacet)
	{
		auto h = pFacet->facet_begin();
		do
		{
			tris.push_back(h->vertex()->id());
		}
		while(++h!=pFacet->facet_begin());
	}
	
	PolyhedronPtr sp(new Polyhedron);
	polyhedron_builder<HalfedgeDS> builder(coords,tris);
	sp->delegate(builder);
	delete p;
	return sp;
}
