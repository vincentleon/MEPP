#include "SegmentController.h"
#include "softICPController.h"
#include <set>

#include <viewer.hxx>
#include <Polyhedron/polyhedron.h>
#include <QGLViewer/vec.h>
#include <QGLViewer/quaternion.h>
#include <boost/graph/graph_concepts.hpp>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

#include <CGAL/Polygon_mesh_processing/self_intersections.h>


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

double SegmentController::energyBorder(Halfedge_handle prev, Halfedge_handle curr, Halfedge_handle next, double alpha, double beta)
{
	
	Point3d p,pN,pP;
	p =  curr->vertex()->point();
	pP = prev->vertex()->point();
	pN = next->vertex()->point(); 

	Vector vi(p.x(),p.y(),p.z());
	Vector viP(pP.x(),pP.y(),pP.z()); 
	Vector viN(pN.x(),pN.y(),pN.z());
	
	Vector distanceTerm(prev->vertex()->point(),curr->vertex()->point());
	Vector curvatureTerm = viN -2*vi + viP;
	
	double eInt = alpha * sqrt(distanceTerm.squared_length()) + beta * sqrt(curvatureTerm.squared_length());
	
	std::vector<double > & prevDescr = prev->vertex()->getSemantic();
	std::vector<double > & currDescr = curr->vertex()->getSemantic();
	
	double eExt = L2Dist(prevDescr,currDescr);
	//return eInt + eExt;
	//return eInt;
	return eExt;
}


void SegmentController::moveBorder(vector< Halfedge_handle >& snake)
{
	for(unsigned i=0;i<snake.size();++i)
	{
		// Circular configuration
		Halfedge_handle prev,next;
		if(i==0) { prev = snake.back();}
		else{prev = snake[i-1];}
		if(i==snake.size()-1) {next = snake[0];}
		else{ next = snake[i+1];}
		
		// Previous energy
		double eMin = energyBorder(prev,snake[i],next);
		Halfedge_handle minCandidate = snake[i];
		// for each neighbor vertex of snake[i]
		Halfedge_around_vertex_circulator h = snake[i]->vertex()->vertex_begin();
		do{
			Halfedge_handle candidate = h->opposite();
			double energy = energyBorder(prev,candidate,next);
			if(energy < eMin)
			{
				eMin = energy;
				minCandidate = candidate;
			}
		}
		while(h!=snake[i]->vertex()->vertex_begin());
		
		// Move snake[i] to location minCandidate whose energy is smallest
		//if(minCandidate == 
		snake[i] = minCandidate;
	}
	
	for(auto pVertex = m_polyhedron->vertices_begin(); pVertex!=m_polyhedron->vertices_end(); ++pVertex)
	{
		pVertex->color(1,1,1);
	}
	for(unsigned b = 0; b < m_borders.size(); ++b)
	{
		for(unsigned v = 0; v < m_borders[b].size(); ++v)
		{
			m_borders[b][v]->vertex()->color(1,0,0);
		}
	}
	
}

snaxel * SegmentController::createSnake(std::vector<Halfedge_handle> & border)
{
	snaxel * snakeBegin = new snaxel();
	snakeBegin->he = border[0];
	snaxel * snake = snakeBegin;
	
	for(unsigned i=1;i<border.size();++i)
	{
		
		snake->next = new snaxel();
		snake->next->he = border[i];
		snake->next->previous = snake;
		snake = snake->next;
	}
	
	// close the loop
	snake->next = snakeBegin;
	snakeBegin->previous = snake;
	return snakeBegin;
}


void SegmentController::moveSnake(snaxel * snakeBegin)
{
	snaxel * snake = snakeBegin;
	do{
		//previous energy
		double eMin = energyBorder(snake->previous->he,snake->he,snake->next->he);
		Halfedge_handle minCandidate = snake->he;
		
		Halfedge_around_vertex_circulator h = snake->he->vertex()->vertex_begin();
		do
		{
			Halfedge_handle candidate = h->opposite();
			double energy = energyBorder(snake->previous->he,candidate,snake->next->he); // Pas sur de cette partie, quid de la connectivité du snake ?
			if(energy < eMin)
			{
				eMin = energy;
				minCandidate = candidate;
			}	
		}
		while(h!=snake->he->vertex()->vertex_begin());
		// replace snake.current with candidate and insert new snaxels if necessary
		// insérer nvx snaxels ici ?
		snake->he = minCandidate;
		snake = snake->next;
	}
	while(snake!=snakeBegin);
	
	// insérer nvx snaxels ici ?
		
}

void SegmentController::replaceSnaxel(snaxel * curr, Halfedge_handle candidate)
{

}




void SegmentController::getBorder()
{
	std::map<Facet_handle,bool> isFaceSelected;
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
		
		if(nbBelongs>=2)
		{
			do
			{
				isFaceSelected[pFacet]=true;
			}
			while(++h!=pFacet->facet_begin());	
		}
	}
	
	std::map<Halfedge_handle,bool> isHEVisited;
	//std::vector< std::vector< Halfedge_handle> > borders;
	for(auto pHalfedge = m_polyhedron->halfedges_begin(); pHalfedge!=m_polyhedron->halfedges_end();
	    ++pHalfedge)
	{
		if(isHEVisited[pHalfedge] || isHEVisited[pHalfedge->opposite()]){continue;}
		if(isFaceSelected[pHalfedge->facet()] && !isFaceSelected[pHalfedge->opposite()->facet()])
		{
			Halfedge_handle borderStart = pHalfedge;
			Halfedge_handle borderCirculator = pHalfedge;
			
			m_borders.resize(m_borders.size()+1);
			m_borders.back().push_back(borderStart);
			
			do{
				Halfedge_around_vertex_circulator h = borderCirculator->vertex()->vertex_begin();
				do{
					isHEVisited[h] = true;
					if(h->opposite() == borderCirculator){continue;}
					
					if(isFaceSelected[h->opposite()->facet()] && !isFaceSelected[h->facet()])
					{
						borderCirculator = h->opposite();
						isHEVisited[borderCirculator] = true;
						m_borders.back().push_back(borderCirculator);
						break;
					}
				}
				while(++h!=pHalfedge->vertex()->vertex_begin());
				
			}
			while(borderCirculator!= borderStart);
		}
	}
	
	
	for(auto pVertex = m_polyhedron->vertices_begin(); pVertex!=m_polyhedron->vertices_end(); ++pVertex)
	{
		pVertex->color(1,1,1);
		
	}
	
	std::cout << "# of borders : " << m_borders.size() << std::endl;
	for(unsigned b = 0; b < m_borders.size(); ++b)
	{
		std::cout << "size of border #" << b << " : " << m_borders[b].size() << std::endl;
		for(unsigned v = 0; v < m_borders[b].size(); ++v)
		{
			m_borders[b][v]->vertex()->color(1,0,0);
		}
		m_snakes.push_back(createSnake(m_borders[b]));
	}
	
}

void SegmentController::optimizeBorders()
{
	for(unsigned b = 0; b < m_snakes.size(); ++b)
	{
		//moveBorder(m_borders[b]);
		snaxel * snakeBegin = m_snakes[b];
		snaxel * snake = snakeBegin;
		do{
			snake->he->vertex()->color(1,1,1);
			snake = snake->next;
		}
		while(snake!=snakeBegin);
		
		moveSnake(m_snakes[b]);
		snake = snakeBegin;
		do{
			snake->he->vertex()->color(0,1,0);
			snake = snake->next;
		}
		while(snake!=snakeBegin);
	}
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
		/*if(it->second >0 )
		{
			//it->first->color(1,1,0);
		}
		else if (it->second<0)
		{
			//it->first->color(0,1,1);
		}*/
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
		//p->keep_largest_connected_components(1);
		p->compute_bounding_box();
		p->compute_normals();
		p->calc_nb_components();
	}
}

PolyhedronPtr SegmentController::fillHoles(PolyhedronPtr p)
{
	p->keep_largest_connected_components(1);
	simplePolyhedron  * poly = convertToSimplePolyhedron(p);
	/*
	unsigned nb_holes = 0;
	for(auto h = poly->halfedges_begin(); h!=poly->halfedges_end(); ++h)
	{
		if(h->is_border())
		{
		std::vector<simplePolyhedron::Facet_handle> patch_facets;
		std::vector<simplePolyhedron::Vertex_handle> patch_vertices;
		bool success = CGAL::cpp11::get<0>(
			CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(
			*poly,
			h,
			std::back_inserter(patch_facets),
			std::back_inserter(patch_vertices),
			CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point,*poly)).geom_traits(simpleKernel())));
			
			++nb_holes;
		}
	}*/
	PolyhedronPtr r = convertToEnrichedPolyhedron(poly);
	return r;
}

void SegmentController::fuseMeshes(Viewer* v, PolyhedronPtr target, PolyhedronPtr model)
{
	this->unionSegments(v);
	
	std::vector<double> coords;
	std::vector<int> tris;
	
	collectVertsAndFaces(target,coords,tris);
	collectVertsAndFaces(model,coords,tris);
	
	PolyhedronPtr p(new Polyhedron);
	polyhedron_builder<HalfedgeDS> builder(coords,tris);
	p->delegate(builder);
	
	// stitch_borders(p)
	simplePolyhedron * sp = convertToSimplePolyhedron(p);
	CGAL::Polygon_mesh_processing::stitch_borders(*sp);
	
	// triangulate_refine_and_fair_hole(p)
	for(simplePolyhedron::Halfedge_iterator h = sp->halfedges_begin();
	    h!=sp->halfedges_end();
		++h)
	{
		if(h->is_border())
		{
			std::vector<simplePolyhedron::Facet_handle> patch_facets;
		std::vector<simplePolyhedron::Vertex_handle> patch_vertices;
		bool success = CGAL::cpp11::get<0>(
			CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(
			*sp,
			h,
			std::back_inserter(patch_facets),
			std::back_inserter(patch_vertices),
			CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point,*sp)).geom_traits(simpleKernel())));
		}
	}
	
	
	PolyhedronPtr newPoly = convertToEnrichedPolyhedron(sp);
	newPoly->compute_normals();
	v->addFrame();
	v->getScenePtr()->add_polyhedron(newPoly);
}

void SegmentController::fuseMeshes2(Viewer* v, PolyhedronPtr target, PolyhedronPtr model, int idHETarget, int idHEModel)
{
	this->unionSegments(v);
	
	simplePolyhedron * t = convertToSimplePolyhedron(target);
	simplePolyhedron * m = convertToSimplePolyhedron(model);
	
	simplePolyhedron::Halfedge_iterator ht = t->halfedges_begin();
	simplePolyhedron::Halfedge_iterator hm = m->halfedges_begin();
	
	for(unsigned i = 0; i < idHETarget; i++)
	{
		ht++;
	}
	for(unsigned i = 0; i < idHEModel; i++)
	{
		hm++;
	}
	
	if(ht->is_border())
	{
		std::vector<simplePolyhedron::Facet_handle> patch_facets;
	std::vector<simplePolyhedron::Vertex_handle> patch_vertices;
	bool success = CGAL::cpp11::get<0>(
		CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(
		*t,
		ht,
		std::back_inserter(patch_facets),
		std::back_inserter(patch_vertices),
		CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point,*t)).geom_traits(simpleKernel())));
	}

	
	
	// triangulate_refine_and_fair_hole(p)
	
	if(hm->is_border())
	{
		std::vector<simplePolyhedron::Facet_handle> patch_facets;
	std::vector<simplePolyhedron::Vertex_handle> patch_vertices;
	bool success = CGAL::cpp11::get<0>(
		CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(
		*m,
		hm,
		std::back_inserter(patch_facets),
		std::back_inserter(patch_vertices),
		CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point,*m)).geom_traits(simpleKernel())));
	}
	
	std::vector<double> coords;
	std::vector<int> tris;
	
	model = convertToEnrichedPolyhedron(m);
	target = convertToEnrichedPolyhedron(t);
	
	collectVertsAndFaces(target,coords,tris);
	collectVertsAndFaces(model,coords,tris);
	
	/*PolyhedronPtr p(new Polyhedron);
	polyhedron_builder<HalfedgeDS> builder(coords,tris);
	p->delegate(builder);*/
	
	// stitch_borders(p)
	//simplePolyhedron * sp = convertToSimplePolyhedron(p);
	//CGAL::Polygon_mesh_processing::stitch_borders(*sp);
	
	
	
	
	//PolyhedronPtr newPoly = convertToEnrichedPolyhedron(sp);
	model->compute_normals();
	target->compute_normals();
	
	v->addFrame();
	v->getScenePtr()->add_polyhedron(model);
	v->addFrame();
	v->getScenePtr()->add_polyhedron(target);
}


void SegmentController::unionSegments(Viewer* v)
{
	float x, y, z;
	double a, b, c, w;
	ScenePtr S = v->getScenePtr();
	
	//Put everything in world coordinates
	for(unsigned i=0; i<v->getScenePtr()->get_nb_polyhedrons();i++)
	{
		Vertex_iterator pVertex = NULL;
		PolyhedronPtr P = S->get_polyhedron(i);
		//P->keep_largest_connected_components(1);
		v->frame(i)->getPosition(x,y,z);
		v->frame(i)->getOrientation(a,b,c,w);
		qglviewer::Vec T(x,y,z);
		qglviewer::Quaternion Q(a,b,c,w);
		
		for(pVertex = P->vertices_begin(); pVertex != P->vertices_end();++pVertex)
		{
			qglviewer::Vec V = Q * qglviewer::Vec(pVertex->point().x(),
							      pVertex->point().y(),
							      pVertex->point().z()) + T;
			pVertex->point() = CGAL::ORIGIN + Vector(V[0],V[1],V[2]);
		}
		v->frame(i)->setPosition(0,0,0);
		v->frame(i)->setOrientation(0,0,0,1);
	}
	v->show();
	v->recreateListsAndUpdateGL();
}


void SegmentController::joinSegments(Viewer * v)
{
	std::vector<double> coords;
	std::vector<int> tris;
	
	std::vector<PolyhedronPtr> newMeshes;
	std::vector<qglviewer::Quaternion> Orientation;
	std::vector<qglviewer::Vec> Position;
	
	this->unionSegments(v);
	
	for(unsigned partID = 0; partID < v->getScenePtr()->get_nb_polyhedrons(); ++partID)
	{
		PolyhedronPtr p = v->getScenePtr()->get_polyhedron(partID);
		p->set_index_vertices();
 		PolyhedronPtr n = fillHoles(p);
		n->compute_normals();
		newMeshes.push_back(n);
		
		/*v->frame(partID)->getOrientation(a,b,c,w); qglviewer::Quaternion Q(a,b,c,w);
		v->frame(partID)->getPosition(x,y,z); qglviewer::Vec P(x,y,z);*/
		
		/*Orientation.push_back(Q);
		Position.push_back(P);*/
		
	}
	for(unsigned partID = 0; partID < v->getScenePtr()->get_nb_polyhedrons(); ++partID)
	{
		v->getScenePtr()->delete_polyhedron(partID);
	}
	
	
	for(unsigned nm = 0; nm<newMeshes.size();++nm)
	{
		v->getScenePtr()->add_polyhedron(newMeshes[nm]);
		//v->addFrame();
		//v->frame(v->get_nb_frames()-1)->rotate(Orientation[nm]);
		//v->frame(v->get_nb_frames()-1)->translate(Position[nm]);
		v->getScenePtr()->setVisible(nm,true);
	}
	v->getScenePtr()->todoIfModeSpace(v,0.0);
	v->recreateListsAndUpdateGL();
	
	/*for(unsigned partID =0; partID <  v->getScenePtr()->get_nb_polyhedrons(); ++partID)
	{	
		int offset = coords.size() / 3;
		PolyhedronPtr p = v->getScenePtr()->get_polyhedron(partID);
		p->set_index_vertices();
		fillHoles(p);
		p->set_index_vertices();
		for(auto pVertex = p->vertices_begin();
		pVertex != p->vertices_end();++pVertex)
		{
			Point3d vertexInWorld = toWorld(v,partID,pVertex->point());
			//Point3d vertexInWorld = pVertex->point();
			coords.push_back(vertexInWorld.x());
			coords.push_back(vertexInWorld.y());
			coords.push_back(vertexInWorld.z());
		}// collect all vertices for builder
		
		for(auto pFacet = p->facets_begin();
		pFacet!=p->facets_end();++pFacet)
		{
			auto h = pFacet->facet_begin();
			do
			{
				tris.push_back(h->vertex()->tag()+offset);
			}
			while(++h!=pFacet->facet_begin());
		}// collect all faces for builder
	}
	
	PolyhedronPtr sp(new Polyhedron);
	polyhedron_builder<HalfedgeDS> builder(coords,tris);
	sp->delegate(builder);
	v->addFrame();
	v->getScenePtr()->add_polyhedron(sp);
	sp->compute_normals();
	sp->calc_nb_components();
	sp->calc_nb_boundaries();*/
	
	//fillHoles(sp);
	/*
	 * Then extract all borders
	 * Then, contour following correspondence (~todo)
	 * Then convert to simplePolyhedron (ok)
	 * Then use CGAL deformation with appropriate ROI (see example)
	 */
	
	/////////////////////////////////////////////////////////
	/// Extract all borders from sp
	/////////////////////////////////////////////////////////
	/*std::map<Halfedge_handle,bool> isVisited;
	std::vector<std::vector<Halfedge_handle> > borders; // all the borders belonging to part "p"
	
	sp->normalize_border();
	for(auto he = sp->border_halfedges_begin(); he!=sp->halfedges_end();++he)
	{
		if(! (isVisited[he] && isVisited[he->opposite()]) )  //halfedge has not been visited, it's a new border loop
		{
			visitBorder(he,isVisited,borders.back()); //Visit the neighbours, if they are border hE, add them to the border
			borders.resize(borders.size()+1); //Add a new, empty borders
		}
	}
	
	for(auto b=borders.begin();b!=borders.end();++b)
	{
		if(b->size() == 0)
		{
			b = borders.erase(b);
			if(b == borders.end()) { break;}
		}
	}
	/////////////////////////////////////////////////////////
	/// Border correspondence
	/////////////////////////////////////////////////////////
	Halfedge_handle f1;
	Halfedge_handle f2;
	Halfedge_handle p1;
	Halfedge_handle p2;
	Halfedge_handle c1;
	Halfedge_handle c2;
	
	c1 = borders[0][0]; // first sewing halfedge
	double thresh = 100.0;
	int borderIndex;
	
	double distMin = std::numeric_limits<double>::max();
	for(unsigned b=1;b<borders.size();++b) // find the corresponding halfedge on other borders
	{
		for(unsigned h =0; h < borders[b].size(); ++h)
		{
			double dist = CGAL::squared_distance(c1->vertex()->point(),borders[b][h]->vertex()->point());
			if(dist<distMin)// && dist < thresh)
			{
				distMin = dist;
				c2 = borders[b][h];
				borderIndex = b;
			}
		}
	}
	std::vector<Halfedge_handle> cBord1; cBord1.push_back(c1);
	std::vector<Halfedge_handle> cBord2; cBord2.push_back(c2);
	
	Vector V1 = c1->vertex()->point() - c1->opposite()->vertex()->point(); 
	Vector V2 = c2->vertex()->point() - c2->opposite()->vertex()->point();
	if( V1*V2 < 0) {c2 = c2->opposite();}
	
	
	
	p1 = c1; f1 = c1;
	p2 = c2; f2 = c2;
	
	bool endB1 = false;
	bool endB2 = false;
	
	
	
	do
	{
		// Find which halfedge to increment (c1 or c2)
		Halfedge_around_vertex_circulator p1N = p1->vertex()->vertex_begin();
		Halfedge_around_vertex_circulator p2N = p2->vertex()->vertex_begin();
		double distMinB1 = std::numeric_limits<double>::max();
		double distMinB2 = std::numeric_limits<double>::max();
		
		// Which hE on Border1 is closest to p2
		Halfedge_handle candidateC1;
		do{
			if(!p1N->opposite()->is_border_edge()){continue;}
			if(p1N->opposite() == p1->opposite()){continue;} // make sure we really increment, and don't go backward
			double dist = CGAL::squared_distance(p2->vertex()->point(),p1N->opposite()->vertex()->point());
			//std::cout << dist << std::endl;
			if(dist < distMinB1)// && dist < thresh)
			{
				distMinB1 = dist;
				candidateC1 = p1N->opposite();
			}
		}
		while(++p1N != p1->vertex()->vertex_begin());
		
		// Which hE on Border2 is closest to p1
		Halfedge_handle candidateC2;
		do{
			if(!p2N->opposite()->is_border_edge()){continue;}
			if(p2N->opposite() == p2->opposite()){continue;} // make sure we really increment, and don't go backward
			double dist = CGAL::squared_distance(p1->vertex()->point(),p2N->opposite()->vertex()->point());
			//std::cout << dist << std::endl;
			if(dist < distMinB2)// && dist < thresh)
			{
				distMinB2 = dist;
				candidateC2 = p2N->opposite();
			}
		}
		while(++p2N != p2->vertex()->vertex_begin());
		
		// Choose the less costly incrementation
		//std::cout << distMinB1 << " " << distMinB2 << std::endl;
		if(distMinB1 < distMinB2)
		{
			//std::cout << "Increment C1 " << std::endl;
			c1 = candidateC1;
		}
		else
		{
			//std::cout << "Increment C2 " << std::endl;
			c2 = candidateC2;
		}
		cBord1.push_back(c1); p1=c1;
		cBord2.push_back(c2); p2=c2;
		if(c1 == f1){endB1 = true;}
		if(endB1 && (c2==f2))
		{break;}
		if(c2 == f2){endB2 = true;}
		if(endB2 && (c1==f1))
		{
			break;
		}
	}
	while( true ); 
	
	//For now, add faces
	for(unsigned b = 0; b < cBord1.size();++b)
	{
		//auto hnew = sp->split_edge(cBord1[b]); hnew->vertex()->point() = cBord2[b]->vertex()->point();
		//hnew = sp->split_edge(cBord2[b]); hnew->vertex()->point() = cBord1[b]->vertex()->point();
		Point3d mem = cBord1[b]->vertex()->point();
		cBord1[b]->vertex()->point() = cBord2[b]->vertex()->point();
		cBord2[b]->vertex()->point() = mem;
		
	}*/
	//sp->compute_normals();
	
	/////////////////////////////////////////////////////////
	/// Non rigid deformation
	/////////////////////////////////////////////////////////
	
	
}

void SegmentController::glueSegments(Viewer * v)
{
	/////////////////////////////////////////////////////////
	/// Extract all borders for all parts in the scene
	/////////////////////////////////////////////////////////
	/*int nbSeg = m_parts.size() + m_mainPart.size();
	
	// Get all borders from all parts m[id_of_part][id_of_border][id_of_halfedge]
	std::vector<std::vector<std::vector<Halfedge_handle> > > m_borderHE(nbSeg);
	
	for(unsigned p=0; p< nbSeg;++p)
	{
		std::map<Halfedge_handle,bool> isVisited;
		std::vector<std::vector<Halfedge_handle> > & borders = m_borderHE[p]; // all the borders belonging to part "p"
		PolyhedronPtr part;
		if(p < m_parts.size())
		{
			part = m_parts[p];
		}
		else
		{
			part  = m_mainPart[p - m_parts.size()];
		}
		part->normalize_border();
		for(auto he = part->border_halfedges_begin(); he!=part->halfedges_end();++he)
		{
			if(! (isVisited[he] && isVisited[he->opposite()]) )  //halfedge has not been visited, it's a new border loop
			{
				visitBorder(he,isVisited,borders.back()); //Visit the neighbours, if they are border hE, add them to the border
				borders.resize(borders.size()+1); //Add a new, empty borders
			}
		}
		
		for(auto b=borders.begin();b!=borders.end();++b)
		{
			if(b->size() == 0)
			{
				b = borders.erase(b);
				if(b == borders.end()) { break;}
			}
		}		
	}
	
	
	
	std::map<Halfedge_handle,Halfedge_handle> corresHE;
	std::map<Halfedge_handle,bool> hasCorres;
	std::map<Halfedge_handle,Point3d> wcCorres;
	
	for(unsigned p = 0; p< m_parts.size(); ++p)
	{
		PolyhedronPtr & part = m_parts[p];
		//for one border of p
		for(unsigned b=0;b<m_borderHE[p].size();++b)
		{
			Halfedge_handle fHE;
			double distMinOverBorder = std::numeric_limits<double>::max();
			for(unsigned h = 0; h< m_borderHE[p][b].size();++h)
			{	
				Halfedge_handle bh = m_borderHE[p][b][h];
				Point3d wcp = toWorld(v,p,bh->vertex()->point());
				
				double distMin = std::numeric_limits<double>::max();
				for(unsigned q = 0; q < m_mainPart.size(); ++q)
				{
					// For now take all mesh as ROI, maybe make it smaller later
					for(Halfedge_iterator hE = m_mainPart[q]->halfedges_begin(); hE!= m_mainPart[q]->halfedges_end();++hE)
					{
						Point3d wcmp = toWorld(v,m_parts.size()+q,hE->vertex()->point());
						
						if(hasCorres[hE]){continue;}
						double dist = CGAL::squared_distance(wcmp,wcp);
						if(dist < distMin)
						{
							distMin = dist;
							corresHE[bh] = hE;
							wcCorres[bh] = wcmp;
							hasCorres[hE] = true;
						}
					}
				}
				if(wcCorres.count(bh)==0){continue;}
				if(distMin<distMinOverBorder)
				{
					fHE = bh;
					distMinOverBorder = distMin;
				}
			}
			// Sew the first edge
			int q = partID;
			part->split_edge(fHE); // create new geometry using euler operators
			fHE->next()->vertex()->point() = toPartFrame(v,p,wcCorres[fHE]); // displace new vertex to coorresponding point
			Halfedge_handle current1;
			Halfedge_handle previous1;
			Halfedge_handle previous2 = corresHE[previous1];
			current1 = fHE;
			do
			{
				Halfedge_around_vertex_circulator hC = current1->vertex()->vertex_begin();
				do
				{
					Halfedge_handle hh = hC->opposite();
					if(hh->is_border() && hh!= current1->opposite())
					{	
						previous1 = current1;
						current1 = hh;
					}
				}
				while( hC != current1->vertex()->vertex_begin());
				
				
				///////// find the correspondence of current1 among the neighbors of 
				Point3d neighcorres;
				Halfedge_handle current2;
				double distMin = std::numeric_limits<double>::max();
				Halfedge_around_vertex_circulator hC2 = previous2->vertex()->vertex_begin();
				do
				{
					Halfedge_handle hh2 = hC->opposite();
					Point3d twCandidate = toWorld(v,q,hh2);
					Point3d twNeigh = toWorld(v,p,current1);
					double dist = CGAL::squared_distance(twNeigh,twCandidate);
					if(dist < distMin)
					{
						distMin = dist;
						neighcorres = twCandidate;
						current2 = hh2;
						corresHE[current1] = current2;
						wcCorres[current1] = toWorld(v,q,current2);
					}
					
				}
				while( hC2 != previous2->vertex()->vertex_begin());
				
				part->split_edge(current1);
				current1->next()->vertex()->point() = toPartFrame(v,p,wcCorres[current1]);
			}
			while( current1 != fHE );
		}
	}*/
}



void SegmentController::sewSegments(Viewer * v, PolyhedronPtr target, PolyhedronPtr model)
{
	
	
	//std::cout << " sew Segments " << std::endl;
	//target->keep_largest_connected_components(1);
	//model->keep_largest_connected_components(1);
	//std::cout << " after largestconnectedCC " << std::endl;
	
	/*unionSegments(v);
	std::cout << " after unionSegments " << std::endl;
	
	v->recreateListsAndUpdateGL();
	
	target->normalize_border();
	//model->normalize_border();
	
	std::map<Halfedge_handle,bool> isVisited;
	std::vector<std::vector<Halfedge_handle> > borders;
	
	/////////////////////////////////////////////////////////
	/// Border correspondence
	/////////////////////////////////////////////////////////
	Halfedge_handle f1;
	Halfedge_handle f2;
	Halfedge_handle p1;
	Halfedge_handle p2;
	Halfedge_handle c1;
	Halfedge_handle c2;
	
	c1 = target->border_halfedges_begin();
	
	//c1->vertex()->color(1,0,0);
	
	double thresh = 100.0;
	int borderIndex;
	
	double distMin = std::numeric_limits<double>::max();
	PolyhedronPtr modelTest = model;
	
	ScenePtr scene_ptr = v->getScenePtr();
	for(unsigned i =0;i<scene_ptr->get_nb_polyhedrons();++i)
	{
		modelTest = scene_ptr->get_polyhedron(i);
		if(modelTest == target){continue;}
		modelTest->normalize_border();
		for(auto h = modelTest->halfedges_begin(); h!=modelTest->halfedges_end();++h)
		{
				if(!h->is_border()){continue;}
				double dist = CGAL::squared_distance(c1->vertex()->point(),h->vertex()->point());
				if(dist<distMin && dist < thresh)
				{
					distMin = dist;
					c2 = h;
					model = modelTest;
				}
		}
	}*/
	/*
	std::vector<Halfedge_handle> cBord1; cBord1.push_back(c1);
	std::vector<Halfedge_handle> cBord2; cBord2.push_back(c2);
	//
	//
	Halfedge_iterator ex1 = c1;
	Halfedge_iterator ex2 = c2;
	for(unsigned i=0;i<10;++i)
	{
		++ex1;
		++ex2;
	}
	Vector V1 = c1->vertex()->point() - ex1->vertex()->point(); 
	Vector V2 = c2->vertex()->point() - ex2->vertex()->point();
	if( V1*V2 < 0) {c2 = c2->opposite();}
	c2->vertex()->color(1,0,0);
	
	
	p1 = c1; f1 = c1;
	p2 = c2; f2 = c2;
	
	bool endB1 = false;
	bool endB2 = false;
	
	do
	{
		// Find which halfedge to increment (c1 or c2)
		Halfedge_around_vertex_circulator p1N = p1->vertex()->vertex_begin();
		Halfedge_around_vertex_circulator p2N = p2->vertex()->vertex_begin();
		double distMinB1 = std::numeric_limits<double>::max();
		double distMinB2 = std::numeric_limits<double>::max();
		
		// Which hE on Border1 is closest to p2
		Halfedge_handle candidateC1;
		do{
			if(!p1N->opposite()->is_border_edge()){continue;}
			if(p1N->opposite() == p1->opposite()){continue;} // make sure we really increment, and don't go backward
			double dist = CGAL::squared_distance(p2->vertex()->point(),p1N->opposite()->vertex()->point());
			//std::cout << dist << std::endl;
			if(dist < distMinB1)// && dist < thresh)
			{
				distMinB1 = dist;
				candidateC1 = p1N->opposite();
			}
		}
		while(++p1N != p1->vertex()->vertex_begin());
		
		// Which hE on Border2 is closest to p1
		Halfedge_handle candidateC2;
		do{
			if(!p2N->opposite()->is_border_edge()){continue;}
			if(p2N->opposite() == p2->opposite()){continue;} // make sure we really increment, and don't go backward
			double dist = CGAL::squared_distance(p1->vertex()->point(),p2N->opposite()->vertex()->point());
			//std::cout << dist << std::endl;
			if(dist < distMinB2)// && dist < thresh)
			{
				distMinB2 = dist;
				candidateC2 = p2N->opposite();
			}
		}
		while(++p2N != p2->vertex()->vertex_begin());
		
		// Choose the less costly incrementation
		//std::cout << distMinB1 << " " << distMinB2 << std::endl;
		if(distMinB1 < distMinB2)
		{
			//std::cout << "Increment C1 " << std::endl;
			c1 = candidateC1;
		}
		else
		{
			//std::cout << "Increment C2 " << std::endl;
			c2 = candidateC2;
		}
		cBord1.push_back(c1); p1=c1;
		cBord2.push_back(c2); p2=c2;
		if(c1 == f1){endB1 = true;}
		if(endB1 && (c2==f2))
		{break;}
		if(c2 == f2){endB2 = true;}
		if(endB2 && (c1==f1))
		{
			break;
		}
	}
	while( true ); 
	std::cout << " after border correspondence " << std::endl;*/
	//for(unsigned b = 0; b < cBord1.size();++b)
	/*{
		Point3d mem = cBord1[b]->vertex()->point();
		auto hnew = target->split_edge(cBord1[b]); hnew->vertex()->point() = cBord2[b]->vertex()->point();
		auto hhnew = model->split_edge(cBord2[b]); hhnew->vertex()->point() = mem;
		deformation.insert_control_vertex(cBord1[b]->vertex());
		Point3d mem = cBord1[b]->vertex()->point();
		cBord1[b]->vertex()->point() = cBord2[b]->vertex()->point();
		cBord2[b]->vertex()->point() = mem;
		
	}*/
	/*for(unsigned b = 0; b < cBord1.size();++b)
	{
		if(b<cBord1.size()-1)
		{
			if(cBord1[b] == cBord1[b+1])
			{
				auto hnew = target->split_edge(cBord1[b]); hnew->vertex()->point() = cBord2[b]->vertex()->point();
			}
			else
			{
				cBord1[b]->vertex()->point() = cBord2[b]->vertex()->point();
			}
		}
		else
		{
			if(cBord1[b]== cBord1[0])
			{
				auto hnew = target->split_edge(cBord1[b]); hnew->vertex()->point() = cBord2[b]->vertex()->point();
			}
		}
	}
	std::cout << " after displacement " << std::endl;
	*/
	model->compute_normals();
	target->compute_normals();
	
	softICPController icp(model,target);
	//icp.buildTreeStructure(2);
	icp.snapRegions(0.2,1);
	
	
	//test_softICP_SVD();
	v->recreateListsAndUpdateGL();
	
}


void SegmentController::fitSegments(Viewer * v, PolyhedronPtr target, PolyhedronPtr model, int & idBordTarget, int & idBordModel)
{
	unionSegments(v);
	
	simplePolyhedron * t = convertToSimplePolyhedron(target);
	simplePolyhedron * m = convertToSimplePolyhedron(model);
	
	t->normalize_border();
	m->normalize_border();
	
	std::map<simplePolyhedron::Halfedge_handle,bool> isVisited;
	std::vector<std::vector<simplePolyhedron::Halfedge_handle> > borders;
	
	/*for(auto he = m->border_halfedges_begin(); he!=m->halfedges_end();++he)
	{
		if(! (isVisited[he] && isVisited[he->opposite()]) )  //halfedge has not been visited, it's a new border loop
		{
			visitBorder(he,isVisited,borders.back()); //Visit the neighbours, if they are border hE, add them to the border
			borders.resize(borders.size()+1); //Add a new, empty borders
		}
	}
	
	for(auto b=borders.begin();b!=borders.end();++b)
	{
		if(b->size() == 0)
		{
			b = borders.erase(b);
			if(b == borders.end()) { break;}
		}
	}*/
	/////////////////////////////////////////////////////////
	/// Border correspondence
	/////////////////////////////////////////////////////////
	simplePolyhedron::Halfedge_handle f1;
	simplePolyhedron::Halfedge_handle f2;
	simplePolyhedron::Halfedge_handle p1;
	simplePolyhedron::Halfedge_handle p2;
	simplePolyhedron::Halfedge_handle c1;
	simplePolyhedron::Halfedge_handle c2;
	
	c1 = t->border_halfedges_begin();
	
	//c1 = borders[0][0]; // first sewing halfedge
	double thresh = 100.0;
	int borderIndex;
	
	double distMin = std::numeric_limits<double>::max();
	/*for(unsigned b=1;b<borders.size();++b) // find the corresponding halfedge on other borders
	{
		for(unsigned h =0; h < borders[b].size(); ++h)
		{*/
	for(auto h = m->border_halfedges_begin(); h!=m->halfedges_end();++h)
	{
			//double dist = CGAL::squared_distance(c1->vertex()->point(),borders[b][h]->vertex()->point());
			double dist = CGAL::squared_distance(c1->vertex()->point(),h->vertex()->point());
			if(dist<distMin)// && dist < thresh)
			{
				distMin = dist;
				c2 = h;
			//	borderIndex = b;
			}
		//}
	}
	std::vector<simplePolyhedron::Halfedge_handle> cBord1; cBord1.push_back(c1);
	std::vector<simplePolyhedron::Halfedge_handle> cBord2; cBord2.push_back(c2);
	
	Vector V1 = c1->vertex()->point() - c1->opposite()->vertex()->point(); 
	Vector V2 = c2->vertex()->point() - c2->opposite()->vertex()->point();
	if( V1*V2 < 0) {c2 = c2->opposite();}
	
	
	
	p1 = c1; f1 = c1;
	p2 = c2; f2 = c2;
	
	bool endB1 = false;
	bool endB2 = false;
	
	do
	{
		// Find which halfedge to increment (c1 or c2)
		simplePolyhedron::Halfedge_around_vertex_circulator p1N = p1->vertex()->vertex_begin();
		simplePolyhedron::Halfedge_around_vertex_circulator p2N = p2->vertex()->vertex_begin();
		double distMinB1 = std::numeric_limits<double>::max();
		double distMinB2 = std::numeric_limits<double>::max();
		
		// Which hE on Border1 is closest to p2
		simplePolyhedron::Halfedge_handle candidateC1;
		do{
			if(!p1N->opposite()->is_border_edge()){continue;}
			if(p1N->opposite() == p1->opposite()){continue;} // make sure we really increment, and don't go backward
			double dist = CGAL::squared_distance(p2->vertex()->point(),p1N->opposite()->vertex()->point());
			//std::cout << dist << std::endl;
			if(dist < distMinB1)// && dist < thresh)
			{
				distMinB1 = dist;
				candidateC1 = p1N->opposite();
			}
		}
		while(++p1N != p1->vertex()->vertex_begin());
		
		// Which hE on Border2 is closest to p1
		simplePolyhedron::Halfedge_handle candidateC2;
		do{
			if(!p2N->opposite()->is_border_edge()){continue;}
			if(p2N->opposite() == p2->opposite()){continue;} // make sure we really increment, and don't go backward
			double dist = CGAL::squared_distance(p1->vertex()->point(),p2N->opposite()->vertex()->point());
			//std::cout << dist << std::endl;
			if(dist < distMinB2)// && dist < thresh)
			{
				distMinB2 = dist;
				candidateC2 = p2N->opposite();
			}
		}
		while(++p2N != p2->vertex()->vertex_begin());
		
		// Choose the less costly incrementation
		//std::cout << distMinB1 << " " << distMinB2 << std::endl;
		if(distMinB1 < distMinB2)
		{
			//std::cout << "Increment C1 " << std::endl;
			c1 = candidateC1;
		}
		else
		{
			//std::cout << "Increment C2 " << std::endl;
			c2 = candidateC2;
		}
		cBord1.push_back(c1); p1=c1;
		cBord2.push_back(c2); p2=c2;
		if(c1 == f1){endB1 = true;}
		if(endB1 && (c2==f2))
		{break;}
		if(c2 == f2){endB2 = true;}
		if(endB2 && (c1==f1))
		{
			break;
		}
	}
	while( true ); 
	
	surface_mesh_deformation deformation(*t); // Test w/ custom weights
	surface_mesh_deformation deformModel(*m);
	
	// Define the region of interest (here, the whole part)
	int count = 0;
	for(auto pHE = t->halfedges_begin(); pHE != t->halfedges_end(); ++pHE)
	{
		pHE->id() = count;
		count++;
	}
	count = 0;
	for(auto pVertex = t->vertices_begin(); pVertex != t->vertices_end(); ++pVertex)
	{
		pVertex->id() = count;
		deformation.insert_roi_vertex(pVertex);
		count++;
	}
	count = 0;
	for(auto pHE = m->halfedges_begin(); pHE != m->halfedges_end(); ++pHE)
	{
		pHE->id() = count;
		count++;
	}
	count = 0;
	for(auto pVertex = m->vertices_begin(); pVertex != m->vertices_end(); ++pVertex)
	{
		pVertex->id() = count;
		deformModel.insert_roi_vertex(pVertex);
		count++;
	}
	
	for(unsigned b = 0; b < cBord1.size();++b)
	{
		/*Point3d mem = cBord1[b]->vertex()->point();
		auto hnew = t->split_edge(cBord1[b]); hnew->vertex()->point() = cBord2[b]->vertex()->point();*/
		//hnew = m->split_edge(cBord2[b]); hnew->vertex()->point() = mem;
		Point3d mem = cBord1[b]->vertex()->point();
		deformation.insert_control_vertex(cBord1[b]->vertex());
		//cBord1[b]->vertex()->point() = cBord2[b]->vertex()->point();
		//cBord2[b]->vertex()->point() = mem;
		mem = cBord2[b]->vertex()->point();
		deformModel.insert_control_vertex(cBord2[b]->vertex());
		
	}
	
	// The definition of the ROI and the control vertices is done, call preprocess
	bool is_matrix_factorization_OK = deformation.preprocess();
	if(!is_matrix_factorization_OK)
	{
		std::cerr << "Error in preprocessing, check documentation of preprocess()" << std::endl;
		return;
	}
	is_matrix_factorization_OK = deformModel.preprocess();
	if(!is_matrix_factorization_OK)
	{
		std::cerr << "Error in preprocessing, check documentation of preprocess()" << std::endl;
		return;
	}
	
	for(unsigned c = 0; c < cBord1.size(); ++c)
	{
		deformation.set_target_position(cBord1[c]->vertex(),cBord2[c]->vertex()->point());
		deformModel.set_target_position(cBord2[c]->vertex(),cBord1[c]->vertex()->point());
	}
	
	/*std::cout <<  "Before deformation :";
	  bool intersecting = CGAL::Polygon_mesh_processing::does_self_intersect(*t,CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, *t)));
  std::cout
    << (intersecting ? "There are self-intersections." : "There is no self-intersection.")
    << std::endl;
    
    std::cout <<  "Before deformation :";
	  intersecting = CGAL::Polygon_mesh_processing::does_self_intersect(*m,CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, *m)));
  std::cout
    << (intersecting ? "There are self-intersections." : "There is no self-intersection.")
    << std::endl;*/
	
  	deformation.deform(10,0.01);
	deformModel.deform(10,0.01);
/*	
	std::cout <<  "After deformation :";
	  intersecting = CGAL::Polygon_mesh_processing::does_self_intersect(*t,CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, *t)));
  std::cout
    << (intersecting ? "There are self-intersections." : "There is no self-intersection.")
    << std::endl;
    
    std::cout <<  "After deformation :";
	  intersecting = CGAL::Polygon_mesh_processing::does_self_intersect(*m,CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, *m)));
  std::cout
    << (intersecting ? "There are self-intersections." : "There is no self-intersection.")
    << std::endl;
	
	
	simplePolyhedron::Halfedge_handle b1 = cBord1[0];
	if(b1->is_border_edge() && !(b1->is_border()))
	{
		b1 = b1->opposite();
	}
	
	simplePolyhedron::Halfedge_handle b2 = cBord2[0];
	if(b2->is_border_edge() && !(b2->is_border()))
	{
		b2 = b2->opposite();
	}
	
	if(b1->is_border()){
		std::vector<simplePolyhedron::Facet_handle> patch_facets;
		std::vector<simplePolyhedron::Vertex_handle> patch_vertices;
		bool success = CGAL::cpp11::get<0>(
			CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(
			*t,
			b1,
			std::back_inserter(patch_facets),
			std::back_inserter(patch_vertices),
			CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point,*t)).geom_traits(simpleKernel())));
	}
	if(b2->is_border()){
		std::vector<simplePolyhedron::Facet_handle> patch_facets;
		std::vector<simplePolyhedron::Vertex_handle> patch_vertices;
		bool success = CGAL::cpp11::get<0>(
			CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(
			*m,
			b2,
			std::back_inserter(patch_facets),
			std::back_inserter(patch_vertices),
			CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point,*m)).geom_traits(simpleKernel())));
	}*/
	
	
	PolyhedronPtr pol = convertToEnrichedPolyhedron(t);
	PolyhedronPtr polM = convertToEnrichedPolyhedron(m);
	
	pol->compute_normals();
	polM->compute_normals();
	v->addFrame();
	v->getScenePtr()->add_polyhedron(pol);
	v->addFrame();
	v->getScenePtr()->add_polyhedron(polM);
}

void SegmentController::alignSegments(Viewer* v, PolyhedronPtr s, PolyhedronPtr t,int sourceFrameID, int targetFrameID)
{
	//fillHoles(s);
	//fillHoles(t);
	
	//s->keep_largest_connected_components(1);
	//t->keep_largest_connected_components(1);
	unsigned sizeS = s->size_of_vertices();
	unsigned sizeT = t->size_of_vertices();
	colorMesh(s,1,1,0);
	colorMesh(t,0,1,1);
	
	double* S = new double[3*sizeS];
	double* T = new double[3*sizeT];
	
	auto sVertex = s->vertices_begin();
	for(unsigned vS=0; vS < sizeS; ++vS,++sVertex)
	{
		Point3d wcp = toWorld(v,sourceFrameID,sVertex->point());
		S[vS*3] = wcp.x();
		S[vS*3+1] = wcp.y();
		S[vS*3+2] = wcp.z();
		++sVertex;
	}
	
	auto tVertex = t->vertices_begin();
	for(unsigned vT=0; vT < sizeT; ++vT,++tVertex)
	{
		Point3d wcp = toWorld(v,targetFrameID,tVertex->point());
		T[vT*3] = wcp.x();
		T[vT*3+1] = wcp.y();
		T[vT*3+2] = wcp.z();
		++tVertex;
	}
	
	//Initialization
	Matrix R = Matrix::eye(3);
	Matrix P(3,1);
	IcpPointToPoint icp2(S,sizeS,3);
	//IcpPointToPlane icp(S,sizeS,3);
	//icp.fit(T,sizeT,R,P,-1);
	icp2.fit(T,sizeT,R,P,-1);
	
	std::cout << endl << "Transformation results:" << endl;
	std::cout << "R:" << endl << R << endl << endl;
	std::cout << "t:" << endl << P << endl << endl;
	
	qglviewer::Quaternion Q;
	double rData[3][3];
	for(unsigned i=0;i<3;++i)
	{
		for(unsigned j=0;j<3;++j)
		{
			rData[i][j] = R.val[i][j];
		}
	}
	Q.setFromRotationMatrix(rData);

	//v->frame(targetFrameID)->translate(P.val[0][0],P.val[1][0],P.val[2][0]);
	//v->frame(targetFrameID)->rotate(Q);
	
	v->frame(targetFrameID)->setTranslation(P.val[0][0],P.val[1][0],P.val[2][0]);
	v->frame(targetFrameID)->setRotation(Q);
	
	v->update();
	v->recreateListsAndUpdateGL();
	
}

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

void visitBorder(Halfedge_handle h, std::map<Halfedge_handle,bool> & isVisited, std::vector<Halfedge_handle> & border)
{
	Vertex_iterator hVertex = h->vertex();
	
	if(isVisited[h] || isVisited[h->opposite()])// halfedge has been visited
	{
		return;
	}
	else if ( h->is_border() )
	{
		border.push_back(h->opposite());
		isVisited[h] = true;
		isVisited[h->opposite()] = true;
		Halfedge_around_vertex_circulator hNext = hVertex->vertex_begin();
		
		do
		{
			Halfedge_handle hNextOp = hNext->opposite();
			visitBorder(hNextOp,isVisited,border);
		}
		while(++hNext!=hVertex->vertex_begin());
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
	unsigned m = 0;
	
	std::vector<double> coords;
	for(auto pVertex = p->vertices_begin();
	    pVertex != p->vertices_end();++pVertex)
	{
		coords.push_back(pVertex->point().x());
		coords.push_back(pVertex->point().y());
		coords.push_back(pVertex->point().z());
		pVertex->id() = m;
		++m;
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
	//delete p;
	return sp;
}

// Compute world coordinates using QGL operations, returns a CGAL::Point3d
Point3d toWorld(Viewer * v, int p, Point3d lcp)
{
	Point3d wcp;
	float x,y,z;
	double a,b,c,w;
	v->frame(p)->getPosition(x,y,z); qglviewer::Vec T(x,y,z);
	v->frame(p)->getOrientation(a,b,c,w); qglviewer::Quaternion Q(a,b,c,w);
	qglviewer::Vec V = Q * qglviewer::Vec(lcp.x(),lcp.y(),lcp.z()) + T;
	wcp = Point3d(V[0],V[1],V[2]);
	return wcp;
}

// Compute local coordinates for a point in world coordinate, in the frame #p
Point3d toPartFrame(Viewer* v, int p, Point3d wcp)
{
	Point3d lcp;
	float x,y,z;
	double a,b,c,w;
	v->frame(p)->getTranslation(x,y,z); qglviewer::Vec T(x,y,z);
	v->frame(p)->getRotation(a,b,c,w); qglviewer::Quaternion Q(a,b,c,w);
	qglviewer::Vec V = Q.inverseRotate(qglviewer::Vec(wcp.x(),wcp.y(),wcp.z()) - T);
	lcp = CGAL::ORIGIN + Vector(V[0],V[1],V[2]);
	return lcp;
}

void collectVertsAndFaces(PolyhedronPtr p, std::vector<double> & coords, std::vector<int> & faces)
{
	unsigned vertexOffset = coords.size() / 3.0;
	
	for(auto pVertex = p->vertices_begin();
	    pVertex != p->vertices_end();++pVertex)
	{
		coords.push_back(pVertex->point().x());
		coords.push_back(pVertex->point().y());
		coords.push_back(pVertex->point().z());
	}// collect all vertices for builder
	
	for(auto pFacet = p->facets_begin();
	pFacet!=p->facets_end();++pFacet)
	{
		auto h = pFacet->facet_begin();
		do
		{
			faces.push_back(h->vertex()->tag()+vertexOffset);
		}
		while(++h!=pFacet->facet_begin());
	}
}

void test_softICP_SVD()
{
	
	// Create dummy triangle data
	std::vector<double> coords;
	std::vector<int> faces;
	
	coords.push_back(0.0); coords.push_back(0.0); coords.push_back(0.0); // P1
	coords.push_back(2.0); coords.push_back(0.0); coords.push_back(0.0); // P2
	coords.push_back(1.0); coords.push_back(1.5); coords.push_back(0.0); // P3
	faces.push_back(0); faces.push_back(1); faces.push_back(2);
	
	
	PolyhedronPtr model(new Polyhedron);
	polyhedron_builder<HalfedgeDS> builder(coords,faces);
	model->delegate(builder);
	
	coords.clear();
	faces.clear();
	
	coords.push_back(1.0); coords.push_back(2.0); coords.push_back(-0.5); // Q1
	coords.push_back(1.0); coords.push_back(4.0); coords.push_back(-0.5); // Q2
	coords.push_back(-0.5); coords.push_back(3.0); coords.push_back(-0.5); // Q3

	faces.push_back(0); faces.push_back(1); faces.push_back(2);

	PolyhedronPtr target(new Polyhedron);
	polyhedron_builder<HalfedgeDS> buildert(coords,faces);
	target->delegate(buildert);
	
	softICPController softicp(model,target);
	
	std::map<Vertex_handle,Vertex_handle> phi;
	auto qVertex = target->vertices_begin();
	for(auto pVertex = model->vertices_begin(); 
	    pVertex!= model->vertices_end();
		++pVertex,++qVertex)
	{
		phi[pVertex] = qVertex;
	}
		
		
	std::vector<pointTransformation> transf;
	for(auto pVertex = model->vertices_begin();
	    pVertex!= model->vertices_end();
		++pVertex)
	{    
		std::vector<Vertex_handle> N;
		std::vector<Vertex_handle> phiN;
		for(auto pVertex = model->vertices_begin();
		    pVertex != model->vertices_end();
			++pVertex)
		{
			N.push_back(pVertex);
			phiN.push_back(phi[pVertex]);
		}
		pointTransformation ti = softicp.computeTransformation(N,phiN);
		transf.push_back(ti);
		std::cout<<std::endl;
	}
	
	int i = 0;
	for(auto pVertex = model->vertices_begin();
	    pVertex!= model->vertices_end();
		++pVertex)
	{
		Point3d p = pVertex->point();
		std::cout << p.x() << " " << p.y() << " " << p.z() << std::endl;
		softicp.applyTransformation(pVertex,transf[i],1,1);
		p = pVertex->point();
		std::cout << p.x() << " " << p.y() << " " << p.z() << std::endl;
		Point3d q = phi[pVertex]->point();
		std::cout << q.x() << " " << q.y() << " " << q.z() << std::endl;
		std::cout << std::endl;
		i++;
	}
	
}
