#ifndef Correspondence_ITEMS_H
#define Correspondence_ITEMS_H

/*!
 * \file Correspondence_Items.h
 * \brief Enrichement of the Facets, Halfedges and vertices of the polyhedra
 * \author Vincent Leon, CRIStAL
   \date 2015
 */

#include <mepp_config.h>
#ifdef BUILD_component_Correspondence

#include "../../../../mepp/Polyhedron/polyhedron_shared_items.h"

/**
 \class	Correspondence_Facet

 \brief	Enriches the Facets of a Polyhedra

 */
template <class Refs, class T, class P, class Norm, class Plane>
class Correspondence_Facet : virtual public MEPP_Common_Facet<Refs, T, Norm>
{
	public:

		Correspondence_Facet() {}
};


/*!
 * \class Correspondence_Halfedge
 * \brief Enriches the Halfedges of a Polyhedra
 */
template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class Correspondence_Halfedge : virtual public MEPP_Common_Halfedge<Refs,Tprev,Tvertex,Tface,Norm>
{
	public:

		Correspondence_Halfedge() {}
};


/*!
 * \class Correspondence_Vertex
 * \brief Enriches the Vertices of a Polyhedra
 */
template <class Refs, class T, class P, class Norm>
class Correspondence_Vertex : virtual public MEPP_Common_Vertex<Refs, T, P, Norm>
{
	public:

		Correspondence_Vertex() {}
		
		std::vector<double> &  getSemantic()
		{
			return m_semantic;
		}
		void setSemantic( std::vector<double> & sem)
		{
			m_semantic = sem;
		}
		
		void pushSemantic( double semEl)
		{
			m_semantic.push_back(semEl);
		}
		
		void setIndex( int index)
		{
			m_index = index;
		}
		
		int getIndex() 
		{
			return m_index;
		}
		
	private :
		std::vector<double> m_semantic;
		int m_index;
		
};
/*!
 * \class Correspondence_Polyhedron
 * \brief Enriches the Polyhedron
 */
template <class kernel, class items>
class Correspondence_Polyhedron : virtual public MEPP_Common_Polyhedron<kernel,items>
{
	public:

		Correspondence_Polyhedron() {}
};

#endif

#endif // Correspondence_ITEMS_H
