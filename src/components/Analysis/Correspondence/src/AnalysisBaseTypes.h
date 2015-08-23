#ifndef __AnalysisBaseTypes_h
#define __AnalysisBaseTypes_h
#include <iostream>

#include "DataStructures/CGAL/EnrichedPolyhedron/CGAL_EnrichedPolyhedron.hxx"
#include "Kernel/MeshFilter.hxx"

namespace Analysis
{
///////////////////////////////////////////////////////////////////////////////
	typedef DataStructures::MyOpenMesh_Mesh		OpenMesh_Mesh_Wraped_Type;
	typedef OpenMesh_Mesh_Wraped_Type::MeshType	OpenMesh_Polyhedron ;

	typedef DataStructures::MyCGAL_Mesh			CGAL_Mesh_Wraped_Type;
	typedef CGAL_Mesh_Wraped_Type::MeshType			CGAL_Polyhedron ;
///////////////////////////////////////////////////////////////////////////////
	// GCC does not inherit exported types from templated base class
	// => therefore we need a non-templated base class!
	class BasePolyhedronTraitsCGAL
	{
		public :
        typedef CGAL_Mesh_Wraped_Type::MeshType					MeshType;
	typedef CGAL_Mesh_Wraped_Type::Vector					Vector;
	typedef CGAL_Mesh_Wraped_Type::Vertex					Vertex;
	typedef CGAL_Mesh_Wraped_Type::Point					Point3d;
	typedef CGAL_Mesh_Wraped_Type::Vertex_handle				Vertex_handle;
	typedef CGAL_Mesh_Wraped_Type::Vertex_iterator				Vertex_iterator;
	typedef CGAL_Mesh_Wraped_Type::Facet_iterator				Facet_iterator;
	typedef CGAL_Mesh_Wraped_Type::Facet_handle				Facet_handle;
	typedef CGAL_Mesh_Wraped_Type::Halfedge_around_facet_circulator		Halfedge_around_facet_circulator;
	typedef CGAL_Mesh_Wraped_Type::Halfedge_around_vertex_circulator	Halfedge_around_vertex_circulator;
	typedef CGAL_Mesh_Wraped_Type::Halfedge					Halfedge;
	typedef CGAL_Mesh_Wraped_Type::Halfedge_handle				Halfedge_handle;
	typedef CGAL_Mesh_Wraped_Type::Halfedge_iterator			Halfedge_iterator;
    };

    class BasePolyhedronTraitsOpenMesh
    {
        public:
        typedef OpenMesh_Mesh_Wraped_Type::MeshType MeshType;
        typedef MeshType::Vertex Vertex;
        typedef MeshType::Point Point3d;
        typedef MeshType::VertexHandle Vertex_handle;
        typedef MeshType::VertexIter Vertex_iterator;
    };
}

#endif
