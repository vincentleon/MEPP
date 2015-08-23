/*!
	\file Correspondence_Component.h
 */


#ifndef Correspondence_COMPONENT_H
#define Correspondence_COMPONENT_H

#include <mepp_config.h>
#ifdef BUILD_component_Correspondence

#include "../../../../mepp/mepp_component.h"

#include "Correspondence_Polyhedron.h"
#include "Shape.h"
#include "geodesic/geodesic_algorithm_exact.h"

#include <vector>

/*!
 * \class Correspondence_Component
 * \brief Correspondence computation on polyhedra
 */
class Correspondence_Component :
  public mepp_component
{
	public:

		/*!
		* \brief Constructor
		*/
		Correspondence_Component(Viewer* v, PolyhedronPtr p);

		/*!
		* \brief Destructor
		*/
		~Correspondence_Component() {}
		
		
		
	private : 
		
		int m_nbLabel;
	
		//std::vector<double> m_centreDescr;
		Vertex_handle m_centreSelection;
		double m_isolineValue;
		
		std::vector<double> m_maxVector;
		
		//std::vector<double> m_ellipse;
		
		bool m_colorCompare2Source;
		bool m_learningMode;
		bool m_saveFeatures;
		bool m_selectionMode;
		bool m_readSelection;
		
		std::vector<Vertex_handle> & m_selection;
		std::map<Vertex_handle,int> m_tag;
		Analysis::Shape m_Shape;
		geodesic::Mesh m_gmesh;
		geodesic::GeodesicAlgorithmExact * m_geoAlg;
		
		void initGeodesicMesh(PolyhedronPtr p);
		
		void saveDescriptor(PolyhedronPtr p);
		void readDescriptor(PolyhedronPtr p);
		void initMaxVector(PolyhedronPtr p);
		void normalize(std::vector<double> & descr);
		void computeDescriptorAllVertices(PolyhedronPtr p);
		void computeGeodesicDistancesAllVertices(PolyhedronPtr p);
		
		std::vector<double> getClosetVertexDescriptor(PolyhedronPtr p, Point3d pt);
		
		void compareToDescrEllipse(PolyhedronPtr p, std::vector<double> & ellipse);
	
		void initializeEllipsoid(PolyhedronPtr p);
		void readSelectionBasedOnColor(PolyyhedronPtr p);
		
		Vertex_handle getSelectionCenter();
		Vertex_handle getFurtherFromSelectionCenter();
		
		
		void tagSelectionVertices(PolyhedronPtr p);
};

double L2Dist(std::vector<double> & descr1, std::vector<double> & descr2);

#endif

#endif // Correspondence_COMPONENT_H
