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

#include <nlopt.hpp>

#include "CGAL/Linear_algebraCd.h"

typedef CGAL::Linear_algebraCd<double>::Matrix myMatrix;
typedef CGAL::Linear_algebraCd<double>::Vector myVector;


/*!
 * \class Correspondence_Component
 * \brief Correspondence computation on polyhedra
 */
class Correspondence_Component :
  public mepp_component
{
	public:

		Correspondence_Component(Viewer* v, PolyhedronPtr p);
		
		~Correspondence_Component() {}
		
		void initParameters(int nbLabel, int meshId);
		
		void learnDescriptor(PolyhedronPtr p);
		
		void readDescriptor(PolyhedronPtr p);
		
		void showDescriptor(PolyhedronPtr p, int dim);
		
		void selectPoint(PolyhedronPtr p);
		
		void paintRegion(PolyhedronPtr p);
		
		void readSelectionBasedOnColor(PolyhedronPtr p);
		
		void initializeEllipsoid(PolyhedronPtr p);
		
		void compareDescriptorToEllipse(PolyhedronPtr p);
		
		void compareDescriptorToGaussian(PolyhedronPtr p);
		
		void computeEllipseParameters(PolyhedronPtr p);
		
		void computeGaussianParameters(PolyhedronPtr p);
		
		double computeEnergy(const std::vector<double> & ellipse);
		
		std::vector<double> getCentreDescr() const;
		
		void setCentreDescriptor(const std::vector<double> & centreDescr);
		
		std::vector<double> getEllipse() const;
		
		void setEllipse(const std::vector<double> & ellipse);
		
	private : 
		
		int m_nbLabel;
	
		Vertex_handle m_centreSelection;
		
		std::vector<double> m_centreDescriptor;
		
		double m_isolineValue;
		
		std::vector<double> m_maxVector;	
		
		std::vector<double> m_ellipse;
		
		myMatrix m_gaussianMatrix;
		
		bool m_colorCompare2Source;
		bool m_learningMode;
		bool m_saveFeatures;
		bool m_selectionMode;
		bool m_readSelection;
		
		std::vector<Vertex_handle> m_selection;
		std::map<Vertex_handle,int> m_tag;
		Analysis::Shape m_Shape;
		geodesic::Mesh m_gmesh;
		geodesic::GeodesicAlgorithmExact * m_geoAlg;
		
		void initGeodesicMesh(PolyhedronPtr p);
		
		void saveDescriptor(PolyhedronPtr p);
		
		void initMaxVector(PolyhedronPtr p);
		void normalize(std::vector<double> & descr);
		void computeDescriptorAllVertices(PolyhedronPtr p);
		void computeGeodesicDistancesAllVertices(PolyhedronPtr p);
		
		std::vector<double> & getClosetVertexDescriptor(PolyhedronPtr p, Point3d pt);
		
		Vertex_handle getSelectionCenter();
		Vertex_handle getFurtherFromSelectionCenter();
		
		void tagSelectionVertices(PolyhedronPtr p);
		

};

double L2Dist(std::vector<double> & descr1, std::vector<double> & descr2);

double objectiveFun(const std::vector<double> & ellipse, std::vector<double> & grad, void *data);


#endif

#endif // Correspondence_COMPONENT_H
