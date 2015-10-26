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

#include "svm.h"

#include "SegmentController.h"

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
		
		void initParameters(int nbLabel, int meshId, std::string meshDir);
		
		void learnDescriptor(PolyhedronPtr p, std::string meshDir);
		
		void readDescriptor(PolyhedronPtr p,std::string meshDir);
		
		void showDescriptor(PolyhedronPtr p, int dim);
		
		void readSelectionBasedOnColor(PolyhedronPtr p);
		
		void initializeEllipsoid(PolyhedronPtr p);
		
		void compareDescriptorToEllipse(PolyhedronPtr p);
		
		void compareDescriptorToEllipseRotation(PolyhedronPtr p);
		
		void compareDescriptorToGaussian(PolyhedronPtr p);
		
		void compareDescriptorToSVM(PolyhedronPtr p);
		
		void computeEllipseParameters(PolyhedronPtr p);
		
		void computeEllipseParametersRotation(PolyhedronPtr p);
		
		void computeGaussianParameters(PolyhedronPtr p);
		
		void learnSVMClassifier(PolyhedronPtr p);
		
		double computeEnergy(const std::vector<double> & ellipse);
		
		double computeEnergyRotation(const std::vector<double> & ellipse);
		
		double computeEnergyGaussian(const std::vector<double> & threshold);
		
		std::vector<double> getCentreDescr() const;
		
		void setCentreDescriptor(const std::vector<double> & centreDescr);
		
		std::vector<double> getEllipse() const;
		
		void setEllipse(const std::vector<double> & ellipse);
		
		myMatrix getMatrix() const;
		myVector getVector() const;
		myMatrix getInverseMatrix() const;
		double getDeterminant() const;
		double getThreshold() const;
		svm_model * getSVM() const;
		
		void setMatrix(const myMatrix & m);
		void setVector(const myVector & v);
		void setInverseMatrix(const myMatrix & m);
		void setDeterminant(double det);
		void setThreshold(double thresh);
		void setSVM(svm_model * svm);	
		
		SegmentController m_segCtr;
		
		void scaleMesh(Polyhedron::Iso_cuboid bbox, PolyhedronPtr p);
	private : 
		
		int m_nbLabel;
	
		Vertex_handle m_centreSelection;
		
		std::vector<double> m_centreDescriptor;
		
		double m_isolineValue;
		
		std::vector<double> m_maxVector;	
		
		std::vector<double> m_ellipse;
		
		myMatrix m_gaussianMatrix;
		myMatrix m_inverseGaussianMatrix;
		double m_gaussianDeterminant;
		double m_threshold;
		myVector m_mu;
		
		bool m_colorCompare2Source;
		bool m_learningMode;
		bool m_saveFeatures;
		bool m_selectionMode;
		bool m_readSelection;
		
		svm_model * m_svmModel;
		
		std::vector<Vertex_handle> m_selection;
		std::map<Vertex_handle,int> m_tag;
		Analysis::Shape m_Shape;
		geodesic::Mesh m_gmesh;
		geodesic::GeodesicAlgorithmExact * m_geoAlg;
		
		
		void initGeodesicMesh(PolyhedronPtr p);
		
		void saveDescriptor(PolyhedronPtr p,std::string meshDir);
		
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

double objectiveFunRotation(const std::vector<double> & ellipse, std::vector<double> & grad, void *data);

double objectiveFunGaussian(const std::vector<double> & threshold, std::vector<double> & grad, void *data);



#endif

#endif // Correspondence_COMPONENT_H
