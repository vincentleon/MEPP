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

#include "../components/Analysis/Curvature/src/Curvature_Component.h"

#include "../components/Analysis/Correspondence/src/SegmentController.h"

#include "svm.h"

#include <eigen3/Eigen/Eigenvalues>


typedef CGAL::Linear_algebraCd<double>::Matrix myMatrix;
typedef CGAL::Linear_algebraCd<double>::Vector myVector;


struct VertDistOrdering {
	bool operator ()(std::pair<Vertex_handle, double> const& a, std::pair<Vertex_handle, double> const& b)
	{
		return (a.second) < (b.second);
	}
};

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
		
		bool readDescriptor(PolyhedronPtr p,std::string meshDir, bool normalize=true);
		
		void showDescriptor(PolyhedronPtr p, int dim);
		
		void readSelectionBasedOnColor(PolyhedronPtr p);
		
		void initializeEllipsoid(PolyhedronPtr p);
		
		void compareDescriptorToEllipse(PolyhedronPtr p);
		
		void compareDescriptorWithSVM(PolyhedronPtr p, unsigned SVM_mode=0);
		
		void compareDescriptorToEllipseRotation(PolyhedronPtr p);
		
		void compareDescriptorToGaussian(PolyhedronPtr p);
		
		void compareDescriptorToSVM(PolyhedronPtr p);
		
		void computeEllipseParameters(PolyhedronPtr p);
		
		void computeEllipseParametersRotation(PolyhedronPtr p);
		
		void computeGaussianParameters(PolyhedronPtr p);
		
		void learnSVMClassifier(PolyhedronPtr p);
		
		void learnSVMPatch(PolyhedronPtr p, unsigned SVM_mode=0);
		
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
		double getRadius() const;
		int getNbCandidates() const;
		
		
		void setMatrix(const myMatrix & m);
		void setVector(const myVector & v);
		void setInverseMatrix(const myMatrix & m);
		void setDeterminant(double det);
		void setThreshold(double thresh);
		void setSVM(svm_model * svm);	
		void setRadius(double radius);
		void setNbCandidates(int nbCandidates);
		
		SegmentController m_segCtr;
		
		void scaleMesh(Polyhedron::Iso_cuboid bbox, PolyhedronPtr p);
		
		
		void saveDescriptor(PolyhedronPtr p,std::string meshDir);
		
		Analysis::Shape & getShape(){return m_Shape;}
		
		std::vector<double> featureMeans;
		std::vector<double> featureSTD;
		
	private : 
		
		int m_nbLabel;
	
		Vertex_handle m_centreSelection;
		
		std::vector<double> m_centreDescriptor;
		
		double m_isolineValue;
		
		double m_patchRadius;
		int m_nbCandidates;
		
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

		void initMaxVector(PolyhedronPtr p);
		void normalize(std::vector<double> & descr);
		void computeDescriptorAllVertices(PolyhedronPtr p);
		void computeGeodesicDistancesAllVertices(PolyhedronPtr p);
		
		std::vector<double> & getClosetVertexDescriptor(PolyhedronPtr p, Point3d pt);
		
		void getInRadius(PolyhedronPtr p, Vertex_handle c, double radius, std::set<Vertex_handle> & vertices);
		
		Vertex_handle getSelectionCenter();
		Vertex_handle getFurtherFromSelectionCenter();
		
		void tagSelectionVertices(PolyhedronPtr p);
		
		

};

double L2Dist(std::vector<double> & descr1, std::vector<double> & descr2);

double objectiveFun(const std::vector<double> & ellipse, std::vector<double> & grad, void *data);

double objectiveFunRotation(const std::vector<double> & ellipse, std::vector<double> & grad, void *data);

double objectiveFunGaussian(const std::vector<double> & threshold, std::vector<double> & grad, void *data);

bool sphere_clip_vector(Point3d &O, double r,const Point3d &P, Vector &V);

void vector_times_transpose_multi(long double pVector[3],long double ppMatrix[3][3], double coeff);

void addM(long double pMatrix[3][3], long double pMatrixSum[3][3]);

double fix_sin(double sine);

double areaFace(Facet_handle &f);

void principal_curv(Vertex_handle pVertex, long double ppMatrix_sum[3][3], int size);

void computePatchBasis(std::vector<Vertex_handle> & sel, Enriched_kernel::Vector_3 & U, Enriched_kernel::Vector_3 & V, Enriched_kernel::Vector_3 & N);

Halfedge_handle create_center_vertex_with_descriptor( PolyhedronPtr p, Facet_iterator f, int nbLabel);

#endif

#endif // Correspondence_COMPONENT_H
