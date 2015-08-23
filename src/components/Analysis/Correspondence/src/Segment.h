#ifndef SEGMENT_H
#define SEGMENT_H
#include <vector>
#include <CGAL/Point_3.h>


namespace Analysis
{

  
  struct vec3
  {
    double m_x;
    double m_y;
    double m_z;
    
    bool operator==(vec3 const&v) const
    {
      return (this->m_x == v.m_x)
      && (this->m_y == v.m_y)
      && (this->m_z == v.m_z);
    }
  };
  
  struct face
  {
    vec3 m_f;
    vec3 m_s;
    vec3 m_t;
    
    
  };
  

  struct find_vec : std::unary_function<vec3, bool> {
    vec3 m_v;
    find_vec(vec3 v):m_v(v){}
    bool operator()(vec3 const& v) const 
    {
        return ( (v.m_x==m_v.m_x) && (v.m_y==m_v.m_y) && (v.m_z==m_v.m_z) );
    }  
};

  struct find_face : std::unary_function<face, bool> {
    face m_face;
    find_face(face f):m_face(f){}
    bool operator()(face const& f) const 
    { 
        return ((f.m_f == m_face.m_f) && (f.m_s == m_face.m_s) && (f.m_t == m_face.m_t))
	||     ((f.m_f == m_face.m_f) && (f.m_s == m_face.m_t) && (f.m_t == m_face.m_s))
	||     ((f.m_f == m_face.m_s) && (f.m_s == m_face.m_t) && (f.m_t == m_face.m_f))
	||     ((f.m_f == m_face.m_s) && (f.m_s == m_face.m_f) && (f.m_t == m_face.m_t))
	||     ((f.m_f == m_face.m_t) && (f.m_s == m_face.m_f) && (f.m_t == m_face.m_s)) 
	||     ((f.m_f == m_face.m_t) && (f.m_s == m_face.m_s) && (f.m_t == m_face.m_f)) ;
    }  
};

class Segment
{
    public:
	
        Segment() :counter(0), m_center(), m_label(-1), m_ipHist(), m_faceMap(), m_sampledHarmonics(), m_harmonics(), m_faces(), m_verts()
        {}

        virtual ~Segment()
        {}

        int getLabel() const
        {
          return m_label;
        }

        void setLabel(const int label)
        {
          m_label = label;
        }

        vec3 & getCenter()
        {
          return m_center;
        }

        void setCenter(double x, double y, double z)
        {
          m_center.m_x = x;
          m_center.m_y = y;
          m_center.m_z = z;
        }

        void addVertice(double x, double y, double z)
        {
          vec3 v;
          v.m_x = x;v.m_y = y;v.m_z = z;
          m_verts.push_back(v);
        }
        
        int addFace(vec3 v[3])
	{
	  face f;
	  f.m_f = v[0];
	  f.m_s = v[1];
	  f.m_t = v[2];
	  m_faces.push_back(f);
	  m_faceMap.push_back(0);
	  /*m_verts.push_back(v[0]);
	  m_verts.push_back(v[1]);
	  m_verts.push_back(v[2]);*/
	}

        void initHarmonics(int nbLabels)
        {
          for(int i=0;i<nbLabels;++i)
          {
            m_harmonics.push_back(new std::vector<double>());
            for(int j=0;j<m_verts.size();++j)
            {
              m_harmonics.back()->push_back(0.0);
            }
          }
        }
        
        void initSampledHarmonics(int nbLabels, int nbSamples)
	{
	  for(unsigned l=0;l<nbLabels;++l)
          {
            m_sampledHarmonics.push_back(new std::vector<double>());
            for(int s=0;s<nbSamples;++s)
            {
              m_sampledHarmonics.back()->push_back(0.0);
            }
          }
	}

        void addHarmonicValue(int label, double value, double x, double y, double z)
        {
          vec3 v; v.m_x = x; v.m_y = y; v.m_z = z;
          int vertID = std::distance(m_verts.begin(),std::find_if(m_verts.begin(),m_verts.end(),find_vec(v)));
          if(vertID == m_verts.size())
          {
            return;
          }
          m_harmonics[label]->at(vertID) = value;
	  return;	
        }
        
        void addHarmonicSample(int label, int sample, double value, double x, double y, double z)
	{
	  vec3 v; v.m_x = x; v.m_y = y; v.m_z = z;
	  int vertID = std::distance(m_verts.begin(),std::find_if(m_verts.begin(),m_verts.end(),find_vec(v)));
	  if(vertID == m_verts.size())
          {
            return;
          }
	  m_sampledHarmonics[label]->at(sample) = value;
	  return;
	}
	
	bool isOnSegment(double x, double y, double z)
	{
	  vec3 v; v.m_x = x; v.m_y = y; v.m_z = z;
	  int vertID = std::distance(m_verts.begin(),std::find_if(m_verts.begin(),m_verts.end(),find_vec(v)));
	  return (vertID != m_verts.size());
	}
	
	bool isOnSegment(face f)
	{
	  int faceID = std::distance(m_faces.begin(),std::find_if(m_faces.begin(),m_faces.end(),find_face(f)));
	  return (faceID != m_faces.size());
	}
	
        double getHarmonicValue(int label, vec3 v)
        {
          int vertID = std::distance(m_verts.begin(),std::find_if(m_verts.begin(),m_verts.end(),find_vec(v)));
          return m_harmonics[label]->at(vertID);
        }
        
        int findFace(face f)
	{
	  
	  int faceID = std::distance(m_faces.begin(),std::find_if(m_faces.begin(),m_faces.end(),find_face(f)));
	  if(faceID == m_faces.size()){return -1;}
	  return faceID;
	}
	

        std::vector<vec3> m_verts;
	std::vector<face> m_faces;
	
        std::vector<std::vector<double>* > m_harmonics;
	std::vector<std::vector<double>* > m_sampledHarmonics;
	std::vector<unsigned> m_faceMap;
	
	int m_segID;
	
	std::vector<long> m_ipHist; // Shape Index Representation
	
    private:
        int m_label;
        vec3 m_center;
        int counter;
	

};
}

#endif // SEGMENT_H
