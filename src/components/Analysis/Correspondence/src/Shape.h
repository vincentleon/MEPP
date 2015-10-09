#ifndef SHAPE_H
#define SHAPE_H

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/concept_check.hpp>
#include <boost/graph/graph_concepts.hpp>

#include "Segment.h"

namespace Analysis
{

class Shape
{
  public:
      Shape() : m_faceMap(), m_vertices(), m_faces(), m_segments()
      {}

      virtual ~Shape()
      {}

      void addSegment(Segment & segment)
      {
        m_segments.push_back(&segment);
      }

      int getSizeSegment() const
      {
        return m_segments.size();
      }

      Segment & getSegment(int i)
      {
        return *m_segments[i];

      }

      void initHarmonics(int nblabels)
      {
        for(int i=0;i<m_segments.size();++i)
        {
          m_segments[i]->initHarmonics(nblabels);
        }
      }

      void addHarmonicValue(int seg, int label, double value, double x, double y, double z)
      {
        for(int i = 0; i<m_segments.size(); ++i) // test for all segs
        {
          m_segments[i]->addHarmonicValue(label,value,x,y,z);
        }
      }
      
      int isOnSegment(face f)
      {
	  for(unsigned s=0;s<m_segments.size();++s)
	  {
	    if(m_segments[s]->isOnSegment(f))
	    {
	      return s;
	    }
	  }
	  return -1;
      }
       
      int isOnSegment(double x, double y, double z)
      {
	for(unsigned s=0;s<m_segments.size();++s)
	{
	  if(m_segments[s]->isOnSegment(x,y,z))
	  {
	    return s;
	  }
	}
	return -1;
      }
       
      void addFace(vec3 faces[3])
      {
	face f;
	f.m_f = faces[0];
	f.m_s = faces[1];
	f.m_t = faces[2];
	m_faces.push_back(f);
      }
      
      void addVertice(vec3 v)
      {
        m_vertices.push_back(v);
      }
      
      void buildCorrespondance()
      {
	for(unsigned f=0; f<m_faces.size(); ++f)
	{
	  int s = isOnSegment(m_faces[f]);
	  int fs = m_segments[s]->findFace(m_faces[f]);
	  if(s == -1)
	  {continue;}
	  m_segments[s]->m_faceMap[fs] = f;
	}
      }
      
      void initFaceLabelsAndSegments()
      {
	int meshID = m_meshID;
	std::ifstream labels;
	std::stringstream ssLabels;
	ssLabels<<"/home/leon/quadruClean/result/"<<meshID<<".labels";
	labels.open(ssLabels.str().c_str());
	if(labels)
	{
		std::string line;
		while(getline(labels,line))
		{	
			m_faceLabels.push_back(atoi(line.c_str()));
		}
	}
	else
	{
		std::cout<<"\tUnable to open "<<ssLabels.str()<<std::endl;
	}
	labels.close();

	std::ifstream parts;
	std::stringstream ssParts;
	ssParts<<"/home/leon/quadruClean/result/"<<meshID<<".parts";
	parts.open(ssParts.str().c_str());
	if(parts)
	{
		std::string line;
		while(getline(labels,line))
		{
			m_faceSegments.push_back(atoi(line.c_str()));
		}
	}
	else
	{
		std::cout<<"\tUnable to open "<<ssParts.str()<<std::endl;
	}
	parts.close();
      }

      double getHarmonicValue(int seg, int label, vec3 v)
      {
        return m_segments[seg]->getHarmonicValue(label,v);
      }
      
      void pushSemantics(std::vector<double> & sem)
      {
	m_Semantics.push_back(sem);
      }
      
      void setNormalizer(std::vector<double> & normalizer)
      {
	m_Normalizer = normalizer;      
      }
      
      std::vector<double> & getNormalizer()
      {
	return m_Normalizer;      
      }   
      
      std::vector<double> getSemantic(int i) const
      {
	return m_Semantics[i];
      }
      
      std::vector<unsigned> m_faceMap;
      
      std::vector<vec3> m_vertices;
      std::vector<face> m_faces;
      std::vector<int> m_faceLabels;
      std::vector<int> m_faceSegments;
      
      int m_meshID;
      int m_categoryLabel;
      
      std::vector<std::vector<double> *> m_featureVectors;

  private:
      std::vector<Segment *> m_segments;
      std::vector<std::vector<double> > m_Semantics;
      std::vector<double> m_Normalizer;
};
}
#endif // SHAPE_H
