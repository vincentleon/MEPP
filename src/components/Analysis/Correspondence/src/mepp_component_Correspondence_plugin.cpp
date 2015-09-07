#include <mepp_config.h>
#ifdef BUILD_component_Correspondence

#include "mepp_component_Correspondence_plugin.hxx"

#include "dialSettings.hxx"

#include <QObject>
#include <QAction>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>
#include <QMdiSubWindow>
#include <QGLFramebufferObject>

#include "Correspondence_Component.h"

typedef boost::shared_ptr<Correspondence_Component> Correspondence_ComponentPtr;

void mepp_component_Correspondence_plugin::post_draw()
{
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
	}
	
	
	
}


void mepp_component_Correspondence_plugin::pre_draw()
{
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
		//polyhedron_ptr->set_index_vertices();
		
		/*glClearStencil(0);
		glEnable(GL_STENCIL_TEST);
		glStencilOp(GL_KEEP,GL_KEEP,GL_REPLACE);
		glPushMatrix();
		
		::glBegin(GL_POINTS);
		::glColor4f(0.0,0.0,0.0,0.0);
		for (Vertex_iterator pVertex = polyhedron_ptr->vertices_begin();
			pVertex !=	polyhedron_ptr->vertices_end();	pVertex++)
		{
			glStencilFunc(GL_ALWAYS,pVertex->tag()+1,0);
			std::cout << pVertex->tag() << std::endl;
			const Point3d& p1 = pVertex->point();
			::glVertex3f(p1[0],p1[1],p1[2]);
		}
		::glEnd();
		
		glPopMatrix();
		glDisable(GL_STENCIL_TEST);*/
		this->PaintStart();
		
	}
}


void mepp_component_Correspondence_plugin::OnMouseLeftDown(QMouseEvent *event)
{
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
		
		if (!doesExistComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewer, polyhedron_ptr))
		{
			int x,y,h,w;
			x = event->x();
			y = event->y();
			h = viewer->height();
			w = viewer->width();
			
			unsigned char color[3];
			//m_fbo->bind();
			glReadPixels(x,h-y,1,1,GL_RGB,GL_UNSIGNED_BYTE,color);
			std::cout << "x,y :" << x << " " << y << std::endl;
			std::cout << "color : " << (int)color[0] << " " << (int)color[1] << " " << (int)color[2]<< std::endl;
			unsigned index = (color[0] << 16) | (color[1] << 8) | color[2];
			if(index ==  0xFFFFFF){index = -1;}
			std::cout <<"index : "<< index << "\n";
			//m_fbo->release();
		}
 	}
}

void mepp_component_Correspondence_plugin::OnMouseMotion(QMouseEvent* event)
{
	if (mw->activeMdiChild() != 0)
	{
		
	}
}

void mepp_component_Correspondence_plugin::OnCorrespondence()
{
	if (mw->activeMdiChild() != 0)
	{
		m_hasNotBeenPainted = true;
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewer, polyhedron_ptr);
		{
			int nbLabel = 8+3;
			int meshID = 10;

			SettingsDialog dial;
			if (dial.exec() == QDialog::Accepted)
			{
				QApplication::setOverrideCursor(Qt::WaitCursor);

				mw->statusBar()->showMessage(tr("Correspondence..."));
				
				//component_ptr->readSelectionBasedOnColor(polyhedron_ptr);
				
				component_ptr->initParameters(nbLabel,meshID);
				
				component_ptr->readDescriptor(polyhedron_ptr);
				
				component_ptr->initializeEllipsoid(polyhedron_ptr);
				
				component_ptr->computeEllipseParameters(polyhedron_ptr);
				
				component_ptr->compareDescriptorToEllipse(polyhedron_ptr);
				
				mw->statusBar()->showMessage(tr("Correspondence is done"));

				component_ptr->set_init(2);
				viewer->recreateListsAndUpdateGL();
				
			}
		}
	}
	QApplication::restoreOverrideCursor();
	
}

void mepp_component_Correspondence_plugin::PaintStart()
{
	if (mw->activeMdiChild() != 0)
	{
		
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
	
		
		//m_fbo = new QGLFramebufferObject(viewer,viewer,QGLFramebufferObject::Depth);
		//m_fbo->bind();
		glDisable(GL_LIGHTING);
		glClearColor(0.0,0.0,0.0,1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		::glBegin(GL_TRIANGLES);
		int i = 0;
		for (Facet_iterator pFacet = polyhedron_ptr->facets_begin();
		     pFacet!= polyhedron_ptr->facets_end(); pFacet++)
		     {
			Halfedge_around_facet_circulator hE = pFacet->facet_begin();
			
			do{
				Vertex_handle pVertex = hE->opposite()->vertex();
				//unsigned id = pVertex->tag();
				unsigned id = i;
				
				unsigned char color[3];
				color[0] = (id>>16) & 0xFF;
				color[1] = (id>>8) & 0xFF; 
				color[2] = id & 0xFF;
				glColor3ubv(color);
				Point3d p = pVertex->point();
				glVertex3f(p.x(),p.y(),p.z());
				i++;
				
			}while(++hE!=pFacet->facet_begin());
			
			
		}
		::glEnd();
		glPopMatrix();
		//m_fbo->release();
	}
	
}



#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(mepp_component_Correspondence_plugin, mepp_component_Correspondence_plugin);
#endif

#endif
