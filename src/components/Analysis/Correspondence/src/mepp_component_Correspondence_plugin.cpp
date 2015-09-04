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
	}
}


void mepp_component_Correspondence_plugin::OnMouseLeftDown(QMouseEvent *event)
{
	/*if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		if (doesExistComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewer, polyhedron_ptr)) // important !!!
		{
			std::cout<< " PAUL " << std::endl;
			if(m_hasNotBeenPainted)
			{
				int x,y,h,w;
				x = event->x();
				y = event->y();
				h = viewer->height();
				w = viewer->width();
				unsigned index;
				glReadPixels(x,y,1,1,GL_STENCIL_INDEX,GL_UNSIGNED_INT,&index);
				
				std::cout << index -1 << std::endl;
			}
		}
 	}*/
}

void mepp_component_Correspondence_plugin::OnMouseMotion(QMouseEvent* event)
{
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
		
		std::cout << "Paul ";
		
		if (doesExistComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewer, polyhedron_ptr)) // important !!!
		{
			std::cout << "George ";
			if(true)
			{
				/*std::cout << "Ringo ";
				int x,y,h,w;
				x = event->x();
				y = event->y();
				h = viewer->height();
				w = viewer->width();
				unsigned index;
				glReadPixels(x,y,1,1,GL_STENCIL_INDEX,GL_UNSIGNED_INT,&index);
				std::cout << index -1 << " ";*/
			}
		}
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
		/*Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
		m_fbo = new QGLFramebufferObject(viewer->width(),viewer->height());
		m_fbo->bind();
		
		glPushMatrix();
		
		::glBegin(GL_TRIANGLES);
		for (Facet_iterator pFacet = polyhedron_ptr->facets_begin();
		     pFacet!= polyhedron_ptr->facets_end(); pFacet++)
		     {
			Halfedge_around_facet_circulator hE = pFacet->vertices_begin();
			
			do{
				Vertex_handle pVertex = hE->opposite()->vertex();
				unsigned id = pVertex->tag();
				unsigned char color[3];
				color[0] = (id>>16) & 0xFF;
				color[1] = (id>>8) & 0xFF; 
				color[2] = id & 0xFF;
				glColor3ubv(color);
				
			}while(++hE!=pFacet->vertices_begin());
			
			
		}
		::glEnd();
		
		glPopMatrix();
		
		m_fbo->release();
		*/
	}
	
}



#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(mepp_component_Correspondence_plugin, mepp_component_Correspondence_plugin);
#endif

#endif
