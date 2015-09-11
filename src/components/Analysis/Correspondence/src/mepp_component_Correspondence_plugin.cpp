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
		
	}
}


void mepp_component_Correspondence_plugin::OnMouseLeftDown(QMouseEvent *event)
{
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
		
 	}
}

void mepp_component_Correspondence_plugin::OnMouseLeftUp(QMouseEvent *event)
{
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
		if(doesExistComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewer, polyhedron_ptr))
		{
			m_hasNotBeenPainted=false;
			std::cout << "HAS BEEN PAINTED " << std::endl;
		}
 	}
}

void mepp_component_Correspondence_plugin::OnMouseMotion(QMouseEvent* event)
{
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
		
 		if (doesExistComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewer, polyhedron_ptr))
		{
			Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewer, polyhedron_ptr);
			
			if(component_ptr->get_init() != 2)
			{
				int x,y,h,w;
				x = event->x();
				y = event->y();
				h = viewer->height();
				w = viewer->width();

				unsigned char color[3];
				m_fbo->bind();
				glReadPixels(x,h-y,1,1,GL_RGB,GL_UNSIGNED_BYTE,color);
				unsigned index = (color[0] << 16) | (color[1] << 8) | color[2];
				if(index ==  0xFFFFFF){index = -1;}
				else{
				Facet_iterator pFacet = m_facets[index];
				Halfedge_around_facet_circulator hE = pFacet->facet_begin();
				do{
					Vertex_handle pVertex = hE->opposite()->vertex();
					pVertex->color(1,0,0);
				}while(++hE!=pFacet->facet_begin());
				}
				m_fbo->release();
				//std::cout <<"index : "<< index << "\n";
				//m_paintedFacets.insert(index);
			}
		}
		viewer->recreateListsAndUpdateGL();
	}
}

void mepp_component_Correspondence_plugin::OnCorrespondence()
{
	
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
		
		Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewer, polyhedron_ptr);
		{
			int nbLabel = 8;
			int meshID = 10;

			SettingsDialog dial;
			if (dial.exec() == QDialog::Accepted)
			{
				if(m_hasNotBeenPainted)
				{
					std::cout << "HAS NOT BEEN PAINTED " << std::endl;
					//PolyhedronPtr polyhedron_ptr = view->getScenePtr()->get_polyhedron();
		
					for(Facet_iterator pFacet = polyhedron_ptr->facets_begin();pFacet!=polyhedron_ptr->facets_end();++pFacet)
					{
						m_facets.push_back(pFacet);
					}
					PaintStart(viewer);
				}
				else{
					std::cout << "Correspondence ...  " << std::endl;
					QApplication::setOverrideCursor(Qt::WaitCursor);

					mw->statusBar()->showMessage(tr("Correspondence..."));
					
					component_ptr->readSelectionBasedOnColor(polyhedron_ptr);
					
					component_ptr->initParameters(nbLabel,meshID);
					
					component_ptr->readDescriptor(polyhedron_ptr);
					
					component_ptr->initializeEllipsoid(polyhedron_ptr);
					viewer->recreateListsAndUpdateGL();
					
					
					
					component_ptr->computeEllipseParameters(polyhedron_ptr);
					
					component_ptr->compareDescriptorToEllipse(polyhedron_ptr);
					
					
					
					mw->statusBar()->showMessage(tr("Correspondence is done"));

					component_ptr->set_init(2);
					viewer->recreateListsAndUpdateGL();
					
					compareToDataset(component_ptr,meshID);
				}
			}
		}
	}
	QApplication::restoreOverrideCursor();
		
}

void mepp_component_Correspondence_plugin::PaintStart(Viewer * view)
{
	std::cout<< "inside PaintStart"; 
	
		PolyhedronPtr polyhedron_ptr = view->getScenePtr()->get_polyhedron();
		
		/*for(Facet_iterator pFacet = polyhedron_ptr->facets_begin();pFacet!=polyhedron_ptr->facets_end();++pFacet)
		{
			m_facets.push_back(pFacet);
		}*/
		
		view->makeCurrent();
		
		QGLFramebufferObjectFormat fboFormat;
		fboFormat.setAttachment(QGLFramebufferObject::CombinedDepthStencil);
		fboFormat.setInternalTextureFormat(GL_RGB);
		fboFormat.setMipmap(false);
		m_fbo = new QGLFramebufferObject(view->width(),view->height(),fboFormat);
		
		m_fbo->bind();
		glDisable(GL_LIGHTING);
		glShadeModel(GL_FLAT);
		glClearColor(0.0,0.0,0.0,1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glBegin(GL_TRIANGLES);
		int i = 0;
		for (Facet_iterator pFacet = polyhedron_ptr->facets_begin();
		     pFacet!= polyhedron_ptr->facets_end(); pFacet++)
		     {
			Halfedge_around_facet_circulator hE = pFacet->facet_begin();
			
			do{
				Vertex_handle pVertex = hE->opposite()->vertex();
				unsigned id = i;
				
				unsigned char color[3];
				color[0] = (id>>16) & 0xFF;
				color[1] = (id>>8) & 0xFF; 
				color[2] = id & 0xFF;
				glColor3ubv(color);
				Point3d p = pVertex->point();
				glVertex3f(p.x(),p.y(),p.z());
			}while(++hE!=pFacet->facet_begin());
			i++;
		}
		glEnd();
		glPopMatrix();
		m_fbo->release();
		glEnable(GL_LIGHTING);
		glShadeModel(GL_SMOOTH);
	
}

void mepp_component_Correspondence_plugin::compareToDataset(Correspondence_ComponentPtr sourceCorrespondence, int sourceID)
{
	Viewer* viewerI = NULL;
	
	PolyhedronPtr polyhedron_ptr_out;
	
	for(int m = 1; m <=20; ++m)
	{
		if(m == sourceID) {continue;}
		emit(mw->get_actionNewEmpty()->trigger());
		
		for(int i=0; i<lwindow.size();i++)
		{
			viewerI = (Viewer*)qobject_cast<QWidget *>(lwindow[i]->widget());
			if(viewerI->getScenePtr()->get_polyhedron()->empty())
			{
				std::stringstream ss;
				ss << "/home/leon/datasetHuman/" << m << ".off";
				viewerI->getScenePtr()->add_mesh(ss.str().c_str(),0,NULL,viewerI);
				
				PolyhedronPtr polyhedron_ptr = viewerI->getScenePtr()->get_polyhedron();	
				Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewerI, polyhedron_ptr);
				component_ptr->initParameters(8,m);
				component_ptr->readDescriptor(polyhedron_ptr);
				component_ptr->setEllipse(sourceCorrespondence->getEllipse());
				component_ptr->setCentreDescriptor(sourceCorrespondence->getCentreDescr());
				component_ptr->compareDescriptorToEllipse(polyhedron_ptr);
				viewerI->recreateListsAndUpdateGL();
			}
		}
	}
}




#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(mepp_component_Correspondence_plugin, mepp_component_Correspondence_plugin);
#endif

#endif
