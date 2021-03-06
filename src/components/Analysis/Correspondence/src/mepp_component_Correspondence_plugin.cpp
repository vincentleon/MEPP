#include <mepp_config.h>
#ifdef BUILD_component_Correspondence

#include "mepp_component_Correspondence_plugin.hxx"

#include "dialSettings.hxx"
#include <Polyhedron/polyhedron_enriched_polyhedron.h>

#include <QObject>
#include <QAction>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>
#include <QMdiSubWindow>
#include <QGLFramebufferObject>

#include <QFileDialog>

#include "../components/Tools/Boolean_Operations/src/BoolPolyhedra.h"

#include "../components/Tools/Boolean_Operations/src/Boolean_Operations_Component.h"

typedef boost::shared_ptr<Boolean_Operations_Component> Boolean_Operations_ComponentPtr;


std::stack<clock_t> timer_stack;

void timer_tic() {
  timer_stack.push(clock());
}

void timer_toc() {
    std::cout << "Time elapsed: "
              << ((double)(clock() - timer_stack.top())) / CLOCKS_PER_SEC
              << std::endl;
    timer_stack.pop();
}

void mepp_component_Correspondence_plugin::post_draw()
{
	if (mw->activeMdiChild() != 0) 
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
		
	}
}

void mepp_component_Correspondence_plugin::post_draw_all_scene()
{
	if(mw->activeMdiChild() !=0)
	{
		Viewer* mainViewer = (Viewer*)qobject_cast<QWidget *>(lwindow[0]->widget());
		PolyhedronPtr mainPoly = mainViewer->getScenePtr()->get_polyhedron(0);
		Correspondence_ComponentPtr mainCorres = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(mainViewer, mainPoly);
		
		
		SegmentController & mainSegCtr = mainCorres->m_segCtr;
		PolyhedronPtr partPoly = mainViewer->getScenePtr()->get_polyhedron(1);
			
			
		/*if(m_icpOk)
		{
			glPushMatrix();
			glDisable(GL_LIGHTING);
			glLineWidth(2);
			
		
			auto phi = mainSegCtr.getPhi();
			
			int id = 0;
			for(auto pVertex = partPoly->vertices_begin(); pVertex!=partPoly->vertices_end();++pVertex)
			{
				unsigned char color[3];
				color[0] = (id*732)%255;
				color[1] = (id*3625 + 20)%255;
				color[2] = (id*4589 + 95)%255;
				glColor3ubv(color);
				Vertex_handle corres = phi[pVertex];
				if(corres == Vertex_handle()){continue;}
				Vec pi(pVertex->point().x(), pVertex->point().y(), pVertex->point().z());

				Vec pj(corres->point().x(), corres->point().y(), corres->point().z());
				draw_link(mainViewer, 0, 1, pi, pj);
				id++;
			}
			glEnable(GL_LIGHTING);
			glPopMatrix();
			//mainViewer->recreateListsAndUpdateGL();
		}*/
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
		//std::cout << "Mouse Left Down" << std::endl;
		
		//std::cout << "correspondence Done : " << m_correspondenceDone << std::endl;
		
		if(m_correspondenceDone)
		{
			m_facets.clear();
			PaintStart(viewer);
			m_PaintingMode = false;
			//std::cout << "Mouse Left Down Correspondence Done" << std::endl;
			Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewer, polyhedron_ptr);
			int tagCC = 0;
			
			int x,y,h,w;
			x = event->x();
			y = event->y();
			h = viewer->height();
			w = viewer->width();
			unsigned char color[3];
			
			m_fbo->bind();
			glReadPixels(x,h-y,1,1,GL_RGB,GL_UNSIGNED_BYTE,color);
			unsigned index = (color[0] << 16) | (color[1] << 8) | color[2];
			
			if( index != 0xFFFFFF )
			{
				Facet_iterator pFacet = m_facets[index];
				Halfedge_around_facet_circulator hE = pFacet->facet_begin();
				Vertex_handle pVertex = hE->vertex();
				tagCC = component_ptr->m_segCtr.m_ccMap[pVertex];
			}
			m_fbo->release();
			delete m_fbo;
			viewer->recreateListsAndUpdateGL();
			viewer->getScenePtr()->set_loadType(1);
			component_ptr->m_segCtr.cutCC(tagCC);
			PolyhedronPtr newPart = component_ptr->m_segCtr.m_parts.back();
			viewer->getScenePtr()->set_loadType(1);
			viewer->getScenePtr()->add_polyhedron(newPart);
			viewer->getScenePtr()->todoIfModeSpace(viewer,0.0);
			component_ptr->m_segCtr.cutComplementCC(tagCC);
			for(unsigned ip=0;ip<component_ptr->m_segCtr.m_mainPart.size();++ip)
			{
				PolyhedronPtr mpart = component_ptr->m_segCtr.m_mainPart[ip];
				viewer->getScenePtr()->add_polyhedron(mpart);
				viewer->getScenePtr()->todoIfModeSpace(viewer,0.0);
			}
			
			viewer->recreateListsAndUpdateGL();
			
			Viewer* mainViewer = (Viewer*)qobject_cast<QWidget *>(lwindow[0]->widget());
			if(viewer != mainViewer && !m_selectedTarget)
			{
				m_selectedTarget = true;
			}
			else	
			{
				m_selectedModel = true;
			}
			
			if(m_selectedTarget && m_selectedModel)
			{
				//std::cout << "both are selected" << std::endl;
				PolyhedronPtr polyhedron_ptr = mainViewer->getScenePtr()->get_polyhedron();
				Correspondence_ComponentPtr mainComponent_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(mainViewer, polyhedron_ptr);
				SegmentController & mainSegCtr = mainComponent_ptr->m_segCtr;
	
				mainViewer->getScenePtr()->add_polyhedron(newPart);
				mainViewer->getScenePtr()->todoIfModeSpace(mainViewer,0.0);
				mainViewer->recreateListsAndUpdateGL();
				
				mainSegCtr.alignSegments(mainViewer,mainSegCtr.m_parts.back(),newPart,mainViewer->get_nb_frames()-2,mainViewer->get_nb_frames()-1);
				//std::cout << "alignSegments ok" << std::endl;
				
				SettingsDialog dial;
				double elasticity = 1/500.0;
				double regionSize = 0.5;
				int itermax = 5;
				/*if (dial.exec() == QDialog::Accepted)
				{*/
					/*elasticity = dial.doubleSpinBox->value();
					regionSize = dial.doubleSpinBox_2->value();
					itermax = dial.spinBox_2->value();*/
					elasticity = 0.001;
					regionSize = 1.1;
					itermax = 2;
					mainSegCtr.softICP(mainViewer,newPart,elasticity,regionSize,itermax);	
					std::list<PolyhedronPtr> & polys = mainSegCtr.getPolyhedronToMerge();
					
					std::cout << "after softICP & remesh" << std::endl;
					
					for(auto p = polys.begin(); p!= polys.end(); ++p)
					{
						mainViewer->getScenePtr()->add_polyhedron(*p);
						//mainViewer->getScenePtr()->todoIfModeSpace(mainViewer,0.0);
					}
					mainViewer->getScenePtr()->todoIfModeSpace(mainViewer,0.0);
					/*
					PolyhedronPtr resultSR = polys.front();
					polys.pop_front();
					
					Boolean_Operations_ComponentPtr boolOp = findOrCreateComponentForViewer<Boolean_Operations_ComponentPtr, Boolean_Operations_Component>(mainViewer, resultSR);
					
					PolyhedronPtr p1 = resultSR;
					PolyhedronPtr p2 = polys.front();
					polys.pop_front();
					PolyhedronPtr res(new Polyhedron);
					
					
					p1->compute_normals();
					p1->compute_bounding_box();
					p2->compute_normals();
					p2->compute_bounding_box();
					
					mainViewer->getScenePtr()->add_polyhedron(p1);
					mainViewer->getScenePtr()->todoIfModeSpace(mainViewer,0.0);
					mainViewer->getScenePtr()->add_polyhedron(p2);
					mainViewer->getScenePtr()->todoIfModeSpace(mainViewer,0.0);
					
					std::cout << "p1 size of facets : " << p1->size_of_facets() << std::endl;
					std::cout << "p1 size of facets : " << p1->size_of_vertices() << std::endl;
					
					std::cout << "p2 size of facets : " << p2->size_of_facets() << std::endl;
					std::cout << "p2 size of facets : " << p2->size_of_vertices() << std::endl;
					
					boolOp->Boolean_Union(p1, p2, res);
					
					res->compute_normals();
					res->compute_bounding_box();
					
					std::cout << "res size of facets : " << res->size_of_facets() << std::endl;
					std::cout << "res size of facets : " << res->size_of_vertices() << std::endl;
					
					std::cout << "first Boolean Union" << std::endl;
					boolOp->cpt_U++;
					while(!polys.empty())
					{
						std::cout << "poly not empty : " << polys.size() << std::endl;
						p1 = res;
						p2 = polys.front();
					
						PolyhedronPtr temp(new Polyhedron);
						std::cout << "size p1 : " << p1->size_of_facets() << std::endl;
						std::cout << "size p2 : " << p2->size_of_facets() << std::endl;
						polys.pop_front();
						boolOp->Boolean_Union(p1, p2, temp);
						boolOp->cpt_U++;
						std::cout << "polys : " << polys.size() << std::endl;
						res = temp;
					}
					mainViewer->getScenePtr()->add_polyhedron(res);
					mainViewer->getScenePtr()->todoIfModeSpace(mainViewer,0.0);*/
				//}
				
				mainViewer->recreateListsAndUpdateGL();
				m_correspondenceDone = false;
				
			}
			
		}
		
 	}
 	
}

void mepp_component_Correspondence_plugin::OnMouseLeftUp(QMouseEvent *event)
{
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
		//m_hasNotBeenPainted=false;
 	}
}

void mepp_component_Correspondence_plugin::OnMouseMotion(QMouseEvent* event)
{
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
		Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewer, polyhedron_ptr);
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
		viewer->recreateListsAndUpdateGL();
	}
}

void mepp_component_Correspondence_plugin::showDescriptor()
{
	Viewer* viewerI = NULL;
	
	int nbLabel = 4;
	for(int i=0; i<lwindow.size();i++)
	{
		viewerI = (Viewer*)qobject_cast<QWidget *>(lwindow[i]->widget());
		if(!viewerI->getScenePtr()->get_polyhedron()->empty())
		{
			for(unsigned j =0 ; j<viewerI->getScenePtr()->get_nb_polyhedrons();++j)
			{
				PolyhedronPtr polyhedron_ptr = viewerI->getScenePtr()->get_polyhedron(j);
				Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewerI, polyhedron_ptr);
			
				/*
				std::string meshIDString = polyhedron_ptr->pName;
				unsigned posB = meshIDString.find_last_of("/");
				unsigned posE = meshIDString.find_last_of(".ply");
					
				int meshID = atoi(meshIDString.substr(posB+1 ,posE).c_str());
				
				component_ptr->initParameters(nbLabel,meshID,meshIDString.substr(0,posB));
							
				component_ptr->readDescriptor(polyhedron_ptr,meshIDString.substr(0,posB));*/
				component_ptr->showDescriptor(polyhedron_ptr,m_currentLabel);
				viewerI->recreateListsAndUpdateGL();
			}
		}
	}
	m_currentLabel = (m_currentLabel + 1) % nbLabel;
}

void mepp_component_Correspondence_plugin::OnPainting()
{
	if(mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
		for(Facet_iterator pFacet = polyhedron_ptr->facets_begin();pFacet!=polyhedron_ptr->facets_end();++pFacet)
		{
			m_facets.push_back(pFacet);
		}
		PaintStart(viewer);
		
	}
	m_hasNotBeenPainted = true;
}

void mepp_component_Correspondence_plugin::OnCorrespondence()
{
	
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
		
		Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewer, polyhedron_ptr);
		{
			int nbLabel = 4;
			std::string meshIDString = polyhedron_ptr->pName;
			unsigned posB = meshIDString.find_last_of("/");
			unsigned posE = meshIDString.find_last_of(".ply");
			if(posE == std::string::npos)
			{
				posE = meshIDString.find_last_of(".off");
			}
			if(posE == std::string::npos)
			{
				posE = meshIDString.find_last_of(".obj");
			}
			
			
			int meshID = atoi(meshIDString.substr(posB+1 ,posE).c_str());

			//SettingsDialog dial;
			/*if (dial.exec() == QDialog::Accepted)
			{*/
				//nbLabel = dial.spinBox->value();
				
				QApplication::setOverrideCursor(Qt::WaitCursor);

				mw->statusBar()->showMessage(tr(" Ellipse Correspondence..."));
				
				component_ptr->readSelectionBasedOnColor(polyhedron_ptr);
				
				component_ptr->initParameters(nbLabel,meshID,meshIDString.substr(0,posB));
				
				saveToPly("/home/leon/Desktop/kring.ply",polyhedron_ptr);
				
				component_ptr->readDescriptor(polyhedron_ptr,meshIDString.substr(0,posB));
				
				component_ptr->initializeEllipsoid(polyhedron_ptr);
				component_ptr->compareDescriptorToEllipse(polyhedron_ptr);
				saveToPly("/home/leon/Desktop/init.ply",polyhedron_ptr);
				
				
				std::cout<< " learn Ellipse " << std::endl;
				timer_tic();
				component_ptr->computeEllipseParameters(polyhedron_ptr);
				timer_toc();
				
				const std::vector<double> ell = component_ptr->getEllipse();
				unsigned b=0; double v=0.0;
				for(unsigned l=0;l<ell.size();++l)
				{
					if(ell[l]>v){v=ell[l];b=l;}
				}
				std::vector<double> largerEll = ell;
				std::vector<double> smallerEll = ell;
				
				component_ptr->setEllipse(smallerEll); // set a smaller Ellipse for the query
				
				std::cout<< " test Ellipse " << std::endl;
				timer_tic();
				component_ptr->compareDescriptorToEllipse(polyhedron_ptr);
				timer_toc();
				
				saveToPly("/home/leon/Desktop/final.ply",polyhedron_ptr);
				
				SegmentController & segCtr = component_ptr->m_segCtr;
				segCtr.tagVerticesAfterCorrespondence();
				
				mw->statusBar()->showMessage(tr("Correspondence is done"));

				component_ptr->set_init(2);
				viewer->showAllScene();
				viewer->recreateListsAndUpdateGL();
			
				//component_ptr->setEllipse(largerEll); // set a larger Ellipse for the segments to add
 				compareToDataset(component_ptr,meshID);
			//}
		}
		m_correspondenceDone = true;
	}
	QApplication::restoreOverrideCursor();
}

void mepp_component_Correspondence_plugin::OnSVMCorrespondance()
{
	if(mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
		Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewer, polyhedron_ptr);
		int nbLabel = 4;
		
		std::string meshIDString = polyhedron_ptr->pName;
		unsigned posB = meshIDString.find_last_of("/");
		unsigned posE = meshIDString.find_last_of(".ply");
		if(posE == std::string::npos)
		{
			posE = meshIDString.find_last_of(".off");
		}
		if(posE == std::string::npos)
		{
			posE = meshIDString.find_last_of(".obj");
		}
		int meshID = atoi(meshIDString.substr(posB+1 ,posE).c_str());
		
		SettingsDialog dial; 
		
		mw->statusBar()->showMessage(tr("SVM Correspondance..."));
		
		component_ptr->readSelectionBasedOnColor(polyhedron_ptr);
		
		component_ptr->initParameters(nbLabel,meshID,meshIDString.substr(0,posB));
		
		component_ptr->readDescriptor(polyhedron_ptr,meshIDString.substr(0,posB));
		
		// Find the point that is the closest to the center of the patch. 
		// Then classify all the facets in its vicinity
		
		unsigned SVM_mode = 0;
		std::cout << "Select SVM Mode : 0 for semDescr, 1 for (d,h,a), 2 for all" << std::endl;
		std::cin >> SVM_mode;
		if(SVM_mode > 2)
		{
			SVM_mode = 0;
		}
		//for(unsigned mode = 0; mode < 3 ; ++mode)
		//{
		//	SVM_mode = mode;
			std::cout << "Learn SVM : " << std::endl;
			timer_tic();
			component_ptr->learnSVMPatch(polyhedron_ptr,SVM_mode);
			timer_toc();
			
			std::cout << "test SVM : " << std::endl;
			timer_tic();
			component_ptr->compareDescriptorWithSVM(polyhedron_ptr,SVM_mode);
			timer_toc();
			
		
			saveToPly("/home/leon/Desktop/svm_mode"+to_string(SVM_mode)+".ply",polyhedron_ptr);
		//}
		
		
		compareToDatasetSVM(component_ptr,meshID,SVM_mode);
		viewer->recreateListsAndUpdateGL();
	}
	
}

void mepp_component_Correspondence_plugin::OnSVM()
{
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
		
		Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewer, polyhedron_ptr);
		{
			int nbLabel = 4;
			std::string meshIDString = polyhedron_ptr->pName;
			unsigned posB = meshIDString.find_last_of("/");
			unsigned posE = meshIDString.find_last_of(".obj");
			
			int meshID = atoi(meshIDString.substr(posB+1 ,posE).c_str());
			
			SettingsDialog dial;
			if (dial.exec() == QDialog::Accepted)
			{
				
				QApplication::setOverrideCursor(Qt::WaitCursor);

				mw->statusBar()->showMessage(tr("SVM Correspondence..."));
				
				component_ptr->readSelectionBasedOnColor(polyhedron_ptr);
				
				component_ptr->initParameters(nbLabel,meshID,meshIDString.substr(0,posB));
				
				component_ptr->readDescriptor(polyhedron_ptr,meshIDString.substr(0,posB));
					
				component_ptr->learnSVMClassifier(polyhedron_ptr);
				
				component_ptr->compareDescriptorWithSVM(polyhedron_ptr);
				
				mw->statusBar()->showMessage(tr("SVM Correspondence is done"));
				
				component_ptr->set_init(2);
				viewer->recreateListsAndUpdateGL();
		
				compareToDatasetSVM(component_ptr,meshID);
			}
		}
	}
	QApplication::restoreOverrideCursor();
}

void mepp_component_Correspondence_plugin::OnMahalanobis()
{
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
		
		Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewer, polyhedron_ptr);
		{
			int nbLabel = 4;
			std::string meshIDString = polyhedron_ptr->pName;
			unsigned posB = meshIDString.find_last_of("/");
			unsigned posE = meshIDString.find_last_of(".ply");
			
			int meshID = atoi(meshIDString.substr(posB+1 ,posE).c_str());
			
			SettingsDialog dial;
			if (dial.exec() == QDialog::Accepted)
			{
				
				QApplication::setOverrideCursor(Qt::WaitCursor);

				mw->statusBar()->showMessage(tr("Mahalanobis..."));
				
				component_ptr->readSelectionBasedOnColor(polyhedron_ptr);
				
				component_ptr->initParameters(nbLabel,meshID,meshIDString.substr(0,posB));
				
				component_ptr->readDescriptor(polyhedron_ptr,meshIDString.substr(0,posB));
				
				component_ptr->computeGaussianParameters(polyhedron_ptr);
				
				component_ptr->compareDescriptorToGaussian(polyhedron_ptr);

				mw->statusBar()->showMessage(tr("Mahalanobis is done"));

				component_ptr->set_init(2);
				viewer->recreateListsAndUpdateGL();
				
				compareToDatasetMahalanobis(component_ptr,meshID);
				
			}
		}
	}
	QApplication::restoreOverrideCursor();
}

void mepp_component_Correspondence_plugin::PaintStart(Viewer * view)
{	
		PolyhedronPtr polyhedron_ptr = view->getScenePtr()->get_polyhedron();
		
		view->makeCurrent();
		
		QGLFramebufferObjectFormat fboFormat;
		fboFormat.setAttachment(QGLFramebufferObject::CombinedDepthStencil);
		fboFormat.setInternalTextureFormat(GL_RGB);
		fboFormat.setMipmap(false);
		m_fbo = new QGLFramebufferObject(view->width(),view->height(),fboFormat);
		
		m_fbo->bind();
		glDisable(GL_LIGHTING);
		glShadeModel(GL_FLAT);
		glClearColor(1.0,1.0,10.0,1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glBegin(GL_TRIANGLES);
		int i = 0;
		for (Facet_iterator pFacet = polyhedron_ptr->facets_begin();
		     pFacet!= polyhedron_ptr->facets_end(); pFacet++)
		     {
			Halfedge_around_facet_circulator hE = pFacet->facet_begin();
			m_facets.push_back(pFacet);
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
	std::vector<std::string> files;
	DIR *dp;
	struct dirent *dirp;
	
	QString path = QFileDialog::getExistingDirectory (mw, tr("Choose directory"),"/home/leon/");
	std::string dir = path.toStdString();
	
	// read all files in directory
	if((dp  = opendir(dir.c_str())) == NULL){ std::cout<< "Unable to open directory "<<dir<<std::endl;}
	while ((dirp = readdir(dp)) != NULL){ files.push_back(std::string(dirp->d_name));}
	closedir(dp);
	
	int nbLabel = 4;
	//SettingsDialog dial;
	/*if (dial.exec() == QDialog::Accepted)
	{
		nbLabel = dial.spinBox->value();
 	}*/
	
	PolyhedronPtr sourcePoly = sourceCorrespondence->get_polyhedron_ptr();
	Polyhedron::Iso_cuboid mbbox = sourcePoly->bbox();

	double elapsedTime = 0.0;
	
	for(unsigned i=0;i<files.size();++i)
	{
		unsigned len = files[i].size();
		bool isHidden = (files[i][0]=='.');
		
		if(!isHidden)
		{
			bool isPLY = (files[i].substr(len-3)=="ply");
			bool isOBJ = (files[i].substr(len-3)=="obj");
			if(isPLY || isOBJ)
			{
				emit(mw->get_actionNewEmpty()->trigger());
				
				for(int j=0; j<lwindow.size();j++)
				{
					viewerI = (Viewer*)qobject_cast<QWidget *>(lwindow[j]->widget());
					if(viewerI->getScenePtr()->get_polyhedron()->empty())
					{
						QSize wSize = lwindow[j]->size();
 						lwindow[j]->resize(wSize.width()/2.0,wSize.height()/2.0);
						std::cout << "File :"<< files[i] << std::endl;
						viewerI->getScenePtr()->add_mesh(path+"/"+files[i].c_str(),0,NULL,viewerI);
						if(viewerI->getScenePtr()->get_polyhedron()->empty())
						{continue;}
						
						PolyhedronPtr polyhedron_ptr = viewerI->getScenePtr()->get_polyhedron();
						
						std::string meshIDString = polyhedron_ptr->pName;
						unsigned posB = meshIDString.find_last_of("/");
						unsigned posE = 0;
						if(isPLY){ posE = meshIDString.find_last_of(".ply"); }
						if(isOBJ){ posE = meshIDString.find_last_of(".obj"); }
						int meshID = atoi(meshIDString.substr(posB+1 ,posE).c_str());
			
						Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewerI, polyhedron_ptr);
						
						component_ptr->scaleMesh(mbbox,polyhedron_ptr);
						
						component_ptr->initParameters(nbLabel,meshID,dir);
						component_ptr->readDescriptor(polyhedron_ptr,dir);
						component_ptr->setEllipse(sourceCorrespondence->getEllipse());
						component_ptr->setCentreDescriptor(sourceCorrespondence->getCentreDescr());
						component_ptr->compareDescriptorToEllipse(polyhedron_ptr);
						component_ptr->m_segCtr.tagVerticesAfterCorrespondence();
						viewerI->recreateListsAndUpdateGL();
					}
				}
			}
		}
	}
}

void mepp_component_Correspondence_plugin::compareToDatasetMahalanobis(Correspondence_ComponentPtr sourceCorrespondence, int sourceID)
{
	Viewer* viewerI = NULL;
	std::vector<std::string> files;
	DIR *dp;
	struct dirent *dirp;
	
	QString path = QFileDialog::getExistingDirectory (mw, tr("Choose directory"),"/home/leon/");
	std::string dir = path.toStdString();
	
	// read all files in directory
	if((dp  = opendir(dir.c_str())) == NULL){ std::cout<< "Unable to open directory "<<dir<<std::endl;}
	while ((dirp = readdir(dp)) != NULL){ files.push_back(std::string(dirp->d_name));}
	closedir(dp);
	
	int nbLabel = 4;
	SettingsDialog dial;
	if (dial.exec() == QDialog::Accepted)
	{
		nbLabel = dial.spinBox->value();
	}
	
	for(unsigned i=0;i<files.size();++i)
	{
		unsigned len = files[i].size();
		bool isHidden = (files[i][0]=='.');
		
		if(!isHidden)
		{
			bool isPLY = (files[i].substr(len-3)=="ply");
			if(isPLY)
			{
				emit(mw->get_actionNewEmpty()->trigger());
				
				for(int j=0; j<lwindow.size();j++)
				{
					viewerI = (Viewer*)qobject_cast<QWidget *>(lwindow[j]->widget());
					if(viewerI->getScenePtr()->get_polyhedron()->empty())
					{
						
						
						viewerI->getScenePtr()->add_mesh(path+"/"+files[i].c_str(),0,NULL,viewerI);
						
						PolyhedronPtr polyhedron_ptr = viewerI->getScenePtr()->get_polyhedron();
						std::string meshIDString = polyhedron_ptr->pName;
						unsigned posB = meshIDString.find_last_of("/");
						unsigned posE = meshIDString.find_last_of(".ply");
						int meshID = atoi(meshIDString.substr(posB+1 ,posE).c_str());
						Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewerI, polyhedron_ptr);
						
						component_ptr->initParameters(nbLabel,meshID,dir);	
						component_ptr->readDescriptor(polyhedron_ptr,dir);
						
						component_ptr->setMatrix(sourceCorrespondence->getMatrix());
						component_ptr->setInverseMatrix(sourceCorrespondence->getInverseMatrix());
						component_ptr->setDeterminant(sourceCorrespondence->getDeterminant());
						component_ptr->setVector(sourceCorrespondence->getVector());
						component_ptr->setThreshold(sourceCorrespondence->getThreshold());
					
						component_ptr->setCentreDescriptor(sourceCorrespondence->getCentreDescr());
						component_ptr->compareDescriptorToGaussian(polyhedron_ptr);
						viewerI->recreateListsAndUpdateGL();
					}
				}
			}
		}
	}
	
}

void mepp_component_Correspondence_plugin::compareToDatasetSVM(Correspondence_ComponentPtr sourceCorrespondence, int sourceID, unsigned SVM_mode)
{	
	Viewer* viewerI = NULL;
	std::vector<std::string> files;
	DIR *dp;
	struct dirent *dirp;
	
	QString path = QFileDialog::getExistingDirectory (mw, tr("Choose directory"),"/home/leon/");
	std::string dir = path.toStdString();
	
	// read all files in directory
	if((dp  = opendir(dir.c_str())) == NULL){ std::cout<< "Unable to open directory "<<dir<<std::endl;}
	while ((dirp = readdir(dp)) != NULL){ files.push_back(std::string(dirp->d_name));}
	closedir(dp);
	
	int nbLabel = 4;
	SettingsDialog dial;
	if (dial.exec() == QDialog::Accepted)
	{
		nbLabel = dial.spinBox->value();
	}
	
	for(unsigned i=0;i<files.size();++i)
	{
		unsigned len = files[i].size();
		bool isHidden = (files[i][0]=='.');
		
		if(!isHidden)
		{
			bool isOBJ = (files[i].substr(len-3)=="obj");
			bool isPLY = (files[i].substr(len-3)=="ply");
			if(true)
			{
				emit(mw->get_actionNewEmpty()->trigger());
				
				for(int j=0; j<lwindow.size();j++)
				{
					viewerI = (Viewer*)qobject_cast<QWidget *>(lwindow[j]->widget());
					if(viewerI->getScenePtr()->get_polyhedron()->empty())
					{
						viewerI->getScenePtr()->add_mesh(path+"/"+files[i].c_str(),0,NULL,viewerI);
						
						PolyhedronPtr polyhedron_ptr = viewerI->getScenePtr()->get_polyhedron();
						std::string meshIDString = polyhedron_ptr->pName;
						unsigned posB = meshIDString.find_last_of("/");
						unsigned posE = meshIDString.find_last_of(".ply");
						if(isPLY){ posE = meshIDString.find_last_of(".ply"); }
						if(isOBJ){ posE = meshIDString.find_last_of(".obj"); }
						int meshID = atoi(meshIDString.substr(posB+1 ,posE).c_str());
						
						if(meshID == sourceID)
						{continue;}
						
						Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewerI, polyhedron_ptr);
						
						component_ptr->initParameters(nbLabel,meshID,dir);	
						bool descrRead = component_ptr->readDescriptor(polyhedron_ptr,meshIDString.substr(0,posB));
						
						if(!descrRead)
						{continue;}
						component_ptr->setSVM(sourceCorrespondence->getSVM());
						component_ptr->setNbCandidates(sourceCorrespondence->getNbCandidates());
						component_ptr->setRadius(sourceCorrespondence->getRadius());
						component_ptr->setCentreDescriptor(sourceCorrespondence->getCentreDescr());
						component_ptr->featureSTD = sourceCorrespondence->featureSTD;
						component_ptr->featureMeans = sourceCorrespondence->featureMeans;
						//component_ptr->compareDescriptorToSVM(polyhedron_ptr);
						
						component_ptr->compareDescriptorWithSVM(polyhedron_ptr,SVM_mode);
				
						viewerI->recreateListsAndUpdateGL();
					}
				}
			}
		}
	}
	
}

void mepp_component_Correspondence_plugin::OnLearn()
{	
	if (mw->activeMdiChild() != 0)
	{	
		
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
		
		int nbLabel = 4;
		SettingsDialog dial;
		if (dial.exec() == QDialog::Accepted)
		{
			nbLabel = dial.spinBox->value();
		}
		
		Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewer, polyhedron_ptr);
		{
			std::string meshIDString = polyhedron_ptr->pName;
			unsigned posB = meshIDString.find_last_of("/");
			unsigned posE = meshIDString.find_last_of(".");
			
			int meshID = atoi(meshIDString.substr(posB+1 ,posE).c_str());
			
			std::cout << "Mesh ID : " << meshID << std::endl;
			component_ptr->initParameters(nbLabel,meshID,meshIDString.substr(0,posB));
			std::cout << "Init " << std::endl;
			component_ptr->learnDescriptor(polyhedron_ptr,meshIDString.substr(0,posB));
			std::cout << "learn descriptor " << std::endl;
			viewer->recreateListsAndUpdateGL();
		}
	}
}

void mepp_component_Correspondence_plugin::OnPrepareData()
{
	Viewer* viewerI = NULL;
	std::vector<std::string> files;
	DIR *dp;
	struct dirent *dirp;
	
	QString path = QFileDialog::getExistingDirectory (mw, tr("Choose directory"),"/home/leon/");	
	std::string dir = path.toStdString();
	
	// read all files in directory
	if((dp  = opendir(dir.c_str())) == NULL){ std::cout<< "Unable to open directory "<<dir<<std::endl;}
	while ((dirp = readdir(dp)) != NULL){ files.push_back(std::string(dirp->d_name));}
	closedir(dp);
	
	int nbLabel = 4;
	SettingsDialog dial;
	if (dial.exec() == QDialog::Accepted)
	{
		nbLabel = dial.spinBox->value();
	}
	
	for(unsigned i=0;i<files.size();++i)
	{
		unsigned len = files[i].size();
		bool isHidden = (files[i][0]=='.');
		
		if(!isHidden)
		{
			bool isPLY = (files[i].substr(len-3)=="obj");
			if(isPLY)
			{
				emit(mw->get_actionNewEmpty()->trigger());
				for(int j=0; j<lwindow.size();j++)
				{
				
					viewerI = (Viewer*)qobject_cast<QWidget *>(lwindow[j]->widget());
					if(viewerI->getScenePtr()->get_polyhedron()->empty())
					{
						viewerI->getScenePtr()->add_mesh(path+"/"+files[i].c_str(),0,NULL,viewerI);
						if(viewerI->getScenePtr()->get_polyhedron()->empty())
						{
							std::cout << files[i] << " not learnt" << std::endl;
							continue;
						}
						
						PolyhedronPtr polyhedron_ptr = viewerI->getScenePtr()->get_polyhedron();
						
						if(polyhedron_ptr->is_pure_triangle())	
						{
							Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewerI, polyhedron_ptr);
							
							std::string meshIDString = files[i];
						        int meshID = atoi(meshIDString.substr(0,len-4).c_str());
							component_ptr->initParameters(nbLabel,meshID,dir);
							component_ptr->learnDescriptor(polyhedron_ptr,dir);
							
							viewerI->setWindowTitle(files[i].c_str());
							viewerI->recreateListsAndUpdateGL();
						}
						else
						{
							lwindow[i]->close();
						}
					}
					
				}
			}
		}
	}	
}

void mepp_component_Correspondence_plugin::OnSelectionBorder()
{
	Viewer* mainViewer = (Viewer*)qobject_cast<QWidget *>(lwindow[0]->widget());
	PolyhedronPtr mainPoly = mainViewer->getScenePtr()->get_polyhedron();
	Correspondence_ComponentPtr mainCorres = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(mainViewer, mainPoly);
	SegmentController & mainSegCtr = mainCorres->m_segCtr;
	
	mainSegCtr.getBorder();
	mainViewer->recreateListsAndUpdateGL();
}

void mepp_component_Correspondence_plugin::OnMoveBorder()
{
	Viewer* mainViewer = (Viewer*)qobject_cast<QWidget *>(lwindow[0]->widget());
	PolyhedronPtr mainPoly = mainViewer->getScenePtr()->get_polyhedron();
	Correspondence_ComponentPtr mainCorres = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(mainViewer, mainPoly);
	SegmentController & mainSegCtr = mainCorres->m_segCtr;
	
	mainSegCtr.optimizeBorders();
	mainViewer->recreateListsAndUpdateGL();
}

void mepp_component_Correspondence_plugin::OnCut()
{	
	//////////////////////////////////////////////////////////////////////////////////////////
	////
	//// Main Viewer : cut the mesh, remove the main polyhedron, add the parts,
	////			setVisible, color the selected ones (red)
	////
	//////////////////////////////////////////////////////////////////////////////////////////
	Viewer* mainViewer = (Viewer*)qobject_cast<QWidget *>(lwindow[0]->widget());
	PolyhedronPtr mainPoly = mainViewer->getScenePtr()->get_polyhedron();
	Correspondence_ComponentPtr mainCorres = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(mainViewer, mainPoly);
	SegmentController & mainSegCtr = mainCorres->m_segCtr;
	
	mainSegCtr.cutSegments(); // Cut the mesh
	mainViewer->getScenePtr()->set_loadType(1); // activate Space mode
	mainViewer->getScenePtr()->todoIfModeSpace(mainViewer,0.0);
	mainViewer->getScenePtr()->delete_polyhedron(0); // Remove the main polyhedron
	for(unsigned p=0; p<mainSegCtr.m_parts.size();++p)
	{
		mainViewer->getScenePtr()->add_polyhedron(mainSegCtr.m_parts[p]);
		mainViewer->getScenePtr()->setVisible(p,true);
		//colorMesh(mainSegCtr.m_parts[p],1,0,0);
	}
	for(unsigned p=0; p<mainSegCtr.m_mainPart.size();++p)
	{
		mainViewer->getScenePtr()->add_polyhedron(mainSegCtr.m_mainPart[p]);
		mainViewer->getScenePtr()->setVisible(p+mainSegCtr.m_parts.size(),true);
	}
	mainViewer->getScenePtr()->todoIfModeSpace(mainViewer,0.0);
	mainViewer->showAllSceneForSpaceMode();
	mainViewer->recreateListsAndUpdateGL();
	
	//////////////////////////////////////////////////////////////////////////////////////////
	////
	//// Secondary Viewers : cut the mesh, remove the main polyhedron, add the parts ,
	////			setVisible, color the selected ones (green)
	////
	//////////////////////////////////////////////////////////////////////////////////////////
	if(mw->activeMdiChild() !=0)
	{
		Viewer* viewerI = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyI = viewerI->getScenePtr()->get_polyhedron();
		Correspondence_ComponentPtr corresI = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewerI, polyI);
		int cMeshI = viewerI->getScenePtr()->get_current_polyhedron();
		SegmentController & segCtrI = corresI->m_segCtr;
		segCtrI.cutSegments(); // Cut the mesh
		viewerI->getScenePtr()->set_loadType(1); // activate Space mode
		viewerI->getScenePtr()->todoIfModeSpace(viewerI,0.0);
		viewerI->getScenePtr()->delete_polyhedron(cMeshI); // Remove the main polyhedron
		for(unsigned p=0; p<segCtrI.m_parts.size();++p)
		{		
			viewerI->getScenePtr()->add_polyhedron(segCtrI.m_parts[p]);
			viewerI->getScenePtr()->setVisible(p,true);
			colorMesh(segCtrI.m_parts[p],0,1,0);
		}
		for(unsigned p=0; p<segCtrI.m_mainPart.size();++p)
		{
			viewerI->getScenePtr()->add_polyhedron(segCtrI.m_mainPart[p]);
			viewerI->getScenePtr()->setVisible(p+segCtrI.m_parts.size(),true);
		}
		viewerI->getScenePtr()->todoIfModeSpace(viewerI,0.0);
		viewerI->showAllSceneForSpaceMode();
		viewerI->recreateListsAndUpdateGL();
	}
}

void mepp_component_Correspondence_plugin::OnSaveParts()
{
	Viewer* mainViewer = (Viewer*)qobject_cast<QWidget *>(lwindow[0]->widget());
	ScenePtr scene = mainViewer->getScenePtr();
	
	// Select directory
	QString path = QFileDialog::getExistingDirectory(mw, tr("Choose directory to save parts"),"/home/leon/");	
	std::string directory = path.toStdString();
	
	// Save all polyhedron data and descriptor to directory
	for(unsigned i=0; i<scene->get_nb_polyhedrons();++i)
	{
		PolyhedronPtr p = scene->get_polyhedron(i);
		//p->keep_largest_connected_components(1);
		//Correspondence_ComponentPtr corres = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(mainViewer,p);
		std::stringstream filename; 
		filename << directory << "/" << i << ".ply";
		//corres->initParameters(4,i,directory);
		//corres->getShape().m_meshID = i;
		//corres->saveDescriptor(p,directory);
		//corres->m_segCtr.unionSegments(mainViewer);
		//p->write_off(filename.str(),true,true);
		std::ofstream file(filename.str().c_str());
		file << "ply\n";
		file << "format ascii 1.0\n";
		file << "element vertex " << p->size_of_vertices() << "\n";
		file << "property float x\n";
		file << "property float y\n";
		file << "property float z\n";
		file << "property uchar red\n";
		file << "property uchar green\n";
		file << "property uchar blue\n";
		file << "element face "<< p->size_of_facets() << "\n";
		file << "property list uchar int vertex_indices\n";
		file << "end_header\n";
		p->set_index_vertices();
		for(auto pVertex = p->vertices_begin(); pVertex != p->vertices_end(); ++pVertex)
		{
			Point3d point = pVertex->point();
			file << point.x() << " " << point.y() << " " << point.z() <<" ";
			file << (int)pVertex->color(0)*255 << " " << (int)pVertex->color(1)*255 << " " << (int)pVertex->color(2)*255 << "\n";
		}
		for(auto pFacet = p->facets_begin(); pFacet != p->facets_end(); ++pFacet)
		{
			Halfedge_around_facet_circulator hE = pFacet->facet_begin();
			file << "3 ";
			do{
				file << hE->vertex()->tag() << " ";
				
			}while(++hE!=pFacet->facet_begin());
			file << "\n";
		}
		file.close();
	}
	
}

void mepp_component_Correspondence_plugin::OnLoadDescriptor()
{
	Viewer* mainViewer = (Viewer*)qobject_cast<QWidget *>(lwindow[0]->widget());
	ScenePtr scene = mainViewer->getScenePtr();
	
	// Select directory
	QString path = QFileDialog::getExistingDirectory (mw, tr("Choose directory to save parts"),"/home/leon/");	
	std::string directory = path.toStdString();
	
	
	// Save all polyhedron data and descriptor to directory
	for(unsigned i=0; i<scene->get_nb_polyhedrons();++i)
	{
		PolyhedronPtr p = scene->get_polyhedron(i);
		
		Correspondence_ComponentPtr corres = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(mainViewer,p);
		
		unsigned posB = p->pName.find_last_of("/");
		unsigned posE = p->pName.find_last_of(".");
		
		unsigned meshID = atoi(p->pName.substr(posB+1 ,posE).c_str());
		corres->initParameters(4,meshID,directory);
		corres->getShape().m_meshID = meshID;
		corres->readDescriptor(p,directory,false);
	}
}

void mepp_component_Correspondence_plugin::OnGlue()
{
	Viewer* mainViewer = (Viewer*)qobject_cast<QWidget *>(lwindow[0]->widget());
	PolyhedronPtr mainPoly = mainViewer->getScenePtr()->get_polyhedron();
	Correspondence_ComponentPtr mainCorres = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(mainViewer, mainPoly);
	int cMesh = mainViewer->getScenePtr()->get_current_polyhedron();
	SegmentController & mainSegCtr = mainCorres->m_segCtr;
	
	mainSegCtr.glueSegments(mainViewer);
	
	mainViewer->recreateListsAndUpdateGL();
	
	/*if(mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
		Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewer, polyhedron_ptr);
		
		component_ptr->m_segCtr.glueSegments();
		for(unsigned p=0;p<component_ptr->m_segCtr.m_parts.size();++p)
		{
			viewer->addFrame();
			viewer->getScenePtr()->add_polyhedron(component_ptr->m_segCtr.m_parts[p]);
		}
		for(unsigned p=0;p<component_ptr->m_segCtr.m_mainPart.size();++p)
		{
			viewer->addFrame();
			viewer->getScenePtr()->add_polyhedron(component_ptr->m_segCtr.m_mainPart[p]);
		}
		viewer->recreateListsAndUpdateGL();
	}*/
}

void mepp_component_Correspondence_plugin::OnUnion()
{
	//////////////////////////////////////////////////////////////
	///
	///	Get the Main viewer and fusion all the parts in it
	///
	//////////////////////////////////////////////////////////////
	SettingsDialog dial;
	double elasticity = 1/500.0;
	double regionSize = 0.5;
	int itermax = 5;
	if (dial.exec() == QDialog::Accepted)
	{
		elasticity = dial.doubleSpinBox->value();
		regionSize = dial.doubleSpinBox_2->value();
		itermax = dial.spinBox_2->value();
		std::cout << "itermax :" << itermax << std::endl;
		
		//std::cout << " Elasticity :" << elasticity << std::endl;
		
		Viewer* mainViewer = (Viewer*)qobject_cast<QWidget *>(lwindow[0]->widget());
		PolyhedronPtr mainPoly = mainViewer->getScenePtr()->get_polyhedron(0);
		Correspondence_ComponentPtr mainCorres = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(mainViewer, mainPoly);

		SegmentController & mainSegCtr = mainCorres->m_segCtr;
		PolyhedronPtr partPoly = mainViewer->getScenePtr()->get_polyhedron(1);
		
		std::cout << "itermax :" << itermax << std::endl;
		
		mainSegCtr.softICP(mainViewer,partPoly,mainPoly,elasticity,regionSize,itermax);		
		
		mainViewer->recreateListsAndUpdateGL();
	}
}

void mepp_component_Correspondence_plugin::OnAddSegment()
{
	//////////////////////////////////////////////////////////////////////////////////////////
	////
	//// Main Viewer : the first one to be open, just init everything
	////
	//////////////////////////////////////////////////////////////////////////////////////////
	Viewer* mainViewer = (Viewer*)qobject_cast<QWidget *>(lwindow[0]->widget());
	int currentMainViewer = mainViewer->getScenePtr()->get_current_polyhedron();
	PolyhedronPtr mainPolyhedron = mainViewer->getScenePtr()->get_polyhedron(); // This is the polyhedron to be replaced
	Correspondence_ComponentPtr mainCorres = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(mainViewer, mainPolyhedron);
	SegmentController & mainSegCtr = mainCorres->m_segCtr;
	
	//std::cout << "After getting mainViewer " << std::endl;
	
	//////////////////////////////////////////////////////////////////////////////////////////
	//// 
	//// Add the selected part in the active viewer to _mainViewer_
	////
	//////////////////////////////////////////////////////////////////////////////////////////
	if(mw->activeMdiChild() !=0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		int currentSecondaryViewer = viewer->getScenePtr()->get_current_polyhedron();
		PolyhedronPtr secondaryPolyhedron = viewer->getScenePtr()->get_polyhedron(); // This is the polyhedron that will replace _mainPolyhedron_
		
		
		
		Correspondence_ComponentPtr corres = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewer,secondaryPolyhedron);
		SegmentController & segCtr = corres->m_segCtr;
		
		mainViewer->getScenePtr()->add_polyhedron(secondaryPolyhedron);
		mainSegCtr.m_parts.push_back(secondaryPolyhedron);
		mainViewer->getScenePtr()->todoIfModeSpace(mainViewer,0.0);
		
		std::cout << "After adding polyhedron " << std::endl;
		mainViewer->recreateListsAndUpdateGL();
		
		mainSegCtr.unionSegments(mainViewer);
		mainSegCtr.alignSegments(mainViewer,mainPolyhedron,secondaryPolyhedron,currentMainViewer,mainViewer->get_nb_frames()-1);
		std::cout << "After alignSegments " << std::endl;
		
		/*mainViewer->getScenePtr()->add_polyhedron(secondaryPolyhedron);
		mainViewer->getScenePtr()->todoIfModeSpace(mainViewer,0.0);*/
		mainViewer->recreateListsAndUpdateGL();	
	}
}

void mepp_component_Correspondence_plugin::OnCompareMethods()
{
	Viewer* viewerI = NULL;
	
	int meshID = 0;
	
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
		
		emit(mw->get_actionNewEmpty()->trigger());
		
		Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewer, polyhedron_ptr);
		{
			int nbLabel = 4;
			std::string meshIDString = polyhedron_ptr->pName;
			unsigned posB = meshIDString.find_last_of("/");
			unsigned posE = meshIDString.find_last_of(".ply");
			
			meshID = atoi(meshIDString.substr(posB+1 ,posE).c_str());

			SettingsDialog dial;
			if (dial.exec() == QDialog::Accepted)
			{
				QApplication::setOverrideCursor(Qt::WaitCursor);
				nbLabel = dial.spinBox->value();
				mw->statusBar()->showMessage(tr(" Ellipse Correspondence..."));
				
				
				component_ptr->readSelectionBasedOnColor(polyhedron_ptr);
				
				component_ptr->initParameters(nbLabel,meshID,meshIDString.substr(0,posB));
				
				component_ptr->readDescriptor(polyhedron_ptr,meshIDString.substr(0,posB));
				std::cout << "Learn SVM : " << std::endl;
				//timer_tic();
				component_ptr->learnSVMClassifier(polyhedron_ptr);
				//timer_toc();
				
				component_ptr->initializeEllipsoid(polyhedron_ptr);
				std::cout << "Ellipse parameters : " << std::endl;
				//timer_tic();
				component_ptr->computeEllipseParameters(polyhedron_ptr);
				//timer_toc();
				std::cout << "compareToEllipse : " << std::endl;
				//timer_tic();
				component_ptr->compareDescriptorToEllipse(polyhedron_ptr);
				//timer_toc();
				mw->statusBar()->showMessage(tr("Correspondence is done"));

				component_ptr->set_init(2);
				viewer->recreateListsAndUpdateGL();
				
				for(int i=0; i<lwindow.size();i++)
				{
					viewerI = (Viewer*)qobject_cast<QWidget *>(lwindow[i]->widget());
					if(viewerI->getScenePtr()->get_polyhedron()->empty())
					{
						
						std::stringstream ss;
						ss << meshIDString.substr(0,posB) << "/"  << meshID << ".ply";
						viewerI->getScenePtr()->add_mesh(ss.str().c_str(),0,NULL,viewerI);
						viewerI->setWindowTitle("SVM");
						PolyhedronPtr polyhedron_copy_ptr = viewerI->getScenePtr()->get_polyhedron();	
						Correspondence_ComponentPtr svmcomponent_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewerI, polyhedron_copy_ptr);
						
						svmcomponent_ptr->readSelectionBasedOnColor(polyhedron_copy_ptr);
						
						svmcomponent_ptr->initParameters(nbLabel,meshID,meshIDString.substr(0,posB));	
						
						svmcomponent_ptr->readDescriptor(polyhedron_copy_ptr,meshIDString.substr(0,posB));
						std::cout << "compareToSVM : " << std::endl;
						//timer_tic();
						svmcomponent_ptr->setSVM(component_ptr->getSVM());
						svmcomponent_ptr->compareDescriptorToSVM(polyhedron_copy_ptr);
						//timer_toc();
						viewerI->recreateListsAndUpdateGL();
					}	
				}
			}
		}
	}
}

void mepp_component_Correspondence_plugin::OnCleanData()
{	
	Viewer* viewerI;
	
	SettingsDialog dial;
	
	int nbLabel = 4;
	
	if (dial.exec() == QDialog::Accepted)
	{
		QApplication::setOverrideCursor(Qt::WaitCursor);
		nbLabel = dial.spinBox->value();
			
		std::vector<std::string> files = getFileList("choose the segmented meshes directory");
		
		QString path = QFileDialog::getExistingDirectory (mw, tr("Choose the non-segmented meshes directory"),"/home/leon/");
		std::string meshName = path.toStdString();
		
		for(unsigned j=0;j<files.size();++j)
		{
			unsigned posE = files[j].find_last_of(".ply"); std::cout << posE << std::endl;
			if(posE == std::string::npos)
			{
				std::cout << "file " << files[j] << " is not .ply, continue;" << std::endl;
				continue;
			}
			std::cout << "file " << files[j] << std::endl;
				
			emit(mw->get_actionNewEmpty()->trigger()); // Create new window for each mesh
			for(int i=0; i<lwindow.size();i++)
			{
				
				viewerI = (Viewer*)qobject_cast<QWidget *>(lwindow[i]->widget());
				ScenePtr scene = viewerI->getScenePtr();
				if(viewerI->getScenePtr()->get_polyhedron()->empty())
				{
					viewerI->getScenePtr()->add_mesh(files[j].c_str(),1,NULL,viewerI);
					PolyhedronPtr segMesh = viewerI->getScenePtr()->get_polyhedron();
					
					if(viewerI->getScenePtr()->get_polyhedron()->empty())
					{
						//std::cout << "couldn't open mesh " << files[j] << std::endl;
						break;
					}
					
					// Load .part and .label files
					Correspondence_ComponentPtr comp_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewerI, segMesh);
					std::string meshIDString = segMesh->pName;
					unsigned posB = meshIDString.find_last_of("/");
					unsigned posE = meshIDString.find_last_of(".obj");
					int meshID = atoi(meshIDString.substr(posB+1 ,posE).c_str());
					std::string meshDir = meshIDString.substr(0,posB);
					comp_ptr->initParameters(nbLabel,meshID,meshDir);
					
					// Add the non-segmented mesh
					std::stringstream nSF;
					nSF<<meshName<<"/"<<meshID<<".obj";
					scene->add_mesh(nSF.str().c_str(),1,NULL,viewerI);
					
					//std::cout << "\t"<<files[j]<<" & "<< nSF.str() << std::endl;
					
					PolyhedronPtr goodMesh = scene->get_polyhedron();
					
					goodMesh->set_index_vertices();
					segMesh->set_index_vertices();
					
					std::cout << segMesh->pName << " & " << goodMesh->pName << std::endl;
					
					
					// Compute facet correspondence
					std::vector<int> faceLabels;
					std::vector<int> faceCC;
					std::vector<Facet_handle> closest;
					int pc = 0;
					Analysis::Shape & shape = comp_ptr->getShape();
					
					std::list<Triangle> triangles;
					int u = 0;
					for(auto pFacet = segMesh->facets_begin(); pFacet != segMesh->facets_end(); pFacet++)
					{
						triangles.push_back(Triangle(pFacet));
						pFacet->tag() = u;
						u++;
					}
					AABB_Tree tree(triangles.begin(),triangles.end()); 
					std::list<AABB_Tree::Primitive_id> primitives;
					
					int v = 0;
					for(auto pFacets = goodMesh->facets_begin();
					pFacets!=goodMesh->facets_end();++pFacets)
					{
						pFacets->tag() = v;
						v++;
						//tree.all_intersected_primitives(Triangle(pFacets),std::back_inserter(primitives));
						Facet_handle inter = tree.closest_point_and_primitive(pFacets->facet_begin()->vertex()->point()).second->facet();
						//Facet_handle inter = primitives.back()->facet();
						primitives.clear();
						closest.push_back(inter);
						pc++;
					}
						
						//int qc = 0;
						/*Point3d p = pFacets->facet_begin()->vertex()->point();
						double distMin = std::numeric_limits<double>::max();
						for(auto qFacets = segMesh->facets_begin();
						qFacets!=segMesh->facets_end();++qFacets)
						{
							Point3d q = qFacets->facet_begin()->vertex()->point();
							double dist = CGAL::squared_distance(p,q);
							if(dist< distMin)
							{
								distMin = dist;
								closest[pc] = qc;
							}
							qc++;
						}
						pc++;
					}*/
					for(unsigned c=0;c<closest.size();++c)
					{
						faceLabels.push_back(shape.m_faceLabels[closest[c]->tag()]);
						faceCC.push_back(shape.m_faceSegments[closest[c]->tag()]);
					}
					
					std::cout << "writing files ... in " << meshName << std::endl;
					
					// Now print the resulting labels and segmentID in files
					std::ofstream fileL;
					std::stringstream ssL;
					ssL<<meshName<<"/"<<shape.m_meshID<<".labelsN";
					fileL.open(ssL.str().c_str());
					for(unsigned c=0;c<faceLabels.size();++c)
					{
						
						fileL<<faceLabels[c]<<"\n";
					}
					fileL.close();
					
					std::ofstream fileCC;
					std::stringstream ssCC;
					ssCC<<meshName<<"/"<<shape.m_meshID<<".partsN";
					fileCC.open(ssCC.str().c_str());
					for(unsigned c=0;c<faceCC.size();++c)
					{
						
						fileCC<<faceCC[c]<<"\n";
					}
					fileCC.close();	
					break;
				}
				
			}
		}
	}
}

void mepp_component_Correspondence_plugin::OnSoftICP()
{
	//////////////////////////////////////////////////////////////
	///
	///	Get the Main viewer and fusion all the parts in it
	///
	//////////////////////////////////////////////////////////////
	SettingsDialog dial;
	double elasticity = 1/500.0;
	double regionSize = 0.5;
	int itermax = 5;

	if (dial.exec() == QDialog::Accepted)
	{
		elasticity = dial.doubleSpinBox->value();
		regionSize = dial.doubleSpinBox_2->value();
		itermax = dial.spinBox_2->value();
	
		Viewer* mainViewer = (Viewer*)qobject_cast<QWidget *>(lwindow[0]->widget());
		PolyhedronPtr mainPoly = mainViewer->getScenePtr()->get_polyhedron(0);
		Correspondence_ComponentPtr mainCorres = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(mainViewer, mainPoly);

		SegmentController & mainSegCtr = mainCorres->m_segCtr;
		PolyhedronPtr partPoly = mainViewer->getScenePtr()->get_polyhedron(1);
		
		mainSegCtr.softICP(mainViewer,partPoly,mainPoly,elasticity,regionSize,itermax);
		mainViewer->recreateListsAndUpdateGL();
		
		m_icpOk = true;
		
		//std::cout << "ICP OK ! " << std::endl;
		/*glLineWidth(2);
		glColor3f(1., 0., 0.);
	
		auto phi = mainSegCtr.getPhi();
		
		for(auto pVertex = mainPoly->vertices_begin(); pVertex!=mainPoly->vertices_end();++pVertex)
		{
			Vertex_handle corres = phi[pVertex];
			if(corres == Vertex_handle()){continue;}
			Vec pi(pVertex->point().x(), pVertex->point().y(), pVertex->point().z());

			Vec pj(corres->point().x(), corres->point().y(), corres->point().z());
			draw_link(mainViewer, 0, 1, pi, pj);
		}
		mainViewer->recreateListsAndUpdateGL();*/
	}
}

void mepp_component_Correspondence_plugin::OnRemesh()
{		
	Viewer* mainViewer = (Viewer*)qobject_cast<QWidget *>(lwindow[0]->widget());
	PolyhedronPtr mainPoly = mainViewer->getScenePtr()->get_polyhedron(0);
	Correspondence_ComponentPtr mainCorres = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(mainViewer, mainPoly);

	SegmentController & mainSegCtr = mainCorres->m_segCtr;
	PolyhedronPtr partPoly = mainViewer->getScenePtr()->get_polyhedron(1);
	
	mainSegCtr.remesh(mainViewer,partPoly,mainPoly);
	
	mainViewer->recreateListsAndUpdateGL();
	
}

vector< string > mepp_component_Correspondence_plugin::getFileList(std::string message)
{
	std::vector<std::string> files;
	DIR *dp;
	struct dirent *dirp;
	
	QString path = QFileDialog::getExistingDirectory (mw, tr(message.c_str()),"/home/leon/");
	std::string dir = path.toStdString();
	
	// read all files in directory
	if((dp  = opendir(dir.c_str())) == NULL){ std::cout<< "Unable to open directory "<<dir<<std::endl;}
	while ((dirp = readdir(dp)) != NULL)
	{
		if(dirp->d_name[0]!='.')
		{
			files.push_back(dir+"/"+std::string(dirp->d_name));
		}
		
	}
	closedir(dp);
	return files;
}

void mepp_component_Correspondence_plugindrawConnections(Viewer* viewer, int frame_i, int frame_j)
{
	PolyhedronPtr pMesh_i = viewer->getScenePtr()->get_polyhedron(frame_i);
	PolyhedronPtr pMesh_j = viewer->getScenePtr()->get_polyhedron(frame_j);
	
	
	
	for(auto pVertex = pMesh_i->vertices_begin(); pVertex!= pMesh_i->vertices_end();++pVertex)
	{
		
		
	}
		
}

void saveToPly(std::string filename, PolyhedronPtr p)
{
	std::ofstream file(filename.c_str());
	file << "ply\n";
	file << "format ascii 1.0\n";
	file << "element vertex " << p->size_of_vertices() << "\n";
	file << "property float x\n";
	file << "property float y\n";
	file << "property float z\n";
	file << "property uchar red\n";
	file << "property uchar green\n";
	file << "property uchar blue\n";
	file << "element face "<< p->size_of_facets() << "\n";
	file << "property list uchar int vertex_indices\n";
	file << "end_header\n";
	p->set_index_vertices();
	for(auto pVertex = p->vertices_begin(); pVertex != p->vertices_end(); ++pVertex)
	{
		Point3d point = pVertex->point();
		file << point.x() << " " << point.y() << " " << point.z() <<" ";
		file << (int)(pVertex->color(0)*255) << " " << (int)(pVertex->color(1)*255) << " " << (int)(pVertex->color(2)*255) << "\n";
	}
	for(auto pFacet = p->facets_begin(); pFacet != p->facets_end(); ++pFacet)
	{
		Halfedge_around_facet_circulator hE = pFacet->facet_begin();
		file << "3 ";
		do{
			file << hE->vertex()->tag() << " ";
			
		}while(++hE!=pFacet->facet_begin());
		file << "\n";
	}
	file.close();
}

#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(mepp_component_Correspondence_plugin, mepp_component_Correspondence_plugin);
#endif

#endif
