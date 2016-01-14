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

#include <stack>
#include <ctime>

#include <dirent.h>

std::stack<clock_t> tictoc_stack;

void tic() {
    tictoc_stack.push(clock());
}

void toc() {
    std::cout << "Time elapsed: "
              << ((double)(clock() - tictoc_stack.top())) / CLOCKS_PER_SEC
              << std::endl;
    tictoc_stack.pop();
}


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
		m_hasNotBeenPainted=false;
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
			if(posE = meshIDString.size())
			{
				posE = meshIDString.find_last_of(".off");
			}
			
			int meshID = atoi(meshIDString.substr(posB+1 ,posE).c_str());

			SettingsDialog dial;
			if (dial.exec() == QDialog::Accepted)
			{
				nbLabel = dial.spinBox->value();
				
				QApplication::setOverrideCursor(Qt::WaitCursor);

				mw->statusBar()->showMessage(tr(" Ellipse Correspondence..."));
				
				component_ptr->readSelectionBasedOnColor(polyhedron_ptr);
				
				component_ptr->initParameters(nbLabel,meshID,meshIDString.substr(0,posB));
				
				component_ptr->readDescriptor(polyhedron_ptr,meshIDString.substr(0,posB));
				
				component_ptr->initializeEllipsoid(polyhedron_ptr);
				
				std::cout << "Compute Ellipse Parameters" << std::endl;
				tic();
				component_ptr->computeEllipseParameters(polyhedron_ptr);
				toc();
				
				component_ptr->compareDescriptorToEllipse(polyhedron_ptr);
				
				mw->statusBar()->showMessage(tr("Correspondence is done"));

				component_ptr->set_init(2);
				viewer->showAllScene();
				viewer->recreateListsAndUpdateGL();
				
				std::cout << "Compare to dataset" << std::endl;
				//tic();
				compareToDataset(component_ptr,meshID);
				//toc();
			}
		}
	}
	QApplication::restoreOverrideCursor();
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
			unsigned posE = meshIDString.find_last_of(".ply");
			
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
				
				component_ptr->compareDescriptorToSVM(polyhedron_ptr);
				
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
	SettingsDialog dial;
	if (dial.exec() == QDialog::Accepted)
	{
		nbLabel = dial.spinBox->value();
 	}
	
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
			if(isPLY)
			{
				emit(mw->get_actionNewEmpty()->trigger());
				
				for(int j=0; j<lwindow.size();j++)
				{
					viewerI = (Viewer*)qobject_cast<QWidget *>(lwindow[j]->widget());
					if(viewerI->getScenePtr()->get_polyhedron()->empty())
					{	
						std::cout << "File :"<< files[i] << std::endl;
						viewerI->getScenePtr()->add_mesh(path+"/"+files[i].c_str(),0,NULL,viewerI);
						if(viewerI->getScenePtr()->get_polyhedron()->empty())
						{continue;}
						
						PolyhedronPtr polyhedron_ptr = viewerI->getScenePtr()->get_polyhedron();
						
						std::string meshIDString = polyhedron_ptr->pName;
						unsigned posB = meshIDString.find_last_of("/");
						unsigned posE = meshIDString.find_last_of(".ply");
						int meshID = atoi(meshIDString.substr(posB+1 ,posE).c_str());
			
						Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewerI, polyhedron_ptr);
						
						component_ptr->scaleMesh(mbbox,polyhedron_ptr);
						
						component_ptr->initParameters(nbLabel,meshID,dir);
						component_ptr->readDescriptor(polyhedron_ptr,dir);
						component_ptr->setEllipse(sourceCorrespondence->getEllipse());
						component_ptr->setCentreDescriptor(sourceCorrespondence->getCentreDescr());
						component_ptr->compareDescriptorToEllipse(polyhedron_ptr);
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

void mepp_component_Correspondence_plugin::compareToDatasetSVM(Correspondence_ComponentPtr sourceCorrespondence, int sourceID)
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
						
						if(meshID == sourceID)
						{continue;}
						
						Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewerI, polyhedron_ptr);
						
						component_ptr->initParameters(nbLabel,meshID,dir);	
						component_ptr->readDescriptor(polyhedron_ptr,meshIDString.substr(0,posB));
						
						component_ptr->setSVM(sourceCorrespondence->getSVM());
						
						component_ptr->compareDescriptorToSVM(polyhedron_ptr);
				
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
			unsigned posE = meshIDString.find_last_of(".ply");
			
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
	QString path = QFileDialog::getExistingDirectory (mw, tr("Choose directory to save parts"),"/home/leon/");	
	std::string directory = path.toStdString();
	
	
	// Save all polyhedron data and descriptor to directory
	for(unsigned i=0; i<scene->get_nb_polyhedrons();++i)
	{
		PolyhedronPtr p = scene->get_polyhedron(i);
		p->keep_largest_connected_components(1);
		Correspondence_ComponentPtr corres = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(mainViewer,p);
		std::stringstream filename; 
		filename << directory << "/" << i << ".off";
		corres->initParameters(4,i,directory);
		corres->getShape().m_meshID = i;
		corres->saveDescriptor(p,directory);
		p->write_off(filename.str(),true,true);
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
	Viewer* mainViewer = (Viewer*)qobject_cast<QWidget *>(lwindow[0]->widget());
	PolyhedronPtr mainPoly = mainViewer->getScenePtr()->get_polyhedron(0);
	Correspondence_ComponentPtr mainCorres = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(mainViewer, mainPoly);
	//int cMesh = mainViewer->getScenePtr()->get_polyhedron(1);
	SegmentController & mainSegCtr = mainCorres->m_segCtr;
	PolyhedronPtr partPoly = mainViewer->getScenePtr()->get_polyhedron(1);
		
	//mainSegCtr.sewSegments(mainViewer,partPoly,mainPoly);
	
	mainSegCtr.softICP(mainViewer,partPoly,mainPoly);
	
	mainViewer->recreateListsAndUpdateGL();
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
				tic();
				component_ptr->learnSVMClassifier(polyhedron_ptr);
				toc();
				
				component_ptr->initializeEllipsoid(polyhedron_ptr);
				std::cout << "Ellipse parameters : " << std::endl;
				tic();
				component_ptr->computeEllipseParameters(polyhedron_ptr);
				toc();
				std::cout << "compareToEllipse : " << std::endl;
				tic();
				component_ptr->compareDescriptorToEllipse(polyhedron_ptr);
				toc();
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
						tic();
						svmcomponent_ptr->setSVM(component_ptr->getSVM());
						svmcomponent_ptr->compareDescriptorToSVM(polyhedron_copy_ptr);
						toc();
						viewerI->recreateListsAndUpdateGL();
					}	
				}
			}
		}
	}
}

#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(mepp_component_Correspondence_plugin, mepp_component_Correspondence_plugin);
#endif

#endif
