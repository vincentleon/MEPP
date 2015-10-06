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
		
 		//if (doesExistComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewer, polyhedron_ptr))
		//{
			Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewer, polyhedron_ptr);
			
			//if(component_ptr->get_init() != 2)
			//{
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
			//}
		//}
		viewer->recreateListsAndUpdateGL();
	}
}

void mepp_component_Correspondence_plugin::showDescriptor()
{
	/*if(mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
		
		Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewer, polyhedron_ptr);
		int nbLabel = 4;
		
		std::string meshIDString = polyhedron_ptr->pName;
		unsigned posB = meshIDString.find_last_of("/");
		unsigned posE = meshIDString.find_last_of(".ply");
			
		int meshID = atoi(meshIDString.substr(posB+1 ,posE).c_str());
		
		component_ptr->initParameters(nbLabel,meshID);
					
		component_ptr->readDescriptor(polyhedron_ptr);
		component_ptr->showDescriptor(polyhedron_ptr,m_currentLabel);
		m_currentLabel = (m_currentLabel + 1) % nbLabel;
		viewer->recreateListsAndUpdateGL();
	}*/
	Viewer* viewerI = NULL;
	
	int nbLabel = 4;
	for(int i=0; i<lwindow.size();i++)
	{
		viewerI = (Viewer*)qobject_cast<QWidget *>(lwindow[i]->widget());
		if(!viewerI->getScenePtr()->get_polyhedron()->empty())
		{
			PolyhedronPtr polyhedron_ptr = viewerI->getScenePtr()->get_polyhedron();
	
			Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewerI, polyhedron_ptr);
			
			
			std::string meshIDString = polyhedron_ptr->pName;
			unsigned posB = meshIDString.find_last_of("/");
			unsigned posE = meshIDString.find_last_of(".ply");
				
			int meshID = atoi(meshIDString.substr(posB+1 ,posE).c_str());
			
			component_ptr->initParameters(nbLabel,meshID);
						
			component_ptr->readDescriptor(polyhedron_ptr);
			component_ptr->showDescriptor(polyhedron_ptr,m_currentLabel);
			viewerI->recreateListsAndUpdateGL();
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
				
				component_ptr->initParameters(nbLabel,meshID);
				
				component_ptr->readDescriptor(polyhedron_ptr);
				
				component_ptr->initializeEllipsoid(polyhedron_ptr);
				
				component_ptr->computeEllipseParameters(polyhedron_ptr);
				
				component_ptr->compareDescriptorToEllipse(polyhedron_ptr);
				
				mw->statusBar()->showMessage(tr("Correspondence is done"));

				component_ptr->set_init(2);
				viewer->recreateListsAndUpdateGL();
				tic();
				compareToDataset(component_ptr,meshID);
				toc();
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
				
				component_ptr->initParameters(nbLabel,meshID);
				
				component_ptr->readDescriptor(polyhedron_ptr);
				
				component_ptr->learnSVMClassifier(polyhedron_ptr);
				
				component_ptr->compareDescriptorToSVM(polyhedron_ptr);
				
				mw->statusBar()->showMessage(tr("SVM Correspondence is done"));
				
				component_ptr->set_init(2);
				viewer->recreateListsAndUpdateGL();
				tic();
				compareToDatasetSVM(component_ptr,meshID);
				toc();
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
			int nbLabel = 8;
			std::string meshIDString = polyhedron_ptr->pName;
			unsigned posB = meshIDString.find_last_of("/");
			unsigned posE = meshIDString.find_last_of(".off");
			
			int meshID = atoi(meshIDString.substr(posB+1 ,posE).c_str());
			
			SettingsDialog dial;
			if (dial.exec() == QDialog::Accepted)
			{
				
				QApplication::setOverrideCursor(Qt::WaitCursor);

				mw->statusBar()->showMessage(tr("Mahalanobis..."));
				
				component_ptr->readSelectionBasedOnColor(polyhedron_ptr);
				
				component_ptr->initParameters(nbLabel,meshID);
				
				component_ptr->readDescriptor(polyhedron_ptr);
				
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
	
	QString path = QFileDialog::getExistingDirectory (mw, tr("Choose directory"),"/home/leon/");
	
	
	for(int m = 1; m <=73; ++m)
	{
		if (m==65 || m == 39) { continue;}
		if(m == sourceID) {continue;}
		emit(mw->get_actionNewEmpty()->trigger());
		
		for(int i=0; i<lwindow.size();i++)
		{
			viewerI = (Viewer*)qobject_cast<QWidget *>(lwindow[i]->widget());
			if(viewerI->getScenePtr()->get_polyhedron()->empty())
			{
				std::stringstream ss;
				ss << path.toStdString() << "/" << m << ".ply";
				
				viewerI->getScenePtr()->add_mesh(ss.str().c_str(),0,NULL,viewerI);
				
				PolyhedronPtr polyhedron_ptr = viewerI->getScenePtr()->get_polyhedron();	
				Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewerI, polyhedron_ptr);
				component_ptr->initParameters(4,m);
				component_ptr->readDescriptor(polyhedron_ptr);
				component_ptr->setEllipse(sourceCorrespondence->getEllipse());
				component_ptr->setCentreDescriptor(sourceCorrespondence->getCentreDescr());
				component_ptr->compareDescriptorToEllipse(polyhedron_ptr);
				viewerI->recreateListsAndUpdateGL();
				
				/*component_ptr->initParameters(8,m);	
				component_ptr->readDescriptor(polyhedron_ptr);
				component_ptr->setMatrix(sourceCorrespondence->getMatrix());
				component_ptr->setCentreDescriptor(sourceCorrespondence->getCentreDescr());
				component_ptr->compareDescriptorToGaussian(polyhedron_ptr);*/
			}
		}
	}
}

void mepp_component_Correspondence_plugin::compareToDatasetMahalanobis(Correspondence_ComponentPtr sourceCorrespondence, int sourceID)
{
	Viewer* viewerI = NULL;
	
	for(int m = 1; m <= 20; ++m)
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

void mepp_component_Correspondence_plugin::compareToDatasetSVM(Correspondence_ComponentPtr sourceCorrespondence, int sourceID)
{
	Viewer* viewerI = NULL;
	
	for(int m = 1; m <= 73; ++m)
	{	
		if (m==65 || m == 39) { continue;}
		if(m == sourceID) {continue;}
		emit(mw->get_actionNewEmpty()->trigger());
		
		for(int i=0; i<lwindow.size();i++)
		{
			viewerI = (Viewer*)qobject_cast<QWidget *>(lwindow[i]->widget());
			if(viewerI->getScenePtr()->get_polyhedron()->empty())
			{
				std::stringstream ss;
				ss << "/home/leon/modelesFTP/result/" << m << ".ply";
				viewerI->getScenePtr()->add_mesh(ss.str().c_str(),0,NULL,viewerI);
				
				PolyhedronPtr polyhedron_ptr = viewerI->getScenePtr()->get_polyhedron();	
				Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewerI, polyhedron_ptr);
				
				component_ptr->initParameters(4,m);	
				component_ptr->readDescriptor(polyhedron_ptr);
				
				component_ptr->setSVM(sourceCorrespondence->getSVM());
				
				component_ptr->compareDescriptorToSVM(polyhedron_ptr);
		
				viewerI->recreateListsAndUpdateGL();
			}
		}
	}
	
}

void mepp_component_Correspondence_plugin::OnLearn()
{
	Viewer* viewerI = NULL;
	for(int m=1;m<=73;++m)
	{
		if (m==65 || m == 39) { continue;}
		emit(mw->get_actionNewEmpty()->trigger());
		for(int i=0; i<lwindow.size();i++)
		{
			
			viewerI = (Viewer*)qobject_cast<QWidget *>(lwindow[i]->widget());
			if(viewerI->getScenePtr()->get_polyhedron()->empty())
			{
				
				std::stringstream ss;
				ss << "/home/leon/modelesFTP/result/" << m << ".ply";
				
				std::cout << "mesh : " << m << std::endl;
				
				viewerI->getScenePtr()->add_mesh(ss.str().c_str(),0,NULL,viewerI);
				
				PolyhedronPtr polyhedron_ptr = viewerI->getScenePtr()->get_polyhedron();
				
				if(polyhedron_ptr->is_pure_triangle())	
				{
					
					Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewerI, polyhedron_ptr);
					
					component_ptr->initParameters(4,m);
					component_ptr->learnDescriptor(polyhedron_ptr);
					
					viewerI->setWindowTitle(ss.str().c_str());
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
						
						if(polyhedron_ptr->is_pure_triangle())
						{
							std::cout << files[i] << " is pure triangle"<<std::endl;
						}
						else if(polyhedron_ptr->is_pure_quad())
						{
							std::cout << files[i] << " is pure quads"<<std::endl;
						}
						else	
						{
							std::cout << files[i] << " is a mix of triangles and quads"<<std::endl;
						}
						std::string meshIDString = files[i];
						
						int meshID = atoi(meshIDString.c_str());
							
						Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewerI, polyhedron_ptr);
						
						component_ptr->initParameters(nbLabel,meshID);
					}
				}
			}
		}
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
			int nbLabel = 8;
			std::string meshIDString = polyhedron_ptr->pName;
			unsigned posB = meshIDString.find_last_of("/");
			unsigned posE = meshIDString.find_last_of(".off");
			
			meshID = atoi(meshIDString.substr(posB+1 ,posE).c_str());

			SettingsDialog dial;
			if (dial.exec() == QDialog::Accepted)
			{
				QApplication::setOverrideCursor(Qt::WaitCursor);
				nbLabel = dial.spinBox->value();
				mw->statusBar()->showMessage(tr(" Ellipse Correspondence..."));
				
				
				component_ptr->readSelectionBasedOnColor(polyhedron_ptr);
				
				component_ptr->initParameters(nbLabel,meshID);
				
				component_ptr->readDescriptor(polyhedron_ptr);
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
						ss << "/home/leon/datasetHuman/" << meshID << ".off";
						viewerI->getScenePtr()->add_mesh(ss.str().c_str(),0,NULL,viewerI);
						viewerI->setWindowTitle("SVM");
						PolyhedronPtr polyhedron_copy_ptr = viewerI->getScenePtr()->get_polyhedron();	
						Correspondence_ComponentPtr svmcomponent_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewerI, polyhedron_copy_ptr);
						
						svmcomponent_ptr->readSelectionBasedOnColor(polyhedron_copy_ptr);
						
						svmcomponent_ptr->initParameters(8,meshID);	
						
						svmcomponent_ptr->readDescriptor(polyhedron_copy_ptr);
						std::cout << "compareToSVM : " << std::endl;
						tic();
						svmcomponent_ptr->setSVM(component_ptr->getSVM());
						toc();
						svmcomponent_ptr->compareDescriptorToSVM(polyhedron_copy_ptr);
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
