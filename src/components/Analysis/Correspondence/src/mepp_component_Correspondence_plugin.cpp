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

#include "Correspondence_Component.h"

typedef boost::shared_ptr<Correspondence_Component> Correspondence_ComponentPtr;

void mepp_component_Correspondence_plugin::post_draw()
{
	// active viewer
	/*if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		if (doesExistComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewer, polyhedron_ptr)) // important !!!
		{
			Correspondence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewer, polyhedron_ptr);
			if (component_ptr->get_init() == 2 && polyhedron_ptr->Correspondence_is_calculated)
			{
				glPushMatrix();
					// here your code
					if (component_ptr->get_displayMinDirections()==true)
					{
						glColor3f(1.0,0.0,0.0);
						::glBegin(GL_LINES);
						for (Vertex_iterator pVertex = polyhedron_ptr->vertices_begin();	pVertex !=	polyhedron_ptr->vertices_end();	pVertex++)
						{
							{

								float RayonMoyen=polyhedron_ptr->average_edge_length_around(pVertex);
								const Point3d& p1 = pVertex->point()-0.4*RayonMoyen*(pVertex->VKminCurv);
								const Point3d& p2 = pVertex->point()+0.4*RayonMoyen*(pVertex->VKminCurv);
								::glVertex3f(p1[0],p1[1],p1[2]);
								::glVertex3f(p2[0],p2[1],p2[2]);
							}


			https://www.google.fr/?gfe_rd=cr&ei=wGPkVcbTOoqEcL79MQ			}
						::glEnd();
					}

					if (component_ptr->get_displayMaxDirections()==true)
					{
						glColor3f(0.0,0.0,1.0);
						::glBegin(GL_LINES);
						for (Vertex_iterator pVertex = polyhedron_ptr->vertices_begin();	pVertex !=	polyhedron_ptr->vertices_end();	pVertex++)
						{
							{

								float RayonMoyen=polyhedron_ptr->average_edge_length_around(pVertex);
								const Point3d& p1 = pVertex->point()-0.4*RayonMoyen*(pVertex->VKmaxCurv);
								const Point3d& p2 = pVertex->point()+0.4*RayonMoyen*(pVertex->VKmaxCurv);
								::glVertex3f(p1[0],p1[1],p1[2]);
								::glVertex3f(p2[0],p2[1],p2[2]);
							}


						}
						::glEnd();
					}
				glPopMatrix();
			}
		}
	}*/
}

void mepp_component_Correspondence_plugin::OnMouseLeftDown(QMouseEvent *event)
{
	if (mw->activeMdiChild() != 0)
	{
		/*Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		if (doesExistComponentForViewer<Correspondence_ComponentPtr, Correspondence_Component>(viewer, polyhedron_ptr)) // important !!!
		{

		}*/
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
				QApplication::setOverrideCursor(Qt::WaitCursor);

				mw->statusBar()->showMessage(tr("Correspondence..."));
				
				component_ptr->readSelectionBasedOnColor(polyhedron_ptr);
				
				component_ptr->selectPoint(polyhedron_ptr);
				
				mw->statusBar()->showMessage(tr("Correspondence is done"));

				component_ptr->set_init(2);
				viewer->recreateListsAndUpdateGL();
				
			}
		}
	}
	QApplication::restoreOverrideCursor();
	
}


#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(mepp_component_Correspondence_plugin, mepp_component_Correspondence_plugin);
#endif

#endif
