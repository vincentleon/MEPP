#ifndef HEADER_MEPP_COMPONENT_CORRESPONDENCE_PLUGIN_INTERFACE_H
#define HEADER_MEPP_COMPONENT_CORRESPONDENCE_PLUGIN_INTERFACE_H

#include <QtGlobal> // important, for QT_VERSION

#include <QObject>

#include <mepp_config.h>

#ifdef BUILD_component_Correspondence

#include "../../../../mepp/mepp_component_plugin_interface.h"

#include <QAction>
#include <QtPlugin>
#include <QGLFramebufferObject>
#include <set>

#include "Correspondence_Component.h"

typedef boost::shared_ptr<Correspondence_Component> Correspondence_ComponentPtr;

/**
 \class	mepp_component_Correspondence_plugin

 \brief	Mepp component Correspondence plugin.

 */
class mepp_component_Correspondence_plugin :
  public QObject,
  public mepp_component_plugin_interface
{
	Q_OBJECT
	Q_INTERFACES(mepp_component_plugin_interface);
#if QT_VERSION >= 0x050000
	Q_PLUGIN_METADATA(IID "mepp_component_Correspondence_plugin")
#endif

	public:
		mepp_component_Correspondence_plugin() : mepp_component_plugin_interface() {}
		~mepp_component_Correspondence_plugin()
		{
			delete actionCorrespondence;
			delete actionShowDescriptor;
			delete actionPainting;
			delete actionMahalanobis;
		}

		void init(mainwindow* mainWindow, QList<QMdiSubWindow *> lw)
		{
			this->mw = mainWindow;
			this->lwindow = lw;
			this->mPluginName = this->metaObject()->className();
			mParentMenu = mainWindow->menuAnalysis_Filtering;
			m_hasNotBeenPainted = true;
			m_currentLabel = 0;
			
			actionCorrespondence = new QAction(tr("Correspondence"), this);
			if(actionCorrespondence)
			{
				connect(actionCorrespondence, SIGNAL(triggered()),this, SLOT(OnCorrespondence()));
			}
			
			actionShowDescriptor = new QAction(tr("Show descriptor"), this);
			if( actionShowDescriptor )
			{
				connect(actionShowDescriptor, SIGNAL(triggered()),this, SLOT(showDescriptor()));
			}
			
			actionPainting = new QAction(tr("Painting"),this);
			if( actionPainting )
			{
				connect(actionPainting, SIGNAL(triggered()),this, SLOT(OnPainting()));
			}
			
			actionMahalanobis = new QAction(tr("Mahalanobis"),this);
			if( actionMahalanobis )
			{
				connect(actionMahalanobis, SIGNAL(triggered()),this, SLOT(OnMahalanobis()));
			}
		}

		QList<QAction*> actions() const
		{
			return QList<QAction*>() << actionPainting
				<< NULL // menu separator
				<< actionCorrespondence
				<< actionShowDescriptor
				<< NULL // menu separator
				<< actionMahalanobis;
		}

		virtual void pre_draw();
		virtual void post_draw();
		virtual void pre_draw_all_scene() {}
		virtual void post_draw_all_scene() {}

		virtual void OnMouseLeftDown(QMouseEvent *event);
		virtual void OnMouseLeftUp(QMouseEvent *event);
		virtual void OnMouseRightDown(QMouseEvent *event) {}
		virtual void OnMouseRightUp(QMouseEvent *event) {}
		virtual void OnMouseMotion(QMouseEvent *event);
		virtual void OnMouseWheel(QWheelEvent *event) {}
		virtual void OnKeyPress(QKeyEvent *event) {}
		virtual void OnKeyRelease(QKeyEvent *event) {}
		
		void compareToDataset(Correspondence_ComponentPtr sourceCorrespondence, int sourceID);
		
		void PaintStart(Viewer * view);
		
		QGLFramebufferObject * m_fbo;
		std::set<int> m_paintedFacets;
		std::vector<Facet_handle> m_facets; // Random access to faces
		
	public slots:

		void OnCorrespondence();
		void showDescriptor();
		void OnPainting();
		void OnMahalanobis();

  private:
    int m_currentLabel;
    QAction *actionCorrespondence;
    QAction *actionShowDescriptor;
    QAction *actionPainting;
    QAction *actionMahalanobis;
    bool m_hasNotBeenPainted;
    bool m_PaintingMode;
};

#endif

#endif
