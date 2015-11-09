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
#include "SegmentController.h"

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
			delete actionSVM;
			delete actionCompare;
			delete actionLearn;
			delete actionPrepareData;
			delete actionGlue;
			delete actionUnion;
		}

		void init(mainwindow* mainWindow, QList<QMdiSubWindow *> lw)
		{
			this->mw = mainWindow;
			this->lwindow = lw;
			this->mPluginName = this->metaObject()->className();
			mParentMenu = mainWindow->menuAnalysis_Filtering;
			m_hasNotBeenPainted = true;
			m_currentLabel = 0;
			
			actionCorrespondence = new QAction(tr("Ellipse Correspondence"), this);
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
			actionSVM = new QAction(tr("SVM correspondence"),this);
			if( actionSVM )
			{
				connect(actionSVM, SIGNAL(triggered()),this, SLOT(OnSVM()));
			}
			actionCompare = new QAction(tr("Ellipse / SVM Comparison"),this);
			if( actionCompare )
			{
				connect(actionCompare, SIGNAL(triggered()),this, SLOT(OnCompareMethods()));
			}
			actionLearn = new QAction(tr("Learn"),this);
			if( actionLearn )
			{
				connect(actionLearn, SIGNAL(triggered()),this, SLOT(OnLearn()));
			}
			actionPrepareData = new QAction(tr("Prepare Data"),this);
			if( actionPrepareData )
			{
				connect(actionPrepareData, SIGNAL(triggered()),this, SLOT(OnPrepareData()));
			}
			actionCut = new QAction(tr("Cut"),this);
			if( actionCut )
			{
				connect(actionCut, SIGNAL(triggered()),this, SLOT(OnCut()));
			}
			actionAddSegment = new QAction(tr("Add segment"),this);
			if( actionAddSegment )
			{
				connect(actionAddSegment, SIGNAL(triggered()),this, SLOT(OnAddSegment()));
			}
			actionGlue = new QAction(tr("Glue"),this);
			if( actionGlue ) 
			{
				connect(actionGlue, SIGNAL(triggered()),this, SLOT(OnGlue()));
			}
			actionUnion = new QAction(tr("Union"),this);
			if( actionUnion ) 
			{
				connect(actionUnion, SIGNAL(triggered()),this, SLOT(OnUnion()));
			}
		}

		QList<QAction*> actions() const
		{
			return QList<QAction*>() << actionPainting
				<< NULL // menu separator
				<< actionCorrespondence
				<< actionShowDescriptor
				<< actionCut
				<< actionAddSegment
				<< actionGlue
				<< NULL // menu separator
				<< actionMahalanobis
				<< NULL // menu separator
				<< actionSVM
				<< NULL // menu separator
				<< actionCompare
				<< NULL // menu separator
				<< actionLearn
				<< NULL // menu separator
				<< actionPrepareData
				<< actionUnion;
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
		void compareToDatasetMahalanobis(Correspondence_ComponentPtr sourceCorrespondence, int sourceID);
		void compareToDatasetSVM(Correspondence_ComponentPtr sourceCorrespondence, int sourceID);
		
		void PaintStart(Viewer * view);
		
		QGLFramebufferObject * m_fbo;
		std::set<int> m_paintedFacets;
		std::vector<Facet_handle> m_facets; // Random access to faces
		
	public slots:

		void OnCorrespondence();
		void showDescriptor();
		void OnPainting();
		void OnMahalanobis();
		void OnSVM();
		void OnCompareMethods();
		void OnLearn();
		void OnPrepareData();
		void OnCut();
		void OnAddSegment();
		void OnGlue();
		void OnUnion();
		
  private:
    int m_currentLabel;
    QAction *actionCorrespondence;
    QAction *actionShowDescriptor;
    QAction *actionPainting;
    QAction *actionMahalanobis;
    QAction *actionSVM;
    QAction *actionCompare;
    QAction *actionLearn;
    QAction *actionPrepareData;
    QAction *actionCut;
    QAction *actionAddSegment;
    QAction *actionGlue;
    QAction *actionUnion;
    bool m_hasNotBeenPainted;
    bool m_PaintingMode;
};

#endif

#endif
