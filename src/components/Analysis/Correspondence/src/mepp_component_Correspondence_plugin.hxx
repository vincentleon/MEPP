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
		{}

		void init(mainwindow* mainWindow, QList<QMdiSubWindow *> lw)
		{
			this->mw = mainWindow;
			this->lwindow = lw;
			this->mPluginName = this->metaObject()->className();
			mParentMenu = mainWindow->menuAnalysis_Filtering;
			m_hasNotBeenPainted = false;
			actionCorrespondence = new QAction(tr("Correspondence"), this);
			if(actionCorrespondence)
			{
				connect(actionCorrespondence, SIGNAL(triggered()),this, SLOT(OnCorrespondence()));
			}
		}

		QList<QAction*> actions() const
		{
			return QList<QAction*>() << actionCorrespondence;
		}

		virtual void pre_draw();
		virtual void post_draw();
		virtual void pre_draw_all_scene() {}
		virtual void post_draw_all_scene() {}

		virtual void OnMouseLeftDown(QMouseEvent *event);
		virtual void OnMouseLeftUp(QMouseEvent *event) {}
		virtual void OnMouseRightDown(QMouseEvent *event) {}
		virtual void OnMouseRightUp(QMouseEvent *event) {}
		virtual void OnMouseMotion(QMouseEvent *event);
		virtual void OnMouseWheel(QWheelEvent *event) {}
		virtual void OnKeyPress(QKeyEvent *event) {}
		virtual void OnKeyRelease(QKeyEvent *event) {}
		
		void PaintStart();
		
		QGLFramebufferObject * m_fbo;
		
	public slots:

		void OnCorrespondence();

  private:
    QAction *actionCorrespondence;
    bool m_hasNotBeenPainted;
};

#endif

#endif
