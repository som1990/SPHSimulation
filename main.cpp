/* Author: Soumitra Goswami
// Date: 5/01/2018
// Description: This is a C++ and QT5 program that simulates a Lagrangian Smooth Particle Hydrodynamics (SPH) simulation
//              Main.cpp is the entry point to the QT widget. 
//				Main.cpp -> sphsimulation.h (QT Widget) -> solver.h (Our actual simulation code)
//              
//              QT works in having signals and slots. "Signals" send out info as messages and slots recieve them.
//              For all the setup of signals and slots check mainwindow.hpp and generation and setup of functionality of UI.
//              The actual visual generation of UI is done via QT Designer.  
//              
*/
#include "mainwindow.hpp"
#include "sphsimulation.h"
#include <QtWidgets/QApplication>
#include <QMainWindow>


/* Brief Overview of workflow with QT5
 ->  A window is generated (mainwindow.h) which contains my UI buttons and an OpenGLWidget. 
 -> This OpenGLWidget is connected to "sphsimulation.h" which contains SLOTS to listen to SIGNALS sent by "mainwindow.h"
 -> "sphsimulation.h" then transfers these incoming SIGNALS to solver.h which contains the bulk of our simulation code.
*/

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);

	QSurfaceFormat format;
	format.setRenderableType(QSurfaceFormat::OpenGL);
	format.setProfile(QSurfaceFormat::CompatibilityProfile);
	format.setVersion(2, 1);
	

	SPHSimulation *widget = new SPHSimulation;
	widget->setFormat(format);
	Mainwindow w;
	w.show();

	return a.exec();
}
