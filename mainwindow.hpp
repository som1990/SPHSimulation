/* mainwindow.hpp
// Author:Soumitra Goswami
// Date: 05/01/2018
// Description: Here we add functionality to our UI setup done via the QT Designer
//              I also translate the SIGNALS from UI and generate my own SIGNALS to be picked up by "sphsimulation.h". 
//
*/
#ifndef MAINWINDOW_HPP
#define MAINWINDOW_HPP
#include <QWidget>
#include "ui_mainwindow.h"

namespace Ui {
	class Mainwindow;

}

class Mainwindow : public QWidget {
	Q_OBJECT

public:
	Mainwindow(QWidget * parent = Q_NULLPTR);
	~Mainwindow();

// Recieve SIGNALs from UI via these SLOTs	
public slots:
	void winQuit();
	void settimestep();
	void setparticleCount();
	void setmass();
	void setradius();

	void setgravity();
	void setpresConst();
	void setpresGamma();
	void setdensity();
	
	void setviscConst();
	void setviscE();
	void setWallSticky();

//	bool pauseSim();

// Transmit these SIGNALs to be picked up by the "sphsimulation.h" SLOTs 
signals:
	void timeStepChanged(float t);
	void pCountChanged(size_t pCount);
	void massChanged(float m);
	void radiusChanged(float r);

	void gravityChanged(float gravY);
	void presConstChanged(float pConst);
	void presGammaChanged(float gamma);
	void densityChanged(float rho);

	void viscConstChanged(float vConst);
	void viscEChanged(float vE);
	void wallStickyChanged(float wSticky);

protected: 
	void keyPressEvent(QKeyEvent *event);

private:
	Ui::Mainwindow *ui;
};

#endif // MAINWINDOW_HPP