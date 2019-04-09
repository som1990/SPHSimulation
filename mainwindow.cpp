#include "mainwindow.hpp"
#include <QDebug>
#include <iostream>
#include <QtWidgets>

#include "sphsimulation.h"
#include "Solver.h"

//Here is where all the connections of Signals and slots takes place from the UI to this class 
Mainwindow::Mainwindow(QWidget * parent) : QWidget(parent),ui(new Ui::Mainwindow) {
	ui->setupUi(this);
//	QPushButton *qButton = ui->quitButton;
	connect(ui->quitButton, SIGNAL(clicked()), this, SLOT(winQuit()));
	connect(ui->pauseButton, SIGNAL(clicked()), ui->openGLWidget, SLOT(pauseSimulation()));
	connect(ui->restartButton, SIGNAL(clicked()), ui->openGLWidget, SLOT(restartSimulation()));

	connect(ui->v_timeStep, SIGNAL(editingFinished()), this, SLOT(settimestep()));
	connect(this, SIGNAL(timeStepChanged(float)), ui->openGLWidget, SLOT(settimestep(float)));

	connect(ui->v_numParticles, SIGNAL(editingFinished()), this, SLOT(setparticleCount()));
	connect(this, SIGNAL(pCountChanged(size_t)), ui->openGLWidget, SLOT(setparticleCount(size_t)));

	connect(ui->v_mass, SIGNAL(editingFinished()), this, SLOT(setmass()));
	connect(this, SIGNAL(massChanged(float)), ui->openGLWidget, SLOT(setmass(float)));

	connect(ui->v_radius, SIGNAL(editingFinished()), this, SLOT(setradius()));
	connect(this, SIGNAL(radiusChanged(float)), ui->openGLWidget, SLOT(setradius(float)));

	connect(ui->v_gravity, SIGNAL(editingFinished()), this, SLOT(setgravity()));
	connect(this, SIGNAL(gravityChanged(float)), ui->openGLWidget, SLOT(setgravity(float)));

	connect(ui->v_pConst, SIGNAL(editingFinished()), this, SLOT(setpresConst()));
	connect(this, SIGNAL(presConstChanged(float)), ui->openGLWidget, SLOT(setpresConst(float)));

	connect(ui->v_pGamma, SIGNAL(editingFinished()), this, SLOT(setpresGamma()));
	connect(this, SIGNAL(presGammaChanged(float)), ui->openGLWidget, SLOT(setpresGamma(float)));

	connect(ui->v_density, SIGNAL(editingFinished()), this, SLOT(setdensity()));
	connect(this, SIGNAL(densityChanged(float)), ui->openGLWidget, SLOT(setdensity(float)));

	connect(ui->v_vConst, SIGNAL(editingFinished()), this, SLOT(setviscConst()));
	connect(this, SIGNAL(viscConstChanged(float)), ui->openGLWidget, SLOT(setviscConst(float)));

	connect(ui->v_vE, SIGNAL(editingFinished()), this, SLOT(setviscE()));
	connect(this, SIGNAL(viscEChanged(float)), ui->openGLWidget, SLOT(setviscE(float)));

	connect(ui->v_wallSticky, SIGNAL(editingFinished()), this, SLOT(setWallSticky()));
	connect(this, SIGNAL(wallStickyChanged(float)), ui->openGLWidget, SLOT(setWallSticky(float)));
}

Mainwindow::~Mainwindow() {
	delete ui;
}

void Mainwindow::winQuit()
{
	qDebug() << "Hi";
	this->close();
}

void Mainwindow::settimestep()
{
	QString s_t = ui->v_timeStep->text();
	float t = s_t.toFloat();
	emit timeStepChanged(t);
}

void Mainwindow::setparticleCount()
{
	QString s_p = ui->v_numParticles->text();
	size_t nparts = s_p.toInt();
	emit pCountChanged(nparts);
}

void Mainwindow::setmass()
{
	QString s_m = ui->v_mass->text();
	float m = s_m.toFloat();
	emit massChanged(m);
}

void Mainwindow::setradius()
{
	QString s_r = ui->v_radius->text();
	float r = s_r.toFloat();
	emit radiusChanged(r);
}

void Mainwindow::setgravity()
{
	QString s_g = ui->v_gravity->text();
	float g = s_g.toFloat();
	emit gravityChanged(g);
}

void Mainwindow::setpresConst()
{
	QString s_pC = ui->v_pConst->text();
	float pC = s_pC.toFloat();
	emit presConstChanged(pC);
}

void Mainwindow::setpresGamma()
{
	QString s_pG = ui->v_pGamma->text();
	float pG = s_pG.toFloat();
	emit presGammaChanged(pG);
}

void Mainwindow::setdensity()
{
	QString s_dens = ui->v_density->text();
	float dens = s_dens.toFloat();
	emit densityChanged(dens);
}

void Mainwindow::setviscConst()
{
	QString s_vC = ui->v_vConst->text();
	float vC = s_vC.toFloat();
	emit viscConstChanged(vC);
}

void Mainwindow::setviscE()
{
	QString s_vE = ui->v_vE->text();
	float vE = s_vE.toFloat();
	emit viscEChanged(vE);
}

void Mainwindow::setWallSticky()
{
	QString s_wallSticky = ui->v_wallSticky->text();
	float wallStick = s_wallSticky.toFloat();
	emit wallStickyChanged(wallStick);
}

void Mainwindow::keyPressEvent(QKeyEvent *event)
{
	if (event->key() == Qt::Key_Escape)
		this->close();
	else
		QWidget::keyPressEvent(event);
}
