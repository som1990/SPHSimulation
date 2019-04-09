#include "sphsimulation.h"
#include <qdebug.h>
#include <qstring.h>
#include <qopenglshaderprogram.h>
#include "Vertex.h"
#include "Particle.h"

/*
static const Vertex sg_vertexes[] = {
	Vertex(QVector3D( 0.00f,  0.75f, 1.00f), QVector3D(1.00f, 0.00f ,0.00f)),
	Vertex(QVector3D( 0.75f, -0.75f, 1.00f), QVector3D(0.00f, 1.00f, 0.00f)),
	Vertex(QVector3D(-0.75f, -0.75f, 1.00f), QVector3D(0.00f, 0.00f, 1.00f))
};
*/

SPHSimulation::SPHSimulation(QWidget *parent):
	QOpenGLWidget(parent)
{
	connect(&timer, SIGNAL(timeout()), this, SLOT(update()));
	//timer.start(16);
}

SPHSimulation::~SPHSimulation()
{
	makeCurrent();
	teardownGL();
}

void printList(std::vector<Vertex> &list1)
{
	for (int i = 0; i < list1.size(); i++)
	{
		qDebug() << qPrintable("Position  ") << i << qPrintable(": ") << list1[i].position() << endl;
		qDebug() << qPrintable("Color  ") << i << qPrintable(": ") << list1[i].color() << endl;
	}
}

void SPHSimulation::initializeGL()
{
	vList.clear();
	pList.clear();
	timer.start(16);
	initializeOpenGLFunctions();
	printContextInformation();
	sph.generateParticleList();
	pList = sph.partList();
//	vList = sph.vList();
//	qDebug() << vList.size() << " Particles extracted";
//	printList(vList);
	//OpenGL 2.0 syntax here
	glClearColor(0.2f, 0.2f, 0.2f, 1.0f);

	//glColor3f(0.0f, 1.0f, 0.0f);
	glPointSize(5.0);
	
// Trying OpenGL 4.0 but couldn't get vertex buffers to work.	
/*	{
		m_program = new QOpenGLShaderProgram();
		m_program->addShaderFromSourceFile(QOpenGLShader::Vertex, "Shaders/simple.vert");
		m_program->addShaderFromSourceFile(QOpenGLShader::Fragment, "Shaders/simple.frag");
		m_program->link();
		m_program->bind();

		m_vertex.create();
		m_vertex.bind();
		m_vertex.setUsagePattern(QOpenGLBuffer::StaticDraw);
//		m_vertex.allocate(sg_vertexes, sizeof(sg_vertexes));
		m_vertex.allocate(&vList, vList.size()*sizeof(Vertex));

		m_object.create();
		m_object.bind();
		m_program->enableAttributeArray(0);
		m_program->enableAttributeArray(1);
		m_program->setAttributeBuffer(0, GL_FLOAT, Vertex::positionOffset(), Vertex::PositionTupleSize, Vertex::stride());
		m_program->setAttributeBuffer(1, GL_FLOAT, Vertex::colorOffset(), Vertex::ColorTupleSize, Vertex::stride());
		

		m_object.release();
		m_vertex.release();
		m_program->release();
	}
*/
}

void SPHSimulation::initializeValues()
{

}

void SPHSimulation::resizeGL(int w, int h)
{
	int side = qMin(w, h);
	glViewport((w-side)/2, (h-side)/2, side, side);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	qreal  aspectratio = qreal(w) / qreal(h);
	glOrtho(-1.0, 1.0, -1.0, 1.0, 4, -2);
	glMatrixMode(GL_MODELVIEW);
}

void SPHSimulation::paintGL()
{
	pList.clear();
	glClear(GL_COLOR_BUFFER_BIT);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_POINT_SMOOTH);
	glBegin(GL_POINTS);
	pList = sph.partList();
	for (std::pair<int,Particle> element : pList)
	{
		int id = element.first;
		Particle p = element.second;
		QVector3D col = p.vertex().color();
		glColor3f(col.x(),col.y(),col.z());
		glVertex2f(p.vertex().position().x(), p.vertex().position().y());

	}
	glEnd();
	//openGL 4.0 stuff 
	
/*	m_program->bind();
	{
		m_object.bind();
		//glEnable(GL_POINT_SMOOTH);
		glDrawArrays(GL_POINTS,0, vList.size());
		qDebug() << vList.size() *  sizeof(Vertex) << " Particles Drawn";
		m_object.release();
	}
	
	m_program->release();
*/
	glDisable(GL_BLEND);
	glDisable(GL_POINT_SMOOTH);
	glFlush();
}

void SPHSimulation::teardownGL()
{
	//m_object.destroy();
	//m_vertex.destroy();
	//delete m_program;
}

void SPHSimulation::update()
{
	if (pause == false && timer.isActive())
	{
		sph.simulate();
	//	timer.stop();
		
	}
	QOpenGLWidget::update();

}



void SPHSimulation::pauseSimulation()
{
	if (pause == false) {
		qDebug() << qPrintable("Pausing Simulation");
		pause = true;
		timer.stop();
		
	} 
	else
	{
		qDebug() << qPrintable("Resuming Simulation");
		pause = false;
		timer.start();
	}
}

void SPHSimulation::restartSimulation()
{
	pList.clear();
	sph.generateParticleList();
	pList = sph.partList();
	repaint();
}

void SPHSimulation::settimestep(float t)
{
	if (sph.dt() != t  )
	{

		sph.settimestep(t);
		qDebug() << "TimeStep Set to: " << double(t)<<endl;
		emit timestepChanged(t);
	}
}

void SPHSimulation::setparticleCount(size_t pCount)
{

	if (pCount != sph.nParticles())
	{
		sph.setparticleCount(pCount);
		qDebug() << "Particle Count Set to: " << int(pCount) << endl;
		emit pCountChanged(pCount);
	}
}

void SPHSimulation::setmass(float m)
{
	if (m != sph.mass())
	{
		sph.setmass(m);
		qDebug() << "Mass Set to: " << double(m) << endl;
		emit massChanged(m);
	}
}

void SPHSimulation::setradius(float r)
{
	if (r != sph.radius())
	{
		sph.setradius(r);
		qDebug() << "Radius set to: " << double(r) << endl;
		emit radiusChanged(r);
	}
}

void SPHSimulation::setgravity(float gravY)
{
	if (gravY != sph.gravity().y())
	{
		sph.setgravity(gravY);
		qDebug() << "Gravity Set to: " << double(gravY) << endl;
		emit gravityChanged(gravY);
	}
}

void SPHSimulation::setpresConst(float pConst)
{
	if (pConst != sph.presConst())
	{
		sph.setpresConst(pConst);
		qDebug() << "Pressure Constant Set to: " << double(pConst) << endl;
		emit presConstChanged(pConst);
	}
}


void SPHSimulation::setpresGamma(float gamma)
{
	if (gamma != sph.presGamma())
	{
		sph.setpresGamma(gamma);
		qDebug() << "presGamma Set to: " << double(gamma) << endl;
		emit presGammaChanged(gamma);
	}
}

void SPHSimulation::setdensity(float rho)
{
	if (rho != sph.initDensity())
	{
		sph.setdensity(rho);
		qDebug() << "initDensity Set to: " << double(rho) << endl;
		emit densityChanged(rho);
	}
}

void SPHSimulation::setviscConst(float vConst)
{
	if (vConst != sph.viscConst())
	{
		sph.setviscConst(vConst);
		qDebug() << "vConst Set to: " << double(vConst) << endl;
		emit viscConstChanged(vConst);
	}
}

void SPHSimulation::setviscE(float vE)
{
	if (vE != sph.viscE())
	{
		sph.setviscE(vE);
		qDebug() << "ViscE Set to: " << double(vE) << endl;
		emit viscEChanged(vE);
	}
}

void SPHSimulation::setWallSticky(float wSticky)
{
	if (wSticky != sph.wallSticky())
	{
		sph.setWallSticky(wSticky);
		qDebug() << "Wall Sticky Set to: " << double(wSticky) << endl;
		emit wallStickyChanged(wSticky);
	}
}

void SPHSimulation::printContextInformation()
{
	QString glType;
	QString glVersion;
	QString glProfile;

	glType = (context()->isOpenGLES()) ? "OpenGL ES" : "OpenGL";
	glVersion = reinterpret_cast<const char*>(glGetString(GL_VERSION));

#define CASE(c) case QSurfaceFormat::c: glProfile = #c; break
	switch (format().profile())
	{
		CASE(NoProfile);
		CASE(CoreProfile);
		CASE(CompatibilityProfile);
	}
#undef CASE

	qDebug() << qPrintable(glType) << qPrintable(glVersion) << "(" << qPrintable(glProfile) << ")";
}
 