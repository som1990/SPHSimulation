#include "Solver.h"
#include <random>
#include <QDebug>
#include <QtMath>

Q_DECL_CONSTEXPR Solver::Solver()
{
}


Q_DECL_CONSTEXPR Solver::Solver(size_t& pNum, float& pMass, float& pRad, QVector3D& fGravity, float& presConst, float& presGamma, float& iDensity, float& vConst, float& vE, float& wSticky, float& t)
{
	num_particles = pNum;
	p_mass = pMass;
	p_radius = pRad;
	force_g = fGravity;
	pres_const = presConst;
	pres_gamma = presGamma;
	init_density = iDensity;
	visc_const = vConst;
	visc_e = vE;
	wall_sticky = wSticky;
	s_dt = t;
}

Q_DECL_CONSTEXPR Solver::~Solver()
{
}


void Solver::settimestep(float &t)
{
	s_dt = t;
}

void Solver::setparticleCount(size_t &pCount)
{
	num_particles = pCount;
}

void Solver::setmass(float &m)
{
	p_mass = m;
}

void Solver::setradius(float &r)
{
	p_radius = r;
}

void Solver::setgravity(float &gravY)
{
	force_g = QVector3D(0, gravY, 0);

}

void Solver::setpresConst(float &pConst)
{
	pres_const = pConst;

}


void Solver::setpresGamma(float &gamma)
{
	pres_gamma = gamma;
}

void Solver::setdensity(float &rho)
{
	init_density = rho;
}

void Solver::setviscConst(float &vConst)
{
	visc_const = vConst;
}

void Solver::setviscE(float &vE)
{
	visc_e = vE;
}

void Solver::setWallSticky(float wSticky)
{
	wall_sticky = wSticky;
}

Q_DECL_CONSTEXPR const float& Solver::dt() const { return s_dt; }

Q_DECL_CONSTEXPR const size_t& Solver::nParticles() const { return num_particles; }

Q_DECL_CONSTEXPR const float& Solver::mass() const { return p_mass; }

Q_DECL_CONSTEXPR const float& Solver::radius() const { return p_radius; }

Q_DECL_CONSTEXPR const QVector3D& Solver::gravity() const { return force_g; }

Q_DECL_CONSTEXPR const float& Solver::presConst() const { return pres_const; }

Q_DECL_CONSTEXPR const float& Solver::presGamma() const { return pres_gamma; }

Q_DECL_CONSTEXPR const float& Solver::initDensity() const { return init_density; }

Q_DECL_CONSTEXPR const float& Solver::viscConst() const { return visc_const; }

Q_DECL_CONSTEXPR const float& Solver::viscE() const { return visc_e; }

Q_DECL_CONSTEXPR const float& Solver::wallSticky() const { return wall_sticky; }

const std::vector<Vertex> Solver::vList() const{ return s_vList; }
std::map<int, Particle>& Solver::partList() { return pList2;  }


void Solver::generateParticleList()
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-0.2, 0.2);
	if (pList2.size() != 0)
	{
		//pList.clear();
		s_vList.clear();
		pList2.clear();
		occupancyVolume.clear();
	}
//#pragma omp parallel for
	for (int i = 0; i < num_particles; i++)
	{
		Particle p;
		p.setmass(p_mass);
		p.setdensity(init_density);
		p.setradius(p_radius);
		p.setvelocity(QVector3D(0, 0, 0));
		p.setacceleration(QVector3D(0, 0, 0));
		p.setid(i);
		QVector3D pos = QVector3D(dis(gen), dis(gen), 1.0f);
		QVector3D col = QVector3D(0.0f, 0.2f, 1.0f);
		Vertex v = Vertex(pos, col);
		p.setvertex(v);
		//pList.push_back(p);
		pList2.insert(std::make_pair(i, p));
		s_vList.push_back(v);
	}
	qDebug() << pList2.size() << qPrintable(" Particles Generated");
	qDebug() << qPrintable("Size of Particle") << sizeof(Particle);
}




void Solver::generateOccupancyVolume()
{
	//Calculating Min and Max to find the boundaries of the volume
	float xMax, xMin, yMax, yMin;
	xMin = 1.0; xMax = -1.0;
	yMin = 1.0; yMax = -1.0;
//#pragma omp parallel for
	for(std::pair<int,Particle> element : pList2)
	{
		int id = element.first;
		Particle pA = element.second;
		QVector3D pos = pA.position();
		float x = pos.x();
		float y = pos.y();
		if (x > xMax) xMax = x;
		else if (x < xMin) xMin = x;

		if (y > yMax) yMax = y;
		else if (y < yMin ) yMin = y;
	}
	//llc - lower left corner from xMin, yMin
	//urc - upper right corner from xMax, yMax
	QVector3D llc = QVector3D(xMin, yMin, 1.0);
	QVector3D urc = QVector3D(xMax, yMax, 1.0);

	// Calculating the number of grids in the Occupancy Volume
	QVector3D bound = urc - llc;
	int Nx = int((bound.x() / p_radius) + 1);
	int Ny = int((bound.y() / p_radius) + 1);
	occVolumeSizeX = Nx;
	occVolumeSizeY = Ny;
	float dx = bound.x() / (Nx - 1);
	float dy = bound.y() / (Ny - 1);
	size_t size = Nx*Ny;
	occupancyVolume.clear();
	occupancyVolume.resize(size);


	//Preparing to populate list of particles 
#pragma omp parallel for
	for (int i = 0; i < occupancyVolume.size(); i++)
	{
		std::vector<size_t> idList; // the scope of this maybe causing issues. May need a better method
		occupancyVolume[i] = idList;
	}

//Populate IDs 
	std::map<int, Particle>::iterator it = pList2.begin();
//#pragma omp parallel for	
	for (std::pair<int, Particle> element : pList2)
	{
		Particle p = element.second;
		QVector3D pPos = p.vertex().position();
		QVector3D temp = pPos - llc;
		if (temp.x() > 0 && temp.y() > 0)
		{
			if (temp.x() < bound.x() && temp.y() < bound.y())
			{
				int ix = int(temp.x() / dx);
				int iy = int(temp.y() / dy);

				int index = ix + Nx*iy;
				occupancyVolume[index].push_back(p.pId());
			}
		}
	}

	//qDebug() << occupancyVolume[0].size() << endl;
}

void Solver::generatevlist()
{
	s_vList.clear();
	for (std::pair<int, Particle> element : pList2)
	{
		int id = element.first;
		std::map<int, Particle>::iterator itA;
		itA = pList2.find(id);
		if (itA != pList2.end())
		{
			s_vList.push_back(itA->second.vertex());
		}
	}
	qDebug() << s_vList.size() << qPrintable(" Particles Detected");
}

std::vector<size_t> Solver::combinedInfList(int i, int j)
{
	int i0j0 = (i - 1) + occVolumeSizeX*(j - 1);
	int i0j = (i - 1) + occVolumeSizeX*(j + 0);
	int i0j1 = (i - 1) + occVolumeSizeX*(j + 1);

	int ij0 = (i + 0) + occVolumeSizeX*(j - 1);
	int ij = (i + 0) + occVolumeSizeX*(j + 0);
	int ij1 = (i + 0) + occVolumeSizeX*(j + 1);

	int i1j0 = (i + 1) + occVolumeSizeX*(j - 1);
	int i1j = (i + 1) + occVolumeSizeX*(j + 0);
	int i1j1 = (i + 1) + occVolumeSizeX*(j + 1);

	std::vector<size_t> listI0J0, listI0J, listI0J1;
	std::vector<size_t> listIJ0, listIJ1, listIJ;
	std::vector<size_t> listI1J0, listI1J, listI1J1;

	if (i > 0 && j > 0)
		listI0J0 = occupancyVolume[i0j0];
	if (i > 0)
		listI0J = occupancyVolume[i0j];
	if (i > 0 && j < (occVolumeSizeY - 1))
		listI0J1 = occupancyVolume[i0j1];
	if (j > 0)
		listIJ0 = occupancyVolume[ij0];
	listIJ = occupancyVolume[ij];
	if (j < (occVolumeSizeY - 1))
		listIJ1 = occupancyVolume[ij1];
	if (i < (occVolumeSizeX - 1) && j > 0)
		listI1J0 = occupancyVolume[i1j0];
	if (i < (occVolumeSizeX - 1))
		listI1J = occupancyVolume[i1j];
	if (i < (occVolumeSizeX - 1) && j < (occVolumeSizeY - 1))
		listI1J1 = occupancyVolume[i1j1];

	std::vector<size_t> combinedList;
	combinedList.reserve(listI0J0.size() + listI0J.size() + listI0J1.size() + listIJ0.size() + listIJ.size() + listIJ1.size() + listI1J0.size() + listI1J.size() + listI1J1.size());
	combinedList.insert(combinedList.end(), listI0J0.begin(), listI0J0.end());
	combinedList.insert(combinedList.end(), listI0J.begin(), listI0J.end());
	combinedList.insert(combinedList.end(), listI0J1.begin(), listI0J1.end());
	combinedList.insert(combinedList.end(), listIJ0.begin(), listIJ0.end());
	combinedList.insert(combinedList.end(), listIJ.begin(), listIJ.end());
	combinedList.insert(combinedList.end(), listIJ1.begin(), listIJ1.end());
	combinedList.insert(combinedList.end(), listI1J0.begin(), listI1J0.end());
	combinedList.insert(combinedList.end(), listI1J.begin(), listI1J.end());
	combinedList.insert(combinedList.end(), listI1J1.begin(), listI1J1.end());

	return combinedList;
}

float Solver::calcWeight(QVector3D posA, QVector3D posB)
{
	float r = (posA - posB).length();
	if (r == 0 || p_radius == 0)
		return 0;
	if ((r / p_radius) >= 1.0)
		return 0;
	float w = pow((1 - (r / p_radius)), 3) * 10 / (M_PI*p_radius*p_radius);
	return w;
}

QVector3D Solver::calcGradWeight(QVector3D posA, QVector3D posB)
{
	QVector3D r = posA - posB;
	float rMag = r.length();
	if( rMag == 0 || p_radius == 0)
		return QVector3D(0, 0, 0);
	if ((rMag / p_radius) >= 1.0)
		return QVector3D(0, 0, 0);
	QVector3D gW = -pow((1 - (rMag / p_radius)), 2) * (30 / (M_PI*pow(p_radius, 3))) * r.normalized();
	return gW;
}



QVector3D colorBlend(QVector3D vec, float minMag, float maxMag)
{
	float mag = fabs(vec.length());
	if (mag < minMag)
		mag = minMag;
	if (mag > maxMag)
		mag = maxMag;
	
	float t = (mag - minMag) / (maxMag - minMag);
	QVector3D res = QVector3D(0.0, 0.2, 1.0)*(1.0-t) + (t)*QVector3D(1.0, 0.2, 0.0);
	return res;
}



void Solver::calcDensity()
{
	for (int j = 0; j < occVolumeSizeY; j++)
	{
		for (int i = 0; i < occVolumeSizeX; i++)
		{
			int ij = i + occVolumeSizeX*j;
			std::vector<size_t>gridParticles = occupancyVolume[ij];
			std::vector<size_t> infList = combinedInfList(i, j);
			std::map<int, Particle>::iterator itA;
			std::map<int, Particle>::iterator itB;
			for(size_t id : gridParticles)
			{ 
				itA = pList2.find(id);
				if (itA != pList2.end())
				{
					Particle pA = itA->second;
					float density = 0;
					for (size_t index : infList)
					{
						itB = pList2.find(index);
						Particle pB;
						if (itB != pList2.end())
						{
							pB = itB->second;
							float wB = calcWeight(pA.position(), pB.position());
							density += pB.mass()*wB;
						}
					}
					itA->second.setdensity(density);
				}
			}
		}
	}
}

float Solver::calcPressure(Particle pA)
{
	float pres = pres_const*(pow((pA.density() / init_density), pres_gamma) - 1);
	return pres;
}


float Solver::calcVisc(Particle pA, Particle pB)
{
	float B = (visc_const * p_radius) / (pA.density() + pB.density());
	QVector3D rAB = pA.position() - pB.position();
	QVector3D vAB = pA.velocity() - pB.velocity();

	float vr = QVector3D::dotProduct(rAB, vAB);
	if (vr >= 0)
		return 0;
	float visc = -B*(vr / (rAB.lengthSquared() + (visc_e * p_radius*p_radius)));
	return visc;
}

QVector3D Solver::calcForces(Particle pA, std::vector<size_t>infList)
{
	return QVector3D(0, 0, 0);
}

void Solver::setForces()
{
	for (int j = 0; j < occVolumeSizeY; j++)
	{
#pragma omp parallel for 
		for (int i = 0; i < occVolumeSizeX; i++)
		{
			int ij = i + occVolumeSizeX*j;
			std::vector<size_t> infList = combinedInfList(i, j);
			std::map<int, Particle>::iterator itA;
			std::vector<size_t>gridParticles = occupancyVolume[ij];
			for(size_t partID: gridParticles)
			{
				itA = pList2.find(partID);
				Particle pA;
				if (itA != pList2.end())
				{
					pA = itA->second;
					
					QVector3D accel(0, 0, 0);
					float presA = calcPressure(pA);
					itA->second.setpressure(presA);
					std::map<int, Particle>::iterator itB;
					for (size_t id : infList)
					{
						itB = pList2.find(id);
						Particle pB;
						if (itB != pList2.end()) {
							pB = itB->second;
							float presB = calcPressure(pB);
							QVector3D gradW = calcGradWeight(pA.vertex().position(), pB.vertex().position());
							float presAB = (presA / (itA->second.density()*itA->second.density())) + (presB / (itB->second.density()*itB->second.density()));
							float viscAB = calcVisc(pA, pB);
							accel += -pB.mass()*(presAB+viscAB)*gradW;
						}
					}
					accel += force_g;
					itA->second.setacceleration(accel);
					//qDebug() << pA.acceleration() << qPrintable(" Acceleration in forces ");
				}
			}
		}
	}
}

void Solver::advectForward(float dtime)
{
	for (std::pair<int, Particle> element : pList2)
	{
		int id = element.first;
		std::map<int, Particle>::iterator itA;
		itA = pList2.find(id);
		if (itA != pList2.end())
		{
			
			QVector3D pos = itA->second.position();
			pos.setX(pos.x() + itA->second.velocity().x()*dtime);
			pos.setY(pos.y() + itA->second.velocity().y()*dtime);
			itA->second.setposition(pos);
		}
	}

}

void Solver::advectVelocity(float dtime)
{
	for (std::pair<int, Particle> element : pList2)
	{
		int id = element.first;
		std::map<int, Particle>::iterator itA;
		itA = pList2.find(id);
		if (itA != pList2.end())
		{
			
			QVector3D vel = itA->second.velocity();
			vel.setX(vel.x() + itA->second.acceleration().x()*dtime*(1 / element.second.mass()));
			vel.setY(vel.y() + itA->second.acceleration().y()*dtime*(1 / element.second.mass()));
			itA->second.setvelocity(vel);
			QVector3D col = colorBlend(vel, 0, 10);
			itA->second.setcolor(col);
		}
	}
}

void Solver::boundaryConditions()
{
	QVector3D llc = QVector3D(-1.0, -1.0, 1.0);
	QVector3D urc = QVector3D(1.0, 1.0, 1.0);

	for (std::pair<int, Particle> element : pList2)
	{
		int id = element.first;
		std::map<int, Particle>::iterator itA;
		itA = pList2.find(id);
		if (itA != pList2.end())
		{
			QVector3D pos = itA->second.position();
			QVector3D vel = itA->second.velocity();
			QVector3D updatedPos = pos;
			if (pos.x() < llc.x())
			{
				updatedPos.setX( llc.x() + wall_sticky*(llc.x() - pos.x()));
				vel.setX(-vel.x()*wall_sticky);
			}
			if (pos.y() < llc.y())
			{
				updatedPos.setY(llc.y() + wall_sticky*(llc.y() - pos.y()));
				vel.setY(-vel.y()*wall_sticky);
			}
			if (pos.x() > urc.x())
			{
				updatedPos.setX(urc.x() + wall_sticky*(urc.x() - pos.x()));
				vel.setX(-vel.x()*wall_sticky);
			}
			if (pos.y() > urc.y())
			{
				updatedPos.setY(urc.y() + wall_sticky*(urc.y() - pos.y()));
				vel.setY(-vel.y()*wall_sticky);
			}
			itA->second.setposition(updatedPos);
			itA->second.setvelocity(vel);
		}
	}
}

void Solver::sixthAdvection(float dtime)
{
	float a = 1 / (4 - powf(4, (1 / 3.0)));
	float b = 1 - 4 * a;

	advectParticles(a*dtime);
	advectParticles(a*dtime);
	advectParticles(b*dtime);
	advectParticles(a*dtime);
	advectParticles(a*dtime);
}

void Solver::advectParticles(float dtime)
{

	//qDebug() << pList2[0].position() << "Position Before";
	advectForward(dtime*0.5);
	boundaryConditions();
	//qDebug() << pList2[0].position() << "Position After";
	generateOccupancyVolume();
	calcDensity();
	setForces();
	//qDebug() << pList2[0].velocity() << "Velocity Before";
	advectVelocity(dtime);
	//qDebug() << pList2[0].velocity() << "Velocity After";
	advectForward(dtime*0.5);
	boundaryConditions();
	
	//qDebug() << qPrintable("Advected Leap Frog");
	//generatevlist();

}

void Solver::simulate()
{
	//	qDebug() << qPrintable("Simulating");
	
	if (advectionType == LEAPFROG)
	{
		advectParticles(s_dt);
	}
	if (advectionType == SIXTH)
	{
		sixthAdvection(s_dt);
	}
}

