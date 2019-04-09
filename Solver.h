/* Solver.h
// Author: Soumitra Goswami
// Date: 05/01/2018
// Description: This solver simulates an Lagrangian Smooth Particle Hydrodynamics (SPH).
//              To understand the implementation in detail Refer to the notes by Bridson (https://www.cs.ubc.ca/~rbridson/fluidsimulation/fluids_notes.pdf)  
//              Since Euler Integration is very unstable, I incorporate two other stable advection
//              schemes:
//               - Leap Frog Advection 
//               - Sixth Advection
//              These schemes give stable performance(especially sixth) even >10000 particles
//              Allows for viscocity, pressure and other forces.
//              Optimized with an occupancy grid so only searches within a certain grid region for points.
//
// 
*/

#ifndef Solver_h__
#define Solver_h__
#include "Particle.h"
#include <vector>
#include <QVector3D>
#include <map>

class Solver
{
public:
	Q_DECL_CONSTEXPR Solver();
	Q_DECL_CONSTEXPR Solver(size_t& pNum, float& pMass, float& pRad, QVector3D& fGravity, float& presConst, float& presGamma, float& iDensity, float& vConst, float& vE, float& wSticky, float& t );
	~Solver();
	void generateParticleList();
	void simulate(); // Entry point of our class
	enum Advection {LEAPFROG, SIXTH }; // allows for two advection types
	Advection advectionType = LEAPFROG;

	
//Gets set by "sphsimulation.cpp" . Doesn't really needs to be SLOTs. 
public slots:
	void settimestep(float &t);
	
	void setparticleCount(size_t &pCount);
	void setmass(float &m);
	void setradius(float &r);
	
	void setgravity(float &gravY);
	void setpresConst(float &pConst);
	void setpresGamma(float &gamma);
	void setdensity(float &rho);

	void setviscConst(float &vConst);
	void setviscE(float &vE);
	void setWallSticky(float wSticky);

//Returns for the variables 
public:
	Q_DECL_CONSTEXPR const float& dt() const;
	Q_DECL_CONSTEXPR const size_t& nParticles() const;
	Q_DECL_CONSTEXPR const float& mass() const;
	Q_DECL_CONSTEXPR const float& radius() const;

	Q_DECL_CONSTEXPR const QVector3D& gravity() const;
	Q_DECL_CONSTEXPR const float& presConst() const;
	Q_DECL_CONSTEXPR const float& presGamma() const;
	Q_DECL_CONSTEXPR const float& initDensity() const;
	Q_DECL_CONSTEXPR const float& viscConst() const;
	Q_DECL_CONSTEXPR const float& viscE() const;
	Q_DECL_CONSTEXPR const float& wallSticky() const;
	const std::vector<Vertex> vList() const;
	std::map<int, Particle> &partList(); //To Do: I need to fix this. This works for a certain low number of particles but fails after a while. Need to allocate pList2 to the heap so it can store more than a few 1000 particles 

private:
	//Not Utilized
	QVector3D calcForces(Particle pA,std::vector<size_t>infList);

	//void setForces()
	//Calculates and sets inter-particulate, pressure and viscocity forces 
	void setForces();

	/* void generateOccupancyVolume()
	// Description: creates the occupancy volume grid. Divides region into grids and stores ids of particles within the each grid.
	*/
	void generateOccupancyVolume();

	/* void generatevlist()
	// Description generates a vertex list based on particle list to be used by opengl to render points.
	*/
	void generatevlist();
	
	/* std::vector<size_t> combinedInfList(int, int )
	// Description: Given grid position. It generates a list of particles influencing the given occupancy grid. This way it limits the number of distance checks for interparticulate forces.
	// Input: @int i - horizontal occupancy grid position
	//        @int j - vertical occupancy grid position
	// Return std::vector<size_t> - returns a list of all the particle ids influencing the particle in the neighboring grids
	*/
	std::vector<size_t> combinedInfList(int i, int j);//Todo: At the moment copying the vector over to the other. Need to change this into an address

	//To Do: This needs to be more efficient. Right now its calculating the density twice since the particles have to calculate distance with eachother twice. 
	/*void calcDensity()
	// Description: Calculates the density of the particle based on number of particles surrounding.
	*/
	void calcDensity();

	
	/* float calcWeight(QVector3D, QVector3D)
	// Description: calculates the attraction weight based on distance between the two particles.
	// Inputs: @QVector3D posA - position of particle A
	//         @QVector3D posB - position of particle B
	// Return: Returns the weight between two particles based on position
	*/
	float calcWeight(QVector3D posA, QVector3D posB);

	/* QVector3D calcGradWeight(QVector3D, QVector3D)
	// Description: calculates the gradient weight based on distance between the two particles.
	// Inputs: @QVector3D posA - position of particle A
	//         @QVector3D posB - position of particle B
	// Return: Returns the gradient weight in each dimension
	*/
	QVector3D calcGradWeight(QVector3D posA, QVector3D posB);

	/* void advectParticles(float)
	// Description: Advects the particles based on the leap frog advection.
	//              Leap frog(LF) advection takes place in these steps:
	//              1) Advects forward half step
	//                 f(x+0.5*dx) = x + u*0.5*dt
	//              2) Calculate the forces/acceleration (a)
	//              3) Advect velocity full step
	//                 f(v+ dv) = u + a*dt
	//              4) Advect forward half step
	//                 f(x+0.5*dx) = x + v*0.5*dt
	//              In our case we do boundary checks twice after each forward advection step
	// Inputs: @float dTime: The default timestep
	*/
	void advectParticles(float dtime);


	/* float calcPressure(Particle pA)
	// Description: Calculates pressure step of Navier Stokes equation
	// Returns : @float - returns the calculated pressure
	*/
	float calcPressure(Particle pA);


	/* float calcVisc(Particle , Particle) 
	// Description: calculates the particulate viscocity between two particles.
	// Inputs: @Particle pA - particle A
	//         @Particle pB - particle B
	// Return : @float - calculated viscocity between the two particles
	*/
	
	float calcVisc(Particle pA, Particle pB);

	/*void advectForward(float)
	// Description: Does a forward advection of position in dtime step
	//              f(x+dx) = x + u*dtime
	// Inputs: @float dtime : The timestep
	*/
	void advectForward(float dtime);

	/* void advectVelocity(float)
	// Description: Does a forward advection of each dimension in dtime step
	//              f(u + du) = u + a*dtime
    // Inputs: @float dtime : The timestep
	*/
	void advectVelocity(float dtime);

	/* void boundaryConditions()
	// Description: finds the edges and calculates a reflected velocity based on the boundary
	//             It also incorporates wall stickiness so that some particles stick to the wall 
	//             and slide down to fake surface tension.
	*/
	void boundaryConditions();

	/* void sixthAdvection(float)
	// Description: Advects the particles using the Sixth advection scheme.
	//              Sixth advection scheme includes using leap frog five times with varying timesteps
	//              It's heavier since it uses leap frog 5 times but it's a lot more stable.
	*/
	void sixthAdvection(float dtime);

	
	//List Variables
	
	std::vector<std::vector<size_t>> occupancyVolume;
	//std::vector<Particle> pList;
	std::map<int,Particle> pList2; // Need to allocate to the heap without the program breaking.
	std::vector<Vertex> s_vList;
	
	
	float force_pres;
	float force_visc;
	int occVolumeSizeX;
	int occVolumeSizeY;

	///////////INPUT VARIABLES///////////////

	float s_dt = 0.005;
	//particle inputs
	size_t num_particles = 1000;
	float p_mass = 1.0;
	float p_radius = 0.15;

	//Gravity Input
	QVector3D force_g = QVector3D(0.0f,-9.8f,0.0f);

	//Pressure Force Input
	float pres_const = 0.05;
	float pres_gamma = 3;
	float init_density= 25.0;
	
	//Viscosity Input
	float visc_const = 3;
	float visc_e = 0.15 ;
	//Boundary Input
	float wall_sticky = 0.5;
};

#endif // Solver_h__