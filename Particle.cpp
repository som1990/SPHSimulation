#include "Particle.h"


Q_DECL_CONSTEXPR const QVector3D& Particle::position() const { return p_vert.position(); }

Q_DECL_CONSTEXPR const QVector3D& Particle::color() const { return p_vert.color(); }

Q_DECL_CONSTEXPR const Vertex& Particle::vertex() const{ return p_vert; }

Q_DECL_CONSTEXPR const float& Particle::mass() const { return p_mass; }

Q_DECL_CONSTEXPR const float& Particle::pressure() const { return p_pressure; }

Q_DECL_CONSTEXPR const float& Particle::density() const { return p_density; }

Q_DECL_CONSTEXPR const float& Particle::radius() const { return p_radius; }

Q_DECL_CONSTEXPR const int& Particle::pId() const { return id; }

void Particle::setvertex(Vertex v) { p_vert = v; }

void Particle::setvelocity(QVector3D vel) { p_vel = vel; }

void Particle::setoldvelocity(QVector3D oldvel){ p_oldVel = oldvel;}

void Particle::setacceleration(QVector3D accel) { p_accel = accel; }

void Particle::setposition(QVector3D& pos) { p_vert.setPosition(pos); }

void Particle::setcolor(QVector3D& col){ p_vert.setColor(col); }

void Particle::setmass(float& m) { p_mass = m; }

void Particle::setpressure(float& p) { p_pressure = p; }

void Particle::setdensity(float& rho) { p_density = rho; }

void Particle::setradius(float& h) { p_radius = h; }

void Particle::setid(int &ind) { id = ind; }

