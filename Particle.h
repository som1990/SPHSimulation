#ifndef Particle_h__
#define Particle_h__

#include "Vertex.h"
#include <QVector3D>

class Particle
{
public:
	Q_DECL_CONSTEXPR Particle();
	Q_DECL_CONSTEXPR Particle(Vertex &v, float &h );
	Q_DECL_CONSTEXPR Particle(Vertex & v, float &h, float& m);
	Q_DECL_CONSTEXPR const QVector3D& velocity() const;
	Q_DECL_CONSTEXPR const QVector3D& oldvelocity() const;
	Q_DECL_CONSTEXPR const QVector3D& acceleration() const;
	
	Q_DECL_CONSTEXPR const QVector3D& position() const;
	Q_DECL_CONSTEXPR const QVector3D& color() const;
	Q_DECL_CONSTEXPR const Vertex& vertex() const;
	Q_DECL_CONSTEXPR const float& mass() const;
	Q_DECL_CONSTEXPR const float& pressure() const;
	Q_DECL_CONSTEXPR const float& density() const;
	Q_DECL_CONSTEXPR const float& radius() const;
	Q_DECL_CONSTEXPR const int& pId() const;
	

	void setvertex(Vertex v);
	void setvelocity(QVector3D vel);
	void setoldvelocity(QVector3D oldvel);
	void setacceleration(QVector3D accel);
	
	void setposition(QVector3D& pos);
	void setcolor(QVector3D& col);
	void setmass(float& m);
	void setpressure(float& p);
	void setdensity(float& rho);
	void setradius(float& h);
	void setid(int &ind);

private:
	Vertex p_vert;
	QVector3D p_vel;
	QVector3D p_oldVel;
	QVector3D p_accel;
	int id;
	float p_mass = 1.0;
	float p_pressure = 0.0 ;
	float p_density = 0.0;
	float p_radius;
};

Q_DECLARE_TYPEINFO(Particle, Q_MOVABLE_TYPE);


Q_DECL_CONSTEXPR inline Particle::Particle() { }
Q_DECL_CONSTEXPR inline Particle::Particle(Vertex &v, float &h) : p_vert(v), p_radius(h) { }
Q_DECL_CONSTEXPR inline Particle::Particle(Vertex &v, float &h, float &m) : p_vert(v), p_radius(h), p_mass(m) { }

Q_DECL_CONSTEXPR inline const QVector3D& Particle::velocity() const { return p_vel; }
Q_DECL_CONSTEXPR inline const QVector3D& Particle::oldvelocity() const { return p_oldVel; }
Q_DECL_CONSTEXPR inline const QVector3D& Particle::acceleration() const { return p_accel; }



#endif // Particle_h__