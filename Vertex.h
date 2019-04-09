#ifndef Vertex_h__
#define Vertex_h__

#include <QVector3D>

class Vertex
{
public:
	Q_DECL_CONSTEXPR Vertex();
	Q_DECL_CONSTEXPR explicit Vertex(const QVector3D &position);
	Q_DECL_CONSTEXPR Vertex(const QVector3D &position, const QVector3D &color);

	Q_DECL_CONSTEXPR const QVector3D& position() const;
	Q_DECL_CONSTEXPR const QVector3D& color() const;
	void setPosition(const QVector3D position);
	void setColor(const QVector3D color);

	static const int PositionTupleSize = 3;
	static const int ColorTupleSize = 3;
	static Q_DECL_CONSTEXPR int positionOffset();
	static Q_DECL_CONSTEXPR int colorOffset();
	static Q_DECL_CONSTEXPR int stride();

private: 
	QVector3D m_position;
	QVector3D m_color;

};
// Note: Q_MOVABLE_TYPE means it can be memcpy'd.
Q_DECLARE_TYPEINFO(Vertex, Q_MOVABLE_TYPE);

Q_DECL_CONSTEXPR inline Vertex::Vertex() {}
Q_DECL_CONSTEXPR inline Vertex::Vertex(const QVector3D &position) : m_position(position) {}
Q_DECL_CONSTEXPR inline Vertex::Vertex(const QVector3D &position, const QVector3D &color) : m_position(position), m_color(color) {}

Q_DECL_CONSTEXPR inline const QVector3D& Vertex::position() const { return m_position; }
Q_DECL_CONSTEXPR inline const QVector3D& Vertex::color() const { return m_color; }
inline void Vertex::setPosition(const QVector3D position) {m_position = position;}
inline void Vertex::setColor(const QVector3D color) { m_color = color; }

Q_DECL_CONSTEXPR inline int Vertex::positionOffset() { return offsetof(Vertex, m_position); }
Q_DECL_CONSTEXPR inline int Vertex::colorOffset() { return offsetof(Vertex, m_color); }
Q_DECL_CONSTEXPR inline int Vertex::stride() { return sizeof(Vertex); }


#endif // Vertex_h__