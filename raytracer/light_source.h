/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		light source classes

***********************************************************/

#include "util.h"

// Base class for a light source.  You could define different types
// of lights here, but point light is sufficient for most scenes you
// might want to render.  Different light sources shade the ray 
// differently.
class LightSource {
public:
	virtual void shade( Ray3D& ) = 0;
	virtual Point3D get_position() const = 0; 
	virtual char get_type() const = 0;
	virtual Vector3D get_u() const = 0;
	virtual Vector3D get_v() const = 0;
	virtual double get_lx() const = 0;
	virtual double get_ly() const = 0;
};

// A point light is defined by its position in world space and its
// colour.
class PointLight : public LightSource {
public:
	PointLight( Point3D pos, Colour col ) : _pos(pos), _col_ambient(col), 
	_col_diffuse(col), _col_specular(col) {}
	PointLight( Point3D pos, Colour ambient, Colour diffuse, Colour specular ) 
	: _pos(pos), _col_ambient(ambient), _col_diffuse(diffuse), 
	_col_specular(specular) {}
	void shade( Ray3D& ray );
	Point3D get_position() const { return _pos; }
	char get_type() const { return 'p'; }
	Vector3D get_u() const { return Vector3D (1, 0, 0); }
	Vector3D get_v() const { return Vector3D (0, 1, 0); }
	double get_lx() const { return 0.0001; }
	double get_ly() const { return 0.0001; }

private:
	Point3D _pos;
	Colour _col_ambient;
	Colour _col_diffuse; 
	Colour _col_specular; 
};

class AreaLight : public LightSource {
public:
	AreaLight( Point3D pos, Colour col, Vector3D u0, Vector3D v0, double x, double y) :  sc(pos),_col_ambient(col),
		_col_diffuse(col), _col_specular(col), u(u0), v(v0), lx(x), ly(y){};
	void shade( Ray3D& ray );
	char get_type() const { return 'a'; }
	Point3D get_position() const { return sc; }
	Vector3D get_u() const { return u; }
	Vector3D get_v() const { return v; }
	double get_lx() const { return lx; }
	double get_ly() const { return ly; }

private:
	Point3D sc;
	Vector3D u,v;
	double lx,ly;

	Colour _col_ambient;
	Colour _col_diffuse; 
	Colour _col_specular; 
};
