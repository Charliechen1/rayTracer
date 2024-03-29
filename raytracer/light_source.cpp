﻿/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include "light_source.h"

void PointLight::shade( Ray3D& ray ) {
	// TODO: implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.  
	Colour amb, dif, spe;

	Vector3D n = ray.intersection.normal;
	n.normalize();

	Vector3D l = _pos - ray.intersection.point;
	l.normalize();

	Vector3D v = -ray.dir;
	v.normalize();

	Vector3D r = 2 * l.dot(n) * n - l;
	r.normalize();
	
	double shiness = ray.intersection.mat -> specular_exp;
	double cos_1 = n.dot(l);
	cos_1 = (cos_1 >= 0.0? cos_1 : 0.0);// only contribution, no penalty

	double cos_2 = r.dot(v);
	cos_2 = (cos_2 >= 0.0? cos_2 : 0.0);// only contribution, no penalty

	amb = ray.intersection.mat ->ambient * _col_ambient;

	dif = cos_1 * ray.intersection.mat ->diffuse * _col_diffuse;

	spe = pow(cos_2, shiness) * ray.intersection.mat -> specular * _col_diffuse;
	
	if(!ray.isShadowed){
		ray.col = ray.col + amb + dif + spe;
	}else{
		ray.col = ray.col + amb; 
	}
	ray.col.clamp();
}

void AreaLight::shade( Ray3D& ray ) {
	// TODO: implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.  
	Colour amb, dif, spe;

	Vector3D n = ray.intersection.normal;
	n.normalize();

	Vector3D l = sc - ray.intersection.point;
	l.normalize();

	Vector3D v = -ray.dir;
	v.normalize();

	Vector3D r = 2 * l.dot(n) * n - l;
	r.normalize();
	
	double shiness = ray.intersection.mat -> specular_exp;
	double cos_1 = n.dot(l);
	cos_1 = (cos_1 >= 0.0? cos_1 : 0.0);// only contribution, no penalty

	double cos_2 = r.dot(v);
	cos_2 = (cos_2 >= 0.0? cos_2 : 0.0);// only contribution, no penalty

	amb = ray.intersection.mat ->ambient * _col_ambient;

	dif = cos_1 * ray.intersection.mat ->diffuse * _col_diffuse;

	spe = pow(cos_2, shiness) * ray.intersection.mat -> specular * _col_diffuse;
	
//	if(ray.ShadowHardness < 0.00000001){
		ray.col = ray.col + 0.01f*(amb + dif + spe);
//	}else{
//		ray.col = ray.col + amb + (1-ray.ShadowHardness) * (dif + spe);
//	}

	ray.col.clamp();	
}