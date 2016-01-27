/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements scene_object.h

***********************************************************/

#include <cmath>
#include <iostream>
#include "scene_object.h"


bool UnitSquare::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSquare, which is
	// defined on the xy-plane, with vertices (0.5, 0.5, 0), 
	// (-0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
	// (0, 0, 1).
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.

	Point3D p0(0, 0, 0); 
	Vector3D n(0, 0, 1); 

	Vector3D dir = worldToModel * ray.dir;
	Point3D O = worldToModel * ray.origin;
	Point3D inter;
	double t;

	if ((p0 - O).dot(n)!= 0 &&  dir.dot(n) != 0){
		t = (p0 - O).dot(n)/dir.dot(n);
		inter = O + t * dir;

		if (((ray.intersection.none || t <= ray.intersection.t_value) && t >= 0.00000001)
			&&(inter[0] < 0.5 && inter[0] > -0.5&& inter[1] < 0.5 && inter[1] > -0.5)){

			ray.intersection.point = modelToWorld * inter;
			ray.intersection.normal = transNorm(worldToModel, n);
			ray.intersection.normal.normalize();
			ray.intersection.obj_index = this ->index;
			ray.intersection.t_value = t;
			ray.intersection.none = false;
			ray.intersection.shape = "squ";
			return true;
		}
	}
	return false;
}


bool UnitSphere::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSphere, which is centred 
	// on the origin.  
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.
	    Ray3D modelRay;
		modelRay.dir = worldToModel * ray.dir;
		modelRay.origin = worldToModel * ray.origin;
		Point3D origin(0, 0, 0);
		double discriminant;

		Vector3D po = modelRay.origin - origin;
		double a = modelRay.dir.dot(modelRay.dir);
		double b = 2 * po.dot(modelRay.dir);
		double c = po.dot(po) - 1;
		double t;
		discriminant = b * b -4 * a * c;

		if (discriminant < 0) {
			return false;
		}
		else if (discriminant == 0 || discriminant <= 0.000001){
			t = (-b + sqrt(discriminant))/ (2 * a);
		}else{
			double t1 = (-b + sqrt(discriminant)) / (2 * a);
			double t2 = (-b - sqrt(discriminant)) / (2 * a);
			t = t2;  // t2 is always the smaller one
			if ( t2 < -0.000001 && t1 > 0.000001) {
				t = t1;
			}else if ( t1 < -0.0000001 ) {
				return false;
			}
		}
		Point3D inter = modelRay.origin + t * modelRay.dir;
        Vector3D normal = inter - origin;
		
        if (!ray.intersection.none && t > ray.intersection.t_value || t < -0.0000001){
			return false;
        }

        ray.intersection.none = false;	
        ray.intersection.point = modelToWorld * inter;
		ray.intersection.normal = transNorm(worldToModel,normal);
        ray.intersection.normal.normalize();
        ray.intersection.t_value = t;
		ray.intersection.obj_index = this -> index;
        ray.intersection.shape = "sph";
		return true;

}


bool UnitCylinder::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
                const Matrix4x4& modelToWorld ) {
      
		Ray3D modelRay;
		Point3D origin(0,0,0),intersection;

		modelRay.dir = worldToModel * ray.dir;
		modelRay.origin = worldToModel * ray.origin;

		double x0=modelRay.origin[0], y0=modelRay.origin[1], z0=modelRay.origin[2];
		double dx=modelRay.dir[0], dy=modelRay.dir[1], dz=modelRay.dir[2];
		
        double t;

        double a = dx*dx + dy*dy;
		double b = 2*(x0*dx + y0*dy);
		double c = x0*x0 + y0*y0 -1;
		double d = b*b - 4*a*c;
		
		Point3D inter;
		Vector3D n;
		  
		double d1 = (-0.5-z0)/dz;
		double d2 = (0.5-z0)/dz;
		double t1 = (-b + sqrt(abs(d)))/(2*a);
		double t2 = (-b - sqrt(abs(d)))/(2*a);

		Point3D D1,D2,T1,T2;
		D1= modelRay.origin + d1*modelRay.dir;
		D2= modelRay.origin + d2*modelRay.dir;
		T1= modelRay.origin + t1*modelRay.dir;
		T2= modelRay.origin + t2*modelRay.dir;

		bool b1 = D1[0] * D1[0] + D1[1] * D1[1] <= 1;
		bool b2 = D2[0] * D2[0] + D2[1] * D2[1] <= 1;
		
	//***************** case 1 *****************
		if (b1 && b2){
			if (d1 > 0 && d2 > 0){
				if (d1 < d2){
					t = d1;
					n = Vector3D(0,0,-1);
				}else{
					t = d2;
					n = Vector3D(0,0,1);
				}
			}else if (d1 > 0 && d2 < 0){
				t = d1;
				n = Vector3D(0,0,-1);
			}else if (d1 < 0 && d2 > 0){
				t = d2;
				n = Vector3D(0,0,1);
			}else{
				return false;
			}
			
			
	//***************** case 2 *****************
		}else if(!b1 && !b2){
			if(d<0 || a==0)
				return false;
			else if(t2>0){                                                                      
				Point3D p = modelRay.origin + t2*modelRay.dir;
				if(abs(p[2]) - 0.5 < -0.0000001){
					t=t2;
					n=Vector3D(p[0],p[1],0);
				}else
					return false;

			}else if(t1>0 && t2<0){
				Point3D p = modelRay.origin + t1*modelRay.dir;
				if(abs(p[2]) - 0.5 < -0.0000001){
					t=t1;
					n=Vector3D(p[0],p[1],0);
				}else
					return false;
			}else 
				return false;

	//***************** case 3 *****************
		}else if (b1 && !b2){
			if (d < 0) return false;
			else{
				Point3D p = modelRay.origin + t1*modelRay.dir;
				if (abs(p[2]) - 0.5 < -0.0000001){
					if (d1 < t1){
						t = d1;
						n = Vector3D(0,0,-1);
					}else{
						t = t1;
						n=Vector3D(p[0],p[1],0);
					}
				}else{
					p = modelRay.origin + t2*modelRay.dir;
					if (d1 < t2){
						t = d1;
						n = Vector3D(0,0,-1);
					}else{
						t = t2;
						n=Vector3D(p[0],p[1],0);
					}
				}
			}
	//***************** case 4 *****************
		}else if (!b1 && b2){
			if (d < 0) return false;
			else{
				Point3D p = modelRay.origin + t1*modelRay.dir;
				if (abs(p[2]) - 0.5 < -0.0000001){
					if (d2 < t1){
						t = d2;
						n = Vector3D(0,0,1);
					}else{
						t = t1;
						n=Vector3D(p[0],p[1],0);
					}
				}else{
					p = modelRay.origin + t2*modelRay.dir;
					if (d2 < t2){
						t = d2;
						n = Vector3D(0,0,-1);
					}else{
						t = t2;
						n=Vector3D(p[0],p[1],0);
					}
				}
			}
		}

		inter = modelRay.origin + t * modelRay.dir;

		if ((!ray.intersection.none && t > ray.intersection.t_value) || t < -0.00000001)
			return false;
        
        ray.intersection.none = false;
        ray.intersection.point = modelToWorld * inter;
		ray.intersection.normal = transNorm(worldToModel,n);
        ray.intersection.normal.normalize();
        ray.intersection.t_value = t;
		ray.intersection.obj_index = this -> index;
		ray.intersection.shape = "cyl";
		return true;
}

bool UnitCone::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
                const Matrix4x4& modelToWorld ) {
        //NOTE: The intersection for a unit cylinder. A unit cylinder here means that
        //      the top and the base are unit circles centered at (0,0,z), where z is
        //      -0.5 and 0.5, respectively. And the height of the cylinder is 1

		Ray3D modelRay;
		Point3D origin(0,0,0),intersection;

		modelRay.dir = worldToModel * ray.dir;
		modelRay.origin = worldToModel * ray.origin;

		double x0=modelRay.origin[0], y0=modelRay.origin[1], z0=modelRay.origin[2];
		double dx=modelRay.dir[0], dy=modelRay.dir[1], dz=modelRay.dir[2];
		
        double t;

        double a = dx*dx + dy*dy - dz*dz;
		double b = 2*(x0*dx + y0*dy - z0*dz);
		double c = x0*x0 + y0*y0 -z0*z0;
		double dis = b*b-4*a*c;

		Point3D inter;
		Vector3D normal;
		  
		double d = (1 - z0)/dz;
		double t1 = (-b + sqrt(dis))/(2*a);
		double t2 = (-b - sqrt(dis))/(2*a);

		Point3D D,T1,T2;
		D = modelRay.origin + d*modelRay.dir;
		T1= modelRay.origin + t1*modelRay.dir;
		T2= modelRay.origin + t2*modelRay.dir;

		bool hitbase = D[0] * D[0] + D[1] * D[1] <= 1;

		if(hitbase){
			if(dis>=0){
				if(T1[2]>0 && T1[2]<1){
					if(d>0 && t1>0){
						if(t1>d)
							t = d;
						 else
							t = t1;
					}else if(d>0 && t1<0){
						t = d;
					}else if(d<0 && t1>0){
						t = t1;
					}else{return false;}

				}else{
					if(d>0 && t2>0){
						if(t2>d)
							t = d;
						 else
							t = t2;
					}else if(d>0 && t2<0){
						t = d;
					}else if(d<0 && t2>0){
						t = t2;
					}else{return false;}
				}
			}else{return false;}

		}else{
			if(dis>=0){
				if(t2>0){
					if(T2[2]>0 && T2[2]<1){
						t=t2;
					}else if(T1[2]>0 && T1[2]<1){ 
						t=t1;
					}else{return false;}
				}else if( t1>0 && t2<0){
					if(T1[2]>0 && T1[2]<1){ 
						t=t1;
					}else{return false;}
				}else{return false;}
			}else{return false;}
		}


		inter = modelRay.origin + t * modelRay.dir;

		if(t==d)
			normal = Vector3D(0,0,1);
		else{
			normal =  Vector3D(inter[0], inter[1], sqrt(inter[0]*inter[0]+ inter[1]*inter[1]));
		}

		if (!ray.intersection.none && t > ray.intersection.t_value)
			return false;
        
        ray.intersection.none = false;
        ray.intersection.point = modelToWorld * inter;
		ray.intersection.normal = transNorm(worldToModel,normal);
        ray.intersection.normal.normalize();
        ray.intersection.t_value = t;
		ray.intersection.obj_index = this -> index;
        ray.intersection.shape = "con";
		return true;
}



