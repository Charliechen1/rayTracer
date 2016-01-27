/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		Implementations of functions in raytracer.h, 
		and the main function which specifies the 
		scene to be rendered.	

***********************************************************/


#include "raytracer.h"
#include "bmp_io.h"
#include <cmath>
#include <iostream>
#include <cstdlib>

bool is_out = true;

Raytracer::Raytracer() : _lightSource(NULL) {
	_root = new SceneDagNode();
}

Raytracer::~Raytracer() {
	delete _root;
}

SceneDagNode* Raytracer::addObject( SceneDagNode* parent, 
		SceneObject* obj, Material* mat ) {
	SceneDagNode* node = new SceneDagNode( obj, mat );
	node->parent = parent;
	node->next = NULL;
	node->child = NULL;
	
	// Add the object to the parent's child list, this means
	// whatever transformation applied to the parent will also
	// be applied to the child.
	if (parent->child == NULL) {
		parent->child = node;
	}
	else {
		parent = parent->child;
		while (parent->next != NULL) {
			parent = parent->next;
		}
		parent->next = node;
	}
	
	return node;;
}

LightListNode* Raytracer::addLightSource( LightSource* light ) {
	LightListNode* tmp = _lightSource;
	_lightSource = new LightListNode( light, tmp );
	return _lightSource;
}

void Raytracer::rotate( SceneDagNode* node, char axis, double angle ) {
	Matrix4x4 rotation;
	double toRadian = 2*M_PI/360.0;
	int i;
	
	for (i = 0; i < 2; i++) {
		switch(axis) {
			case 'x':
				rotation[0][0] = 1;
				rotation[1][1] = cos(angle*toRadian);
				rotation[1][2] = -sin(angle*toRadian);
				rotation[2][1] = sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'y':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][2] = sin(angle*toRadian);
				rotation[1][1] = 1;
				rotation[2][0] = -sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'z':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][1] = -sin(angle*toRadian);
				rotation[1][0] = sin(angle*toRadian);
				rotation[1][1] = cos(angle*toRadian);
				rotation[2][2] = 1;
				rotation[3][3] = 1;
			break;
		}
		if (i == 0) {
		    node->trans = node->trans*rotation; 	
			angle = -angle;
		} 
		else {
			node->invtrans = rotation*node->invtrans; 
		}	
	}
}

void Raytracer::translate( SceneDagNode* node, Vector3D trans ) {
	Matrix4x4 translation;
	
	translation[0][3] = trans[0];
	translation[1][3] = trans[1];
	translation[2][3] = trans[2];
	node->trans = node->trans*translation; 	
	translation[0][3] = -trans[0];
	translation[1][3] = -trans[1];
	translation[2][3] = -trans[2];
	node->invtrans = translation*node->invtrans; 
}

void Raytracer::scale( SceneDagNode* node, Point3D origin, double factor[3] ) {
	Matrix4x4 scale;
	
	scale[0][0] = factor[0];
	scale[0][3] = origin[0] - factor[0] * origin[0];
	scale[1][1] = factor[1];
	scale[1][3] = origin[1] - factor[1] * origin[1];
	scale[2][2] = factor[2];
	scale[2][3] = origin[2] - factor[2] * origin[2];
	node->trans = node->trans*scale; 	
	scale[0][0] = 1/factor[0];
	scale[0][3] = origin[0] - 1/factor[0] * origin[0];
	scale[1][1] = 1/factor[1];
	scale[1][3] = origin[1] - 1/factor[1] * origin[1];
	scale[2][2] = 1/factor[2];
	scale[2][3] = origin[2] - 1/factor[2] * origin[2];
	node->invtrans = scale*node->invtrans; 
}

Matrix4x4 Raytracer::initInvViewMatrix( Point3D eye, Vector3D view, 
		Vector3D up ) {
	Matrix4x4 mat; 
	Vector3D w;
	view.normalize();
	up = up - up.dot(view)*view;
	up.normalize();
	w = view.cross(up);

	mat[0][0] = w[0];
	mat[1][0] = w[1];
	mat[2][0] = w[2];
	mat[0][1] = up[0];
	mat[1][1] = up[1];
	mat[2][1] = up[2];
	mat[0][2] = -view[0];
	mat[1][2] = -view[1];
	mat[2][2] = -view[2];
	mat[0][3] = eye[0];
	mat[1][3] = eye[1];
	mat[2][3] = eye[2];

	return mat; 
}

void Raytracer::traverseScene( SceneDagNode* node, Ray3D& ray ) {
	SceneDagNode *childPtr;

	// Applies transformation of the current node to the global
	// transformation matrices.
	_modelToWorld = _modelToWorld*node->trans;
	_worldToModel = node->invtrans*_worldToModel; 
	if (node->obj) {
		// Perform intersection.
		if (node->obj->intersect(ray, _worldToModel, _modelToWorld)) {
			ray.intersection.mat = node->mat;
		}
	}
	// Traverse the children.
	childPtr = node->child;
	while (childPtr != NULL) {
		traverseScene(childPtr, ray);
		childPtr = childPtr->next;
	}

	// Removes transformation of the current node from the global
	// transformation matrices.
	_worldToModel = node->trans*_worldToModel;
	_modelToWorld = _modelToWorld*node->invtrans;
}

void Raytracer::computeShading( Ray3D& ray ) {
	LightListNode* curLight = _lightSource;
	for (;;) {
		if (curLight == NULL) break;
		// Each lightSource provides its own shading function.

		// Implement shadows here if needed.
		Ray3D shaRay;
		char type = (curLight -> light) ->get_type();
		if (type == 'p'){
			
			Point3D l_Pos= (curLight->light)->get_position();
			shaRay.origin = l_Pos;
			shaRay.dir = ray.intersection.point - l_Pos;
			traverseScene(_root, shaRay);
		
			if (!shaRay.intersection.none && shaRay.intersection.obj_index != ray.intersection.obj_index){
				ray.isShadowed = true;	
			}
			
			curLight->light->shade(ray);
		}
		if (type == 'a'){
			
			int shaCoun = 0;
			Point3D sc = (curLight -> light) -> get_position();
			Vector3D u = (curLight -> light) -> get_u();
			Vector3D v = (curLight -> light) -> get_v();
			double lx = (curLight -> light) ->get_lx();
			double ly = (curLight -> light) ->get_ly();

			for (double x = 0; x < 1.0f; x += 0.1f)
				for(double y = 0; y < 1.0f; y += 0.1f){

					Point3D si = sc + (0.5 - x) * lx * u + (0.5 - y) * ly * v;
					shaRay.origin = ray.intersection.point;
					shaRay.dir = si - ray.intersection.point;
					traverseScene(_root, shaRay);

		//			printf ("(%f, %f, %f) \n",shaRay.origin[0],shaRay.origin[1],shaRay.origin[2]);
				
					if (!shaRay.intersection.none && shaRay.intersection.obj_index != ray.intersection.obj_index){
						
					}else{
						curLight->light->shade(ray);
					}
			}
		}
		curLight = curLight->next;
	}
}

void Raytracer::initPixelBuffer() {
	int numbytes = _scrWidth * _scrHeight * sizeof(unsigned char);
	_rbuffer = new unsigned char[numbytes];
	_gbuffer = new unsigned char[numbytes];
	_bbuffer = new unsigned char[numbytes];
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {
			_rbuffer[i*_scrWidth+j] = 0;
			_gbuffer[i*_scrWidth+j] = 0;
			_bbuffer[i*_scrWidth+j] = 0;
		}
	}
}

void Raytracer::loading_texture(char *filename, Texture tex){


}
void Raytracer::flushPixelBuffer( const char *file_name ) {
	bmp_write( file_name, _scrWidth, _scrHeight, _rbuffer, _gbuffer, _bbuffer );
	delete _rbuffer;
	delete _gbuffer;
	delete _bbuffer;
}



Colour Raytracer::shadeRay( Ray3D& ray,int reflection_d ) {
	Colour col(0.0, 0.0, 0.0); 
	traverseScene(_root, ray); 
	// Don't bother shading if the ray didn't hit 
	// anything.
	if (!ray.intersection.none) {
		computeShading(ray); 
		col = ray.col;  
		Vector3D l = ray.dir;
		Vector3D n = ray.intersection.normal;
		l.normalize();
		n.normalize();

		Colour reflCol (0.0, 0.0, 0.0);
		Colour refrCol (0.0, 0.0, 0.0);

		if (reflection_d > 0){
			if (ray.intersection.mat -> refl_coe > 0.00000001 ){
				if (!ray.is_Out){
					n = -n;
				}
				Ray3D reflRay, refrRay;
				Vector3D r = l - 2 * n.dot(l) * n;
				r.normalize();
				// implementing the refration effect
				if (ray.intersection.mat -> transparency > 0.000000001){
					double coe = 1 / ray.intersection.mat -> index_of_refr;
					l.normalize();
					if (!ray.is_Out){
						coe = 1/coe;
					}
					double cos_cri = sqrt(1-1/(coe*coe)); 
					double c = -n.dot(l);

					if (1 + coe*coe*(c*c-1)>0){ //not full reflection

						Vector3D T = coe * l + (coe * c - sqrt(1 + coe*coe*(c*c - 1)))*n;
						refrRay.dir = T;
						refrRay.origin = ray.intersection.point;
						refrRay.is_Out = !ray.is_Out;
						traverseScene(_root, refrRay);
						col = col + ray.intersection.mat -> transparency * shadeRay(refrRay,reflection_d - (1-refrRay.is_Out));
						col.clamp();
					}
				}
				//refraction over
				reflRay = Ray3D(ray.intersection.point,r);
				reflRay.is_Out = ray.is_Out;
				traverseScene(_root, reflRay); 

				if (GLOSSY_MODE){
			
					Vector3D u = n.cross(r);		
					u.normalize();
					Vector3D v = r.cross(u);
					v.normalize();
					double coe = ray.intersection.mat -> glossy_coe;

					if (coe > 0.0000001){
						reflCol = 0.1 * reflCol;
						for (double x = 0.0f; x < 1.0f; x += (1.0f / 3))
							for (double y = 0.0f; y < 1.0f; y += (1.0f / 3)){
								Vector3D _r = r + (0.5 - x) * coe * u + (0.5 - y) * coe * v;
								_r.normalize();
								Ray3D f_reflRay = Ray3D (ray.intersection.point, _r);
								traverseScene(_root, f_reflRay); 
								col = col + 0.1 * ray.intersection.mat->refl_coe * shadeRay(f_reflRay,reflection_d - 1);	
								col.clamp();
							}
					}
				}else{
					col = col + ray.intersection.mat->refl_coe * shadeRay(reflRay,reflection_d - 1);
					col.clamp();
				}	
			}
		}
	}

	// You'll want to call shadeRay recursively (with a different ray, 
	// of course) here to implement reflection/refraction effects.  
	return col; 
}	

void Raytracer::render( int width, int height, Point3D eye, Vector3D view, 
		Vector3D up, double fov, const char* fileName ) {
	Matrix4x4 viewToWorld;
	_scrWidth = width;
	_scrHeight = height;
	double factor = (double(height)/2)/tan(fov*M_PI/360.0);

	initPixelBuffer();
	viewToWorld = initInvViewMatrix(eye, view, up);

	// Construct a ray for each pixel.
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {
			// Sets up ray origin and direction in view space, 
			// image plane is at z = -1.
			float prec = (float)(i * _scrWidth + j)/(_scrHeight * _scrWidth) * 100;
	//		printf("%f % \n", prec);
			if(ANTI_ALIASING_MODE == false && DEPTH_MODE == false){
					Point3D origin(0, 0, 0);
					Point3D imagePlane;
					imagePlane[0] = (-double(width)/2 + 0.5 + j)/factor;
					imagePlane[1] = (-double(height)/2 + 0.5 + i)/factor;
					imagePlane[2] = -1;

					// TODO: Convert ray to world space and call 
					// shadeRay(ray) to generate pixel colour. 	
			
					Ray3D ray;
			
					ray.origin = viewToWorld * origin;
					ray.dir = viewToWorld * imagePlane - ray.origin;
					ray.dir.normalize();
					Colour col = shadeRay(ray, 2); 
					_rbuffer[i*width+j] += int(col[0] * 255 );
					_gbuffer[i*width+j] += int(col[1] * 255 );
					_bbuffer[i*width+j] += int(col[2] * 255 );
			}
			if(ANTI_ALIASING_MODE){
				// Setup 4 rays from {i, j, -1},{i+1, j, -1},{i, j+1, -1},{i+1, j+1, -1}
				// for each pixel to do the antiliasing job
				for (int a = i; a < i+2; a++) {
					for (int b = j; b < j+2; b++){
						Point3D origin(0, 0, 0);
						Point3D imagePlane;
						imagePlane[0] = (-double(width)/2 + 0.5 + b)/factor;
						imagePlane[1] = (-double(height)/2 + 0.5 + a)/factor;
						imagePlane[2] = -1;

						// TODO: Convert ray to world space and call 
						// shadeRay(ray) to generate pixel colour. 	
			
						Ray3D ray;
			
						ray.origin = viewToWorld * origin;
						ray.dir = viewToWorld * imagePlane - ray.origin;
						ray.dir.normalize();
						Colour col = shadeRay(ray, 2); 
						_rbuffer[i*width+j] += int(col[0] * 255 * 1/4);
						_gbuffer[i*width+j] += int(col[1] * 255 * 1/4);
						_bbuffer[i*width+j] += int(col[2] * 255 * 1/4);
					}
				}
			}
			if (DEPTH_MODE){

				Point3D origin(0, 0, 0);
				Point3D imagePlane;

				imagePlane[0] = (-double(width)/2 + 0.5 + j)/factor;
				imagePlane[1] = (-double(height)/2 + 0.5 + i)/factor;
				imagePlane[2] = -1;

				Vector3D dir = imagePlane - origin;
				dir.normalize();
				
				Point3D focal_point = viewToWorld * (origin + FOCAL_DEPTH * dir);

					// TODO: Convert ray to world space and call 
					// shadeRay(ray) to generate pixel colour. 	
			
				for (int a = -1; a < 2; a++) {
					for (int b = -1; b < 2; b++){
					Ray3D ray;
					Point3D this_origin = viewToWorld * origin + Vector3D (a/5.0f, b/5.0f, 0);
					Vector3D this_dir = focal_point - this_origin;
					this_dir.normalize();

					ray.origin = this_origin;
					ray.dir = this_dir;

					Colour col = shadeRay(ray, 2); 
					_rbuffer[i*width+j] += int(col[0] * 255 / 9.0f);
					_gbuffer[i*width+j] += int(col[1] * 255 / 9.0f);
					_bbuffer[i*width+j] += int(col[2] * 255 / 9.0f);

					}
				}
			}
		}
	}

	flushPixelBuffer(fileName);
}



int main(int argc, char* argv[])
{	
	// Build your scene and setup your camera here, by calling 
	// functions from Raytracer.  The code here sets up an example
	// scene and renders it from two different view points, DO NOT
	// change this if you're just implementing part one of the 
	// assignment.  
	Raytracer raytracer;
	int width = 320; 
	int height = 240; 

	if (argc == 3) {
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}

	// Camera parameters.

	Point3D eye(4, 2, 1);
	Vector3D view(-4, -2, -6);
	Vector3D up(0, 1, 0);
	double fov = 60;

	// Defines a texture for shading
	Texture metal;
	


	// Defines a material for shading.
	Material gold( Colour(0.3, 0.3, 0.3), Colour(0.75164, 0.60648, 0.22648), 
			Colour(0.628281, 0.555802, 0.366065), 
			51.2, 0.5, 0, 0.2);

	Material jade( Colour(0, 0, 0), Colour(0.54, 0.89, 0.63), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8, 0.3, 0, 0.1 );
	
	Material silver( Colour(0.19225, 0.19225, 0.19225), Colour(0.50754, 0.50754, 0.50754), 
			Colour(0.508273, 0.508273, 0.508273), 
			51.2, 0.5, 0, 0.2 );

	Material pearl( Colour(0.25, 0.20725, 0.20725), Colour(1, 0.829, 0.829), 
			Colour(0.296648, 0.296648, 0.296648), 
			11.264, 0.1, 0, 0 );
	
	Material mirror( Colour(0.05, 0.05, 0.05), Colour(0.05, 0.05, 0.05),
             Colour(0, 0, 0),
             90, 0.75, 0, 0 );
	Material glass  (Colour(0.0, 0.0, 0.0),Colour(0.588235, 0.670588, 0.729412),
		     Colour(0.9, 0.9, 0.9),
             96, 0.5, 1.55, 0.8, 0);
	// Defines a point light source.
	
	raytracer.addLightSource( new PointLight(Point3D(0,0,5), Colour(0.9,0.9,0.9)) );

	// Defines an area light source
	//raytracer.addLightSource( new AreaLight(Point3D (0, 0, 5), 
	//	Colour (0.9, 0.9, 0.9), Vector3D (1, 0, 0), Vector3D (0, 1, 0), 3, 3) );

	
	// Add a unit square into the scene with material mat.

	SceneDagNode* sphere = raytracer.addObject( new UnitSphere(0), &silver );
	SceneDagNode* cone = raytracer.addObject( new UnitCone(1), &gold );
	SceneDagNode* cylinder = raytracer.addObject( new UnitCylinder(2), &pearl );
	SceneDagNode* plane = raytracer.addObject( new UnitSquare(3), &jade );
	// Apply some transformations to the unit square.
	double factor1[3] = { 1.0, 2.0, 1.0 };
	double factor2[3] = { 8.0, 8.0, 8.0 };
	double factor3[3] = { 0.8, 1.8, 1.2 };
	double factor4[3] = { 0.9, 0.9, 3.0 };
	double factor5[3] = { 0.8, 0.8, 1.6 };

	raytracer.translate(sphere, Vector3D(0, 0, -7));	
	raytracer.rotate(sphere, 'x', -45); 
	raytracer.rotate(sphere, 'z', 45); 
	raytracer.scale(sphere, Point3D(0, 0, 0), factor1);

	raytracer.translate(cone, Vector3D(-1, 0, -5));	
	raytracer.scale(cone, Point3D(0, 0, 0), factor4);
	
	raytracer.translate(cylinder, Vector3D(1.8, 1.8, -4.5));	
	raytracer.scale(cylinder, Point3D(0, 0, 0), factor5);

	raytracer.translate(plane, Vector3D(0, 0, -7));	
	raytracer.rotate(plane, 'z', 45); 
	raytracer.scale(plane, Point3D(0, 0, 0), factor2);

	// Render the scene, feel free to make the image smaller for
	// testing purposes.	
	raytracer.render(width, height, eye, view, up, fov, "view1.bmp");
	
	// Render it from a different point of view.
	
	Point3D eye2(-10, 0, -5);
	Vector3D view2(10, 0, 0);
	raytracer.render(width, height, eye2, view2, up, fov, "view2.bmp");

	return 0;
}
