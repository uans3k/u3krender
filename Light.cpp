#include "stdafx.h"
#include "Light.h"

void Light::init(Point4D *pos, Vector4D *dir, Color ambient, Color diffuse, Color specular, float kc, float kl, float kq, float spotInner, float, float pf, int attr)
{	
	this->ambient = ambient;
	this->diffuse = diffuse;
	this->specular = specular;
	this->kc = kc;	
	this->kl = kl;
	this->kq = kq;
	this->spotInner = spotInner;
	this->spotOutter = spotOutter;
	this->pf = pf;
	this->attr = attr;


	if (pos) {
		Point4D_Init(&oPos,pos);
		Point4D_Init(&tPos, &oPos);
	}

	if (dir) {
		Vector4D_Init(&oDir, dir);
		Vector4D_Normalize(&oDir);
		Vector4D_Init(&tDir,&oDir);
	}
}

void Light::initAmbientLight(Color ambient)
{	
	this->ambient = ambient;
	this->attr = LIGHT_ATTR_AMBIENT;
	open();
}

void Light::initInfinityLight(Vector4D *dir,Color diffuse, Color specular)
{

	this->diffuse = diffuse;
	this->specular = specular;


	if (dir) {
		Vector4D_Init(&oDir, dir);
		Vector4D_Normalize(&oDir);
		Vector4D_Init(&tDir, &oDir);
	}
	this->attr = LIGHT_ATTR_INFINITY;
	open();
}
 void Light::initPointLight(Point4D *pos,
	 Color diffuse, Color specular,
	float kc, float kl, float kq)
{

	this->diffuse = diffuse;
	this->specular = specular;
	this->kc = kc;
	this->kl = kl;
	this->kq = kq;


	if (pos) {
		Point4D_Init(&oPos, pos);
		Point4D_Init(&tPos, &oPos);
	}
	this->attr = LIGHT_ATTR_POINT;
	open();
}

void Light::initSpotLight(Point4D * pos, Vector4D * dir, Color diffuse, Color specular, float kc, float kl, float kq,float pf)
{
	
	this->diffuse = diffuse;
	this->specular = specular;
	this->kc = kc;
	this->kl = kl;
	this->kq = kq;
	this->pf = pf;



	if (pos) {
		Point4D_Init(&oPos, pos);
		Point4D_Init(&tPos, &oPos);
	}

	if (dir) {
		Vector4D_Init(&oDir, dir);
		Vector4D_Normalize(&oDir);
		Vector4D_Init(&tDir, &oDir);
	}

	this->attr = LIGHT_ATTR_SPOT;
	open();
}
