#pragma once
#include "GraphicsStruct.h"
#include "Math3D.h"

#define LIGHT_STATE_OFF 0
#define LIGHT_STATE_ON 1

#define LIGHT_ATTR_AMBIENT 0
#define LIGHT_ATTR_INFINITY 1
#define LIGHT_ATTR_POINT 2
#define LIGHT_ATTR_SPOT 3

struct Light
{	
	int state;
	int attr;

	Color ambient;
	Color diffuse;
	Color specular;

	
	Point4D oPos;
	Point4D tPos;
	Vector4D oDir;
	Vector4D tDir;

	float kc, kl, kq;
	float spotInner;  
	float spotOutter;  
	float pf;
	void init(Point4D *pos, Vector4D *dir,
			  Color ambient,Color diffuse,Color specular,
			  float kc,float kl,float kq,
		      float spotInner,float spotOutter,float pf,
		      int attr);
	void initAmbientLight(Color ambient);
	void initInfinityLight(Vector4D *dir,Color diffuse, Color specular);
    void initPointLight(Point4D *pos,
		 Color diffuse, Color specular,
		float kc, float kl, float kq);
	void initSpotLight(Point4D *pos, Vector4D *dir,
		 Color diffuse, Color specular,
		float kc, float kl, float kq,float pf);
	inline void open() {
		state = LIGHT_STATE_ON;
	}
	inline void close() {
		state = LIGHT_STATE_OFF;
	}
};