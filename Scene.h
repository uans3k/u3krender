#pragma once
#include "Object.h"
#include "Light.h"
#include <vector>
#include "Bitmap.h"

class Scene3D
{
public:
	int addObject(Object *object);
	int addLight(Light *light);
	std::vector<Object *> getObjectVector();
	std::vector<Light *>  getLightVector();
	int setCamera(Camera *camera);
	int removeCamera();
	Camera* getCamera();
	void removeAll();
	Scene3D();
	~Scene3D();
private:
	Camera *camera=NULL;
	int hasCamera = 0;
	std::vector<Object *> objectVector;
	std::vector<Light *> lightVector;

};


class Scene2D
{
public:
	std::vector<Object2D *> objectVector;
	int addObject(Object2D *object);
	void removeAll();
	Scene2D();
	~Scene2D();
};
