#include "stdafx.h"
#include "Scene.h"


int Scene3D::addObject(Object *object)
{
	objectVector.push_back(object);
	return 0;
}

int Scene3D::addLight(Light * light)
{
	lightVector.push_back(light);
	return 0;
}

std::vector<Object*> Scene3D::getObjectVector()
{
	return objectVector;
}

std::vector<Light*> Scene3D::getLightVector()
{
	return lightVector;
}

int Scene3D::setCamera(Camera * camera)
{
	hasCamera = 1;
	this->camera = camera;
	return 1;
}

int Scene3D::removeCamera()
{
	this->camera == NULL;
	hasCamera = 0;
	return 1;
}

Camera * Scene3D::getCamera()
{
	if (hasCamera) {
		return camera;
	}
	return NULL;
}

void Scene3D::removeAll()
{
	objectVector.clear();
	lightVector.clear();
}

Scene3D::Scene3D():objectVector(200),lightVector(50)
{
	objectVector.clear();
	lightVector.clear();
}


Scene3D::~Scene3D()
{
}



int Scene2D::addObject(Object2D * object)
{
	objectVector.push_back(object);
	return 1;
}

void Scene2D::removeAll()
{
	objectVector.clear();
}

Scene2D::Scene2D():objectVector(200)
{
	objectVector.clear();
}

Scene2D::~Scene2D()
{
}
