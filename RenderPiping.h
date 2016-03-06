#pragma once
#include <vector>

#include "Math3D.h"
#include "Object.h"
#include "GraphicsStruct.h"
#include "Light.h"
#include "Canvas.h"


int  cmpByAvgZ(const void *arg1, const void *arg2);

int  cmpByMaxZ(const void *arg1, const void *arg2);

int  cmpByMinZ(const void *arg1, const void *arg2);



#define RENDLIST_SORT_AVGZ 0
#define RENDLIST_SORT_MAXZ 1
#define RENDLIST_SORT_MINZ 2


#define RENDLIST_CLIP_X   0x0001
#define RENDLIST_CLIP_Y   0x0002
#define RENDLIST_CLIP_Z   0x0004
#define RENDLIST_CLIP_XYZ 0x0007

class RenderList
{
public:
	int num;
	int max;

	PolySelfTriangle4D *polyDataList;
	PolySelfTriangle4D **polyPtrList;
	RenderList(int maxItems=200);
	~RenderList();
	void clear();
	int insertObject(Object *obj);
	int insertPoly4D(PolyTriangle4D *poly);
	int insertPoly4D(PolySelfTriangle4D *poly);
	void sortByZ(int sortMode=RENDLIST_SORT_AVGZ);

private:

	

};

class RenderPiping {

#define RENDERPIPING_CULL_NONE 0x0
#define RENDERPIPING_CULL_X 0x01
#define RENDERPIPING_CULL_Y 0x02
#define RENDERPIPING_CULL_Z 0x04
#define RENDERPIPING_CULL_XYZ 0x07

#define RENDERPIPING_CLIP_NONE 0x00
#define RENDERPIPING_CLIP_X   0x01
#define RENDERPIPING_CLIP_Y   0x02
#define RENDERPIPING_CLIP_Z   0x04
#define RENDERPIPING_CLIP_XYZ 0x07


public:
	RenderPiping();
	~RenderPiping();

	void render(Canvas *canvas);
	void clear();
	void setCamera(Camera *camera);
	int  sendObject(Object *object);
	void addLight(Light *light);
	void setSortZEnable(int enable);
private:
	//object
	void worldTransfrom(Object *object);
	void cull(Object *object);
	void removeBack(Object *object);
	void cameraTransform(Object *object);
	void cameraTransform(Light *light);


	//render list
	void clipPoly();
	void lightPoly();
	void transformPerspective();
	void transformScreen();
	void transformPerspectiveScreen();

	int cullMode= RENDERPIPING_CULL_XYZ;
	int clipMode= RENDERPIPING_CLIP_XYZ;
	
	int isSortZ=1;
	RenderList renderList=RenderList(1000);
	std::vector<Light*> lightVector= std::vector<Light*>(100);
	Camera *camera;
};

