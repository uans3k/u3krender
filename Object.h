#pragma once
#include "Math3D.h"
#include "Camera.h"
#include "GraphicsStruct.h"
#include "Bitmap.h"

class Object2D {
public:
	Bitmap* bitmap;
	RECT rect;
	Point2D pos;
	RECT selfRect;
	Point2D pinPos;
	int width;
	int height;
	Object2D();
	Object2D(float x, float y, int width, int height,Bitmap *bitmap);
	Object2D(int left,int top,int  right, int bottom, Bitmap *bitmap);
	void translate(float x, float y);
	void scale(float sx, float sy);
	~Object2D();
};


class Object
{
public:

#define OBJECT_CULL_NONE 0x0
#define OBJECT_CULL_X 0x01
#define OBJECT_CULL_Y 0x02
#define OBJECT_CULL_Z 0x04
#define OBJECT_CULL_ALL 0x07

#define OBJECT_TRANS_LOCAL            0
#define OBJECT_TRANS_TRANS            1
#define OBJECT_TRANS_LOCAL2TRANS      2

#define OBJECT_STATE_NULL             0x0000
#define OBJECT_STATE_UNACTIVE           0x0001
#define OBJECT_STATE_UNVISIBLE          0x0002 
#define OBJECT_STATE_CULLED           0x0004

#define OBJECT_PRECOMPUTE_POLYNORMAL 0x01
#define OBJECT_PRECOMPUTE_VERTEXTNORMAL 0x02
#define OBJECT_PRECOMPUTE_ALL 0x03

	friend class RenderPiping;
		
	int   id;           // numeric id of this object
	char  name[64];     // ASCII name of object just for kicks
	int   state;        // state of object
	int   attr;         // attributes of object
	int   mati;         // material index overide (-1) - no material (new)

	Point4D worldPosition;  // position of object in world

	Vector4D dir;       // rotation angles of object in local
						// cords or unit direction vector user defined???

	Vector4D ux, uy, uz;  // local axes to track full orientation
						  // this is updated automatically during
						  // rotation calls

	//vertex

	Vertex4D *originVList;//current frame origin vertices
	Vertex4D *transVList; //current frame transformed vertices
	float *avgRadius;  // average radius of object used for collision detection
	float *maxRadius;   // maximum radius of object
	int verticesNum;   // number of vertices per frame of this object
	Vertex4D *headOriginVList;
	Vertex4D *headTransVList;				 
	int framesNum;     // number of frames
	int currFrame;     // current animation frame (0) if single frame
	int totalVerticesNum; // total number of vertices framesNum* verticesNum

	//texture
	Bitmap *texture;
	Point2D *uvList;       
	int uvNum;


	//poly
	int polygonsNum;           // number of polygons in object mesh
	PolyTriangle4D *polyList;      // polygons 
	Vector4D *originNList;
	Vector4D *transNList;		//polyonsNum 
	Vector4D *headOriginNList; //polyonsNum * frame
	Vector4D *headTransNList;

	void translate(Vector3D *v);

	void rotate(float thX,float thY,float thZ,int allFrame = 0);

	void scale(float sX, float sY, float sZ, int allFrame = 0);

	void initPreComputeNormal(int mode = OBJECT_PRECOMPUTE_ALL);

	void setFrame(int frame);

	Object();
	~Object();
private:

};

