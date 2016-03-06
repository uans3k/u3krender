#pragma 
#include "Math3D.h"

class Camera
{
public:
	
	Point4D pos;
	float viewDist;
	float fov;
	float clipNearZ;     
	float clipFarZ;
	
	Plane3D clipLeftPlane, clipRightPlane, clipUpPlane, clipDownPlane;


	float aspectRatio; //   viewPlaneWidth / viewPlaneHeight

	float viewPlaneWidth;   
	float viewPlaneHeight;
	float viewPortWidth;    
	float viewPortHeight;
	float viewPortCenterX;  
	float viewPortCenterY;

	Matrix4x4 cameraMatrix;
	Matrix4x4 perspectiveMatrix;
	Matrix4x4 screenMatrix;

	Camera(Point4D *pos,float viewPortWidth,float viewPortHeight,float clipNearZ,float clipFarZ,float fov=90.0f);
	~Camera();

	void translate(float x, float y, float z);

private:
};

class CameraUVN :public Camera {
public:
	Vector4D u;
	Vector4D v;
	Vector4D n;
	Point4D targetPoint;
	CameraUVN(Point4D *pos, Point4D *targetP, float viewPortWidth, float viewPortHeight, float clipNearZ, float clipFarZ, float fov = 90.0f);
	~CameraUVN();
};
