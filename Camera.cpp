#include "Camera.h"
#include <math.h>

Camera::Camera(Point4D * pos, float viewPortWidth, float viewPortHeight, float clipNearZ, float clipFarZ, float fov)
{
	Point4D_Copy(&this->pos, pos);
	this->viewPortWidth = viewPortWidth;
	this->viewPortHeight = viewPortHeight;
	this->clipFarZ = clipFarZ;
	
	this->fov = fov;

	aspectRatio = viewPortWidth / viewPortHeight;



	viewPortCenterX = 0.5f*viewPortWidth - 0.5f;
	viewPortCenterY = 0.5f*viewPortWidth - 0.5f;

	viewPlaneWidth = 2.0f; //1.for perspective matrix 2.for make model size absolutely (let eyes same)
	viewPlaneHeight = 2.0f / aspectRatio;

	float tanTh = tan(DEG_TO_RAD(fov / 2));
	viewDist = 0.5*viewPlaneWidth / tanTh;
	this->clipNearZ = MAX(viewDist, clipNearZ);

	Matrix_Init_4x4(&cameraMatrix, 1, 0, 0, 0,//
		0, 1, 0, 0,//
		0, 0, 1, 0,
		-this->pos.x, -this->pos.y, -this->pos.z, 1);
}

Camera::~Camera()
{

}

void Camera::translate(float x, float y, float z)
{
	pos.x += x;
	pos.y += y;
	pos.z += z;
	Matrix4x4 translateM,tmpM;
	Matrix_Build_Translation_4x4(&translateM, -x, -y, -z);
	Matrix_Mul_4x4(&cameraMatrix,&translateM,&tmpM);
	MATRIX_COPY_4X4(&tmpM, &cameraMatrix);
}

CameraUVN::CameraUVN(Point4D *pos, Point4D *targetP, float viewPortWidth, float viewPortHeight, float clipNearZ, float clipFarZ, float fov) :Camera(pos, viewPortWidth, viewPortHeight, clipNearZ, clipFarZ, fov)
{
	Vector4D vY;
	Matrix4x4 tmpM, uvnM;
	Point4D_Init(&this->targetPoint, targetP);
	Vector4D_Build(pos, targetP, &n);
	Vector4D_Init(&vY, 0.0f, 1.0f, 0.0f);
	Vector4D_Cross(&vY, &n, &u);
	Vector4D_Cross(&n, &u, &v);

	Vector4D_Normalize(&u);
	Vector4D_Normalize(&v);
	Vector4D_Normalize(&n);

	Matrix_Init_4x4(&uvnM, u.x, v.x, n.x, 0.0f,
							u.y, v.y, n.y, 0.0f,
							u.z, v.z, n.z, 0.0f,
							0, 0, 0, 1.0f);
	Matrix_Mul_4x4(&cameraMatrix, &uvnM, &tmpM);
	MATRIX_COPY_4X4(&tmpM, &cameraMatrix);
}

CameraUVN::~CameraUVN()
{
}
