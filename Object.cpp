#include "stdafx.h"
#include "Object.h"


void Object::translate(Vector3D * v)
{
	Point4D_Move(&worldPosition, v);
}

void Object::rotate(float thX, float thY, float thZ, int allFrame)
{
	Matrix4x4 rotate;
	Point4D tmpP;
	Vector4D tmpN;
	int vNum, nNum;
	Vector4D *nList;
	Vertex4D *vList;
	Matrix_Build_Rotation_4x4(&rotate, thX, thY, thZ);

	if (allFrame) {
		vList = headOriginVList;
		nList = headOriginNList;
		vNum = totalVerticesNum;
		nNum = framesNum*polygonsNum;
	}
	else
	{
		vList = originVList;
		nList = originNList;
		vNum = verticesNum;
		nNum = polygonsNum;
	}

	for (int i = 0;i < vNum;i++) {
		Matrix_Mul_Vector4D_4x4(&vList[i].p, &rotate, &tmpP);
		Point4D_Copy(&vList[i].p, &tmpP);
		if (vList[i].attr&VERTEX4D_ATTR_NORMAL) {
			Matrix_Mul_Vector4D_4x4(&vList[i].n, &rotate, &tmpP);
			Point4D_Copy(&vList[i].n, &tmpP);
		}
	}
	for (int i = 0;i < nNum;i++) {
		Matrix_Mul_Vector4D_4x4(&nList[i], &rotate, &tmpN);
		Point4D_Copy(&nList[i], &tmpN);
	}

	//axes rotate
	Matrix_Mul_Vector4D_4x4(&ux, &rotate, &tmpN);
	Vector4D_Copy(&ux, &tmpN);

	Matrix_Mul_Vector4D_4x4(&uy, &rotate, &tmpN);
	Vector4D_Copy(&uy, &tmpN);

	Matrix_Mul_Vector4D_4x4(&uz, &rotate, &tmpN);
	Vector4D_Copy(&uz, &tmpN);


}

void Object::scale(float sX, float sY, float sZ, int allFrame)
{
	float maxScale;
	int vNum, nNum;
	Vertex4D *vList;
	Vector4D *nList;

	if (allFrame) {
		vList = headOriginVList;
		nList = headOriginNList;
		vNum = totalVerticesNum;
		nNum = framesNum*polygonsNum;
	}
	else
	{
		vList = originVList;
		nList = originNList;
		vNum = verticesNum;
		nNum = polygonsNum;
	}

	for (int i = 0;i < vNum;i++) {
		vList[i].p.x *= sX;
		vList[i].p.y *= sY;
		vList[i].p.z *= sZ;
	}
	for (int i = 0;i < nNum;i++) {
		nList[i].x *= sX;
		nList[i].y *= sY;
		nList[i].z *= sZ;
	}

	maxScale = MAX(sY, sX);
	maxScale = MAX(maxScale, sZ);
	if (!allFrame) {

		avgRadius[currFrame] *= maxScale;
		maxRadius[currFrame] *= maxScale;
	}
	else
	{
		for (int i = 0;i < framesNum;i++) {
			avgRadius[i] *= maxScale;
			maxRadius[i] *= maxScale;
		}
	}
}

void Object::initPreComputeNormal(int mode)
{
	int indexP0, indexP1, indexP2,indexN;
	Vector4D tmpV1, tmpV2;
	Vector4D *normal;
	Point4D *p0, *p1, *p2;
	Vector4D *n0, *n1, *n2;
	int totalVerticesNum;

	if (mode == OBJECT_PRECOMPUTE_ALL) {
		//init 
		totalVerticesNum = verticesNum*framesNum;
		int* polyTouchNum = new int[totalVerticesNum];
		for(int vidx = 0;vidx < totalVerticesNum;vidx++) {
			Vector4D_Zero(&headOriginVList[vidx].n);
			polyTouchNum[vidx] = 0;
			SET_BIT(headOriginVList[vidx].attr, VERTEX4D_ATTR_NORMAL);
			SET_BIT(headTransVList[vidx].attr, VERTEX4D_ATTR_NORMAL);
		}

		for (int fidx = 0;fidx < framesNum;fidx++) {
			for (int pidx = 0;pidx < polygonsNum;pidx++) {
				indexP0 = polyList[pidx].vertextIndex[0];
				indexP1 = polyList[pidx].vertextIndex[1];
				indexP2 = polyList[pidx].vertextIndex[2];
				indexN = polyList[pidx].normorIndex;
				normal = &headOriginNList[fidx*polygonsNum + indexN];
				p0 = &headOriginVList[fidx*verticesNum + indexP0].p;
				n0 = &headOriginVList[fidx*verticesNum + indexP0].n;
				p1 = &headOriginVList[fidx*verticesNum + indexP1].p;
				n1 = &headOriginVList[fidx*verticesNum + indexP1].n;
				p2 = &headOriginVList[fidx*verticesNum + indexP2].p;
				n2 = &headOriginVList[fidx*verticesNum + indexP2].n;

				Vector4D_Build(p0, p1, &tmpV1);
				Vector4D_Build(p0, p2, &tmpV2);
				Vector4D_Cross(&tmpV1, &tmpV2, normal);

				Vector4D_Add(n0, normal);
				Vector4D_Add(n1, normal);
				Vector4D_Add(n2, normal);

				polyTouchNum[fidx*verticesNum + indexP0]++;
				polyTouchNum[fidx*verticesNum + indexP1]++;
				polyTouchNum[fidx*verticesNum + indexP2]++;

				Vector4D_Normalize(normal);
			}
		}

		for (int vidx = 0;vidx < totalVerticesNum;vidx++) {
			if (polyTouchNum[vidx] >= 1) {
				headOriginVList[vidx].nx /= polyTouchNum[vidx];
				headOriginVList[vidx].ny /= polyTouchNum[vidx];
				headOriginVList[vidx].nz /= polyTouchNum[vidx];
				Vector4D_Normalize(&headOriginVList[vidx].n);
			}
		}
		delete[] polyTouchNum;
	}
	else if (mode == OBJECT_PRECOMPUTE_POLYNORMAL)
	{
		for (int fidx = 0;fidx < framesNum;fidx++) {
			for (int pidx = 0;pidx < polygonsNum;pidx++) {
				indexP0 = polyList[pidx].vertextIndex[0];
				indexP1 = polyList[pidx].vertextIndex[1];
				indexP2 = polyList[pidx].vertextIndex[2];
				indexN = polyList[pidx].normorIndex;
				normal = &headOriginNList[fidx*polygonsNum + indexN];
				p0 = &headOriginVList[fidx*verticesNum + indexP0].p;
				p1 = &headOriginVList[fidx*verticesNum + indexP1].p;
				p2 = &headOriginVList[fidx*verticesNum + indexP2].p;

				Vector4D_Build(p0, p1, &tmpV1);
				Vector4D_Build(p0, p2, &tmpV2);
				Vector4D_Cross(&tmpV1, &tmpV2, normal);

				Vector4D_Normalize(normal);
			}
		}
	}


}

void Object::setFrame(int frame)
{
	originVList = &headOriginVList[frame*verticesNum];
	transVList  = &headTransVList[frame*verticesNum];
	originNList = &headOriginNList[frame*polygonsNum];
	transNList = &headTransNList[frame*polygonsNum];
}

Object::Object()
{

}


Object::~Object()
{
}


Object2D::Object2D()
{
}

Object2D::Object2D(float x, float y, int width, int height, Bitmap * bitmap)
{

	pos.x = x;
	pos.y = y;
	this->width = width;
	this->height = height;
	rect.left = x;
	rect.top = y;
	rect.right = rect.left+width-1;
	rect.bottom = rect.top+height-1;
	this->bitmap = bitmap;
}


Object2D::Object2D(int left, int top, int right, int bottom, Bitmap * bitmap)
{	
	pos.x = left;
	pos.y = top;
	rect.left = left;
	rect.top = top;
	rect.right = right;
	rect.bottom = bottom;
	width = rect.right - rect.left+1;
	height = rect.bottom - rect.top + 1;
	this->bitmap = bitmap;
}

void Object2D::translate(float x, float y)
{
	pos.x += x;
	pos.y += y;
	rect.left += x;
	rect.right += x;
	rect.top += y;
	rect.bottom += y;
}

void Object2D::scale(float sx, float sy)
{
	width *= sx;
	height *=sy;
	rect.right = rect.left + width - 1;
	rect.bottom = rect.top + height - 1;
}

Object2D::~Object2D()
{
}
