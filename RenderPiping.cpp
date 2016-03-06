#include "stdafx.h"
#include "RenderPiping.h"


RenderList::RenderList(int maxItems)
{
	max = maxItems;
	polyDataList = new  PolySelfTriangle4D[max];
	polyPtrList = new PolySelfTriangle4D*[max];
	num = 0;
}

RenderList::~RenderList()
{
	if (polyDataList) {
		delete[] polyDataList;
	}
	if (polyPtrList) {
		delete[] polyPtrList;
	}
}

void RenderList::clear()
{
	num = 0;
}

int RenderList::insertObject(Object * obj)
{
	if (IS_BIT(obj->state, OBJECT_STATE_UNACTIVE) || IS_BIT(obj->state, OBJECT_STATE_UNVISIBLE) || IS_BIT(obj->state, OBJECT_STATE_CULLED)) {
		return 1;
	}

	for (int i = 0; i < obj->polygonsNum; i++)
	{
		if (!insertPoly4D(&obj->polyList[i])) {
			return 0;
		}
	}

	return 1;
}

int RenderList::insertPoly4D(PolyTriangle4D * poly)
{
	if (num == max) {
		return 0;
	}
	if (IS_BIT(poly->state, POLY4D_STATE_UNACTIVE) || IS_BIT(poly->state, POLY4D_STATE_BACKFACE) || IS_BIT(poly->state, POLY4D_STATE_CLIPPED)) {
		return 1;
	}
	polyDataList[num].initFromTrans(poly);
	polyPtrList[num] = &polyDataList[num];
	num++;
	return 1;
}

int RenderList::insertPoly4D(PolySelfTriangle4D * poly)
{
	if (num == max) {
		return 0;
	}
	if (IS_BIT(poly->state, POLY4D_STATE_UNACTIVE) || IS_BIT(poly->state, POLY4D_STATE_BACKFACE) || IS_BIT(poly->state, POLY4D_STATE_CLIPPED)) {
		return 1;
	}
	polyDataList[num].init(poly);
	polyPtrList[num] = &polyDataList[num];
	num++;
	return 1;
}


void RenderList::sortByZ(int mode)
{
	switch (mode)
	{
	case RENDLIST_SORT_MAXZ:
		qsort((void*)polyPtrList, num, sizeof(PolySelfTriangle4D*), cmpByMaxZ);
		break;
	case RENDLIST_SORT_MINZ:
		qsort((void*)polyPtrList, num, sizeof(PolySelfTriangle4D*), cmpByMinZ);
		break;

	case RENDLIST_SORT_AVGZ:
		qsort((void*)polyPtrList, num, sizeof(PolySelfTriangle4D*), cmpByAvgZ);
	default:
		break;
	}
}




int  cmpByAvgZ(const void * arg1, const void * arg2)
{
	PolySelfTriangle4D *poly1 = *((PolySelfTriangle4D **)arg1);
	PolySelfTriangle4D *poly2 = *((PolySelfTriangle4D **)arg2);
	float avgZ1, avgZ2;
	avgZ1 = 0.33333f*(poly1->transVList[0].z + poly1->transVList[1].z + poly1->transVList[2].z);
	avgZ2 = 0.33333f*(poly2->transVList[0].z + poly2->transVList[1].z + poly2->transVList[2].z);
	if (avgZ1 < avgZ2) {
		return 1;
	}
	else if (avgZ1>avgZ2) {
		return -1;
	}
	else
	{
		return 0;
	}
}

int cmpByMaxZ(const void * arg1, const void * arg2)
{
	PolySelfTriangle4D *poly1 = *((PolySelfTriangle4D **)arg1);
	PolySelfTriangle4D *poly2 = *((PolySelfTriangle4D **)arg2);
	float avgZ1, avgZ2;

	avgZ1 = poly1->transVList[0].z;
	avgZ1 = MAX(avgZ1, poly1->transVList[1].z);
	avgZ1 = MAX(avgZ1, poly1->transVList[2].z);

	avgZ2 = poly2->transVList[0].z;
	avgZ2 = MAX(avgZ2, poly2->transVList[1].z);
	avgZ2 = MAX(avgZ2, poly2->transVList[2].z);

	if (avgZ1 < avgZ2) {
		return 1;
	}
	else if (avgZ1>avgZ2) {
		return -1;
	}
	else
	{
		return 0;
	}
}

int cmpByMinZ(const void * arg1, const void * arg2)
{
	PolySelfTriangle4D *poly1 = *((PolySelfTriangle4D **)arg1);
	PolySelfTriangle4D *poly2 = *((PolySelfTriangle4D **)arg2);
	float avgZ1, avgZ2;

	avgZ1 = poly1->transVList[0].z;
	avgZ1 = MIN(avgZ1, poly1->transVList[1].z);
	avgZ1 = MIN(avgZ1, poly1->transVList[2].z);

	avgZ2 = poly2->transVList[0].z;
	avgZ2 = MIN(avgZ2, poly2->transVList[1].z);
	avgZ2 = MIN(avgZ2, poly2->transVList[2].z);

	if (avgZ1 < avgZ2) {
		return 1;
	}
	else if (avgZ1>avgZ2) {
		return -1;
	}
	else
	{
		return 0;
	}
}



void RenderPiping::worldTransfrom(Object * object)
{
	//vertices first transfor
	for (int i = 0;i < object->verticesNum;i++) {
		Point4D_Move(&object->originVList[i].p, &object->transVList[i].p, &object->worldPosition);
		Vector4D_Copy(&object->transVList[i].n, &object->originVList[i].n);
	}
}

void RenderPiping::cameraTransform(Object * object)
{
	Point4D tmpP;
	Vector4D tmpN;
	for (int i = 0;i < object->verticesNum;i++) {
		Matrix_Mul_Point4D_4x4(&object->transVList[i].p, &camera->cameraMatrix, &tmpP);
		Point4D_Copy(&object->transVList[i].p, &tmpP);
		if (IS_BIT(object->transVList[i].attr, VERTEX4D_ATTR_NORMAL)) {
			Matrix_Mul_Vector4D_4x4(&object->transVList[i].n, &camera->cameraMatrix, &tmpN);
			Vector4D_Copy(&object->transVList[i].n, &tmpN);
		}
	}
	//poly normal first trans
	for (int i = 0;i < object->polygonsNum;i++) {
		Matrix_Mul_Vector4D_4x4(&object->originNList[i], &camera->cameraMatrix, &object->transNList[i]);
	}
}

void RenderPiping::cameraTransform(Light * light)
{	//light first transform
	Matrix_Mul_Vector4D_4x4(&light->oPos, &camera->cameraMatrix, &light->tPos);
	Matrix_Mul_Vector4D_4x4(&light->oDir, &camera->cameraMatrix, &light->tDir);
}

void RenderPiping::cull(Object * object)
{
	Point4D posInCamera;
	float test;
	Matrix_Mul_Vector4D_4x4(&object->worldPosition, &camera->cameraMatrix, &posInCamera);
	if (IS_BIT(cullMode, RENDERPIPING_CULL_Z)) {

		if (posInCamera.z + object->maxRadius[object->currFrame]<camera->clipNearZ || posInCamera.z - object->maxRadius[object->currFrame]>camera->clipFarZ) {
			SET_BIT(object->state, OBJECT_STATE_CULLED);
			return;
		}
	}
	if (IS_BIT(cullMode, RENDERPIPING_CULL_X)) {
		//test=tan(fov/2)*z
		//tan(fov/2)=0.5*viewWidth/viewDist;
		test = posInCamera.z*(0.5f)*camera->viewPlaneWidth / camera->viewDist;
		if ((posInCamera.x - object->maxRadius[object->currFrame] > test) || (posInCamera.x + object->maxRadius[object->currFrame] < -test)) {
			SET_BIT(object->state, OBJECT_STATE_CULLED);
		}
	}

	if (IS_BIT(cullMode, RENDERPIPING_CULL_Y)){
		//test=tan(fov/2)*z
		//tan(fov/2)=0.5*viewHeight/viewDist;
		test = posInCamera.z * (0.5f)*camera->viewPlaneHeight / camera->viewDist;
		if ((posInCamera.y - object->maxRadius[object->currFrame]> test) || (posInCamera.y + object->maxRadius[object->currFrame]< -test)) {
			SET_BIT(object->state, OBJECT_STATE_CULLED);
		}
	}
}

void RenderPiping::removeBack(Object * object)
{
	//left hand and time cycle stand normal of polygon
	//
	PolyTriangle4D *poly;
	int indexP0, indexP1, indexP2;
	Vector4D tmpVLook;
	float test;


	if (IS_BIT(object->state,OBJECT_STATE_CULLED)) {
		return;
	}

	for (int i = 0;i < object->polygonsNum;i++) {
		poly = &object->polyList[i];
		if ((IS_BIT(poly->state, POLY4D_STATE_UNACTIVE)) || IS_BIT(poly->state, POLY4D_STATE_BACKFACE)|| IS_BIT(poly->attr, POLY4D_ATTR_2SIDED)) {
			continue;
		}
		indexP0 = poly->vertextIndex[0];
		indexP1 = poly->vertextIndex[1];
		indexP2 = poly->vertextIndex[2];

	

		Vector4D_Build(&object->transVList[indexP0].p, &camera->pos, &tmpVLook);
		test = Vector4D_Dot(&tmpVLook, &poly->originNList[poly->normorIndex]);
		if (test<0.0f) {
			SET_BIT(poly->state, POLY4D_STATE_BACKFACE);
		}
	}
}

void RenderPiping::clipPoly()
{
#define CLIP_CODE_L 1
#define CLIP_CODE_G 2
#define CLIP_CODE_I 3

	int clipCode[3];//for v0,v1,v2
	int inNum;
	int polyNum;
	int v0, v1, v2;//index of vertext ,assure vertext order (assure normal same)w
	float tanTh, testL, testG;
	float leftFraction, rightFraction;;
	Point4D insertRightP, insertLeftP;
	Point2D insertRightT, insertLeftT;
	Vector4D insertRightN, insertLeftN;
	PolySelfTriangle4D* currPoly, tmpPoly;//current poly

	Point2D_Zero(&insertRightT);
	Point2D_Zero(&insertLeftT);
	Point4D_Zero(&insertRightP);
	Point4D_Zero(&insertLeftP);

	polyNum = renderList.num;
	for (int i = 0;i < polyNum;i++) {
		currPoly = renderList.polyPtrList[i];
		//clip x plane
		if (IS_BIT(clipMode, RENDERPIPING_CLIP_X)) {
			inNum = 0;
			tanTh = camera->viewPlaneWidth / camera->viewDist;
			for (int vindex = 0;vindex<3;vindex++) {
				testG = tanTh*currPoly->transVList[vindex].z;
				testL = -testG;
				if (currPoly->transVList[vindex].x>testG) {
					clipCode[vindex] = CLIP_CODE_G;
				}
				else if (currPoly->transVList[vindex].x < testL)
				{
					clipCode[vindex] = CLIP_CODE_L;
				}
				else
				{
					clipCode[vindex] = CLIP_CODE_I;
					inNum++;
				}
			}

			if (inNum == 0) {
				SET_BIT(currPoly->state, POLY4D_STATE_CLIPPED);
				continue;
			}

		}
		//clip y plane
		if (IS_BIT(clipMode, RENDERPIPING_CLIP_Y)) {
			inNum = 0;
			tanTh = camera->viewPlaneHeight / camera->viewDist;
			for (int vindex = 0;vindex <3;vindex++) {
				testG = tanTh*currPoly->transVList[vindex].z;
				testL = -testG;
				if (currPoly->transVList[vindex].y>testG) {
					clipCode[vindex] = CLIP_CODE_G;
				}
				else if (currPoly->transVList[vindex].y < testL)
				{
					clipCode[vindex] = CLIP_CODE_L;
				}
				else
				{
					clipCode[vindex] = CLIP_CODE_I;
					inNum++;
				}
			}

			if (inNum == 0) {
				SET_BIT(currPoly->state, POLY4D_STATE_CLIPPED);
				continue;
			}

		}

		if (IS_BIT(clipMode, RENDERPIPING_CLIP_Z)) {
			inNum = 0;
			testL = camera->clipNearZ;
			testG = camera->clipFarZ;
			for (int vindex = 0;vindex<3;vindex++) {
				if (currPoly->transVList[vindex].z>testG) {
					clipCode[vindex] = CLIP_CODE_G;
				}
				else if (currPoly->transVList[vindex].z < testL)
				{
					clipCode[vindex] = CLIP_CODE_L;
				}
				else
				{
					clipCode[vindex] = CLIP_CODE_I;
					inNum++;
				}
			}

			if (inNum == 0) {
				SET_BIT(currPoly->state, POLY4D_STATE_CLIPPED);
				continue;
			}
			else if ((clipCode[0] == CLIP_CODE_L) || (clipCode[1] == CLIP_CODE_L) || (clipCode[2] == CLIP_CODE_L))
			{
				if (inNum == 1) {

					if (clipCode[0] == CLIP_CODE_I) {
						v0 = 0;
						v1 = 1;
						v2 = 2;
					}
					else if (clipCode[1] == CLIP_CODE_I)
					{
						v0 = 1;
						v1 = 2;
						v2 = 0;
					}
					else {
						v0 = 2;
						v1 = 0;
						v2 = 1;
					}

					if (clipCode[v1] == CLIP_CODE_L&&clipCode[v2]== CLIP_CODE_L) {
						//right
						rightFraction = (camera->clipNearZ - currPoly->transVList[v0].z) / (currPoly->transVList[v1].z - currPoly->transVList[v0].z);
						currPoly->transVList[v1].x = currPoly->transVList[v0].x + (currPoly->transVList[v1].x - currPoly->transVList[v0].x)*rightFraction;
						currPoly->transVList[v1].y = currPoly->transVList[v0].y + (currPoly->transVList[v1].y - currPoly->transVList[v0].y)*rightFraction;
						currPoly->transVList[v1].z = camera->clipNearZ;

						//left
						leftFraction = (camera->clipNearZ - currPoly->transVList[v0].z) / (currPoly->transVList[v2].z - currPoly->transVList[v0].z);
						currPoly->transVList[v2].x = currPoly->transVList[v0].x + (currPoly->transVList[v2].x - currPoly->transVList[v0].x) *leftFraction;
						currPoly->transVList[v2].y = currPoly->transVList[v0].y + (currPoly->transVList[v2].y - currPoly->transVList[v0].y)*leftFraction;
						currPoly->transVList[v2].z = camera->clipNearZ;

						if (IS_BIT(currPoly->attr, POLY4D_ATTR_SHADE_MODE_GOURAUD)) {
							//phong  may be so slow when many clip traingle?
							//right 
							currPoly->transVList[v1].nx = currPoly->transVList[v0].nx+ (currPoly->transVList[v1].nx - currPoly->transVList[v0].nx)*rightFraction;
							currPoly->transVList[v1].ny = currPoly->transVList[v0].ny + (currPoly->transVList[v1].ny - currPoly->transVList[v0].ny)*rightFraction;
							currPoly->transVList[v1].nz = currPoly->transVList[v0].nz + (currPoly->transVList[v1].nz - currPoly->transVList[v0].nz)*rightFraction;
							Vector4D_Normalize(&currPoly->transVList[v1].n);

							//left
							currPoly->transVList[v2].nx = currPoly->transVList[v0].nx + (currPoly->transVList[v2].nx - currPoly->transVList[v0].nx)*leftFraction;
							currPoly->transVList[v2].ny = currPoly->transVList[v0].ny + (currPoly->transVList[v2].ny - currPoly->transVList[v0].ny)*leftFraction;
							currPoly->transVList[v2].nz = currPoly->transVList[v0].nz + (currPoly->transVList[v2].nz - currPoly->transVList[v0].nz)*leftFraction;
							Vector4D_Normalize(&currPoly->transVList[v2].n);
						}

						if (IS_BIT(currPoly->attr, POLY4D_ATTR_SHADE_MODE_TEXTURE)) {
							//right
							currPoly->transVList[v1].u = currPoly->transVList[v0].u + (currPoly->transVList[v1].u - currPoly->transVList[v0].u)*rightFraction;
							currPoly->transVList[v1].v = currPoly->transVList[v0].v + (currPoly->transVList[v1].v - currPoly->transVList[v0].v)*rightFraction;
							//left
							currPoly->transVList[v2].u = currPoly->transVList[v0].u + (currPoly->transVList[v2].u - currPoly->transVList[v0].u)*leftFraction;
							currPoly->transVList[v2].v = currPoly->transVList[v0].v + (currPoly->transVList[v2].v - currPoly->transVList[v0].v)*leftFraction;
						}
					}
					else
					{


						tmpPoly.init(currPoly);
						//right

						rightFraction = (camera->clipNearZ - currPoly->transVList[v1].z) / (currPoly->transVList[v2].z - currPoly->transVList[v1].z);
						insertRightP.x = currPoly->transVList[v1].x + (currPoly->transVList[v2].x - currPoly->transVList[v1].x)*rightFraction;
						insertRightP.y = currPoly->transVList[v1].y + (currPoly->transVList[v2].y - currPoly->transVList[v1].y)*rightFraction;
						insertRightP.z = camera->clipNearZ;
						//left
						leftFraction = (camera->clipNearZ - currPoly->transVList[v0].z) / (currPoly->transVList[v2].z - currPoly->transVList[v0].z);
						insertLeftP.x = currPoly->transVList[v0].x + (currPoly->transVList[v2].x - currPoly->transVList[v0].x)*leftFraction;
						insertLeftP.y = currPoly->transVList[v0].y + (currPoly->transVList[v2].y - currPoly->transVList[v0].y)*leftFraction;
						insertLeftP.z = camera->clipNearZ;

						currPoly->transVList[v2].setPosition(&insertRightP);
						tmpPoly.transVList[v1].setPosition(&insertRightP);
						tmpPoly.transVList[v2].setPosition(&insertLeftP);


						if (IS_BIT(currPoly->state, POLY4D_ATTR_SHADE_MODE_TEXTURE)) {
							//right
							insertRightT.x = currPoly->transVList[v1].u + (currPoly->transVList[v2].u - currPoly->transVList[v1].u)*rightFraction;
							insertRightT.y = currPoly->transVList[v1].v + (currPoly->transVList[v2].v - currPoly->transVList[v1].v)*rightFraction;
							//left
							insertLeftT.x = currPoly->transVList[v0].u + (currPoly->transVList[v2].u - currPoly->transVList[v0].u)*leftFraction;
							insertLeftT.y = currPoly->transVList[v0].v + (currPoly->transVList[v2].v - currPoly->transVList[v0].v)*leftFraction;

							currPoly->transVList[v2].setTexutre(&insertRightT);
							tmpPoly.transVList[v1].setTexutre(&insertRightT);
							tmpPoly.transVList[v2].setTexutre(&insertLeftT);
						}
						renderList.insertPoly4D(&tmpPoly);
					}
				}
				else if (inNum == 2)
				{
					if (clipCode[0] == CLIP_CODE_L) {
						v0 = 0;
						v1 = 1;
						v2 = 2;
					}
					else if (clipCode[1] == CLIP_CODE_L)
					{
						v0 = 1;
						v1 = 2;
						v2 = 0;
					}
					else {
						v0 = 2;
						v1 = 0;
						v2 = 1;
					}


					tmpPoly.init(currPoly);
					//right	
					rightFraction = (camera->clipNearZ - currPoly->transVList[v0].z) / (currPoly->transVList[v2].z - currPoly->transVList[v0].z);
					insertRightP.x = currPoly->transVList[v0].x + (currPoly->transVList[v2].x - currPoly->transVList[v0].x)*rightFraction;
					insertRightP.y = currPoly->transVList[v0].y + (currPoly->transVList[v2].y - currPoly->transVList[v0].y)*rightFraction;
					insertRightP.z = camera->clipNearZ;

					//left
					leftFraction = (camera->clipNearZ - currPoly->transVList[v0].z) / (currPoly->transVList[v1].z - currPoly->transVList[v0].z);
					insertLeftP.x = currPoly->transVList[v0].x + (currPoly->transVList[v1].x - currPoly->transVList[v0].x)*leftFraction;
					insertLeftP.y = currPoly->transVList[v0].y + (currPoly->transVList[v1].y - currPoly->transVList[v0].y)*leftFraction;
					insertLeftP.z = camera->clipNearZ;

					currPoly->transVList[v0].setPosition(&insertRightP);
					tmpPoly.transVList[v2].setPosition(&insertRightP);
					tmpPoly.transVList[v0].setPosition(&insertLeftP);

					if (IS_BIT(currPoly->state, POLY4D_ATTR_SHADE_MODE_GOURAUD)) {
						//phong argorithm , may be so slow when many clip traingle?
						//right
						insertRightN.x = currPoly->transVList[v0].nx + (currPoly->transVList[v2].nx - currPoly->transVList[v0].nx)*rightFraction;
						insertRightN.y= currPoly->transVList[v0].ny + (currPoly->transVList[v2].ny - currPoly->transVList[v0].ny)*rightFraction;
						insertRightN.z = currPoly->transVList[v0].nz + (currPoly->transVList[v2].nz - currPoly->transVList[v0].nz)*rightFraction;

						//right
						insertLeftN.x = currPoly->transVList[v0].nx + (currPoly->transVList[v1].nx - currPoly->transVList[v0].nx)*leftFraction;
						insertLeftN.y = currPoly->transVList[v0].ny + (currPoly->transVList[v1].ny - currPoly->transVList[v0].ny)*leftFraction;
						insertLeftN.z = currPoly->transVList[v0].nz + (currPoly->transVList[v1].nz - currPoly->transVList[v0].nz)*leftFraction;
						
						Vector4D_Normalize(&insertRightN); //this cost a lot of cycles maybe;
						Vector4D_Normalize(&insertLeftN);

						currPoly->transVList[v0].setNormal(&insertRightN);
						tmpPoly.transVList[v2].setNormal(&insertRightN);
						tmpPoly.transVList[v0].setNormal(&insertLeftN);
					}

					if (IS_BIT(currPoly->state, POLY4D_ATTR_SHADE_MODE_TEXTURE)) {
						//right
						insertRightT.x = currPoly->transVList[v0].u + (currPoly->transVList[v2].u - currPoly->transVList[v0].u)*rightFraction;
						insertRightT.y = currPoly->transVList[v0].v + (currPoly->transVList[v2].v - currPoly->transVList[v0].v)*rightFraction;
						//left
						insertLeftT.x = currPoly->transVList[v0].u + (currPoly->transVList[v1].u - currPoly->transVList[v0].u)*leftFraction;
						insertLeftT.y = currPoly->transVList[v0].v + (currPoly->transVList[v1].v - currPoly->transVList[v0].v)*leftFraction;

						currPoly->transVList[v0].setTexutre(&insertRightT);
						tmpPoly.transVList[v2].setTexutre(&insertRightT);
						tmpPoly.transVList[v0].setTexutre(&insertLeftT);
					}
					renderList.insertPoly4D(&tmpPoly);
				}
			}
		}

	}
}

void RenderPiping::lightPoly()
{
	uint32 rBase, gBase, bBase;
	uint32 rTmp, gTmp, bTmp,
		   rSum, gSum, bSum,
		   rVSum[3], gVSum[3], bVSum[3];
	int intensity;
	Vector4D l, v, r;//light vector,view vector,reflection vector
	float lLength, vLength, rLength, decay;
	float dotLN, dotVR, dotLDir, cosThNx;





	PolySelfTriangle4D *poly;
	for (int pi = 0;pi < renderList.num;pi++) {
		poly = renderList.polyPtrList[pi];
		if (IS_BIT(poly->state, POLY4D_STATE_UNACTIVE) || IS_BIT(poly->state, POLY4D_STATE_BACKFACE) || IS_BIT(poly->state, POLY4D_STATE_CLIPPED) || IS_BIT(poly->state, POLY4D_STATE_LIT)) {
			continue;
		}

		SET_BIT(poly->state, POLY4D_STATE_LIT);




		rSum = 0;
		gSum = 0;
		bSum = 0;

		rVSum[0] = 0;
		gVSum[0] = 0;
		bVSum[0] = 0;


		rVSum[1] = 0;
		gVSum[1] = 0;
		bVSum[1] = 0;


		rVSum[2] = 0;
		gVSum[2] = 0;
		bVSum[2] = 0;


		if (!IS_BIT(poly->attr, POLY4D_ATTR_SHADE_MODE_TEXTURE)) {
			
			rBase = poly->color.r;
			gBase = poly->color.g;
			bBase = poly->color.b;
		}
		else
		{	//as white 
			rBase = 255;
			gBase = 255;
			bBase = 255;
		}

		if (IS_BIT(poly->attr, POLY4D_ATTR_SHADE_MODE_FLAT)) {


			for (Light * light : lightVector) {
				//Iambient=RSambient *ISambient=CSambient
				//Ilambertian=RSdifusee*IS*MAX(l*n,0);
				//Ispecular=RSspecular*IS*MAX(r*v,0)*  (n*l>1?1:0)
				//IS=(I0/kc+kl*d+kq*d2)*c     I0*c=Color  
				//r=a-2(a*n)n;  a=-l  ==> r=(l*n)n-l;


				if (light->attr == LIGHT_ATTR_AMBIENT) {
					rSum += (rBase*light->ambient.r) >> 8;
					gSum += (gBase*light->ambient.g) >> 8;
					bSum += (bBase*light->ambient.b) >> 8;

				}
				else if (light->attr == LIGHT_ATTR_INFINITY) {
					//no decay

					Vector4D_Copy(&l, &light->tDir);
					dotLN = Vector4D_Dot(&l, &poly->transN);
					
					
					if (dotLN > 0) {
						//lambertian
						intensity = (128.0f * dotLN); //for translate integer and avoid float argorithm ,assure 2~7;
						rSum += (rBase*light->diffuse.r*intensity) >> 15;
						gSum += (gBase*light->diffuse.g*intensity) >> 15;
						bSum += (bBase*light->diffuse.b*intensity) >> 15;

						//scpecular   infinity too strong,so ignore?

						Vector4D_Copy(&r, &poly->transN);
						Vector4D_Mul_Scalar(2 * dotLN, &r);
						Vector4D_Sub(&r, &l);//get r no normalize
						Vector4D_Build(&poly->transVList[0].p, &camera->pos, &v);//get v no normalize
						dotVR = Vector4D_Dot(&v, &r);
						
						if (dotVR > 0) {
							rLength = Vector4D_Length(&r); 
							vLength = Vector4D_Length(&v); 
							intensity = (128.0f * dotVR) / (vLength*rLength);
							rSum += (rBase*light->specular.r*intensity) >>15;
							gSum += (gBase*light->specular.g*intensity) >> 15;
							bSum += (bBase*light->specular.b*intensity) >> 15;

						}
					}
				}
				else if (light->attr == LIGHT_ATTR_POINT) {
					//decay 1/kc+kl*d+kq*d2

					Vector4D_Build(&poly->transVList[0].p, &light->tPos, &l);
					dotLN = Vector4D_Dot(&l, &poly->transN);
					if (dotLN>0) {
						//lambertian
						lLength = Vector4D_Length(&l);
						decay = light->kc + light->kl*lLength + light->kq*lLength*lLength;
						intensity = (128 * dotLN) / (decay*lLength);//for translate integer and avoid float argorithm
						rSum += (rBase*light->diffuse.r*intensity)>>15;
						gSum += (gBase*light->diffuse.g*intensity)>>15;
						bSum += (bBase*light->diffuse.b*intensity)>>15;
						//specular

						Vector4D_Copy(&r, &poly->transN);
						Vector4D_Mul_Scalar(2 * dotLN, &r);
						Vector4D_Sub(&r, &l);//get r no normalize
						Vector4D_Build(&poly->transVList[0].p, &camera->pos, &v);//get v no normalize
						dotVR = Vector4D_Dot(&v, &r);
						if (dotVR > 0) {
							rLength = Vector4D_Length(&r); //real time argorithm 8%  error!
							vLength = Vector4D_Length(&v); //real time
							intensity = (128.0f * dotVR) / (decay*vLength*rLength);
							rSum += (rBase*light->specular.r*intensity) >> 15;
							gSum += (gBase*light->specular.g*intensity) >> 15;
							bSum += (bBase*light->specular.b*intensity) >> 15;
						}
					}
				}
				else if (light->attr == LIGHT_ATTR_SPOT)
				{
					//decay (1/kc+kl*d+kq*d2 )*cosTh  cotTh=dotLDir/lLength
					Vector4D_Build(&light->tPos, &poly->transVList[0].p, &l);//init l to point to vertex 
					dotLDir = Vector4D_Dot(&l, &light->tDir);//for get cosTh; 
					lLength = Vector4D_Length(&l);
					cosThNx = dotLDir / lLength;
					for (int pfN = 1;pfN < light->pf;pfN++) {
						cosThNx *= cosThNx;
					}
					if (dotLDir > 0) {
						Vector4D_Mul_Scalar(-1, &l); //reverse l to point to light source
						dotLN = Vector4D_Dot(&l, &poly->transN);


						if (dotLN > 0) {
							
							decay = light->kc + light->kl*lLength + light->kq*lLength*lLength;
							intensity = (128 * cosThNx* dotLN) / (decay*lLength);
							rSum += (rBase*light->diffuse.r*intensity) / (128 * 255);
							gSum += (gBase*light->diffuse.g*intensity) / (128 * 255);
							bSum += (bBase*light->diffuse.b*intensity) / (128 * 255);
						}

						Vector4D_Copy(&r, &poly->transN);
						Vector4D_Mul_Scalar(2 * dotLN, &r);
						Vector4D_Sub(&r, &l);//get r no normalize
						Vector4D_Build(&poly->transVList[0].p, &camera->pos, &v);//get v no normalize
						dotVR = Vector4D_Dot(&v, &r);

						if (dotVR > 0) {
							rLength = Vector4D_Length(&r); //real time argorithm 8%  error!
							vLength = Vector4D_Length(&v); //real time
							intensity = (128 * cosThNx * dotVR )/(decay*vLength*rLength);
							rSum += (rBase*light->specular.r*intensity) >> 15;
							gSum += (gBase*light->specular.g*intensity) >> 15;
							bSum += (bBase*light->specular.b*intensity) >> 15;
						}
					}
				}
			}

			rSum = MIN(rSum, 255);
			gSum = MIN(gSum, 255);
			bSum = MIN(bSum, 255);
		
			poly->lightColor[0].r = rSum;
			poly->lightColor[0].g = gSum;
			poly->lightColor[0].b = bSum;

		}
		else if (IS_BIT(poly->attr, POLY4D_ATTR_SHADE_MODE_GOURAUD)) {
			for (Light * light : lightVector) {
				//Iambient=RSambient *ISambient=CSambient
				//Ilambertian=RSdifusee*IS*MAX(l*n,0);
				//Ispecular=RSspecular*IS*MAX(r*v,0)*  (n*l>1?1:0)
				//IS=(I0/kc+kl*d+kq*d2)*c     I0*c=Color  
				//r=a-2(a*n)n;  a=-l  ==> r=(l*n)n-l;


				if (light->attr == LIGHT_ATTR_AMBIENT) {
					rTmp = (rBase*light->ambient.r)>>8;
					gTmp = (gBase*light->ambient.g)>>8;
					bTmp = (bBase*light->ambient.b)>>8;

					rVSum[0] += rTmp;
					gVSum[0] += gTmp;
					bVSum[0] += bTmp;

					rVSum[1] += rTmp;
					gVSum[1] += gTmp;
					bVSum[1] += bTmp;


					rVSum[2] += rTmp;
					gVSum[2] += gTmp;
					bVSum[2] += bTmp;

				}
				else if (light->attr == LIGHT_ATTR_INFINITY) {
					//no decay

					Vector4D_Copy(&l, &light->tDir);

					for (int vi = 0;vi< 3;vi++) {

						dotLN = Vector4D_Dot(&l, &poly->transVList[vi].n);
						if (dotLN > 0) {
							//lambertian
							intensity = (128 * dotLN); //for translate integer and avoid float argorithm ,assure 2~7;



							rVSum[vi] += (rBase*light->diffuse.r*intensity) >> 15;
							gVSum[vi] += (gBase*light->diffuse.g*intensity) >> 15;
							bVSum[vi] += (bBase*light->diffuse.b*intensity) >> 15;

							//scpecular   infinity too strong,so ignore?

							Vector4D_Copy(&r, &poly->transVList[vi].n);
							Vector4D_Mul_Scalar(2 * dotLN, &r);
							Vector4D_Sub(&r, &l);//get r no normalize
							Vector4D_Build(&poly->transVList[vi].p, &camera->pos, &v);//get v no normalize
							dotVR = Vector4D_Dot(&v, &r);
							if (dotVR > 0) {
								rLength = Vector4D_Length(&r); //real time argorithm 8%  error!
								vLength = Vector4D_Length(&v); //real time
								intensity = (128 * dotVR) / (vLength*rLength);
								rVSum[vi] += (rBase*light->specular.r*intensity) >> 15;
								gVSum[vi] += (gBase*light->specular.g*intensity) >> 15;
								bVSum[vi] += (bBase*light->specular.b*intensity) >> 15;
							}
						}
					}
				}
				else if (light->attr == LIGHT_ATTR_POINT) {
					//decay 1/kc+kl*d+kq*d2
					for (int vi = 0;vi < 3;vi++) {
						Vector4D_Build(&poly->transVList[vi].p, &light->tPos, &l);
						dotLN = Vector4D_Dot(&l, &poly->transVList[vi].n);
						if (dotLN > 0) {
							//lambertian
							lLength = Vector4D_Length(&l);
							decay = light->kc + light->kl*lLength + light->kq*lLength*lLength;
							intensity = (128.0f * dotLN) / (decay*lLength);//for translate integer and avoid float argorithm
							rVSum[vi] += (rBase*light->diffuse.r*intensity) >> 15;
							gVSum[vi] += (gBase*light->diffuse.g*intensity) >> 15;
							bVSum[vi] += (bBase*light->diffuse.b*intensity) >> 15;

							//specular
							Vector4D_Copy(&r, &poly->transVList[vi].n);
							Vector4D_Mul_Scalar(2 * dotLN, &r);
							Vector4D_Sub(&r, &l);//get r no normalize
							Vector4D_Build(&poly->transVList[vi].p, &camera->pos, &v);//get v no normalize
							dotVR = Vector4D_Dot(&v, &r);
							if (dotVR > 0) {
								rLength = Vector4D_Length(&r); //real time argorithm 8%  error!
								vLength = Vector4D_Length(&v); //real time
								intensity = (128.0f * dotVR)/ (decay*vLength*rLength);
								rVSum[vi] += (rBase*light->specular.r*intensity) >> 15;
								gVSum[vi] += (gBase*light->specular.g*intensity) >> 15;
								bVSum[vi] += (bBase*light->specular.b*intensity) >> 15;
							}
						}
					}
				}
				else if (light->attr == LIGHT_ATTR_SPOT)
				{

					for (int vi = 0;vi < 3;vi++) {
						//decay (1/kc+kl*d+kq*d2 )*cosTh  cotTh=dotLDir/lLength
						Vector4D_Build(&light->tPos, &poly->transVList[vi].p, &l);//init l to point to vertex 
						dotLDir = Vector4D_Dot(&l, &light->tDir);//for get cosTh; 
						lLength = Vector4D_Length(&l);
						cosThNx = dotLDir / lLength;
						for (int pfN = 1;pfN < light->pf;pfN++) {
							cosThNx *= cosThNx;
						}

						if (dotLDir > 0) {
							Vector4D_Mul_Scalar(-1, &l); //reverse l to point to light source
							dotLN = Vector4D_Dot(&l, &poly->transVList[vi].n);
							if (dotLN > 0) {
								
								decay = light->kc + light->kl*lLength + light->kq*lLength*lLength;
								intensity = (128.0f * cosThNx* dotLN) / (decay*lLength);
								rVSum[vi] += (rBase*light->diffuse.r*intensity) >> 15;
								gVSum[vi] += (gBase*light->diffuse.g*intensity) >> 15;
								bVSum[vi] += (bBase*light->diffuse.b*intensity) >> 15;
							}

							Vector4D_Copy(&r, &poly->transVList[vi].n);
							Vector4D_Mul_Scalar(2 * dotLN, &r);
							Vector4D_Sub(&r, &l);//get r no normalize
							Vector4D_Build(&poly->transVList[vi].p, &camera->pos, &v);//get v no normalize
							dotVR = Vector4D_Dot(&v, &r);
							if (dotVR > 0) {
								rLength = Vector4D_Length(&r); //real time argorithm 8%  error! 
								vLength = Vector4D_Length(&v); //real time
								intensity = (128.0f * cosThNx* dotVR) /(decay*vLength*rLength);
								rVSum[vi] += (rBase*light->specular.r*intensity) >> 15;
								gVSum[vi] += (gBase*light->specular.g*intensity) >> 15;
								bVSum[vi] += (bBase*light->specular.b*intensity) >> 15;
							}
						}
					}
				}
			}
			for (int vi = 0;vi<3;vi++) {
				rVSum[vi] = MIN(rVSum[vi], 255);
				gVSum[vi] = MIN(gVSum[vi], 255);
				bVSum[vi] = MIN(bVSum[vi], 255);

				poly->lightColor[vi].r = rVSum[vi];
				poly->lightColor[vi].g = gVSum[vi];
				poly->lightColor[vi].b = bVSum[vi];
			}
		}
		else {
			//constant shade
			poly->lightColor[0].r = rBase;
			poly->lightColor[0].g = gBase;
			poly->lightColor[0].b = bBase;
		}

	}
}

void RenderPiping::transformPerspective()
{
	Vector4D tmpV;
	float divZ;
	for (int i = 0;i < renderList.num;i++) {
		PolySelfTriangle4D *poly = renderList.polyPtrList[i];
		if (IS_BIT(poly->state, POLY4D_STATE_UNACTIVE) || IS_BIT(poly->state, POLY4D_STATE_BACKFACE) || IS_BIT(poly->state, POLY4D_STATE_CLIPPED)) {
			continue;
		}
		for (int j = 0;j < 3;j++) {
			divZ = camera->viewDist / poly->transVList[j].p.z;
			poly->transVList[j].p.x *= divZ;
			poly->transVList[j].p.y = poly->transVList[j].p.y*camera->aspectRatio*divZ;
		}

	}
}

void RenderPiping::transformScreen()
{
	Vector4D tmpV;
	float divZ;
	float a, b;
	a = 0.5f*camera->viewPortWidth - 0.5f;//width and  port centerX
	b = 0.5f*camera->viewPortHeight - 0.5f;//height and  port centerY
	for (int i = 0;i < renderList.num;i++) {
		PolySelfTriangle4D *poly = &renderList.polyDataList[i];
		if (IS_BIT(poly->state, POLY4D_STATE_UNACTIVE) || IS_BIT(poly->state, POLY4D_STATE_CLIPPED) || IS_BIT(poly->state, POLY4D_STATE_BACKFACE)) {
			continue;
		}
		for (int j = 0;j < 3;j++) {
			poly->transVList[j].p.x = a + a*poly->transVList[j].p.x;
			poly->transVList[j].p.y = b - b*poly->transVList[j].p.y;
		}
	}
}

void RenderPiping::transformPerspectiveScreen()
{
	Vector4D tmpV;
	float divZ;
	float a, b;
	a = 0.5f*camera->viewPortWidth - 0.5f;
	b = 0.5f*camera->viewPortHeight - 0.5f;
	for (int i = 0;i < renderList.num;i++) {
		PolySelfTriangle4D *poly = &renderList.polyDataList[i];
		if (IS_BIT(poly->state, POLY4D_STATE_UNACTIVE) || IS_BIT(poly->state, POLY4D_STATE_CLIPPED) || IS_BIT(poly->state, POLY4D_STATE_BACKFACE)) {
			continue;
		}
		for (int j = 0;j < 3;j++) {
			divZ = camera->viewDist / poly->transVList[j].p.z;
			poly->transVList[j].p.x = a + a*poly->transVList[j].p.x*divZ;
			poly->transVList[j].p.y = b - a*poly->transVList[j].p.y*divZ;
		}

	}
}

RenderPiping::RenderPiping()
{
	lightVector.clear();
}

RenderPiping::~RenderPiping()
{
}

void RenderPiping::render(Canvas * canvas)
{

	clipPoly();
	lightPoly();
	transformPerspective();
	transformScreen();

	if (isSortZ) {
		renderList.sortByZ();
	}

	PolySelfTriangle4D *poly;
	for (int i = 0;i < renderList.num;i++) {
		poly = renderList.polyPtrList[i];
		if (IS_BIT(poly->state, POLY4D_STATE_UNACTIVE) || IS_BIT(poly->state, POLY4D_STATE_BACKFACE) || IS_BIT(poly->state, POLY4D_STATE_CLIPPED)) {
			continue;
		}
		canvas->drawPolyTriangle(poly);
	}
}

void RenderPiping::clear()
{
	renderList.clear();
	lightVector.clear();
}

void RenderPiping::setCamera(Camera * camera)
{
	this->camera = camera;
}

int RenderPiping::sendObject(Object *object)
{
	CLEAR_BIT(object->state, OBJECT_STATE_CULLED);
	for (int i = 0;i < object->polygonsNum;i++) {
		CLEAR_BIT(object->polyList[i].state, POLY4D_STATE_LIT | POLY4D_STATE_BACKFACE | POLY4D_STATE_CLIPPED);
	}
	
	worldTransfrom(object);
	cull(object);
	removeBack(object);
	cameraTransform(object);
	if (!renderList.insertObject(object)) {
		return 0;
	}
	return 1;
}

void RenderPiping::addLight(Light * light)
{	
	cameraTransform(light);
	lightVector.push_back(light);
}

void RenderPiping::setSortZEnable(int enable)
{
	isSortZ = enable;
}




