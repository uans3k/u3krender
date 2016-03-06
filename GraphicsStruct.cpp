#include "stdafx.h"
#include "GraphicsStruct.h"


void PolyTriangle4D::initColorPoly(Vertex4D * originVList, Vertex4D * transVList, int vertextIndex[], Vector4D  * originNList, Vector4D *transNList,int normalIndex,Color color,int shadeMode)
{
	state = 0;
	this->attr = POLY4D_ATTR_SHADE_MODE_GOURAUD;
	this->color = color;

	this->vertextIndex[0] = vertextIndex[0];
	this->vertextIndex[1] = vertextIndex[1];
	this->vertextIndex[2] = vertextIndex[2];

	this->normorIndex = normalIndex;

	this->originVList = originVList;
	this->transVList = transVList;

	this->originNList = originNList;
	this->transNList = transNList;

}

void PolyTriangle4D::initTexturePoly(Vertex4D * originVList, Vertex4D * transVList, int vertextIndex[], Vector4D * originNList, Vector4D * transNList, int normalIndex, Bitmap * texture, Point2D * uvList, int textureIndex[], int shadeMode)
{
	state = 0;
	attr = 0;
	SET_BIT(attr, shadeMode | POLY4D_ATTR_SHADE_MODE_TEXTURE);

	this->originVList = originVList;
	this->transVList = transVList;
	this->vertextIndex[0] = vertextIndex[0];
	this->vertextIndex[1] = vertextIndex[1];
	this->vertextIndex[2] = vertextIndex[2];



	this->normorIndex = normalIndex;
	this->originNList = originNList;
	this->transNList = transNList;

	this->texture = texture;
	this->uvList = uvList;
	this->textureIndex[0] = textureIndex[0];
	this->textureIndex[1] = textureIndex[1];
	this->textureIndex[2] = textureIndex[2];

}

void PolyTriangle4D::reset()
{
	CLEAR_BIT(state, POLY4D_STATE_BACKFACE | POLY4D_STATE_CLIPPED);
}


void PolySelfTriangle4D::initFromTrans(const PolyTriangle4D * poly4D)
{
	state = poly4D->state;
	attr = poly4D->attr;
	color = poly4D->color;
	lightColor[0] = poly4D->lightColor[0];
	lightColor[1] = poly4D->lightColor[1];
	lightColor[2] = poly4D->lightColor[2];

	material = poly4D->material;

	originVList[0]=poly4D->transVList[poly4D->vertextIndex[0]];
	originVList[1]=poly4D->transVList[poly4D->vertextIndex[1]];
	originVList[2] =poly4D->transVList[poly4D->vertextIndex[2]];

	transVList[0] = poly4D->transVList[poly4D->vertextIndex[0]];
	transVList[1] = poly4D->transVList[poly4D->vertextIndex[1]];
	transVList[2] = poly4D->transVList[poly4D->vertextIndex[2]];

	

	if (IS_BIT(poly4D->attr, POLY4D_ATTR_SHADE_MODE_TEXTURE)&&poly4D->texture&&poly4D->uvList) {
		this->texure = poly4D->texture;
		Point2D_Copy(&originVList[0].t, &poly4D->uvList[poly4D->textureIndex[0]]);
		Point2D_Copy(&originVList[1].t, &poly4D->uvList[poly4D->textureIndex[1]]);
		Point2D_Copy(&originVList[2].t, &poly4D->uvList[poly4D->textureIndex[2]]);
		Point2D_Copy(&transVList[0].t, &poly4D->uvList[poly4D->textureIndex[0]]);
		Point2D_Copy(&transVList[1].t, &poly4D->uvList[poly4D->textureIndex[1]]);
		Point2D_Copy(&transVList[2].t, &poly4D->uvList[poly4D->textureIndex[2]]);
		SET_BIT(originVList[0].attr, VERTEX4D_ATTR_TEXTURE);
		SET_BIT(originVList[1].attr, VERTEX4D_ATTR_TEXTURE);
		SET_BIT(originVList[2].attr, VERTEX4D_ATTR_TEXTURE);
		SET_BIT(transVList[0].attr, VERTEX4D_ATTR_TEXTURE);
		SET_BIT(transVList[0].attr, VERTEX4D_ATTR_TEXTURE);
		SET_BIT(transVList[0].attr, VERTEX4D_ATTR_TEXTURE);

	}
	if (poly4D->transNList) {
		Vector4D_Copy(&originN, &poly4D->transNList[poly4D->normorIndex]);
		Vector4D_Copy(&transN, &poly4D->transNList[poly4D->normorIndex]);
	}

}

void PolySelfTriangle4D::init(const PolySelfTriangle4D * poly)
{
	memcpy(this, poly, sizeof(PolySelfTriangle4D));
}

uint16 Color::getRGB566()
{
	int ir = r;
	int ig = g;
	int ib = b;
	return R_G_B_TO_RGB565(ir,ig,ib);
}


