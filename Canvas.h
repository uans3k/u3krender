#pragma once
#include "ddraw.h"
#include "CustomDefine.h"
#include "Math3D.h"
#include "GraphicsStruct.h"
#include "Bitmap.h"

struct ZBuffer
{
	int pitch; //line byte;
	int width; //width pix
	int height;//height pix
	int size;//pitch*height;
	uint32 *zBuffer=NULL;
};

class Canvas {
#define CANVAS_TRA_TYPE_GENERAL 0x00
#define CANVAS_TRA_TYPE_BOTTOM 0x01
#define CANVAS_TRA_TYPE_TOP 0x02
#define CANVAS_TRA_TYPE_FLAT 0x03

#define CANVAS_SAMPLE_POINT 0 
#define CANVAS_SAMPLE_BILINEAR 1

public:
	Canvas();
	~Canvas();
	int reset();
	int setBuffer(uint8 *buffer, int width, int height, int pitch, int bpp);
	int setBufferAndZbuffer(uint8 *buffer, int width, int height, int pitch, int bpp, ZBuffer *zBuffer);


	int setClipper(int left, int top, int right, int bottom);
	int clearClipper();


	int drawTest16(uint16 color);

	/*	
		use 1/z buffer and perspecte adjust 
		wo assume that z<2~12 4012 (it's enough always be clipped) and z>0.5 that mean tan (th/2)<2 ,view angle th<126.86бу

		use pre compute u*width,v*height for more flexible later //wo can do this offline too but now wo do real time;
		use fixp19 for x,y,z ,because of  wo use left-top to avoid 1/z buffer write one pix (dzdx interpolation cause 2 triangle collision,that's too bad),wo need more accuracy;
		use fixp20 for r/z,g/z,b/z  so integer  2~11 2048 ,is't enough and 2~8 color accuracy ,if dc<2~8  per pix ignore (no change),that mean 1 bit change but dx or dy >2`8 256 (and z it' too long),think it don't change reasonablely. 
		use fixp20 for u/z,v/z  u,v max 2`11 2048,wo think enough,and 2~8 uv accuracy ,if du<2~8  per pix ignore (no change),that mean 1 bit change but dx or dy >2`8 256 (and z it' too long),think it don't change reasonablely. 
		use fixp28 for invz  mean that invz<8      z>0.125 ,for normal pespective ,tan (th/2)<2 ,assume z>0.5,view angle th<126.86бу
				and when wo fixp28 >>8 if z=2`12 ,it's also has 2`8 accuraracy,so all 2`8 accuracy at least(z too long) ,well! 
		 
	*/
	int drawPolyTriangle(PolySelfTriangle4D *poly);

	int drawTriangleTextureGouraudInvZ16(float x0, float y0, float z0, float u0, float v0, Color c0,
		float x1, float y1, float z1, float u1, float v1, Color c1,
		float x2, float y2, float z2, float u2, float v2, Color c2,
		uint16 *texture, int width, int height, int textPitch);

	int drawTopTriangleTextureGouraudInvZ16(float x0, float y0, float z0, float u0, float v0,Color c0,
		float x1, float y1, float z1, float u1, float v1, Color c1,
		float x2, float y2, float z2, float u2, float v2, Color c2,
		uint16 *texture, int width, int height, int textPitch);

	int drawBottomTriangleTextureGouraudInvZ16(float x0, float y0, float z0, float u0, float v0, Color c0,
		float x1, float y1, float z1, float u1, float v1, Color c1,
		float x2, float y2, float z2, float u2, float v2, Color c2,
		uint16 *texture, int width, int height, int textPitch);


	int drawTriangleTextureFlatInvZ16(float x0, float y0, float z0, float u0, float v0,
		float x1, float y1, float z1, float u1, float v1,
		float x2, float y2, float z2, float u2, float v2,
		Color c,
		uint16 *texture, int width, int height, int textPitch);

	int drawTopTriangleTextureFlatInvZ16(float x0, float y0, float z0, float u0, float v0,
								      float x1, float y1, float z1, float u1, float v1,
									  float x2, float y2, float z2, float u2, float v2,
									  Color c,
									  uint16 *texture, int width, int height, int textPitch);


	int drawBottomTriangleTextureFlatInvZ16(float x0, float y0, float z0, float u0, float v0,
		float x1, float y1, float z1, float u1, float v1,
		float x2, float y2, float z2, float u2, float v2,
		Color c,
		uint16 *texture, int width, int height, int textPitch);


	int drawTriangleTextureInvZ16(float x0, float y0, float z0,float u0,float v0,
								  float x1, float y1, float z1, float u1,float v1,
								  float x2, float y2, float z2, float u2, float v2,
								  uint16 *texture,int width,int height,int textPitch);

	int drawTopTriangleTextureInvZ16(float x0, float y0, float z0, float u0, float v0,
		float x1, float y1, float z1, float u1, float v1,
		float x2, float y2, float z2, float u2, float v2,
		uint16 *texture, int width, int height, int textPitch);

	int drawBottomTriangleTextureInvZ16(float x0, float y0, float z0, float u0, float v0,
		float x1, float y1, float z1, float u1, float v1,
		float x2, float y2, float z2, float u2, float v2,
		uint16 *texture, int width, int height, int textPitch);


	int drawTriangleGouraudInvZ16(float x0, float y0, float z0, Color c0,
									float x1, float y1, float z1, Color c1,
									float x2, float y2, float z2, Color c2);

	int drawTopTriangleGouraudInvZ16(float x0, float y0, float z0, Color c0,
									 float x1, float y1, float z1, Color c1,
									 float x2, float y2, float z2, Color c2);

	int drawBottomTriangleGouraudInvZ16(float x0, float y0, float z0, Color c0,
		float x1, float y1, float z1, Color c1,
		float x2, float y2, float z2, Color c2);


	int drawPolyTriangleSolidInvZ16(PolySelfTriangle4D *poly);

	int drawTriangleSolidInvZ16(float x0, float y0, float z0,
								float x1, float y1, float z1,
								float x2, float y2, float z2,
								Color c);

	int drawTopTriangleSolidInvZ16(float x0, float y0, float z0, 
								 float x1, float y1, float z1,
								 float x2, float y2,float z2,
								 Color c);
	int drawBottomTriangleSolidInvZ(float x0, float y0, float z0,
								   float x1, float y1, float z1,
								   float x2, float y2, float z2,
								   Color c);

	/*	
		for accelorate , color texture uv .etc, use fixp16 to avoid float to int

		texture size must be 2`n remember , that also make pitch non fill (because of 16*2^n =32*2`(n-1)   not need fill)
	*/	
	int drawPolyTriangleGouraudDirect16(PolySelfTriangle4D *poly);
	int drawTriangleGouraudDirect16(float x0, float y0, float x1, float y1, float x2, float y2, Color c0, Color c1, Color c2);



	int drawPolyTriangleTextureGouraud16(PolySelfTriangle4D *poly);
	int drawTriangleTextureGouraud16(float x0, float y0, float u0, float v0, Color c0,
		float x1, float y1, float u1, float v1, Color c1,
		float x2, float y2, float u2, float v2, Color c2,
		uint16* textBuffer, int width, int height, int pitch);

	int drawTopTriangleTextureGouraud16(float x0, float y0, float u0, float v0,Color c0,
		float x1, float y1, float u1, float v1, Color c1,
		float x2, float y2, float u2, float v2, Color c2,
		uint16* textBuffer, int width, int height, int pitch);
	int drawBottomTriangleTextureGouraud16(float x0, float y0, float u0, float v0, Color c0,
		float x1, float y1, float u1, float v1, Color c1,
		float x2, float y2, float u2, float v2, Color c2,
		uint16* textBuffer, int width, int height, int pitch);


	int drawPolyTriangleTextureFlat16(PolySelfTriangle4D *poly);
	int drawTriangleTextureFlat16(float x0, float y0, float u0, float v0,
								float x1, float y1, float u1, float v1,
								float x2, float y2, float u2, float v2,
								Color c,
								uint16* textBuffer, int width, int height, int pitch);

	int drawTopTriangleTextureFlat16(float x0, float y0, float u0, float v0, 
									float x1, float y1, float u1, float v1, 
									float x2, float y2, float u2, float v2,
									Color c,
									uint16* textBuffer, int width, int height, int pitch);
	int drawBottomTriangleTextureFlat16(float x0, float y0, float u0, float v0,
										float x1, float y1, float u1, float v1,
										float x2, float y2, float u2, float v2,
										Color c,
										uint16* textBuffer, int width, int height, int pitch);

	int drawPolyTriangleTexture16(PolySelfTriangle4D *poly);
	int drawTriangleTexture16(float x0, float y0, float u0, float v0, float x1, float y1, float u1, float v1, float x2, float y2, float u2, float v2, uint16* textBuffer, int width, int height, int pitch);
	int drawTopTriangleTexture16(float x0, float y0, float u0, float v0, float x1, float y1, float u1, float v1, float x2, float y2, float u2, float v2, uint16* textBuffer, int width, int height, int pitch);
	int drawBottomTriangleTexture16(float x0, float y0, float u0, float v0, float x1, float y1, float u1, float v1, float x2, float y2, float u2, float v2, uint16* textBuffer, int width, int height, int pitch);


	int drawPolyTriangleGouraud16(PolySelfTriangle4D *poly);
	int drawTriangleGouraud16(float x0, float y0, float x1, float y1, float x2, float y2, Color c0, Color c1, Color c2);
	int drawTopTriangleGouraud16(float x0, float y0, float x1, float y1, float x2, float y2, Color c0, Color c1, Color c2);
	int drawBottomTriangleGouraud16(float x0, float y0, float x1, float y1, float x2, float y2, Color c0, Color c1, Color c2);

	/*
		Float base draw,left top priciple,for 3d
	*/
	int drawPolyTriangleSolidF16(PolySelfTriangle4D *poly);
	int drawTriangleSolidF16(float x0, float y0, float x1, float y1, float x2, float y2, Color color);
	int drawTopTriangleF16(float x0, float y0, float x1, float y1, float x2, float y2,Color color);
	int drawBottomTriangleF16(float x0, float y0, float x1, float y1, float x2, float y2, Color color);

	/*
		Integer base draw,for 2d
	*/
	int drawHLine16(int x1, int x2, int y, Color c);
	int drawVLine16(int y1, int y2, int x, Color c);
	int drawLine16(int x1, int y1, int x2, int y2, Color c);
	int drawPolyTriangleSolidI16(PolySelfTriangle4D *poly);
	int drawTriangleWireI16(int x1, int y1, int x2, int y2, int x3, int y3, Color c);
	int drawTriangleSolidI16(int x1, int y1, int x2, int y2, int x3, int y3, Color c);
	
	//x1,y1 for  horn point
	//x2,y2,x3,y3 for  flat point
	
	int drawTopTriangleI16(int x1, int y1, int x2, int y2, int x3, int y3, Color c);
	int drawBottomTriangleI16(int x1, int y1, int x2, int y2, int x3, int y3, Color c);

	//draw bitmap
	int drawBitmap(float left, float top,float right,float bottom, Bitmap *bitmap,int sampleMode=CANVAS_SAMPLE_POINT);
	int drawBitmapPoint16(float left, float top, float right, float bottom,uint16 *buffer, int width, int height, int pitch);


private:
	int hasBuffer=0;
	int enableZbuffer = 0;
	long pitch = 0;
	int width = 0;
	int height = 0;
	int bpp = 0;
	RECT cliper = { 0,0,0,0 };
	ZBuffer *zBuffer=NULL;
	uint8 *drawBuffer = NULL;

	/*
	0 for trivially reject
	1 for trivailly accept or clip
	*/
	int clipLine(int& x1, int& y1, int& x2, int& y2);

};
