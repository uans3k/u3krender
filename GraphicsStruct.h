#pragma once
#include "Math3D.h"
#include "Bitmap.h"

#define R_G_B_TO_RGB565_FIXP16(r,g,b) (  ((r) >> (FIXP16_SHIFT+3) << 11) + ((g)>> (FIXP16_SHIFT + 2)<<5) +( (b)>> (FIXP16_SHIFT + 3) )  )
#define R_G_B_TO_RGB565(r,g,b) ( ( (r) >> 3 << 11 ) + ( (g)>>2<<5 ) + ( (b)>> 3 )  )
#define R5_G6_B5_TO_RGB565(r,g,b) ( ( (r)<< 11 ) + ( (g)<<5 ) + (b)  )
#define RGB565_TO_R5_G6_B5(rgb,r,g,b) {(r)=((rgb)>>11);(g)=(((rgb)>>5)&0x3f);(b)=((rgb)&0x1f);}

struct Color
{	
	union
	{
		uint32 argb;                    
		uint8 argbM[4];             
		struct { uint8 b, g, r,a; };
	}; 
	uint16 getRGB566();
};

/*
left- hand , so  p0 p1 p2,left-hand cycle  ==> v1=p1-p0,v2=p2-p0 ==>   v1 x v2 left-hand  cyclejsd s
*/
#define VERTEX4D_ATTR_NULL             0x0000  
#define VERTEX4D_ATTR_POINT            0x0001
#define VERTEX4D_ATTR_NORMAL           0x0002
#define VERTEX4D_ATTR_TEXTURE          0x0004

#define VERTEX_FLAGS_OVERRIDE_MASK          0xf000 
#define VERTEX_FLAGS_OVERRIDE_CONSTANT      0x1000
#define VERTEX_FLAGS_OVERRIDE_FLAT          0x2000
#define VERTEX_FLAGS_OVERRIDE_GOURAUD       0x4000
#define VERTEX_FLAGS_OVERRIDE_TEXTURE       0x8000

struct Vertex4D
{
	union
	{
		float M[12];


		struct
		{
			float x, y, z, w;     // point
			float nx, ny, nz, nw; // normal 
			float u, v;       // texture 

			float intensity;           // final vIndex intensity after lighting
			int   attr;        // attributes/ extra texture coordinates
		};


		struct
		{
			Point4D  p;       // the point 
			Vector4D n;       // the normal
			Point2D  t;       // texture coordinates
		};

	};

	inline void init(Vertex4D *v) {
		memcpy(this, v, sizeof(Vertex4D));
	}

	inline void setPosition(Point4D *p) {
		Point4D_Copy(&this->p, p);
	}
	inline void setTexutre(Point2D *t) {
		Point2D_Copy(&this->t, t);
	}
	inline void setNormal(Vector4D *normal) {
		Vector4D_Copy(&this->n, normal);
	}
};


#define POLY4D_STATE_NONE             0x0000
#define POLY4D_STATE_UNACTIVE             0x0001  
#define POLY4D_STATE_CLIPPED            0x0002  
#define POLY4D_STATE_BACKFACE           0x0004  
#define POLY4D_STATE_LIT                0x0008

#define POLY4D_ATTR_2SIDED                0x0001
#define POLY4D_ATTR_TRANSPARENT           0x0002
#define POLY4D_ATTR_SHADE_MODE_CONSTANT   0x0004 
#define POLY4D_ATTR_SHADE_MODE_FLAT       0x0008
#define POLY4D_ATTR_SHADE_MODE_GOURAUD    0x0010
#define POLY4D_ATTR_SHADE_MODE_TEXTURE    0x0020 
#define POLY4D_ATTR_ENABLE_MATERIAL       0x0040 
#define POLY4D_ATTR_DISABLE_MATERIAL      0x0080 




struct PolyTriangle4D
{
	int state;           // state information
	int attr;            // physical attributes of polygon
	Color color;           // color of polygon
	Color lightColor[3];    // colors after lighting

	int      material;            // material index (-1) no material 

	//vertex
	Vertex4D *originVList; //origin vertext list
	Vertex4D *transVList; // transform vertext list 
	int vertextIndex[3];           // the indices into the vIndex list
	
	//normal
	Vector4D *originNList;
	Vector4D *transNList;
	int normorIndex;

	//texture
	Bitmap *texture;
	Point2D *uvList;
	int textureIndex[3];

	void initColorPoly(Vertex4D *originVList, Vertex4D *transVList, int vertextIndex[], Vector4D  * originNList, Vector4D *transNList, int normalIndex, Color color,int shadeMode= POLY4D_ATTR_SHADE_MODE_FLAT);
	void initTexturePoly(Vertex4D *originVList, Vertex4D *transVList, int vertextIndex[],
						 Vector4D  * originNList, Vector4D *transNList, int normalIndex,
						 Bitmap *texture,Point2D *unList,int textureIndex[],int shadeMode = POLY4D_ATTR_SHADE_MODE_GOURAUD);
	
	void reset();
};

struct PolySelfTriangle4D
{
	int      state;           // state information
	int      attr;            // physical attributes of polygon
	Color    color;           // color of polygon
	Color    lightColor[3];    // colors after lighting, 0 for flat shading
	int      material;    // material index (-1) for no material 

	float    avgZ;   // average z of vertices, used for simple sorting
	Vertex4D originVList[3];  // the vertices of this triangle 
	Vertex4D transVList[3];
	Vector4D originN;
	Vector4D transN;

	Bitmap *texure;

	void initFromTrans(const PolyTriangle4D *poly);
	void init(const PolySelfTriangle4D *poly);
};


