#pragma once
#include "CustomDefine.h"
#include "Canvas.h"
#include "ddraw.h"
#include "Scene.h"
#include "RenderPiping.h"




class Surface
{


public:

	/**
	Pix Format
	*/
#define PIX_FMT_8 8
#define PIX_FMT_555 15
#define PIX_FMT_565 16
#define PIX_FMT_888 24
#define PIX_FMT_8888 32

#define SURFACE_MODE_3D 0
#define SURFACE_MODE_2D 1


	Surface(int mode = SURFACE_MODE_3D);
	~Surface();
	int attachWindow(HWND windowHandle, int bpp);
	void setScene(Scene3D *scene);
	void setScene(Scene2D *scene);
	Camera* createCamera();
	int createZBuffer();
	int setZBufferEnable(int enable);
	void clearZBuffer();
	void removeZBuffer();
	void setBackgroud(uint16 color);
	int render();


private:
	int mode;

	Scene3D *scene3d;
	Scene2D *scene2d;
	Camera *initCamera;
	uint16 backgroudColor=0xffff;
    RenderPiping renderPiping;
	ZBuffer zBuffer;

	int isAttach;
	int isLocked;

	int enableZBuffer=0;
	int hasZBuffer=0;

	LPDIRECTDRAW7 ddraw;
	DDSURFACEDESC2 ddSurfaceDesc;
	LPDIRECTDRAWSURFACE7 ddSurfacePrimary;
	LPDIRECTDRAWSURFACE7 ddSurfaceBack;
	Canvas *canvas;
	
	uint8 *primaryBuffer;
	uint8 *backBuffer;

	RECT windowRect;
	long pitch;
	int windowWidth;
	int windowHeight;
	int bpp;
	int ddPixFMT;
	int windowed;

	Canvas* lock();
	int unlock();
	int clear(uint16 color);


	inline void initSurfeceDesc() {
		memset(&ddSurfaceDesc, 0, sizeof(ddSurfaceDesc));
		ddSurfaceDesc.dwSize = sizeof(ddSurfaceDesc);
	};
	int flip();
	int clearSurface(LPDIRECTDRAWSURFACE7 surface, uint16 color, RECT *region);
	int lockSurface(LPDIRECTDRAWSURFACE7 surface, uint8 **buffer, long * pitch);
	int unlockSurface(LPDIRECTDRAWSURFACE7 surface);

};
