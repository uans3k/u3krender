#include "Surface.h"
#include <Windows.h>
#include <windowsx.h>
#include "Light.h"
Surface::Surface(int mode)
{	
	this->mode = mode;
	isAttach = 0;
	isLocked = 0;
	ddraw = NULL;
	primaryBuffer = NULL;
	backBuffer = NULL;
	canvas = new Canvas();
	SinCosTables_Build();
	renderPiping.setSortZEnable(1);
}

Surface::~Surface()
{
	if (ddSurfaceBack) {
		ddSurfaceBack->Release();
	}
	if (ddSurfacePrimary) {
		ddSurfacePrimary->Release();
	}
	if (ddraw) {
		ddraw->Release();
	}
	if (canvas) {
		delete canvas;
	}
	if (!initCamera) {
		delete initCamera;
	}
}


int Surface::attachWindow(HWND windowHandle, int bpp) {

	DWORD style=GetWindowStyle(windowHandle);
	this->windowed = !(style&WS_POPUP);

	//params init
	GetWindowRect(windowHandle, &windowRect);
	windowWidth = windowRect.right - windowRect.left;
	windowHeight = windowRect.bottom - windowRect.top;
	this->bpp = bpp;


	if (FAILED(DirectDrawCreateEx(NULL, (void **)&ddraw, IID_IDirectDraw7, NULL))) {
		return -1;
	};

	// set cooperation level to fullscreen mode 
	if (FAILED(ddraw->SetCooperativeLevel(windowHandle,
		DDSCL_ALLOWMODEX | DDSCL_FULLSCREEN |
		DDSCL_EXCLUSIVE | DDSCL_ALLOWREBOOT | DDSCL_MULTITHREADED))) {
		return -1;
	}


	if (FAILED(ddraw->SetDisplayMode(windowWidth, windowHeight, bpp, 0, 0))) {
		return -1;
	}


	initSurfeceDesc();


	ddSurfaceDesc.dwFlags = DDSD_CAPS | DDSD_BACKBUFFERCOUNT;
	ddSurfaceDesc.ddsCaps.dwCaps = DDSCAPS_PRIMARYSURFACE | DDSCAPS_FLIP | DDSCAPS_COMPLEX;
	ddSurfaceDesc.dwBackBufferCount = 1;


	ddraw->CreateSurface(&ddSurfaceDesc, &ddSurfacePrimary, NULL);

	DDPIXELFORMAT pixFMT;
	memset(&pixFMT, 0, sizeof(pixFMT));
	ddSurfacePrimary->GetPixelFormat(&pixFMT);
	ddPixFMT = bpp;

	ddSurfaceDesc.ddsCaps.dwCaps = DDSCAPS_BACKBUFFER;

	if (FAILED(ddSurfacePrimary->GetAttachedSurface(&ddSurfaceDesc.ddsCaps, &ddSurfaceBack))) {
		return -1;
	}

	clearSurface(ddSurfaceBack, 0, &windowRect);
	clearSurface(ddSurfacePrimary, 0, &windowRect);



	initCamera = createCamera();
	isAttach = 1;

	return 0;
}

void Surface::setScene(Scene3D *scene)
{	
	this->scene3d = scene;
}

void Surface::setScene(Scene2D * scene)
{
	this->scene2d = scene;
}



Camera * Surface::createCamera()
{
	Point4D pos, target;
	Point4D_Init(&pos, 0, 0, 0);
	Point4D_Init(&target, 0, 0, 1);
	return new CameraUVN(&pos, &target, windowWidth, windowHeight, 3, 200);
}

int Surface::createZBuffer()
{	
	if (!isAttach) {
		return 0;
	}


	if (hasZBuffer) {
		return 1;
	}
	else
	{
		hasZBuffer = 1;
		zBuffer.height = windowHeight;
		zBuffer.width = windowWidth;
		zBuffer.pitch = windowWidth * 4;
		zBuffer.size = zBuffer.pitch*zBuffer.height;
		zBuffer.zBuffer = new uint32[zBuffer.width*zBuffer.height];
		if (!zBuffer.zBuffer) {
			clearZBuffer();
			return 0;
		}
	}
	return 1;
}

int Surface::setZBufferEnable(int enable)
{
	if (hasZBuffer) {
		enableZBuffer = enable;
		renderPiping.setSortZEnable(enable);
	}
	else
	{
		enableZBuffer = 0;
		return 0;
	}
	return 1;
}




void Surface::clearZBuffer()
{
	if (hasZBuffer&&zBuffer.zBuffer) {
		memset((void*)zBuffer.zBuffer, 0, zBuffer.pitch*zBuffer.height);
	}
}

void Surface::removeZBuffer()
{
	if (hasZBuffer) {
		hasZBuffer = 0;
		zBuffer.height = 0;
		zBuffer.width = 0;
		zBuffer.pitch = 0;
		zBuffer.size = 0;

		if (zBuffer.zBuffer) {
			delete[] zBuffer.zBuffer;
		}
	}
}

void Surface::setBackgroud(uint16 color)
{
	backgroudColor = color;
}



int Surface::render()
{	
	clear(backgroudColor);
	Canvas *c=lock();
	//canvas->drawTriangleSolid16(400, 300, 300, 200, 500, 300, 0x2222);
	//canvas->setClipper(200,100, 600, 500);

	if (mode == SURFACE_MODE_2D) {
		if (scene2d) {
			for (Object2D *obj : scene2d->objectVector) {
				canvas->drawBitmap(obj->rect.left, obj->rect.top, obj->rect.right, obj->rect.bottom, obj->bitmap);
			}
		}
	}
	else if(mode==SURFACE_MODE_3D)
	{
		//3d render piping
		renderPiping.clear();
		if (scene3d->getCamera()) {
			renderPiping.setCamera(scene3d->getCamera());
		}
		else
		{
			renderPiping.setCamera(initCamera);
		}
		for (Light *light : scene3d->getLightVector()) {
			renderPiping.addLight(light);
		}
		for (Object *object : scene3d->getObjectVector()) {
			if (!renderPiping.sendObject(object)) {
				break;
			}
		}

		renderPiping.render(c);
	}

	

	unlock();

	return 0;
}




//public 


Canvas* Surface::lock()
{	
	if (!isLocked&&isAttach) {
		uint8 *buffer;
		if (windowed) {
		}
		else
		{	
			lockSurface(ddSurfaceBack, &buffer, &pitch);
			if (enableZBuffer) {
				clearZBuffer();
				canvas->setBufferAndZbuffer(buffer, windowWidth, windowHeight, pitch, bpp,&zBuffer);
			}
			else
			{
				canvas->setBuffer(buffer, windowWidth, windowHeight, pitch, bpp);
			}
		}
		isLocked = 1;
		return canvas;
	}
	return NULL;
}


int Surface::unlock()
{
	if (isLocked) {
		if (windowed) {

		}
		else
		{
			unlockSurface(ddSurfaceBack);
			flip();
			canvas->reset();
		}
		isLocked = 0;
		return 1;
	}
	return 0;
}



int Surface::clear(uint16 color)
{
	return clearSurface(ddSurfaceBack, color, &windowRect);
}


//private

int Surface::flip()
{
	if (windowed) {
	}
	else
	{
		while (FAILED(ddSurfacePrimary->Flip(NULL, DDFLIP_WAIT)));
	}
	return 0;
}

int Surface::clearSurface(LPDIRECTDRAWSURFACE7 surface, uint16 color, RECT * region)
{
	DDBLTFX ddBltFx;

	memset(&ddBltFx, 0, sizeof(ddBltFx));
	ddBltFx.dwSize = sizeof(ddBltFx);


	ddBltFx.dwFillColor = color;


	surface->Blt(region,     // ptr to dest rectangle
		NULL,       // ptr to source surface, NA            
		NULL,       // ptr to source rectangle, NA
		DDBLT_COLORFILL | DDBLT_WAIT,   // fill and wait                   
		&ddBltFx);  // ptr to DDBLTFX structure

	return 0;
}

int  Surface::lockSurface(LPDIRECTDRAWSURFACE7 surface, uint8 **buffer, long * pitch)
{
	if (!isAttach || ddraw == NULL) {
		return -1;
	}
	initSurfeceDesc();
	surface->Lock(NULL, &ddSurfaceDesc, DDLOCK_WAIT | DDLOCK_SURFACEMEMORYPTR, NULL);
	if (pitch) {
		*pitch = ddSurfaceDesc.lPitch;
	}
	*buffer = (uint8 *)ddSurfaceDesc.lpSurface;
	return 1;
}

int Surface::unlockSurface(LPDIRECTDRAWSURFACE7 surface)
{
	if (!isAttach || ddraw == NULL) {
		return -1;
	}
	surface->Unlock(NULL);

	return 0;
}





