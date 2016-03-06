// SoftRender.cpp : 定义应用程序的入口点。
//

#include "stdafx.h"
#include "SoftRender.h"
#include "Surface.h"
#include "Timer.h"
#include "Math3D.h"
#include "Bitmap.h"

#define MAX_LOADSTRING 100
#define WINDOW_APP 0
#define WINDOW_WIDTH 1024
#define WINDOW_HEIGHT 768
#define WINDOW_BPP 16 

#define KEY_DOWN(vk_code) ((GetAsyncKeyState(vk_code) & 0x8000) ? 1 : 0)
#define KEY_UP(vk_code)   ((GetAsyncKeyState(vk_code) & 0x8000) ? 0 : 1)



// 全局变量: 
HINSTANCE hInst;                                // 当前实例
HWND mainWindowHandle;	//
WCHAR szTitle[MAX_LOADSTRING];                  // 标题栏文本
WCHAR szWindowClass[MAX_LOADSTRING];            // 主窗口类名

Scene3D *scene3d;
Scene2D *scene2d;
Object *testObject;
Object2D *testObject2d;
Surface *surface;
Camera *camera;
Light *testLight;
Bitmap *testBitmap;

Timer *timer; 


// 此代码模块中包含的函数的前向声明: 
ATOM                MyRegisterClass(HINSTANCE hInstance);
BOOL                InitInstance(HINSTANCE, int);
LRESULT CALLBACK    WndProc(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK    About(HWND, UINT, WPARAM, LPARAM);
int Game_Main(void *params = NULL);
int Game_Init(void *params = NULL);
int Game_ShutDown(void *params = NULL);

int APIENTRY wWinMain(_In_ HINSTANCE hInstance,
                     _In_opt_ HINSTANCE hPrevInstance,
                     _In_ LPWSTR    lpCmdLine,
                     _In_ int       nCmdShow)
{
    UNREFERENCED_PARAMETER(hPrevInstance);
    UNREFERENCED_PARAMETER(lpCmdLine);

    // TODO: 在此放置代码


    // 初始化全局字符串
    LoadStringW(hInstance, IDS_APP_TITLE, szTitle, MAX_LOADSTRING);
    LoadStringW(hInstance, IDC_SOFTRENDER, szWindowClass, MAX_LOADSTRING);
    MyRegisterClass(hInstance);

    // 执行应用程序初始化: 
    if (!InitInstance (hInstance, nCmdShow))
    {
        return FALSE;
    }

	
	Game_Init();
	


    HACCEL hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_SOFTRENDER));
    MSG msg;

    // 主消息循环: 
	for (;;)
    {	
		if (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE)) {
			if (msg.message == WM_QUIT)break;

			if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg))
			{
				TranslateMessage(&msg);
				DispatchMessage(&msg);
			}
		}
		Game_Main();
    }

    return (int) msg.wParam;
}



//
//  函数: MyRegisterClass()
//
//  目的: 注册窗口类。
//
ATOM MyRegisterClass(HINSTANCE hInstance)
{
    WNDCLASSEXW wcex;

    wcex.cbSize = sizeof(WNDCLASSEX);

    wcex.style          = CS_DBLCLKS|CS_HREDRAW | CS_OWNDC | CS_VREDRAW;
    wcex.lpfnWndProc    = WndProc;
    wcex.cbClsExtra     = 0;
    wcex.cbWndExtra     = 0;
    wcex.hInstance      = hInstance;
    wcex.hIcon          = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_SOFTRENDER));
    wcex.hCursor        = LoadCursor(nullptr, IDC_ARROW);
    wcex.hbrBackground  = (HBRUSH)GetStockObject(BLACK_BRUSH);
    wcex.lpszMenuName   = MAKEINTRESOURCEW(IDC_SOFTRENDER);
    wcex.lpszClassName  = szWindowClass;
    wcex.hIconSm        = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

    return RegisterClassExW(&wcex);
}

//
//   函数: InitInstance(HINSTANCE, int)
//
//   目的: 保存实例句柄并创建主窗口
//
//   注释: 
//
//        在此函数中，我们在全局变量中保存实例句柄并
//        创建和显示主程序窗口。
//
BOOL InitInstance(HINSTANCE hInstance, int nCmdShow)
{
	hInst = hInstance; // 将实例句柄存储在全局变量中

    mainWindowHandle = CreateWindowW(szWindowClass, szTitle, WINDOW_APP ?(WS_OVERLAPPED|WS_SYSMENU|WS_VISIBLE):(WS_POPUP|WS_VISIBLE),
      CW_USEDEFAULT, CW_USEDEFAULT, WINDOW_WIDTH,WINDOW_HEIGHT, nullptr, nullptr, hInstance, nullptr);

   if (!mainWindowHandle)
   {
      return FALSE;
   }



   //if(WINDOW_APP)
   //{
	  // RECT windowrRect = { 0,0,WINDOW_WIDTH,WINDOW_HEIGHT };
	  //

	  // AdjustWindowRectEx(&windowrRect,
		 //  GetWindowStyle(mainWindowHandle),
		 //  GetMenu(mainWindowHandle) != NULL,
		 //  GetWindowExStyle(mainWindowHandle));
	  // 
	  //
	  // MoveWindow(mainWindowHandle,
		 //  10, // x position
		 //  10, // y position
		 //  windowrRect.right - windowrRect.left, // width
		 //  windowrRect.bottom - windowrRect.top, // height
		 //  TRUE);
	  //
   //}
   

   ShowWindow(mainWindowHandle, nCmdShow);

   return TRUE;
}

//
//  函数: WndProc(HWND, UINT, WPARAM, LPARAM)
//
//  目的:    处理主窗口的消息。
//
//  WM_COMMAND  - 处理应用程序菜单
//  WM_PAINT    - 绘制主窗口
//  WM_DESTROY  - 发送退出消息并返回
//
//
LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    switch (message)
    {
    case WM_COMMAND:
        {
            int wmId = LOWORD(wParam);
            // 分析菜单选择: 
            switch (wmId)
            {
            case IDM_ABOUT:
                DialogBox(hInst, MAKEINTRESOURCE(IDD_ABOUTBOX), hWnd, About);
                break;
            case IDM_EXIT:
                DestroyWindow(hWnd);
                break;
            default:
                return DefWindowProc(hWnd, message, wParam, lParam);
            }
        }
        break;
    case WM_PAINT:
        {
            PAINTSTRUCT ps;
            HDC hdc = BeginPaint(hWnd, &ps);
            // TODO: 在此处添加使用 hdc 的任何绘图代码...
            EndPaint(hWnd, &ps);
        }
        break;
    case WM_DESTROY:
        PostQuitMessage(0);
        break;
    default:
        return DefWindowProc(hWnd, message, wParam, lParam);
    }
    return 0;
}

// “关于”框的消息处理程序。
INT_PTR CALLBACK About(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
    UNREFERENCED_PARAMETER(lParam);
    switch (message)
    {
    case WM_INITDIALOG:
        return (INT_PTR)TRUE;

    case WM_COMMAND:
        if (LOWORD(wParam) == IDOK || LOWORD(wParam) == IDCANCEL)
        {
            EndDialog(hDlg, LOWORD(wParam));
            return (INT_PTR)TRUE;
        }
        break;
    }
    return (INT_PTR)FALSE;
}


int pause=0;
float stepY=0.1f,stepX = 0.1f,stepZ=0.1;


int Game_Main(void *params) {

	
	timer->start();
	
	if (!pause) {
	/* if (testObject->worldPosition.x > 5 || testObject->worldPosition.x < -5) {
			stepX = -stepX;
		}

		if (testObject->worldPosition.y > 4 || testObject->worldPosition.y < -4) {
			stepY = -stepY;
		}
		if (testObject->worldPosition.z < 5 || testObject->worldPosition.z > 15) {
			stepZ = -stepZ;
		}

		testObject->worldPosition.x += stepX;
		testObject->worldPosition.y += stepY;
		testObject->worldPosition.z += stepZ;

		testObject->rotate(0, 4, 0);*/
		surface->render();

	}
	/*surface->clear(0xffff);*/
	//Canvas* canvas=surface->lock();

	////canvas->setClipper(400,200,800,600);

	////draw;
	///*surface->drawTest16(0xaaaa);*/
	//canvas->drawHLine16(100, 700, 300, 0);
	//canvas->drawVLine16(100, 500, 400, 0);
	//canvas->drawLine16(100, 300, 700, 500, 0);
	//canvas->drawLine16(100, 300, 200, 50, 0);
	//canvas->drawLine16(700, 300, 100, 100, 0);
	//canvas->drawLine16(700, 300, 600, 550,0);
	//canvas->drawTriangleWire16(300, 400, 500, 100, 700, 400,0);

	////surface->drawTopTriangle(400, 500, 200, 100, 600, 100,0);
	////surface->drawBottomTriangle16(400, 100, 200, 500, 600, 500, 0);

	//static int x1 = 100;
	//canvas->drawTriangleSolid16(x1, 100, 200, 300, 600, 500, 0x6666);

	//x1+=x1==500?-400:5;
	//surface->unlock();

	timer->wait();
//for camera
#define VK_W 87
#define VK_S 83
#define VK_A 65
#define VK_D 68

#define VK_Q 81
#define VK_E 69
#define VK_Z 90
#define VK_C 67

//for light
#define VK_T 84
#define VK_G 71
#define VK_F 70
#define VK_H 72

#define VK_R 82
#define VK_Y 89
#define VK_V 86
#define VK_N 78

//for object
#define VK_U 85
#define VK_O 79
#define VK_I 73 
#define VK_K 75
#define VK_J 74
#define VK_L 76



	if (KEY_DOWN(VK_ESCAPE)) {
		PostMessage(mainWindowHandle, WM_DESTROY, 0, 0);
		return -1;
	}
	if (KEY_DOWN(VK_SPACE)) {
		pause = !pause;
	}

	//for camera 
	if (KEY_DOWN(VK_W)) {
		camera->translate(0, 0, 0.5f);
	}
	if (KEY_DOWN(VK_S)) {
		camera->translate(0, 0, -0.5f);
	}
	if (KEY_DOWN(VK_A)) {
		camera->translate(-0.5f, 0, 0);
	}
	if (KEY_DOWN(VK_D)) {
		camera->translate(0.5f, 0, 0);
	}
	if (KEY_DOWN(VK_Z)) {
		camera->translate(0, 0.5f, 0);
	}
	if (KEY_DOWN(VK_C)) {
		camera->translate(0, -0.5f, 0);
	}

	//for light
	if (KEY_DOWN(VK_R)) {
		if (testLight->attr==LIGHT_ATTR_AMBIENT){
			testLight->ambient.r = MIN(testLight->ambient.r + 3, 255);
			testLight->ambient.g = MIN(testLight->ambient.g + 3, 255);
			testLight->ambient.b = MIN(testLight->ambient.b + 3, 255);
		}
		else
		{
			testLight->diffuse.r = MIN(testLight->diffuse.r + 3, 255);
			testLight->diffuse.g = MIN(testLight->diffuse.g + 3, 255);
			testLight->diffuse.b = MIN(testLight->diffuse.b + 3, 255);
		}
	}
	if (KEY_DOWN(VK_Y)) {
		if (testLight->attr==LIGHT_ATTR_AMBIENT) {
			testLight->ambient.r = MAX(testLight->ambient.r - 3, 0);
			testLight->ambient.g = MAX(testLight->ambient.r - 3, 0);
			testLight->ambient.b = MAX(testLight->ambient.r - 3, 0);
		}
		else
		{
			testLight->diffuse.r = MIN(testLight->diffuse.r - 3, 255);
			testLight->diffuse.g = MIN(testLight->diffuse.g - 3, 255);
			testLight->diffuse.b = MIN(testLight->diffuse.b - 3, 255);
		}
	}

	//for object
	if(KEY_DOWN(VK_U)){
		testObject->rotate(0, -3, 0);
		testObject2d->scale(0.9f, 1);
	}
	if (KEY_DOWN(VK_O)) {
		testObject->rotate(0, 3, 0);
		testObject2d->scale(1, 0.9f);
	}
	if (KEY_DOWN(VK_I)) {
		testObject->worldPosition.z +=0.3;
		testObject2d->translate(0, -5);
	}
	if (KEY_DOWN(VK_K)) {
		testObject->worldPosition.z -=0.3;
		testObject2d->translate(0, 5);
	}
	if (KEY_DOWN(VK_J)) {
		testObject->worldPosition.x -=0.3;
		testObject2d->translate(-5,0);
	}
	if (KEY_DOWN(VK_L)) {
		testObject->worldPosition.x +=0.3;
		testObject2d->translate(5,0);
	}

	return 0;
}


int Game_Init(void *params) {
	timer = new Timer(30);
	scene3d = new Scene3D();
	scene2d = new Scene2D();
	testObject = new Object();
	testLight = new Light();

	surface = new Surface();


	testBitmap = Bitmap::loadBitmap("D:\\texture\\tech01.bmp");
	testObject2d = new Object2D(100,100,500,500,testBitmap);
	scene2d->addObject(testObject2d);

	surface->attachWindow(mainWindowHandle, WINDOW_BPP);
	
	surface->createZBuffer();
	surface->setZBufferEnable(1);

	surface->setScene(scene3d);
	surface->setScene(scene2d);
	
	surface->setBackgroud(0xffff);
	
	

	//vertex
	Vertex4D *originVList=new Vertex4D[8];
	Vertex4D *transVList = new Vertex4D[8];
	Point4D_Init(&originVList[0].p, -1, 1, -1);
	Point4D_Init(&originVList[1].p, 1, 1, -1);
	Point4D_Init(&originVList[2].p, 1, -1, -1);
	Point4D_Init(&originVList[3].p, -1, -1, -1);

	Point4D_Init(&originVList[4].p, -1, 1, 1);
	Point4D_Init(&originVList[5].p, 1, 1, 1);
	Point4D_Init(&originVList[6].p, 1, -1, 1);
	Point4D_Init(&originVList[7].p, -1, -1, 1);

	//texure
	Bitmap* texture = testBitmap;
	Point2D *uvList = new Point2D[4];
	int uvNum = 4;
	Point2D_Init(&uvList[0], 0, 0);
	Point2D_Init(&uvList[1], 1, 0);
	Point2D_Init(&uvList[2], 1, 1);
	Point2D_Init(&uvList[3], 0, 1);
	int textureIndex[3];


	PolyTriangle4D *polyList= new PolyTriangle4D[12];
	Vector4D *originNList = new Vector4D[12];
	Vector4D *transNList = new Vector4D[12];
	int vidx[3];
	int nidx;
	Color color;

	//0
	vidx[0] = 0;
	vidx[1] = 1;
	vidx[2] = 2;
	nidx = 0;
	color.argb = 0xffaa00aa;
	textureIndex[0] = 0;
	textureIndex[1] = 1;
	textureIndex[2] = 2;
	//polyList[0].initColorPoly(originVList, transVList, vidx, originNList,transNList,nidx,color);
	polyList[0].initTexturePoly(originVList, transVList, vidx, originNList, transNList, nidx, texture, uvList, textureIndex);

	//1
	vidx[0] = 0;
	vidx[1] = 2;
	vidx[2] = 3;
	nidx = 1;
	color.argb = 0xffaaaa00;
	textureIndex[0] = 0;
	textureIndex[1] = 2;
	textureIndex[2] = 3;
	//polyList[1].initColorPoly(originVList, transVList, vidx, originNList, transNList, nidx, color);
	polyList[1].initTexturePoly(originVList, transVList, vidx, originNList, transNList, nidx, texture, uvList, textureIndex);

	//2
	vidx[0] = 0;
	vidx[1] = 4;
	vidx[2] = 1;
	nidx = 2;
	textureIndex[0] = 3;
	textureIndex[1] = 0;
	textureIndex[2] = 2;
	color.argb = 0xff002200;
	//polyList[2].initColorPoly(originVList, transVList, vidx, originNList, transNList, nidx, color);
	polyList[2].initTexturePoly(originVList, transVList, vidx, originNList, transNList, nidx, texture, uvList, textureIndex);
	//3
	vidx[0] = 4;
	vidx[1] = 5;
	vidx[2] = 1;
	nidx = 3;
	color.argb = 0xff000022;
	textureIndex[0] = 0;
	textureIndex[1] = 1;
	textureIndex[2] = 2;
	//polyList[3].initColorPoly(originVList, transVList, vidx, originNList, transNList, nidx, color);
	polyList[3].initTexturePoly(originVList, transVList, vidx, originNList, transNList, nidx, texture, uvList, textureIndex);

	//4
	vidx[0] = 3;
	vidx[1] = 2;
	vidx[2] = 7;
	nidx = 4;
	color.argb = 0xff440000;
	polyList[4].initColorPoly(originVList, transVList, vidx, originNList, transNList, nidx, color);

	//5
	vidx[0] = 2;
	vidx[1] = 6;
	vidx[2] = 7;
	nidx = 5;
	color.argb = 0xff004400;
	polyList[5].initColorPoly(originVList, transVList, vidx, originNList, transNList, nidx, color);

	//6
	vidx[0] = 1;
	vidx[1] = 6;
	vidx[2] = 2;
	nidx = 6;
	color.argb = 0xff000044;
	textureIndex[0] = 0;
	textureIndex[1] = 2;
	textureIndex[2] = 3;
	//polyList[6].initColorPoly(originVList, transVList, vidx, originNList, transNList, nidx, color);
	polyList[6].initTexturePoly(originVList, transVList, vidx, originNList, transNList, nidx, texture, uvList, textureIndex);
	//7
	vidx[0] = 1;
	vidx[1] = 5;
	vidx[2] = 6;
	nidx = 7;
	color.argb = 0xff880000;
	textureIndex[0] = 0;
	textureIndex[1] = 1;
	textureIndex[2] = 2;
	//polyList[7].initColorPoly(originVList, transVList, vidx, originNList, transNList, nidx, color);
	polyList[7].initTexturePoly(originVList, transVList, vidx, originNList, transNList, nidx, texture, uvList, textureIndex);

	//8
	vidx[0] = 0;
	vidx[1] = 3;
	vidx[2] = 7;
	nidx = 8;
	color.argb = 0xff008800;
	polyList[8].initColorPoly(originVList, transVList, vidx, originNList, transNList, nidx, color);

	//9
	vidx[0] = 4;
	vidx[1] = 0;
	vidx[2] = 7;
	nidx = 9;
	color.argb = 0xff000088;
	polyList[9].initColorPoly(originVList, transVList, vidx, originNList, transNList, nidx, color);
	
	//10
	vidx[0] = 6;
	vidx[1] = 5;
	vidx[2] = 4;
	nidx = 10;
	color.argb = 0xff888800;
	polyList[10].initColorPoly(originVList, transVList, vidx, originNList, transNList, nidx, color);

	//11
	vidx[0] = 7;
	vidx[1] = 6;
	vidx[2] = 4;
	nidx = 11;
	color.argb = 0xff008888;
	polyList[11].initColorPoly(originVList, transVList, vidx, originNList, transNList, nidx, color);


	testObject->framesNum = 1;
	testObject->currFrame = 0;
	testObject->avgRadius = new float[1];
	testObject->maxRadius = new float[1]; 
	testObject->avgRadius[0]= 1.7320f;
	testObject->maxRadius[0] = 1.7320f;
	Point4D_Init(&testObject->worldPosition, 0, 0, 9);
	testObject->verticesNum = 8;
	testObject->totalVerticesNum = 8;

	//vertex
	testObject->originVList = originVList;
	testObject->transVList = transVList;
	testObject->originNList= originNList;
	testObject->transNList = transNList;

	//normal
	testObject->headOriginVList= originVList;
	testObject->headTransVList = transVList;
	testObject->headOriginNList = originNList;
	testObject->headTransNList = originNList;

	//texture
	testObject->texture = texture;
	testObject->uvList = uvList;
	testObject->uvNum = uvNum;

	//poly
	testObject->polygonsNum = 12;
	testObject->polyList = polyList;

	testObject->initPreComputeNormal();
	scene3d->addObject(testObject);

	/*testObject->rotate(0, -45, 0);*/

	//light
	//Color ambient;
	//ambient.argb = 0xff888888;
	//testLight->initAmbientLight(ambient);
	

	//
	//Color diffuse;
	//diffuse.argb = 0xffffffff;
	//Color specular;
	//specular.argb = 0xff000000;
	//Vector4D dir;
	//Vector4D_Init(&dir, 0, 0, -1);
	//testLight->initInfinityLight(&dir,diffuse,specular);
	//

	Color diffuse;
	diffuse.argb = 0xffffffff;
	Color specular;
	specular.argb = 0xff000000;
	Point4D pos;
	Point4D_Init(&pos, 0, 0, 3);
	float kc = 1;
	float kl = 0;
	float kq = 0;
	testLight->initPointLight(&pos, diffuse, specular,kc,kl,kq);

	/*Color diffuse;
	diffuse.argb = 0xffffffff;
	Color specular;
	specular.argb = 0xff000000;
	Point4D pos;
	Vector4D dir;
	Point4D_Init(&dir, 0, 0, 1);
	Point4D_Init(&pos, 0, 0, 3);
	float kc = 1;
	float kl = 0.00;
	float kq = 0;
	int pf = 1;
	testLight->initSpotLight(&pos,&dir ,diffuse, specular, kc, kl, kq,pf);*/


	
	scene3d->addLight(testLight);
	return 0;

}

int Game_ShutDown(void *params){
	if (testObject) {
		delete testObject;
	}
	if (scene3d){
		delete scene3d;
	}
	if (surface) {
		delete surface;
	}
	return 0;
}