#pragma once
#include "CustomDefine.h"

#pragma pack(push, 1)
//14byte not align 4byte
struct  BitmapFileHeader {
	uint16 bfType;
	uint32 bfSize;
	uint16 bfReserved1;
	uint16 bfReserved2;
	uint32 bfOffBits;
};


//40byte align 4byte
struct BitmapInfoHeader {
	uint32 biSize;
	uint32 biWidth;
	int32 biHeight;
	uint16 biPlanes;
	uint16 biBitCount;
#define BI_NONE        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L
	uint32 biCompression;
	uint32 biSizeImage;
	uint32 biXPelsPerMeter;
	uint32 biYPelsPerMeter;
	uint32 biClrUsed;
	uint32 biClrImportant;
};

struct BitmapRGB{
	uint8 rgbBlue;
	uint8 rgbGreen;
	uint8 rgbRed;
	uint8 rgbAlpha;
};



#pragma pack(pop)

class Bitmap
{
public:
#define BITMAP_MAG 0x4D42
#define BITMAP_PFMT_RGB565 1
#define BITMAP_PFMT_RGB888  2
#define BITMAP_PFMT_ARGB8888 3
	
	static Bitmap* createBitmap(int width, int height, int pixFormat);
	static Bitmap* loadBitmap(char *fileName,int pixFormat= BITMAP_PFMT_RGB565);
	int flipBitmap();
	int getPitch();
	int getWidth();
	int getHeight();
	uint8* getBuffer();
	int getBitCount();
private:

	BitmapFileHeader fileHeader;
	BitmapInfoHeader infoHeader;
	BitmapRGB Palette[256];
	uint8 *buffer;
	int pixFormat;
	long pitch;//byte per line

	Bitmap();
	~Bitmap();
};

