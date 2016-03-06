#include "Bitmap.h"
#include <iostream>
#include "GraphicsStruct.h"
#include <math.h>






Bitmap * Bitmap::createBitmap(int w, int h, int pixFormat)
{	
	Bitmap* bitmap = new Bitmap();
	int bpp;
	int size;
	//init set 0
	memset(&bitmap->fileHeader, 0, sizeof(BitmapFileHeader));
	memset(&bitmap->infoHeader, 0, sizeof(BitmapInfoHeader));

	bitmap->pixFormat = pixFormat;

	if (pixFormat==BITMAP_PFMT_RGB888||pixFormat==BITMAP_PFMT_ARGB8888) {
		if (pixFormat == BITMAP_PFMT_RGB888) {
			bpp = 24;
		}
		else
		{
			bpp=32;
		}
		//file header
		size = w * abs(h) * (bpp >> 3);//byte
		bitmap->fileHeader.bfType = BITMAP_MAG;
		bitmap->fileHeader.bfOffBits = sizeof(BitmapFileHeader) + sizeof(BitmapInfoHeader);
		bitmap->fileHeader.bfSize = bitmap->fileHeader.bfOffBits + size;//there not need align 4byte

														
		bitmap->infoHeader.biSize = sizeof(BitmapInfoHeader); //info header
		bitmap->infoHeader.biWidth = w;
		bitmap->infoHeader.biHeight = h;
		bitmap->infoHeader.biBitCount = bpp;
		bitmap->infoHeader.biPlanes = 1;
		bitmap->infoHeader.biCompression = BI_NONE;

		bitmap->pitch = ((w * bpp + 31) / 32) * 8; //align 4byte  bit==> 32*((w*bpp+31)/32) byte=bit/8
		bitmap->infoHeader.biSizeImage = bitmap->pitch*abs(h);

		bitmap->buffer = new uint8[bitmap->infoHeader.biSizeImage];
	}
	return bitmap;
}

Bitmap * Bitmap::loadBitmap(char *fileName, int pixFormat)
{	
	Bitmap* bitmap = new  Bitmap();
	uint8* tmpBuffer;
	int lineLength,pitchTmp; //pix length ,because of 4byte align (mean that  be divided by 2 or 4 exactly; )so just do divide linelength=pitch (8) linelength=pitch>>1 (16)  linelength=pitch>>2 (32)
	FILE * file;
	fopen_s(&file,fileName, "rb");
	int a, r, g, b;
	int height;

	bitmap->pixFormat = pixFormat;

	fread(&bitmap->fileHeader, sizeof(BitmapFileHeader), 1, file);
	if (bitmap->fileHeader.bfType != BITMAP_MAG) {
		delete bitmap;
		return NULL;
	}
	fread(&bitmap->infoHeader, sizeof(BitmapInfoHeader), 1, file);

	tmpBuffer = new uint8[bitmap->infoHeader.biSizeImage];
	fseek(file,bitmap->fileHeader.bfOffBits,SEEK_SET);
	fread(tmpBuffer, bitmap->infoHeader.biSizeImage,1,file);

	if (bitmap->infoHeader.biBitCount == 24) {

		switch (pixFormat)
		{
		case BITMAP_PFMT_RGB565: {
			pitchTmp = bitmap->infoHeader.biSizeImage / bitmap->infoHeader.biHeight;
			bitmap->pitch = ((bitmap->infoHeader.biWidth * 16 + 31) / 32) * 4;//bpp =16;
			lineLength = bitmap->pitch >> 1;//2byte per pix
			height = abs(bitmap->infoHeader.biHeight);
			uint16* trueBuffer = (uint16 *)(new uint8[bitmap->pitch*height]);
			uint16 color;
			for (int yidx = 0;yidx < height;yidx++) {
				for (int xidx = 0;xidx < bitmap->infoHeader.biWidth;xidx++) {
					b = tmpBuffer[yidx*pitchTmp+xidx * 3];
					g = tmpBuffer[yidx*pitchTmp+xidx * 3 + 1];
					r = tmpBuffer[yidx*pitchTmp+xidx * 3 + 2];
					color = R_G_B_TO_RGB565(r, g, b);
					trueBuffer[yidx*lineLength+xidx] = color;
				}
			}
			bitmap->buffer = (uint8 *)trueBuffer;
			bitmap->infoHeader.biSizeImage = bitmap->pitch*height;
			bitmap->infoHeader.biBitCount = 16;
			bitmap->fileHeader.bfSize = bitmap->fileHeader.bfOffBits + bitmap->infoHeader.biSizeImage;
		}
			break;

		case BITMAP_PFMT_RGB888:
		default:
			bitmap->buffer = tmpBuffer;
			bitmap->pitch = bitmap->infoHeader.biSizeImage / bitmap->infoHeader.biHeight;
			break;
		}
	}
	else if(bitmap->infoHeader.biBitCount==32)
	{
		switch (pixFormat)
		{
		case BITMAP_PFMT_RGB565: {
			pitchTmp = bitmap->infoHeader.biSizeImage / bitmap->infoHeader.biHeight;
			bitmap->pitch = ((bitmap->infoHeader.biWidth * 16 + 31) / 32) * 4;//bpp =16;
			lineLength = bitmap->pitch >> 1;//2byte per pix
			height = abs(bitmap->infoHeader.biHeight);
			uint16* trueBuffer = (uint16 *)(new uint8[bitmap->pitch*height]);
			uint16 color;
			for (int yidx = 0;yidx < height;yidx++) {
				for (int xidx = 0;xidx < bitmap->infoHeader.biWidth;xidx++) {
					b = tmpBuffer[yidx*pitchTmp + xidx * 4];
					g = tmpBuffer[yidx*pitchTmp + xidx * 4 + 1];
					r = tmpBuffer[yidx*pitchTmp + xidx * 4 + 2];
					color = R_G_B_TO_RGB565(r, g, b);
					trueBuffer[yidx*lineLength+xidx] = color;
				}
			}
			bitmap->buffer = (uint8 *)trueBuffer;
			bitmap->fileHeader.bfSize = bitmap->fileHeader.bfOffBits + bitmap->pitch*height;
			bitmap->infoHeader.biBitCount = 16;
			bitmap->infoHeader.biSizeImage = bitmap->pitch*height;
		}
			break;

		case BITMAP_PFMT_ARGB8888:
		default:
			bitmap->buffer = tmpBuffer;
			bitmap->pitch = bitmap->infoHeader.biSizeImage / bitmap->infoHeader.biHeight;
			break;
		}
	}
	else
	{
		return NULL;
	}

	if (tmpBuffer) {
		delete[] tmpBuffer;
	}
	fclose(file);

	bitmap->flipBitmap();

	
	return bitmap;
}

int Bitmap::flipBitmap()
{
	uint8 *tmpb = new uint8[infoHeader.biSizeImage];
	if (!tmpb) {
		return -1;
	}
	int height = abs(infoHeader.biHeight);
	memcpy(tmpb, buffer, infoHeader.biSizeImage);

	for (int index = 0; index < height; index++) {
		memcpy(&buffer[((height - 1) - index)*pitch],&tmpb[index*pitch],pitch);
	}
	if (tmpb) {
		delete[] tmpb;
	}
	return 1;
	
}

int Bitmap::getPitch()
{
	return pitch;
}

int Bitmap::getWidth()
{
	return infoHeader.biWidth;
}

int Bitmap::getHeight()
{
	return infoHeader.biHeight;
}

uint8 * Bitmap::getBuffer()
{
	return buffer;
}

int Bitmap::getBitCount()
{
	return infoHeader.biBitCount;
}

Bitmap::Bitmap()
{
	
}

Bitmap::~Bitmap()
{
}
