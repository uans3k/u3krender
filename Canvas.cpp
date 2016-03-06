#include "Canvas.h"
#include "MemHandler.h"
#include <math.h>


Canvas::Canvas()
{
}


Canvas::~Canvas()
{
}

int Canvas::reset()
{
	hasBuffer = 0;
	enableZbuffer = 0;
	drawBuffer = NULL;
	this->width = 0;
	this->height = 0;
	this->pitch = 0;
	this->bpp = 0;
	setClipper(0, 0, 0, 0);
	return 0;
}

int Canvas::setBuffer(uint8 * buffer, int width, int height, int pitch, int bpp)
{
	drawBuffer = buffer;
	hasBuffer = 1;
	this->width = width;
	this->height = height;
	this->pitch = pitch;
	this->bpp = bpp;
	setClipper(0, 0, width-1, height-1);
	return 0;
}

int Canvas::setBufferAndZbuffer(uint8 * buffer, int width, int height, int pitch, int bpp, ZBuffer * zBuffer)
{	
	if(zBuffer->width>=width&&zBuffer->height>=height&&zBuffer->zBuffer){
		enableZbuffer = 1;
		this->zBuffer = zBuffer;
	}

	drawBuffer = buffer;
	hasBuffer = 1;
	this->width = width;
	this->height = height;
	this->pitch = pitch;
	this->bpp = bpp;
	setClipper(0, 0, width - 1, height - 1);
	return 0;
}


int Canvas::drawPolyTriangleSolidI16(PolySelfTriangle4D *poly)
{	

	return drawTriangleSolidI16(poly->transVList[0].x, poly->transVList[0].y, //
								poly->transVList[1].x, poly->transVList[1].y,//
								poly->transVList[2].x, poly->transVList[2].y,//
								poly->color);
}



int Canvas::drawTest16(uint16 color)
{
	uint16 *buffer = (uint16 *)drawBuffer;
	int lineLength = pitch >> 1;
	int x, y = 10;
	for (int i = 0;i < 8;i++) {
		for (int j = 0;j < width;j++) {
			buffer[j + y*lineLength] = 0;
		}
		y += 50;
	}
	return 0;
}

int Canvas::drawPolyTriangle(PolySelfTriangle4D * poly)
{
	if(enableZbuffer){
		if (IS_BIT(poly->attr, POLY4D_ATTR_SHADE_MODE_TEXTURE)) {
			Bitmap *texture;
			texture = poly->texure;
			if (IS_BIT(poly->attr, POLY4D_ATTR_SHADE_MODE_CONSTANT)) {
				drawTriangleTextureInvZ16(poly->transVList[0].x, poly->transVList[0].y, poly->transVList[0].z, poly->transVList[0].u, poly->transVList[0].v,
					poly->transVList[1].x, poly->transVList[1].y, poly->transVList[1].z, poly->transVList[1].u, poly->transVList[1].v,
					poly->transVList[2].x, poly->transVList[2].y, poly->transVList[2].z, poly->transVList[2].u, poly->transVList[2].v,
					(uint16 *)texture->getBuffer(),texture->getWidth(),texture->getHeight(),texture->getPitch());
			}
			else if (IS_BIT(poly->attr, POLY4D_ATTR_SHADE_MODE_FLAT))
			{
				drawTriangleTextureFlatInvZ16(poly->transVList[0].x, poly->transVList[0].y, poly->transVList[0].z, poly->transVList[0].u, poly->transVList[0].v,
					poly->transVList[1].x, poly->transVList[1].y, poly->transVList[1].z, poly->transVList[1].u, poly->transVList[1].v,
					poly->transVList[2].x, poly->transVList[2].y, poly->transVList[2].z, poly->transVList[2].u, poly->transVList[2].v,
					poly->lightColor[0],
					(uint16 *)texture->getBuffer(), texture->getWidth(), texture->getHeight(), texture->getPitch());
			}
			else if (IS_BIT(poly->attr, POLY4D_ATTR_SHADE_MODE_GOURAUD))
			{
				drawTriangleTextureGouraudInvZ16(poly->transVList[0].x, poly->transVList[0].y, poly->transVList[0].z, poly->transVList[0].u, poly->transVList[0].v, poly->lightColor[0],
					poly->transVList[1].x, poly->transVList[1].y, poly->transVList[1].z, poly->transVList[1].u, poly->transVList[1].v, poly->lightColor[1],
					poly->transVList[2].x, poly->transVList[2].y, poly->transVList[2].z, poly->transVList[2].u, poly->transVList[2].v, poly->lightColor[2],
					(uint16 *)texture->getBuffer(), texture->getWidth(), texture->getHeight(), texture->getPitch());
			}
		}
		else
		{
		if (IS_BIT(poly->attr, POLY4D_ATTR_SHADE_MODE_CONSTANT) || IS_BIT(poly->attr, POLY4D_ATTR_SHADE_MODE_FLAT)) {
				drawTriangleSolidInvZ16(poly->transVList[0].x, poly->transVList[0].y, poly->transVList[0].z,
					poly->transVList[1].x, poly->transVList[1].y, poly->transVList[1].z,
					poly->transVList[2].x, poly->transVList[2].y, poly->transVList[2].z,
					poly->lightColor[0]);
			}
			else if (IS_BIT(poly->attr, POLY4D_ATTR_SHADE_MODE_GOURAUD))
			{
				drawTriangleGouraudInvZ16(poly->transVList[0].x, poly->transVList[0].y, poly->transVList[0].z, poly->lightColor[0],
					poly->transVList[1].x, poly->transVList[1].y, poly->transVList[1].z, poly->lightColor[1],
					poly->transVList[2].x, poly->transVList[2].y, poly->transVList[2].z, poly->lightColor[2]);
			}
		}
	}
	else
	{

		if (IS_BIT(poly->attr, POLY4D_ATTR_SHADE_MODE_TEXTURE)) {
			Bitmap *texture;
			texture = poly->texure;
			if (IS_BIT(poly->attr, POLY4D_ATTR_SHADE_MODE_CONSTANT)) {
			   drawTriangleTexture16(poly->transVList[0].x, poly->transVList[0].y, poly->transVList[0].u, poly->transVList[0].v,
					poly->transVList[1].x, poly->transVList[1].y, poly->transVList[1].u, poly->transVList[1].v,
					poly->transVList[2].x, poly->transVList[2].y, poly->transVList[2].u, poly->transVList[2].v,
					(uint16 *)texture->getBuffer(), texture->getWidth(), texture->getHeight(), texture->getPitch());
			}
			else if (IS_BIT(poly->attr, POLY4D_ATTR_SHADE_MODE_FLAT))
			{
				drawTriangleTextureFlat16(poly->transVList[0].x, poly->transVList[0].y, poly->transVList[0].u, poly->transVList[0].v,
					poly->transVList[1].x, poly->transVList[1].y, poly->transVList[1].u, poly->transVList[1].v,
					poly->transVList[2].x, poly->transVList[2].y, poly->transVList[2].u, poly->transVList[2].v,
					poly->lightColor[0],
					(uint16 *)texture->getBuffer(), texture->getWidth(), texture->getHeight(), texture->getPitch());
			}
			else if (IS_BIT(poly->attr, POLY4D_ATTR_SHADE_MODE_GOURAUD))
			{
				drawTriangleTextureGouraud16(poly->transVList[0].x, poly->transVList[0].y, poly->transVList[0].u, poly->transVList[0].v, poly->lightColor[0],
					poly->transVList[1].x, poly->transVList[1].y, poly->transVList[1].u, poly->transVList[1].v, poly->lightColor[1],
					poly->transVList[2].x, poly->transVList[2].y, poly->transVList[2].u, poly->transVList[2].v, poly->lightColor[2],
					(uint16 *)texture->getBuffer(), texture->getWidth(), texture->getHeight(), texture->getPitch());
			}
		}
		else
		{
			if (IS_BIT(poly->attr, POLY4D_ATTR_SHADE_MODE_CONSTANT) || IS_BIT(poly->attr, POLY4D_ATTR_SHADE_MODE_FLAT)) {
				drawTriangleSolidF16(poly->transVList[0].x, poly->transVList[0].y, //
					poly->transVList[1].x, poly->transVList[1].y,//
					poly->transVList[2].x, poly->transVList[2].y,//
					poly->lightColor[0]);
			}
			else if (IS_BIT(poly->attr, POLY4D_ATTR_SHADE_MODE_GOURAUD))
			{
				drawTriangleGouraud16(poly->transVList[0].x, poly->transVList[0].y, //
					poly->transVList[1].x, poly->transVList[1].y,//
					poly->transVList[2].x, poly->transVList[2].y,//
					poly->lightColor[0], poly->lightColor[1], poly->lightColor[2]);
			}
		}
	}

	return 1;
	
}

int Canvas::drawTriangleTextureGouraudInvZ16(float x0, float y0, float z0, float u0, float v0, Color c0, float x1, float y1, float z1, float u1, float v1, Color c1, float x2, float y2, float z2, float u2, float v2, Color c2, uint16 * texture, int width, int height, int textPitch)
{
	float  tmp, interX, interZ, interU, interV, fraction;
	Color tmpC, interC;
	
	//line so return 0 as trivially reject
	if ((FEQUALE1(y0, y1) && FEQUALE1(y1, y2)) || (FEQUALE1(x0, x1) && FEQUALE1(x1, x2))) {
		return 0;
	}

	//select sort for y small -> big
	if (y1 < y0) {
		SWAP(y1, y0, tmp);
		SWAP(x1, x0, tmp);
		SWAP(z1, z0, tmp);
		SWAP(u1, u0, tmp);
		SWAP(v1, v0, tmp);
		SWAP(c1, c0, tmpC);
		
	}
	if (y2 < y0) {
		SWAP(y2, y0, tmp);
		SWAP(x2, x0, tmp);
		SWAP(z2, z0, tmp);
		SWAP(u2, u0, tmp);
		SWAP(v2, v0, tmp);
		SWAP(c2, c0, tmpC);
	}
	if (y2 < y1) {
		SWAP(y1, y2, tmp);
		SWAP(x1, x2, tmp);
		SWAP(z1, z2, tmp);
		SWAP(u1, u2, tmp);
		SWAP(v1, v2, tmp);
		SWAP(c1, c2, tmpC);
	}

	//trivially reject
	if (y0 > cliper.bottom || y2<cliper.top || (x0<cliper.left&&x1<cliper.left&&x2<cliper.left) || (x0>cliper.right&&x1>cliper.right&&x2>cliper.right)) {
		return 0;
	}

	if (FEQUALE1(y0, y1)) {
		drawTopTriangleTextureGouraudInvZ16(x2, y2, z2, u2, v2,c2,
			x0, y0, z0, u0, v0,c0,
			x1, y1, z1, u1, v1,c1,
			texture, width, height, textPitch);
	}
	else if (FEQUALE1(y1, y2))
	{
		drawBottomTriangleTextureGouraudInvZ16(x0, y0, z0, u0, v0,c0,
			x1, y1, z1, u1, v1,c1,
			x2, y2, z2, u2, v2,c2,
			texture, width, height, textPitch);
	}
	else
	{
		fraction = (y1 - y0) / (y2 - y0);
		interX = x0 + (x2 - x0)*fraction;
		interZ = z0 + (z2 - z0)*fraction;
		interU = u0 + (u2 - u0)*fraction;
		interV = v0 + (v2 - v0)*fraction;
		interC.r = c0.r + (c2.r - c0.r)*fraction;
		interC.g = c0.g + (c2.g - c0.g)*fraction;
		interC.b = c0.b + (c2.b - c0.b)*fraction;

		drawBottomTriangleTextureGouraudInvZ16(x0, y0, z0, u0, v0,c0,
			x1, y1, z1, u1, v1,c1,
			interX, y1, interZ, interU, interV,interC,
			texture, width, height, textPitch);

		drawTopTriangleTextureGouraudInvZ16(x2, y2, z2, u2, v2,c2,
			x1, y1, z1, u1, v1,c1,
			interX, y1, interZ, interU, interV,interC,
			texture, width, height, textPitch);
	}

	return 1;
}

int Canvas::drawTopTriangleTextureGouraudInvZ16(float x0, float y0, float z0, float u0, float v0, Color c0, float x1, float y1, float z1, float u1, float v1, Color c1, float x2, float y2, float z2, float u2, float v2, Color c2, uint16 * texture, int width, int height, int textPitch)
{
	fixp19 xl, xr, dxdyl, dxdyr;

	fixp20  invzFixp20,
		ul, ur, vl, vr,
		dudyl, dudyr, dvdyl, dvdyr,
		invzU0, invzV0, invzU1, invzV1, invzU2, invzV2, invzU, invzV,
		dudx,dvdx,
		rl, rr, gl, gr, bl, br,
		drdyl, drdyr, dgdyl, dgdyr, dbdyl, dbdyr,
		invzR0, invzG0, invzB0, invzR1, invzG1, invzB1, invzR2, invzG2, invzB2, invzR, invzG, invzB,
		drdx, dgdx, dbdx;


	fixp28 invz0, invz1, invz2, invz, zl, zr, dzdyl, dzdyr, dzdx;

	int idx, idy,
		iys, iye, ixs, ixe,
		iu, iv,
		ir, ig, ib,
		rDraw,gDraw,bDraw,
		r0,g0,b0,r1,g1,b1,r2,g2,b2,
		tmpi, textShift, destBufferLingLength, zBufferLineLength;

	float dx, dy, tmpf;
	uint16 *destBuffer,drawColor;
	uint32 *destZBuffer;

	destBufferLingLength = pitch >> 1;
	zBufferLineLength = zBuffer->pitch >> 2;//32 bit zbuffer wo fix  zBuffer deepth
	textShift = Log_Base2(width);


	r0 = c0.r;
	r1 = c1.r;
	r2 = c2.r;

	g0 = c0.g;
	g1 = c1.g;
	g2 = c2.g;

	b0 = c0.b;
	b1 = c1.b;
	b2 = c2.b;



	if (x2 < x1) {
		SWAP(x2, x1, tmpf);
		SWAP(y2, y1, tmpf);
		SWAP(z2, z1, tmpf);
		SWAP(u2, u1, tmpf);
		SWAP(v2, v1, tmpf);

		SWAP(r2, r1, tmpi);
		SWAP(g2, g1, tmpi);
		SWAP(b2, b1, tmpi);
	}

	invz0 = FIXP28_FLOAT_TO(1.0f / z0);
	invz1 = FIXP28_FLOAT_TO(1.0f / z1);
	invz2 = FIXP28_FLOAT_TO(1.0f / z2);


	//pre compute pix 
	tmpf = width - 1;
	invzU0 = FIXP20_FLOAT_TO((u0*tmpf) / z0);
	invzU1 = FIXP20_FLOAT_TO((u1*tmpf) / z1);
	invzU2 = FIXP20_FLOAT_TO((u2*tmpf) / z2);

	tmpf = height - 1;
	invzV0 = FIXP20_FLOAT_TO((v0*tmpf) / z0);
	invzV1 = FIXP20_FLOAT_TO((v1*tmpf) / z1);
	invzV2 = FIXP20_FLOAT_TO((v2*tmpf) / z2);

   
	invzR0 = FIXP20_FLOAT_TO(r0 / z0);
	invzG0 = FIXP20_FLOAT_TO(g0 / z0);
	invzB0 = FIXP20_FLOAT_TO(b0 / z0);

	invzR1 = FIXP20_FLOAT_TO(r1 / z1);
	invzG1 = FIXP20_FLOAT_TO(g1 / z1);
	invzB1 = FIXP20_FLOAT_TO(b1 / z1);

	invzR2 = FIXP20_FLOAT_TO(r2 / z2);
	invzG2 = FIXP20_FLOAT_TO(g2 / z2);
	invzB2 = FIXP20_FLOAT_TO(b2 / z2);



	//drawTraingleGouraud16 has ignored  y0 == y1 == y2 so just div
	dy = y0 - y1;


	dxdyl = FIXP19_FLOAT_TO((x0 - x1) / dy);

	dzdyl = (invz0 - invz1) / dy;
	dudyl = (invzU0 - invzU1) / dy;
	dvdyl = (invzV0 - invzV1) / dy;

	drdyl = (invzR0 - invzR1) / dy;
	dgdyl = (invzG0 - invzG1) / dy;
	dbdyl = (invzB0 - invzB1) / dy;

	dy = y0 - y2;
	dxdyr = FIXP19_FLOAT_TO((x0 - x2) / dy);
	dzdyr = (invz0 - invz2) / dy;
	dudyr = (invzU0 - invzU2) / dy;
	dvdyr = (invzV0 - invzV2) / dy;

	drdyr = (invzR0 - invzR2) / dy;
	dgdyr = (invzG0 - invzG2) / dy;
	dbdyr = (invzB0 - invzB2) / dy;

	//clip bottom
	if (y0>cliper.bottom)
	{
		iye = cliper.bottom;
	}
	else
	{
		iye = ceil(y0) - 1;
	}

	//clip top
	if (y1 < cliper.top) {
		iys = cliper.top;

	}
	else
	{
		iys = ceil(y1);
	}

	dy = iys - y1;

	xl = FIXP19_FLOAT_TO(x1) + dy*dxdyl;
	zl = invz1 + dy*dzdyl;
	ul = invzU1 + dy*dudyl;
	vl = invzV1 + dy*dvdyl;
	rl = invzR1 + dy*drdyl;
	gl = invzG1 + dy*dgdyl;
	bl = invzB1 + dy*dbdyl;


	dy = iys - y2;

	xr = FIXP19_FLOAT_TO(x2) + dy*dxdyr;
	zr = invz2 + dy*dzdyr;
	ur = invzU2 + dy*dudyr;
	vr = invzV2 + dy*dvdyr;
	rr = invzR2 + dy*drdyr;
	gr = invzG2 + dy*dgdyr;
	br = invzB2 + dy*dbdyr;






	//draw

	//pre round up

	xl = FIXP19_ROUNDUP(xl);
	xr = FIXP19_ROUNDUP(xr);



	destBuffer = (uint16 *)drawBuffer + iys*destBufferLingLength;
	destZBuffer = zBuffer->zBuffer + iys*zBufferLineLength;

	if ((x0 < cliper.left) || (x0 > cliper.right) ||
		(x1 <cliper.left) || (x1 > cliper.right) ||
		(x2 < cliper.left) || (x2 >cliper.right)) {
		//should clip


		for (int yidx = iys;yidx <= iye;yidx++) {

			ixs = FIXP19_TO_INT(xl);
			ixe = FIXP19_TO_INT(xr);


			if (ixe < cliper.left || ixs>cliper.right) {

				xl += dxdyl;
				xr += dxdyr;

				zl += dzdyl;
				zr += dzdyr;

				ul += dudyl;
				ur += dudyr;

				vl += dvdyl;
				vr += dvdyr;


				rl += drdyl;
				rr += drdyr;

				gl += dgdyl;
				gr += dgdyr;

				bl += dbdyl;
				br += dbdyr;


				destBuffer += destBufferLingLength;
				destZBuffer += zBufferLineLength;
				continue;
			}

			idx = ixe - ixs;
			if (idx>0) {
				dzdx = (zr - zl) / idx;
				dudx = (ur - ul) / idx;
				dvdx = (vr - vl) / idx;

				drdx = (rr - rl) / idx;
				dgdx = (gr - gl) / idx;
				dbdx = (br - bl) / idx;

			}
			else
			{
				dudx = (ur - ul);
				dvdx = (vr - vl);
				dzdx = (zr - zl);

				drdx = (rr - rl);
				dgdx = (gr - gl);
				dbdx = (br - gl);
			}


			invzU = ul;
			invzV = vl;
			invz = zl;

			invzR = rl;
			invzG = gl;
			invzB = bl;


			if (ixs< cliper.left) {
				idx = cliper.left - ixs;

				invzU += idx*dudx;
				invzV += idx*dvdx;
				invz += idx*dzdx;

				invzR += idx*drdx;
				invzG += idx*dgdx;
				invzB += idx*dbdx;

				ixs = cliper.left;
			}
			if (ixe > cliper.right) {
				ixe = cliper.right;
			}


			for (int xidx = ixs;xidx < ixe;xidx++) {
				if (invz > destZBuffer[xidx]) {
					invzFixp20 = FIXP28_TO_FIXP20(invz);

					iu = (invzU / invzFixp20);
					iv = (invzV / invzFixp20) << textShift;

					ir = invzR / invzFixp20;
					ig = invzG / invzFixp20;
					ib = invzB / invzFixp20;
					drawColor= texture[iu + iv];
					RGB565_TO_R5_G6_B5(drawColor, rDraw, gDraw, bDraw);
					rDraw = (rDraw*ir) >> 8;
					gDraw = (gDraw*ig) >> 8;
					bDraw = (bDraw*ib) >> 8;
					drawColor = R5_G6_B5_TO_RGB565(rDraw, gDraw, bDraw);

					destBuffer[xidx] = drawColor;
					destZBuffer[xidx] = invz;
				}

				invz += dzdx;
				invzU += dudx;
				invzV += dvdx;

				invzR += drdx;
				invzG += dgdx;
				invzB += dbdx;
			}

			xl += dxdyl;
			xr += dxdyr;

			zl += dzdyl;
			zr += dzdyr;

			ul += dudyl;
			ur += dudyr;

			vl += dvdyl;
			vr += dvdyr;


			rl += drdyl;
			rr += drdyr;

			gl += dgdyl;
			gr += dgdyr;

			bl += dbdyl;
			br += dbdyr;




			destBuffer += destBufferLingLength;
			destZBuffer += zBufferLineLength;
		}
	}
	else
	{
		for (int yidx = iys;yidx <= iye;yidx++) {

			ixs = FIXP19_TO_INT(xl);
			ixe = FIXP19_TO_INT(xr);

			idx = ixe - ixs;
			if (idx>0) {
				dzdx = (zr - zl) / idx;
				dudx = (ur - ul) / idx;
				dvdx = (vr - vl) / idx;

				drdx = (rr - rl) / idx;
				dgdx = (gr - gl) / idx;
				dbdx = (br - bl) / idx;
			}
			else
			{
				dzdx = (zr - zl);
				dudx = (ur - ul);
				dvdx = (vr - vl);

				drdx = (rr - rl);
				dgdx = (gr - gl);
				dbdx = (br - bl);
			}


			invzU = ul;
			invzV = vl;
			invz = zl;


			invzR = rl;
			invzG = gl;
			invzB = bl;


			for (int xidx = ixs;xidx < ixe;xidx++) {
				if (invz > destZBuffer[xidx]) {
					invzFixp20 = FIXP28_TO_FIXP20(invz);
					iu = (invzU / invzFixp20);
					iv = (invzV / invzFixp20) << textShift;

					ir = invzR / invzFixp20;
					ig = invzG / invzFixp20;
					ib = invzB / invzFixp20;
					drawColor = texture[iu + iv];
					RGB565_TO_R5_G6_B5(drawColor, rDraw, gDraw, bDraw);
					rDraw = (rDraw*ir) >> 8;
					gDraw = (gDraw*ig) >> 8;
					bDraw = (bDraw*ib) >> 8;
					drawColor = R5_G6_B5_TO_RGB565(rDraw, gDraw, bDraw);
					destBuffer[xidx] = drawColor;

					destZBuffer[xidx] = invz;
				}

				invz += dzdx;
				invzU += dudx;
				invzV += dvdx;


				invzR += drdx;
				invzG += dgdx;
				invzB += dbdx;
			}


			xl += dxdyl;
			xr += dxdyr;

			zl += dzdyl;
			zr += dzdyr;

			ul += dudyl;
			ur += dudyr;

			vl += dvdyl;
			vr += dvdyr;


			rl += drdyl;
			rr += drdyr;

			gl += dgdyl;
			gr += dgdyr;

			bl += dbdyl;
			br += dbdyr;



			destBuffer += destBufferLingLength;
			destZBuffer += zBufferLineLength;
		}
	}


	return 1;
}

int Canvas::drawBottomTriangleTextureGouraudInvZ16(float x0, float y0, float z0, float u0, float v0, Color c0, float x1, float y1, float z1, float u1, float v1, Color c1, float x2, float y2, float z2, float u2, float v2, Color c2, uint16 * texture, int width, int height, int textPitch)
{
	fixp19 xl, xr, dxdyl, dxdyr;

	fixp20  invzFixp20,
		ul, ur, vl, vr,
		dudyl, dudyr, dvdyl, dvdyr,
		invzU0, invzV0, invzU1, invzV1, invzU2, invzV2, invzU, invzV,
		dudx, dvdx,
		rl, rr, gl, gr, bl, br,
		drdyl, drdyr, dgdyl, dgdyr, dbdyl, dbdyr,
		invzR0, invzG0, invzB0, invzR1, invzG1, invzB1, invzR2, invzG2, invzB2, invzR, invzG, invzB,
		drdx, dgdx, dbdx;


	fixp28 invz0, invz1, invz2, invz, zl, zr, dzdyl, dzdyr, dzdx;

	int idx, idy,
		iys, iye, ixs, ixe,
		iu, iv,
		ir, ig, ib,
		rDraw, gDraw, bDraw,
		r0, g0, b0, r1, g1, b1, r2, g2, b2,
		tmpi, textShift, destBufferLingLength, zBufferLineLength;

	float dx, dy, tmpf;
	uint16 *destBuffer, drawColor;
	uint32 *destZBuffer;

	destBufferLingLength = pitch >> 1;
	zBufferLineLength = zBuffer->pitch >> 2;//32 bit zbuffer wo fix  zBuffer deepth
	textShift = Log_Base2(width);


	r0 = c0.r;
	r1 = c1.r;
	r2 = c2.r;

	g0 = c0.g;
	g1 = c1.g;
	g2 = c2.g;

	b0 = c0.b;
	b1 = c1.b;
	b2 = c2.b;



	if (x2 < x1) {
		SWAP(x2, x1, tmpf);
		SWAP(y2, y1, tmpf);
		SWAP(z2, z1, tmpf);
		SWAP(u2, u1, tmpf);
		SWAP(v2, v1, tmpf);

		SWAP(r2, r1, tmpi);
		SWAP(g2, g1, tmpi);
		SWAP(b2, b1, tmpi);
	}

	invz0 = FIXP28_FLOAT_TO(1.0f / z0);
	invz1 = FIXP28_FLOAT_TO(1.0f / z1);
	invz2 = FIXP28_FLOAT_TO(1.0f / z2);


	//pre compute pix 
	tmpf = width - 1;
	invzU0 = FIXP20_FLOAT_TO((u0*tmpf) / z0);
	invzU1 = FIXP20_FLOAT_TO((u1*tmpf) / z1);
	invzU2 = FIXP20_FLOAT_TO((u2*tmpf) / z2);

	tmpf = height - 1;
	invzV0 = FIXP20_FLOAT_TO((v0*tmpf) / z0);
	invzV1 = FIXP20_FLOAT_TO((v1*tmpf) / z1);
	invzV2 = FIXP20_FLOAT_TO((v2*tmpf) / z2);


	invzR0 = FIXP20_FLOAT_TO(r0 / z0);
	invzG0 = FIXP20_FLOAT_TO(g0 / z0);
	invzB0 = FIXP20_FLOAT_TO(b0 / z0);

	invzR1 = FIXP20_FLOAT_TO(r1 / z1);
	invzG1 = FIXP20_FLOAT_TO(g1 / z1);
	invzB1 = FIXP20_FLOAT_TO(b1 / z1);

	invzR2 = FIXP20_FLOAT_TO(r2 / z2);
	invzG2 = FIXP20_FLOAT_TO(g2 / z2);
	invzB2 = FIXP20_FLOAT_TO(b2 / z2);



	//drawTraingleGouraud16 has ignored  y0 == y1 == y2 so just div
	dy = y1 - y0;


	dxdyl = FIXP19_FLOAT_TO((x1 - x0) / dy);

	dzdyl = (invz1 - invz0) / dy;
	dudyl = (invzU1 - invzU0) / dy;
	dvdyl = (invzV1 - invzV0) / dy;

	drdyl = (invzR1 - invzR0) / dy;
	dgdyl = (invzG1 - invzG0) / dy;
	dbdyl = (invzB1 - invzB0) / dy;

	dy = y2 - y0;
	dxdyr = FIXP19_FLOAT_TO((x2 - x0) / dy);
	dzdyr = (invz2 - invz0) / dy;
	dudyr = (invzU2 - invzU0) / dy;
	dvdyr = (invzV2 - invzV0) / dy;

	drdyr = (invzR2 - invzR0) / dy;
	dgdyr = (invzG2 - invzG0) / dy;
	dbdyr = (invzB2 - invzB0) / dy;

	//clip bottom
	if (y1>cliper.bottom)
	{
		iye = cliper.bottom;
	}
	else
	{
		iye = ceil(y1) - 1;
	}

	//clip top
	if (y0 < cliper.top) {
		iys = cliper.top;

	}
	else
	{
		iys = ceil(y0);
	}

	dy = iys - y0;

	xl = FIXP19_FLOAT_TO(x0) + dy*dxdyl;
	zl = invz0 + dy*dzdyl;
	ul = invzU0 + dy*dudyl;
	vl = invzV0 + dy*dvdyl;
	rl = invzR0 + dy*drdyl;
	gl = invzG0 + dy*dgdyl;
	bl = invzB0 + dy*dbdyl;



	xr = FIXP19_FLOAT_TO(x0) + dy*dxdyr;
	zr = invz0 + dy*dzdyr;
	ur = invzU0 + dy*dudyr;
	vr = invzV0 + dy*dvdyr;
	rr = invzR0 + dy*drdyr;
	gr = invzG0 + dy*dgdyr;
	br = invzB0 + dy*dbdyr;






	//draw

	//pre round up

	xl = FIXP19_ROUNDUP(xl);
	xr = FIXP19_ROUNDUP(xr);



	destBuffer = (uint16 *)drawBuffer + iys*destBufferLingLength;
	destZBuffer = zBuffer->zBuffer + iys*zBufferLineLength;

	if ((x0 < cliper.left) || (x0 > cliper.right) ||
		(x1 <cliper.left) || (x1 > cliper.right) ||
		(x2 < cliper.left) || (x2 >cliper.right)) {
		//should clip


		for (int yidx = iys;yidx <= iye;yidx++) {

			ixs = FIXP19_TO_INT(xl);
			ixe = FIXP19_TO_INT(xr);


			if (ixe < cliper.left || ixs>cliper.right) {

				xl += dxdyl;
				xr += dxdyr;

				zl += dzdyl;
				zr += dzdyr;

				ul += dudyl;
				ur += dudyr;

				vl += dvdyl;
				vr += dvdyr;


				rl += drdyl;
				rr += drdyr;

				gl += dgdyl;
				gr += dgdyr;

				bl += dbdyl;
				br += dbdyr;


				destBuffer += destBufferLingLength;
				destZBuffer += zBufferLineLength;
				continue;
			}

			idx = ixe - ixs;
			if (idx>0) {
				dzdx = (zr - zl) / idx;
				dudx = (ur - ul) / idx;
				dvdx = (vr - vl) / idx;

				drdx = (rr - rl) / idx;
				dgdx = (gr - gl) / idx;
				dbdx = (br - bl) / idx;

			}
			else
			{
				dudx = (ur - ul);
				dvdx = (vr - vl);
				dzdx = (zr - zl);

				drdx = (rr - rl);
				dgdx = (gr - gl);
				dbdx = (br - gl);
			}


			invzU = ul;
			invzV = vl;
			invz = zl;

			invzR = rl;
			invzG = gl;
			invzB = bl;


			if (ixs< cliper.left) {
				idx = cliper.left - ixs;

				invzU += idx*dudx;
				invzV += idx*dvdx;
				invz += idx*dzdx;

				invzR += idx*drdx;
				invzG += idx*dgdx;
				invzB += idx*dbdx;

				ixs = cliper.left;
			}
			if (ixe > cliper.right) {
				ixe = cliper.right;
			}


			for (int xidx = ixs;xidx < ixe;xidx++) {
				if (invz > destZBuffer[xidx]) {
					invzFixp20 = FIXP28_TO_FIXP20(invz);

					iu = (invzU / invzFixp20);
					iv = (invzV / invzFixp20) << textShift;

					ir = invzR / invzFixp20;
					ig = invzG / invzFixp20;
					ib = invzB / invzFixp20;
					drawColor = texture[iu + iv];
					RGB565_TO_R5_G6_B5(drawColor, rDraw, gDraw, bDraw);
					rDraw = (rDraw*ir) >> 8;
					gDraw = (gDraw*ig) >> 8;
					bDraw = (bDraw*ib) >> 8;
					drawColor = R5_G6_B5_TO_RGB565(rDraw, gDraw, bDraw);

					destBuffer[xidx] = drawColor;
					destZBuffer[xidx] = invz;
				}

				invz += dzdx;
				invzU += dudx;
				invzV += dvdx;

				invzR += drdx;
				invzG += dgdx;
				invzB += dbdx;
			}

			xl += dxdyl;
			xr += dxdyr;

			zl += dzdyl;
			zr += dzdyr;

			ul += dudyl;
			ur += dudyr;

			vl += dvdyl;
			vr += dvdyr;


			rl += drdyl;
			rr += drdyr;

			gl += dgdyl;
			gr += dgdyr;

			bl += dbdyl;
			br += dbdyr;




			destBuffer += destBufferLingLength;
			destZBuffer += zBufferLineLength;
		}
	}
	else
	{
		for (int yidx = iys;yidx <= iye;yidx++) {

			ixs = FIXP19_TO_INT(xl);
			ixe = FIXP19_TO_INT(xr);

			idx = ixe - ixs;
			if (idx>0) {
				dzdx = (zr - zl) / idx;
				dudx = (ur - ul) / idx;
				dvdx = (vr - vl) / idx;

				drdx = (rr - rl) / idx;
				dgdx = (gr - gl) / idx;
				dbdx = (br - bl) / idx;
			}
			else
			{
				dzdx = (zr - zl);
				dudx = (ur - ul);
				dvdx = (vr - vl);

				drdx = (rr - rl);
				dgdx = (gr - gl);
				dbdx = (br - bl);
			}


			invzU = ul;
			invzV = vl;
			invz = zl;


			invzR = rl;
			invzG = gl;
			invzB = bl;


			for (int xidx = ixs;xidx < ixe;xidx++) {
				if (invz > destZBuffer[xidx]) {
					invzFixp20 = FIXP28_TO_FIXP20(invz);
					iu = (invzU / invzFixp20);
					iv = (invzV / invzFixp20) << textShift;

					ir = invzR / invzFixp20;
					ig = invzG / invzFixp20;
					ib = invzB / invzFixp20;
					drawColor = texture[iu + iv];
					RGB565_TO_R5_G6_B5(drawColor, rDraw, gDraw, bDraw);
					rDraw = (rDraw*ir) >> 8;
					gDraw = (gDraw*ig) >> 8;
					bDraw = (bDraw*ib) >> 8;
					drawColor = R5_G6_B5_TO_RGB565(rDraw, gDraw, bDraw);
					destBuffer[xidx] = drawColor;

					destZBuffer[xidx] = invz;
				}

				invz += dzdx;
				invzU += dudx;
				invzV += dvdx;


				invzR += drdx;
				invzG += dgdx;
				invzB += dbdx;
			}


			xl += dxdyl;
			xr += dxdyr;

			zl += dzdyl;
			zr += dzdyr;

			ul += dudyl;
			ur += dudyr;

			vl += dvdyl;
			vr += dvdyr;


			rl += drdyl;
			rr += drdyr;

			gl += dgdyl;
			gr += dgdyr;

			bl += dbdyl;
			br += dbdyr;



			destBuffer += destBufferLingLength;
			destZBuffer += zBufferLineLength;
		}
	}


	return 1;
}

int Canvas::drawTriangleTextureFlatInvZ16(float x0, float y0, float z0, float u0, float v0, float x1, float y1, float z1, float u1, float v1, float x2, float y2, float z2, float u2, float v2, Color c, uint16 * texture, int width, int height, int textPitch)
{
	float  tmp, interX, interZ, interU, interV, fraction;
	//line so return 0 as trivially reject
	if ((FEQUALE1(y0, y1) && FEQUALE1(y1, y2)) || (FEQUALE1(x0, x1) && FEQUALE1(x1, x2))) {
		return 0;
	}

	//select sort for y small -> big
	if (y1 < y0) {
		SWAP(y1, y0, tmp);
		SWAP(x1, x0, tmp);
		SWAP(z1, z0, tmp);
		SWAP(u1, u0, tmp);
		SWAP(v1, v0, tmp);
	}
	if (y2 < y0) {
		SWAP(y2, y0, tmp);
		SWAP(x2, x0, tmp);
		SWAP(z2, z0, tmp);
		SWAP(u2, u0, tmp);
		SWAP(v2, v0, tmp);
	}
	if (y2 < y1) {
		SWAP(y1, y2, tmp);
		SWAP(x1, x2, tmp);
		SWAP(z1, z2, tmp);
		SWAP(u1, u2, tmp);
		SWAP(v1, v2, tmp);
	}

	//trivially reject
	if (y0 > cliper.bottom || y2<cliper.top || (x0<cliper.left&&x1<cliper.left&&x2<cliper.left) || (x0>cliper.right&&x1>cliper.right&&x2>cliper.right)) {
		return 0;
	}

	if (FEQUALE1(y0, y1)) {
		drawTopTriangleTextureFlatInvZ16(x2, y2, z2, u2, v2,
			x0, y0, z0, u0, v0,
			x1, y1, z1, u1, v1,
			c,
			texture, width, height, textPitch);
	}
	else if (FEQUALE1(y1, y2))
	{
		drawBottomTriangleTextureFlatInvZ16(x0, y0, z0, u0, v0,
			x1, y1, z1, u1, v1,
			x2, y2, z2, u2, v2,
			c,
			texture, width, height, textPitch);
	}
	else
	{
		fraction = (y1 - y0) / (y2 - y0);
		interX = x0 + (x2 - x0)*fraction;
		interZ = z0 + (z2 - z0)*fraction;
		interU = u0 + (u2 - u0)*fraction;
		interV = v0 + (v2 - v0)*fraction;
		drawBottomTriangleTextureFlatInvZ16(x0, y0, z0, u0, v0,
			x1, y1, z1, u1, v1,
			interX, y1, interZ, interU, interV,
			c,
			texture, width, height, textPitch);

		drawTopTriangleTextureFlatInvZ16(x2, y2, z2, u2, v2,
			x1, y1, z1, u1, v1,
			interX, y1, interZ, interU, interV,
			c,
			texture, width, height, textPitch);
	}

	return 1;
}

int Canvas::drawTopTriangleTextureFlatInvZ16(float x0, float y0, float z0, float u0, float v0, float x1, float y1, float z1, float u1, float v1, float x2, float y2, float z2, float u2, float v2, Color c, uint16 * texture, int width, int height, int textPitch)
{
	fixp19 xl, xr, dxdyl, dxdyr;

	fixp20  invzFixp20,
		ul, ur, vl, vr,
		dudyl, dudyr, dvdyl, dvdyr,
		invzU0, invzV0, invzU1, invzV1, invzU2, invzV2, invzU, invzV,
		dudx, dvdx;


	fixp28 invz0, invz1, invz2, invz, zl, zr, dzdyl, dzdyr, dzdx;

	int idx, idy,
		iys, iye, ixs, ixe,
		iu, iv,
		rBase,bBase,gBase,drawColor,
		rDraw,gDraw,bDraw,
		tmpi, textShift, destBufferLingLength, zBufferLineLength;

	float dx, dy, tmpf;
	uint16 *destBuffer;
	uint32 *destZBuffer;

	destBufferLingLength = pitch >> 1;
	zBufferLineLength = zBuffer->pitch >> 2;//32 bit zbuffer wo fix  zBuffer deepth
	textShift = Log_Base2(width);

	rBase = c.r;
	gBase = c.g;
	bBase = c.b;

	if (x2 < x1) {
		SWAP(x2, x1, tmpf);
		SWAP(y2, y1, tmpf);
		SWAP(z2, z1, tmpf);
		SWAP(u2, u1, tmpf);
		SWAP(v2, v1, tmpf);
	}

	invz0 = FIXP28_FLOAT_TO(1.0f / z0);
	invz1 = FIXP28_FLOAT_TO(1.0f / z1);
	invz2 = FIXP28_FLOAT_TO(1.0f / z2);

	tmpf = width - 1;
	invzU0 = FIXP20_FLOAT_TO((u0*tmpf) / z0);
	invzU1 = FIXP20_FLOAT_TO((u1*tmpf) / z1);
	invzU2 = FIXP20_FLOAT_TO((u2*tmpf) / z2);

	tmpf = height - 1;
	invzV0 = FIXP20_FLOAT_TO((v0*tmpf) / z0);
	invzV1 = FIXP20_FLOAT_TO((v1*tmpf) / z1);
	invzV2 = FIXP20_FLOAT_TO((v2*tmpf) / z2);





	//drawTraingleGouraud16 has ignored  y0 == y1 == y2 so just div
	dy = y0 - y1;


	dxdyl = FIXP19_FLOAT_TO((x0 - x1) / dy);

	dzdyl = (invz0 - invz1) / dy;
	dudyl = (invzU0 - invzU1) / dy;
	dvdyl = (invzV0 - invzV1) / dy;


	dy = y0 - y2;
	dxdyr = FIXP19_FLOAT_TO((x0 - x2) / dy);
	dzdyr = (invz0 - invz2) / dy;
	dudyr = (invzU0 - invzU2) / dy;
	dvdyr = (invzV0 - invzV2) / dy;


	//clip bottom
	if (y0>cliper.bottom)
	{
		iye = cliper.bottom;
	}
	else
	{
		iye = ceil(y0) - 1;
	}

	//clip top
	if (y1 < cliper.top) {
		iys = cliper.top;

	}
	else
	{
		iys = ceil(y1);
	}

	dy = iys - y1;

	xl = FIXP19_FLOAT_TO(x1) + dy*dxdyl;
	zl = invz1 + dy*dzdyl;
	ul = invzU1 + dy*dudyl;
	vl = invzV1 + dy*dvdyl;

	dy = iys - y2;

	xr = FIXP19_FLOAT_TO(x2) + dy*dxdyr;
	zr = invz2 + dy*dzdyr;
	ur = invzU2 + dy*dudyr;
	vr = invzV2 + dy*dvdyr;





	//draw

	//pre round up

	xl = FIXP19_ROUNDUP(xl);
	xr = FIXP19_ROUNDUP(xr);



	destBuffer = (uint16 *)drawBuffer + iys*destBufferLingLength;
	destZBuffer = zBuffer->zBuffer + iys*zBufferLineLength;

	if ((x0 < cliper.left) || (x0 > cliper.right) ||
		(x1 <cliper.left) || (x1 > cliper.right) ||
		(x2 < cliper.left) || (x2 >cliper.right)) {
		//should clip


		for (int yidx = iys;yidx <= iye;yidx++) {

			ixs = FIXP19_TO_INT(xl);
			ixe = FIXP19_TO_INT(xr);


			if (ixe < cliper.left || ixs>cliper.right) {

				xl += dxdyl;
				xr += dxdyr;

				zl += dzdyl;
				zr += dzdyr;

				ul += dudyl;
				ur += dudyr;

				vl += dvdyl;
				vr += dvdyr;



				destBuffer += destBufferLingLength;
				destZBuffer += zBufferLineLength;
				continue;
			}

			idx = ixe - ixs;
			if (idx>0) {
				dudx = (ur - ul) / idx;
				dvdx = (vr - vl) / idx;

				dzdx = (zr - zl) / idx;
			}
			else
			{
				dudx = (ur - ul);
				dvdx = (vr - vl);
				dzdx = (zr - zl);
			}


			invzU = ul;
			invzV = vl;
			invz = zl;


			if (ixs< cliper.left) {
				idx = cliper.left - ixs;
				invzU += idx*dudx;
				invzV += idx*dvdx;
				invz += idx*dzdx;

				ixs = cliper.left;
			}
			if (ixe > cliper.right) {
				ixe = cliper.right;
			}


			for (int xidx = ixs;xidx < ixe;xidx++) {
				if (invz > destZBuffer[xidx]) {
					invzFixp20 = FIXP28_TO_FIXP20(invz);
					iu = (invzU / invzFixp20);
					iv = (invzV / invzFixp20) << textShift;

					drawColor=texture[iu + iv];
					RGB565_TO_R5_G6_B5(drawColor,rDraw,gDraw,bDraw);
					rDraw = (rDraw*rBase) >> 8;
					gDraw = (gDraw*gBase) >> 8;
					bDraw = (bDraw*bBase) >> 8;
					drawColor=R5_G6_B5_TO_RGB565(rDraw, gDraw, bDraw);
					destBuffer[xidx] = drawColor;

					destZBuffer[xidx] = invz;
				}

				invzU += dudx;
				invzV += dvdx;
				invz += dzdx;
			}

			xl += dxdyl;
			xr += dxdyr;

			zl += dzdyl;
			zr += dzdyr;

			ul += dudyl;
			ur += dudyr;

			vl += dvdyl;
			vr += dvdyr;



			destBuffer += destBufferLingLength;
			destZBuffer += zBufferLineLength;
		}
	}
	else
	{
		for (int yidx = iys;yidx <= iye;yidx++) {

			ixs = FIXP19_TO_INT(xl);
			ixe = FIXP19_TO_INT(xr);

			idx = ixe - ixs;
			if (idx>0) {
				dudx = (ur - ul) / idx;
				dvdx = (vr - vl) / idx;

				dzdx = (zr - zl) / idx;
			}
			else
			{
				dudx = (ur - ul);
				dvdx = (vr - vl);
				dzdx = (zr - zl);
			}


			invzU = ul;
			invzV = vl;
			invz = zl;

			for (int xidx = ixs;xidx < ixe;xidx++) {
				if (invz > destZBuffer[xidx]) {
					invzFixp20 = FIXP28_TO_FIXP20(invz);
					iu = (invzU / invzFixp20);
					iv = (invzV / invzFixp20) << textShift;

					drawColor = texture[iu + iv];
					RGB565_TO_R5_G6_B5(drawColor, rDraw, gDraw, bDraw);
					rDraw = (rDraw*rBase) >> 8;
					gDraw = (gDraw*gBase) >> 8;
					bDraw = (bDraw*bBase) >> 8;
					drawColor = R5_G6_B5_TO_RGB565(rDraw, gDraw, bDraw);
					destBuffer[xidx] = drawColor;

					destZBuffer[xidx] = invz;
				}

				invzU += dudx;
				invzV += dvdx;
				invz += dzdx;
			}

			xl += dxdyl;
			xr += dxdyr;

			zl += dzdyl;
			zr += dzdyr;

			ul += dudyl;
			ur += dudyr;

			vl += dvdyl;
			vr += dvdyr;



			destBuffer += destBufferLingLength;
			destZBuffer += zBufferLineLength;
		}
	}


	return 1;
}

int Canvas::drawBottomTriangleTextureFlatInvZ16(float x0, float y0, float z0, float u0, float v0, float x1, float y1, float z1, float u1, float v1, float x2, float y2, float z2, float u2, float v2, Color c, uint16 * texture, int width, int height, int textPitch)
{
	fixp19 xl, xr, dxdyl, dxdyr;

	fixp20  invzFixp20,
		ul, ur, vl, vr,
		dudyl, dudyr, dvdyl, dvdyr,
		invzU0, invzV0, invzU1, invzV1, invzU2, invzV2, invzU, invzV,
		dudx, dvdx;


	fixp28 invz0, invz1, invz2, invz, zl, zr, dzdyl, dzdyr, dzdx;

	int idx, idy,
		iys, iye, ixs, ixe,
		iu, iv,
		rBase, bBase, gBase, drawColor,
		rDraw, gDraw, bDraw,
		tmpi, textShift, destBufferLingLength, zBufferLineLength;

	float dx, dy, tmpf;
	uint16 *destBuffer;
	uint32 *destZBuffer;

	destBufferLingLength = pitch >> 1;
	zBufferLineLength = zBuffer->pitch >> 2;//32 bit zbuffer wo fix  zBuffer deepth
	textShift = Log_Base2(width);


	rBase = c.r;
	gBase = c.g;
	bBase = c.b;


	if (x2 < x1) {
		SWAP(x2, x1, tmpf);
		SWAP(y2, y1, tmpf);
		SWAP(z2, z1, tmpf);
		SWAP(u2, u1, tmpf);
		SWAP(v2, v1, tmpf);
	}

	invz0 = FIXP28_FLOAT_TO(1.0f / z0);
	invz1 = FIXP28_FLOAT_TO(1.0f / z1);
	invz2 = FIXP28_FLOAT_TO(1.0f / z2);

	tmpf = width - 1;
	invzU0 = FIXP20_FLOAT_TO((u0*tmpf) / z0);
	invzU1 = FIXP20_FLOAT_TO((u1*tmpf) / z1);
	invzU2 = FIXP20_FLOAT_TO((u2*tmpf) / z2);

	tmpf = height - 1;
	invzV0 = FIXP20_FLOAT_TO((v0*tmpf) / z0);
	invzV1 = FIXP20_FLOAT_TO((v1*tmpf) / z1);
	invzV2 = FIXP20_FLOAT_TO((v2*tmpf) / z2);





	//drawTraingleGouraud16 has ignored  y0 == y1 == y2 so just div


	dy = y1 - y0;
	dxdyl = FIXP19_FLOAT_TO((x1 - x0) / dy);
	dzdyl = (invz1 - invz0) / dy;
	dudyl = (invzU1 - invzU0) / dy;
	dvdyl = (invzV1 - invzV0) / dy;



	dy = y2 - y0;
	dxdyr = FIXP19_FLOAT_TO((x2 - x0) / dy);
	dzdyr = (invz2 - invz0) / dy;
	dudyr = (invzU2 - invzU0) / dy;
	dvdyr = (invzV2 - invzV0) / dy;


	//clip bottom
	if (y1>cliper.bottom)
	{
		iye = cliper.bottom;
	}
	else
	{
		iye = ceil(y1) - 1;
	}

	//clip top
	if (y0 < cliper.top) {
		iys = cliper.top;

	}
	else
	{
		iys = ceil(y0);
	}

	dy = iys - y0;

	xl = FIXP19_FLOAT_TO(x0) + dy*dxdyl;
	zl = invz0 + dy*dzdyl;
	ul = invzU0 + dy*dudyl;
	vl = invzV0 + dy*dvdyl;


	xr = FIXP19_FLOAT_TO(x0) + dy*dxdyr;
	zr = invz0 + dy*dzdyr;
	ur = invzU0 + dy*dudyr;
	vr = invzV0 + dy*dvdyr;





	//draw

	//pre round up

	xl = FIXP19_ROUNDUP(xl);
	xr = FIXP19_ROUNDUP(xr);



	destBuffer = (uint16 *)drawBuffer + iys*destBufferLingLength;
	destZBuffer = zBuffer->zBuffer + iys*zBufferLineLength;

	if ((x0 < cliper.left) || (x0 > cliper.right) ||
		(x1 <cliper.left) || (x1 > cliper.right) ||
		(x2 < cliper.left) || (x2 >cliper.right)) {
		//should clip


		for (int yidx = iys;yidx <= iye;yidx++) {

			ixs = FIXP19_TO_INT(xl);
			ixe = FIXP19_TO_INT(xr);


			if (ixe < cliper.left || ixs>cliper.right) {

				xl += dxdyl;
				xr += dxdyr;

				zl += dzdyl;
				zr += dzdyr;

				ul += dudyl;
				ur += dudyr;

				vl += dvdyl;
				vr += dvdyr;



				destBuffer += destBufferLingLength;
				destZBuffer += zBufferLineLength;
				continue;
			}

			idx = ixe - ixs;
			if (idx>0) {
				dudx = (ur - ul) / idx;
				dvdx = (vr - vl) / idx;

				dzdx = (zr - zl) / idx;
			}
			else
			{
				dudx = (ur - ul);
				dvdx = (vr - vl);
				dzdx = (zr - zl);
			}


			invzU = ul;
			invzV = vl;
			invz = zl;


			if (ixs< cliper.left) {
				idx = cliper.left - ixs;
				invzU += idx*dudx;
				invzV += idx*dvdx;
				invz += idx*dzdx;

				ixs = cliper.left;
			}
			if (ixe > cliper.right) {
				ixe = cliper.right;
			}


			for (int xidx = ixs;xidx < ixe;xidx++) {
				if (invz > destZBuffer[xidx]) {
					invzFixp20 = FIXP28_TO_FIXP20(invz);
					iu = (invzU / invzFixp20);
					iv = (invzV / invzFixp20) << textShift;

					drawColor = texture[iu + iv];
					RGB565_TO_R5_G6_B5(drawColor, rDraw, gDraw, bDraw);
					rDraw = (rDraw*rBase) >> 8;
					gDraw = (gDraw*gBase) >> 8;
					bDraw = (bDraw*bBase) >> 8;
					drawColor = R5_G6_B5_TO_RGB565(rDraw, gDraw, bDraw);
					destBuffer[xidx] = drawColor;

					destZBuffer[xidx] = invz;
				}

				invzU += dudx;
				invzV += dvdx;
				invz += dzdx;
			}

			xl += dxdyl;
			xr += dxdyr;

			zl += dzdyl;
			zr += dzdyr;

			ul += dudyl;
			ur += dudyr;

			vl += dvdyl;
			vr += dvdyr;



			destBuffer += destBufferLingLength;
			destZBuffer += zBufferLineLength;
		}
	}
	else
	{
		for (int yidx = iys;yidx <= iye;yidx++) {

			ixs = FIXP19_TO_INT(xl);
			ixe = FIXP19_TO_INT(xr);

			idx = ixe - ixs;
			if (idx>0) {
				dudx = (ur - ul) / idx;
				dvdx = (vr - vl) / idx;

				dzdx = (zr - zl) / idx;
			}
			else
			{
				dudx = (ur - ul);
				dvdx = (vr - vl);
				dzdx = (zr - zl);
			}


			invzU = ul;
			invzV = vl;
			invz = zl;

			for (int xidx = ixs;xidx < ixe;xidx++) {
				if (invz > destZBuffer[xidx]) {
					invzFixp20 = FIXP28_TO_FIXP20(invz);
					iu = (invzU / invzFixp20);
					iv = (invzV / invzFixp20) << textShift;

					drawColor = texture[iu + iv];
					RGB565_TO_R5_G6_B5(drawColor, rDraw, gDraw, bDraw);
					rDraw = (rDraw*rBase) >> 8;
					gDraw = (gDraw*gBase) >> 8;
					bDraw = (bDraw*bBase) >> 8;
					drawColor = R5_G6_B5_TO_RGB565(rDraw, gDraw, bDraw);
					destBuffer[xidx] = drawColor;

					destZBuffer[xidx] = invz;
				}

				invzU += dudx;
				invzV += dvdx;
				invz += dzdx;
			}

			xl += dxdyl;
			xr += dxdyr;

			zl += dzdyl;
			zr += dzdyr;

			ul += dudyl;
			ur += dudyr;

			vl += dvdyl;
			vr += dvdyr;



			destBuffer += destBufferLingLength;
			destZBuffer += zBufferLineLength;
		}
	}


	return 1;
}

int Canvas::drawTriangleTextureInvZ16(float x0, float y0, float z0, float u0, float v0, float x1, float y1, float z1, float u1, float v1, float x2, float y2, float z2, float u2, float v2, uint16 * texture, int width, int height, int textPitch)
{
	float  tmp, interX,interZ, interU, interV, fraction;
	//line so return 0 as trivially reject
	if ((FEQUALE1(y0, y1) && FEQUALE1(y1, y2)) || (FEQUALE1(x0, x1) && FEQUALE1(x1, x2))) {
		return 0;
	}

	//select sort for y small -> big
	if (y1 < y0) {
		SWAP(y1, y0, tmp);
		SWAP(x1, x0, tmp);
		SWAP(z1, z0, tmp);
		SWAP(u1, u0, tmp);
		SWAP(v1, v0, tmp);
	}
	if (y2 < y0) {
		SWAP(y2, y0, tmp);
		SWAP(x2, x0, tmp);
		SWAP(z2, z0, tmp);
		SWAP(u2, u0, tmp);
		SWAP(v2, v0, tmp);
	}
	if (y2 < y1) {
		SWAP(y1, y2, tmp);
		SWAP(x1, x2, tmp);
		SWAP(z1, z2, tmp);
		SWAP(u1, u2, tmp);
		SWAP(v1, v2, tmp);
	}

	//trivially reject
	if (y0 > cliper.bottom || y2<cliper.top || (x0<cliper.left&&x1<cliper.left&&x2<cliper.left) || (x0>cliper.right&&x1>cliper.right&&x2>cliper.right)) {
		return 0;
	}

	if (FEQUALE1(y0, y1)) {
		drawTopTriangleTextureInvZ16(x2, y2, z2,u2, v2,
			x0, y0,z0, u0, v0,
			x1, y1, z1,u1, v1,
			texture, width, height, textPitch);
	}
	else if (FEQUALE1(y1, y2))
	{
		drawBottomTriangleTextureInvZ16(x0, y0,z0, u0, v0,
			x1, y1,z1, u1, v1,
			x2, y2,z2, u2, v2,
			texture, width, height, textPitch);
	}
	else
	{
		fraction = (y1 - y0) / (y2 - y0);
		interX = x0 + (x2 - x0)*fraction;
		interZ = z0 + (z2 - z0)*fraction;
		interU = u0 + (u2 - u0)*fraction;
		interV = v0 + (v2 - v0)*fraction;
		drawBottomTriangleTextureInvZ16(x0, y0,z0 ,u0, v0,
			x1, y1,z1, u1, v1,
			interX, y1,interZ, interU, interV,
			texture, width, height, textPitch);

		drawTopTriangleTextureInvZ16(x2, y2, z2,u2, v2,
			x1, y1, z1,u1, v1,
			interX, y1,interZ, interU, interV,
			texture, width, height, textPitch);
	}

	return 1;
}

int Canvas::drawTopTriangleTextureInvZ16(float x0, float y0, float z0, float u0, float v0, float x1, float y1, float z1, float u1, float v1, float x2, float y2, float z2, float u2, float v2, uint16 * texture, int width, int height, int textPitch)
{
	fixp19 xl, xr, dxdyl, dxdyr;

	fixp20  invzFixp20,
		ul, ur, vl, vr,
		dudyl, dudyr, dvdyl, dvdyr,
		invzU0, invzV0, invzU1, invzV1,invzU2,invzV2, invzU, invzV,
		dudx, dvdx;


	fixp28 invz0, invz1, invz2, invz, zl, zr, dzdyl, dzdyr, dzdx;

	int idx, idy,
		iys, iye, ixs, ixe,
		iu, iv,
		tmpi,textShift, destBufferLingLength, zBufferLineLength;

	float dx, dy, tmpf;
	uint16 *destBuffer;
	uint32 *destZBuffer;

	destBufferLingLength = pitch >> 1;
	zBufferLineLength = zBuffer->pitch >> 2;//32 bit zbuffer wo fix  zBuffer deepth
	textShift = Log_Base2(width);



	if (x2 < x1) {
		SWAP(x2, x1, tmpf);
		SWAP(y2, y1, tmpf);
		SWAP(z2, z1, tmpf);
		SWAP(u2, u1, tmpf);
		SWAP(v2, v1, tmpf);
	}

	invz0 = FIXP28_FLOAT_TO(1.0f / z0);
	invz1 = FIXP28_FLOAT_TO(1.0f / z1);
	invz2 = FIXP28_FLOAT_TO(1.0f / z2);

	tmpf = width - 1;
	invzU0 = FIXP20_FLOAT_TO((u0*tmpf)/ z0);
	invzU1 = FIXP20_FLOAT_TO((u1*tmpf) / z1);
	invzU2 = FIXP20_FLOAT_TO((u2*tmpf) / z2);

	tmpf = height - 1;
	invzV0 = FIXP20_FLOAT_TO((v0*tmpf) / z0);
	invzV1 = FIXP20_FLOAT_TO((v1*tmpf) / z1);
	invzV2 = FIXP20_FLOAT_TO((v2*tmpf) / z2);





	//drawTraingleGouraud16 has ignored  y0 == y1 == y2 so just div
	dy = y0 - y1;

	
	dxdyl = FIXP19_FLOAT_TO((x0 - x1) / dy);

	dzdyl = (invz0 - invz1) / dy;
	dudyl = (invzU0 - invzU1) / dy;
	dvdyl = (invzV0 - invzV1) / dy;


	dy = y0 - y2;
	dxdyr = FIXP19_FLOAT_TO((x0 - x2) / dy);
	dzdyr = (invz0 - invz2) / dy;
	dudyr = (invzU0 - invzU2) / dy;
	dvdyr = (invzV0 - invzV2) / dy;


	//clip bottom
	if (y0>cliper.bottom)
	{
		iye = cliper.bottom;
	}
	else
	{
		iye = ceil(y0) - 1;
	}

	//clip top
	if (y1 < cliper.top) {
		iys = cliper.top;

	}
	else
	{
		iys = ceil(y1);
	}

	dy = iys - y1;

	xl = FIXP19_FLOAT_TO(x1) + dy*dxdyl;
	zl = invz1 + dy*dzdyl;
	ul = invzU1 + dy*dudyl;
	vl = invzV1 + dy*dvdyl;

	dy = iys - y2;

	xr = FIXP19_FLOAT_TO(x2) + dy*dxdyr;
	zr = invz2 + dy*dzdyr;
	ur = invzU2 + dy*dudyr;
	vr = invzV2 + dy*dvdyr;





	//draw

	//pre round up

	xl = FIXP19_ROUNDUP(xl);
	xr = FIXP19_ROUNDUP(xr);



	destBuffer = (uint16 *)drawBuffer + iys*destBufferLingLength;
	destZBuffer = zBuffer->zBuffer + iys*zBufferLineLength;

	if ((x0 < cliper.left) || (x0 > cliper.right) ||
		(x1 <cliper.left) || (x1 > cliper.right) ||
		(x2 < cliper.left) || (x2 >cliper.right)) {
		//should clip


		for (int yidx = iys;yidx <= iye;yidx++) {
		
			ixs = FIXP19_TO_INT(xl);
			ixe = FIXP19_TO_INT(xr);


			if (ixe < cliper.left || ixs>cliper.right) {

				xl += dxdyl;
				xr += dxdyr;

				zl += dzdyl;
				zr += dzdyr;

				ul += dudyl;
				ur += dudyr;

				vl += dvdyl;
				vr += dvdyr;



				destBuffer += destBufferLingLength;
				destZBuffer += zBufferLineLength;
				continue;
			}

			idx = ixe - ixs;
			if (idx>0) {
				dudx = (ur - ul) / idx;
				dvdx = (vr - vl) / idx;
		
				dzdx = (zr - zl) / idx;
			}
			else
			{
				dudx = (ur - ul);
				dvdx = (vr - vl);
				dzdx = (zr - zl);
			}


			invzU = ul;
			invzV = vl;
			invz = zl;


			if (ixs< cliper.left) {
				idx = cliper.left - ixs;
				invzU += idx*dudx;
				invzV += idx*dvdx;
				invz += idx*dzdx;

				ixs = cliper.left;
			}
			if (ixe > cliper.right) {
				ixe = cliper.right;
			}


			for (int xidx = ixs;xidx < ixe;xidx++) {
				if (invz > destZBuffer[xidx]) {
					invzFixp20 = FIXP28_TO_FIXP20(invz);
					iu = (invzU / invzFixp20);
					iv = (invzV / invzFixp20)<<textShift;
			
					destBuffer[xidx] = texture[iu+iv];
					destZBuffer[xidx] = invz;
				}

				invzU += dudx;
				invzV += dvdx;
				invz += dzdx;
			}

			xl += dxdyl;
			xr += dxdyr;

			zl += dzdyl;
			zr += dzdyr;

			ul += dudyl;
			ur += dudyr;

			vl += dvdyl;
			vr += dvdyr;



			destBuffer += destBufferLingLength;
			destZBuffer += zBufferLineLength;
		}
	}
	else
	{
		for (int yidx = iys;yidx <= iye;yidx++) {

			ixs = FIXP19_TO_INT(xl);
			ixe = FIXP19_TO_INT(xr);

			idx = ixe - ixs;
			if (idx>0) {
				dudx = (ur - ul) / idx;
				dvdx = (vr - vl) / idx;

				dzdx = (zr - zl) / idx;
			}
			else
			{
				dudx = (ur - ul);
				dvdx = (vr - vl);
				dzdx = (zr - zl);
			}


			invzU = ul;
			invzV = vl;
			invz = zl;

			for (int xidx = ixs;xidx < ixe;xidx++) {
				if (invz > destZBuffer[xidx]) {
					invzFixp20 = FIXP28_TO_FIXP20(invz);
					iu = (invzU / invzFixp20);
					iv = (invzV / invzFixp20) << textShift;

					destBuffer[xidx] = texture[iu + iv];
					destZBuffer[xidx] = invz;
				}

				invzU += dudx;
				invzV += dvdx;
				invz += dzdx;
			}

			xl += dxdyl;
			xr += dxdyr;

			zl += dzdyl;
			zr += dzdyr;

			ul += dudyl;
			ur += dudyr;

			vl += dvdyl;
			vr += dvdyr;



			destBuffer += destBufferLingLength;
			destZBuffer += zBufferLineLength;
		}
	}


	return 1;
}

int Canvas::drawBottomTriangleTextureInvZ16(float x0, float y0, float z0, float u0, float v0, float x1, float y1, float z1, float u1, float v1, float x2, float y2, float z2, float u2, float v2, uint16 * texture, int width, int height, int textPitch)
{
	fixp19 xl, xr, dxdyl, dxdyr;

	fixp20  invzFixp20,
		ul, ur, vl, vr,
		dudyl, dudyr, dvdyl, dvdyr,
		invzU0, invzV0, invzU1, invzV1, invzU2, invzV2, invzU, invzV,
		dudx, dvdx;


	fixp28 invz0, invz1, invz2, invz, zl, zr, dzdyl, dzdyr, dzdx;

	int idx, idy,
		iys, iye, ixs, ixe,
		iu, iv,
		tmpi, textShift, destBufferLingLength, zBufferLineLength;

	float dx, dy, tmpf;
	uint16 *destBuffer;
	uint32 *destZBuffer;

	destBufferLingLength = pitch >> 1;
	zBufferLineLength = zBuffer->pitch >> 2;//32 bit zbuffer wo fix  zBuffer deepth
	textShift = Log_Base2(width);



	if (x2 < x1) {
		SWAP(x2, x1, tmpf);
		SWAP(y2, y1, tmpf);
		SWAP(z2, z1, tmpf);
		SWAP(u2, u1, tmpf);
		SWAP(v2, v1, tmpf);
	}

	invz0 = FIXP28_FLOAT_TO(1.0f / z0);
	invz1 = FIXP28_FLOAT_TO(1.0f / z1);
	invz2 = FIXP28_FLOAT_TO(1.0f / z2);

	tmpf = width - 1;
	invzU0 = FIXP20_FLOAT_TO((u0*tmpf) / z0);
	invzU1 = FIXP20_FLOAT_TO((u1*tmpf) / z1);
	invzU2 = FIXP20_FLOAT_TO((u2*tmpf) / z2);

	tmpf = height - 1;
	invzV0 = FIXP20_FLOAT_TO((v0*tmpf) / z0);
	invzV1 = FIXP20_FLOAT_TO((v1*tmpf) / z1);
	invzV2 = FIXP20_FLOAT_TO((v2*tmpf) / z2);





	//drawTraingleGouraud16 has ignored  y0 == y1 == y2 so just div


	dy = y1 - y0;
	dxdyl = FIXP19_FLOAT_TO((x1 - x0) / dy);
	dzdyl = (invz1 - invz0) / dy;
	dudyl = (invzU1 - invzU0) / dy;
	dvdyl = (invzV1 - invzV0) / dy;



	dy = y2 - y0;
	dxdyr = FIXP19_FLOAT_TO((x2 - x0) / dy);
	dzdyr = (invz2 - invz0) / dy;
	dudyr = (invzU2 - invzU0) / dy;
	dvdyr = (invzV2 - invzV0) / dy;


	//clip bottom
	if (y1>cliper.bottom)
	{
		iye = cliper.bottom;
	}
	else
	{
		iye = ceil(y1) - 1;
	}

	//clip top
	if (y0 < cliper.top) {
		iys = cliper.top;

	}
	else
	{
		iys = ceil(y0);
	}

	dy = iys - y0;

	xl = FIXP19_FLOAT_TO(x0) + dy*dxdyl;
	zl = invz0 + dy*dzdyl;
	ul = invzU0 + dy*dudyl;
	vl = invzV0 + dy*dvdyl;


	xr = FIXP19_FLOAT_TO(x0) + dy*dxdyr;
	zr = invz0 + dy*dzdyr;
	ur = invzU0 + dy*dudyr;
	vr = invzV0 + dy*dvdyr;





	//draw

	//pre round up

	xl = FIXP19_ROUNDUP(xl);
	xr = FIXP19_ROUNDUP(xr);



	destBuffer = (uint16 *)drawBuffer + iys*destBufferLingLength;
	destZBuffer = zBuffer->zBuffer + iys*zBufferLineLength;

	if ((x0 < cliper.left) || (x0 > cliper.right) ||
		(x1 <cliper.left) || (x1 > cliper.right) ||
		(x2 < cliper.left) || (x2 >cliper.right)) {
		//should clip


		for (int yidx = iys;yidx <= iye;yidx++) {

			ixs = FIXP19_TO_INT(xl);
			ixe = FIXP19_TO_INT(xr);


			if (ixe < cliper.left || ixs>cliper.right) {

				xl += dxdyl;
				xr += dxdyr;

				zl += dzdyl;
				zr += dzdyr;

				ul += dudyl;
				ur += dudyr;

				vl += dvdyl;
				vr += dvdyr;



				destBuffer += destBufferLingLength;
				destZBuffer += zBufferLineLength;
				continue;
			}

			idx = ixe - ixs;
			if (idx>0) {
				dudx = (ur - ul) / idx;
				dvdx = (vr - vl) / idx;

				dzdx = (zr - zl) / idx;
			}
			else
			{
				dudx = (ur - ul);
				dvdx = (vr - vl);
				dzdx = (zr - zl);
			}


			invzU = ul;
			invzV = vl;
			invz = zl;


			if (ixs< cliper.left) {
				idx = cliper.left - ixs;
				invzU += idx*dudx;
				invzV += idx*dvdx;
				invz += idx*dzdx;

				ixs = cliper.left;
			}
			if (ixe > cliper.right) {
				ixe = cliper.right;
			}


			for (int xidx = ixs;xidx < ixe;xidx++) {
				if (invz > destZBuffer[xidx]) {
					invzFixp20 = FIXP28_TO_FIXP20(invz);
					iu = (invzU / invzFixp20);
					iv = (invzV / invzFixp20) << textShift;

					destBuffer[xidx] = texture[iu + iv];
					destZBuffer[xidx] = invz;
				}

				invzU += dudx;
				invzV += dvdx;
				invz += dzdx;
			}

			xl += dxdyl;
			xr += dxdyr;

			zl += dzdyl;
			zr += dzdyr;

			ul += dudyl;
			ur += dudyr;

			vl += dvdyl;
			vr += dvdyr;



			destBuffer += destBufferLingLength;
			destZBuffer += zBufferLineLength;
		}
	}
	else
	{
		for (int yidx = iys;yidx <= iye;yidx++) {

			ixs = FIXP19_TO_INT(xl);
			ixe = FIXP19_TO_INT(xr);

			idx = ixe - ixs;
			if (idx>0) {
				dudx = (ur - ul) / idx;
				dvdx = (vr - vl) / idx;

				dzdx = (zr - zl) / idx;
			}
			else
			{
				dudx = (ur - ul);
				dvdx = (vr - vl);
				dzdx = (zr - zl);
			}


			invzU = ul;
			invzV = vl;
			invz = zl;

			for (int xidx = ixs;xidx < ixe;xidx++) {
				if (invz > destZBuffer[xidx]) {
					invzFixp20 = FIXP28_TO_FIXP20(invz);
					iu = (invzU / invzFixp20);
					iv = (invzV / invzFixp20) << textShift;

					destBuffer[xidx] = texture[iu + iv];
					destZBuffer[xidx] = invz;
				}

				invzU += dudx;
				invzV += dvdx;
				invz += dzdx;
			}

			xl += dxdyl;
			xr += dxdyr;

			zl += dzdyl;
			zr += dzdyr;

			ul += dudyl;
			ur += dudyr;

			vl += dvdyl;
			vr += dvdyr;



			destBuffer += destBufferLingLength;
			destZBuffer += zBufferLineLength;
		}
	}


	return 1;
}

int Canvas::drawTriangleGouraudInvZ16(float x0, float y0, float z0, Color c0, float x1, float y1, float z1, Color c1, float x2, float y2, float z2, Color c2)
{
	float  tmp, interX,interZ, fraction;
	Color tmpC, interC;
	//line so return 0 as trivially reject
	if ((FEQUALE1(y0, y1) && FEQUALE1(y1, y2)) || (FEQUALE1(x0, x1) && FEQUALE1(x1, x2))) {
		return 0;
	}

	//select sort for y small -> big
	if (y1 < y0) {
		SWAP(y1, y0, tmp);
		SWAP(x1, x0, tmp);
		SWAP(z1,z0,tmp)
		SWAP(c1, c0, tmpC);
	}
	if (y2 < y0) {
		SWAP(y2, y0, tmp);
		SWAP(x2, x0, tmp);
		SWAP(z2, z0, tmp)
		SWAP(c2, c0, tmpC);
	}
	if (y2 < y1) {
		SWAP(y1, y2, tmp);
		SWAP(x1, x2, tmp);
		SWAP(z1, z2, tmp)
		SWAP(c1, c2, tmpC);
	}

	//trivially reject
	if (y0 > cliper.bottom || y2<cliper.top || (x0<cliper.left&&x1<cliper.left&&x2<cliper.left) || (x0>cliper.right&&x1>cliper.right&&x2>cliper.right)) {
		return 0;
	}

	if (FEQUALE1(y0, y1)) {
		drawTopTriangleGouraudInvZ16(x2, y2,z2, c2, 
									 x0, y0,z0, c0, 
									 x1, y1,z1, c1);
	}
	else if (FEQUALE1(y1, y2))
	{
		drawBottomTriangleGouraudInvZ16(x0, y0,z0,c0, 
										x1, y1,z1,c1, 
										x2, y2, z2,c2);
	}
	else
	{
		fraction = (y1 - y0) / (y2 - y0);
		interX = x0 + (x2 - x0)*fraction;
		interZ = z0 + (z2 - z0)*fraction;
		interC.r = c0.r + (c2.r - c0.r)*fraction + 0.5;
		interC.g = c0.g + (c2.g - c0.g)*fraction + 0.5;
		interC.b = c0.b + (c2.b - c0.b)*fraction + 0.5;
		drawBottomTriangleGouraudInvZ16(x0, y0, z0,c0,x1, y1,z1 ,c1,interX, y1,interZ , interC);
		drawTopTriangleGouraudInvZ16(x2, y2,z2,c2, x1, y1,z1,c1, interX, y1,interZ,interC);
	}
}

int Canvas::drawTopTriangleGouraudInvZ16(float x0, float y0, float z0, Color c0,
										 float x1, float y1, float z1, Color c1, 
										 float x2, float y2, float z2, Color c2)
{	

	//fixp16 xl, xr, dxdyl, dxdyr;

	fixp19 xl, xr, dxdyl, dxdyr;

	/*float xl, xr, dxdyl, dxdyr;*/

	fixp22  invzFixp22,
			rl, rr, gl, gr, bl, br,
			drdyl, drdyr, dgdyl, dgdyr, dbdyl, dbdyr, 
			invzR0, invzG0, invzB0,invzR1,invzG1,invzB1,invzR2, invzG2, invzB2,invzR,invzG,invzB,
			drdx, dgdx, dbdx;

	
	fixp28 invz0, invz1, invz2,invz,zl,zr ,dzdyl, dzdyr,dzdx;

	int idx, idy, 
		iys, iye, ixs, ixe,
		ir,ig,ib,
		r0, g0, b0, r1, g1, b1, r2, g2, b2, 
		tmpi, destBufferLingLength,zBufferLineLength;

	float dx, dy, tmpf;
	uint16 *destBuffer,drawColor;
	uint32 *destZBuffer;

	destBufferLingLength = pitch >> 1;
	zBufferLineLength = zBuffer->pitch >> 2;//32 bit zbuffer wo fix  zBuffer deepth

	r0 = c0.r;
	r1 = c1.r;
	r2 = c2.r;

	g0 = c0.g;
	g1 = c1.g;
	g2 = c2.g;

	b0 = c0.b;
	b1 = c1.b;
	b2 = c2.b;

	if (x2 < x1) {
		SWAP(x2, x1, tmpf);
		SWAP(y2, y1, tmpf);
		SWAP(z2, z1, tmpf);
		SWAP(r2, r1, tmpi);
		SWAP(g2, g1, tmpi);
		SWAP(b2, b1, tmpi);
	}

	invz0 = FIXP28_FLOAT_TO(1.0f / z0);
	invz1 = FIXP28_FLOAT_TO(1.0f / z1);
	invz2 = FIXP28_FLOAT_TO(1.0f / z2);

	invzR0 = FIXP22_FLOAT_TO(r0 / z0);
	invzG0 = FIXP22_FLOAT_TO(g0 / z0);
	invzB0 = FIXP22_FLOAT_TO(b0 / z0);

	invzR1 = FIXP22_FLOAT_TO(r1 / z1);
	invzG1 = FIXP22_FLOAT_TO(g1 / z1);
	invzB1 = FIXP22_FLOAT_TO(b1 / z1);

	invzR2 = FIXP22_FLOAT_TO(r2 / z2);
	invzG2 = FIXP22_FLOAT_TO(g2 / z2);
	invzB2 = FIXP22_FLOAT_TO(b2 / z2);





	//drawTraingleGouraud16 has ignored  y0 == y1 == y2 so just div
	dy = y0 - y1;

	/*dxdyl = FIXP16_FLOAT_TO((x0 - x1) / dy);*/
	dxdyl = FIXP19_FLOAT_TO((x0 - x1) / dy);
	/*dxdyl = (x0 - x1) / dy;*/
	dzdyl = (invz0 - invz1) / dy;
	drdyl = (invzR0 - invzR1) / dy;
	dgdyl = (invzG0 - invzG1) / dy;
	dbdyl = (invzB0 - invzB1) / dy;


	dy = y0 - y2;
	/*dxdyr = FIXP16_FLOAT_TO((x0 - x2) / dy);*/
	dxdyr = FIXP19_FLOAT_TO((x0 - x2) / dy);
	/*dxdyr = (x0 - x2) / dy;*/
	dzdyr = (invz0 - invz2) / dy;
	drdyr = (invzR0 - invzR2) / dy;
	dgdyr = (invzG0 - invzG2) / dy;
	dbdyr = (invzB0 - invzB2) / dy;





	//clip bottom
	if (y0>cliper.bottom)
	{
		iye = cliper.bottom;
	}
	else
	{
		iye = ceil(y0) - 1;
	}

	//clip top
	if (y1 < cliper.top) {
		iys = cliper.top;

	}
	else
	{
		iys = ceil(y1);
	}

	dy = iys - y1;

	/*xl = FIXP16_FLOAT_TO(x1) + dy*dxdyl;*/
	xl = FIXP19_FLOAT_TO(x1) + dy*dxdyl;
	/*xl = x1 + dy*dxdyl;*/
	zl = invz1 + dy*dzdyl;

	rl = invzR1 + dy*drdyl;
	gl = invzG1 + dy*dgdyl;
	bl = invzB1 + dy*dbdyl;
	
	dy = iys - y2;

	/*xr = FIXP16_FLOAT_TO(x2) + dy*dxdyr;*/
	xr = FIXP19_FLOAT_TO(x2) + dy*dxdyr;
	/*xr = x2 + dy*dxdyr;*/
	zr = invz2  + dy*dzdyr;
	rr = invzR2 + dy*drdyr;
	gr = invzG2 + dy*dgdyr;
	br = invzB2 + dy*dbdyr;




	//draw

	//pre round up
	/*xl = FIXP16_ROUNDUP(xl);
	xr = FIXP16_ROUNDUP(xr);*/
	xl = FIXP19_ROUNDUP(xl);
	xr = FIXP19_ROUNDUP(xr);
	/*xl += 0.5;
	xr += 0.5;*/


	destBuffer = (uint16 *)drawBuffer + iys*destBufferLingLength;
	destZBuffer = zBuffer->zBuffer + iys*zBufferLineLength;

	if ((x0 < cliper.left) || (x0 > cliper.right) ||
		(x1 <cliper.left) || (x1 > cliper.right) ||
		(x2 < cliper.left) || (x2 >cliper.right)) {
		//should clip


		for (int yidx = iys;yidx <= iye;yidx++) {
			/*ixs = FIXP16_TO_INT(xl);
			ixe = FIXP16_TO_INT(xr);*/
			ixs = FIXP19_TO_INT(xl);
			ixe = FIXP19_TO_INT(xr);
		/*	ixs = xl;
			ixe = xr;*/

			if (ixe < cliper.left || ixs>cliper.right) {

				xl += dxdyl;
				xr += dxdyr;

				zl += dzdyl;
				zr += dzdyr;

				rl += drdyl;
				rr += drdyr;

				gl += dgdyl;
				gr += dgdyr;

				bl += dbdyl;
				br += dbdyr;

				

				destBuffer += destBufferLingLength;
				destZBuffer += zBufferLineLength;
				continue;
			}

			idx = ixe - ixs;
			if (idx>0) {
				drdx = (rr - rl) / idx;
				dgdx = (gr - gl) / idx;
				dbdx = (br - bl) / idx;
				dzdx = (zr - zl) / idx;
			}
			else
			{
				drdx = (rr - rl);
				dgdx = (gr - gl);
				dbdx = (br - gl);
				dzdx = (zr - zl);
			}


			invzR = rl;
			invzG = gl;
			invzB = bl;
			invz = zl;


			if (ixs< cliper.left) {
				idx = cliper.left - ixs;
				invzR += idx*drdx;
				invzG += idx*dgdx;
				invzB += idx*dbdx;
				invz += idx*dzdx;
				ixs = cliper.left;
			}
			if (ixe > cliper.right) {
				ixe = cliper.right;
			}


			for (int xidx = ixs;xidx < ixe;xidx++) {
				if (invz > destZBuffer[xidx]) {
					invzFixp22 = FIXP28_TO_FIXP22(invz);
					ir = invzR /invzFixp22;
					ig = invzG / invzFixp22;
					ib = invzB / invzFixp22;
					drawColor = R_G_B_TO_RGB565(ir, ig, ib);

					destBuffer[xidx] = drawColor;

					destZBuffer[xidx] = invz;
				}
				invzR += drdx;
				invzG += dgdx;
				invzB += dbdx;
				invz += dzdx;
			}

			xl += dxdyl;
			xr += dxdyr;

		
			zl += dzdyl;
			zr += dzdyr;

			rl += drdyl;
			rr += drdyr;

			gl += dgdyl;
			gr += dgdyr;

			bl += dbdyl;
			br += dbdyr;



			destBuffer += destBufferLingLength;
			destZBuffer += zBufferLineLength;
		}
	}
	else
	{
		//no clip
		for (int yidx = iys;yidx <= iye;yidx++) {
			/*ixs = FIXP16_TO_INT(xl);
			ixe = FIXP16_TO_INT(xr);*/

			ixs = FIXP19_TO_INT(xl);
			ixe = FIXP19_TO_INT(xr);
			/*ixs = xl;
			ixe = xr;*/


			idx = ixe - ixs;
			if (idx > 0) {
				drdx = (rr - rl) / idx;
				dgdx = (gr - gl) / idx;
				dbdx = (br - bl) / idx;
				dzdx = (zr - zl) / idx;
			}
			else
			{
				drdx = (rr - rl);
				dgdx = (gr - gl);
				dbdx = (br - gl);
				dzdx = (zr - zl);
			}


			invzR = rl;
			invzG = gl;
			invzB = bl;
			invz = zl;

			
			for (int xidx = ixs;xidx < ixe;xidx++) {
		
				if (invz > destZBuffer[xidx]) {
					invzFixp22 = FIXP28_TO_FIXP22(invz);
					ir = invzR / invzFixp22;
					ig = invzG / invzFixp22;
					ib = invzB / invzFixp22;
					drawColor = R_G_B_TO_RGB565(ir, ig, ib);

					destBuffer[xidx] = drawColor;
					destZBuffer[xidx] = invz;
				}
				invzR += drdx;
				invzG += dgdx;
				invzB += dbdx;
				invz += dzdx;
			}

			xl += dxdyl;
			xr += dxdyr;


			zl += dzdyl;
			zr += dzdyr;

			rl += drdyl;
			rr += drdyr;

			gl += dgdyl;
			gr += dgdyr;

			bl += dbdyl;
			br += dbdyr;



			destBuffer += destBufferLingLength;
			destZBuffer += zBufferLineLength;
		}
	}


	return 1; 
}

int Canvas::drawBottomTriangleGouraudInvZ16(float x0, float y0, float z0, Color c0, float x1, float y1, float z1, Color c1, float x2, float y2, float z2, Color c2)
{

	/*fixp16 xl, xr, dxdyl, dxdyr;*/
	fixp19 xl, xr, dxdyl, dxdyr;
	//float xl, xr, dxdyl, dxdyr;

	fixp22  invzFixp22,
		rl, rr, gl, gr, bl, br,
		drdyl, drdyr, dgdyl, dgdyr, dbdyl, dbdyr,
		invzR0, invzG0, invzB0, invzR1, invzG1, invzB1, invzR2, invzG2, invzB2, invzR, invzG, invzB,
		drdx, dgdx, dbdx;


	fixp28 invz0, invz1, invz2, invz, zl, zr, dzdyl, dzdyr, dzdx;

	int idx, idy,
		iys, iye, ixs, ixe,
		ir, ig, ib,
		r0, g0, b0, r1, g1, b1, r2, g2, b2,
		tmpi, destBufferLingLength, zBufferLineLength;

	float dx, dy, tmpf;
	uint16 *destBuffer, drawColor;
	uint32 *destZBuffer;

	destBufferLingLength = pitch >> 1;
	zBufferLineLength = zBuffer->pitch >> 2;//32 bit zbuffer wo fix  zBuffer deepth

	r0 = c0.r;
	r1 = c1.r;
	r2 = c2.r;

	g0 = c0.g;
	g1 = c1.g;
	g2 = c2.g;

	b0 = c0.b;
	b1 = c1.b;
	b2 = c2.b;

	if (x2 < x1) {
		SWAP(x2, x1, tmpf);
		SWAP(y2, y1, tmpf);
		SWAP(z2, z1, tmpf);
		SWAP(r2, r1, tmpi);
		SWAP(g2, g1, tmpi);
		SWAP(b2, b1, tmpi);
	}

	invz0 = FIXP28_FLOAT_TO(1.0f / z0);
	invz1 = FIXP28_FLOAT_TO(1.0f / z1);
	invz2 = FIXP28_FLOAT_TO(1.0f / z2);

	invzR0 = FIXP22_FLOAT_TO(r0 / z0);
	invzG0 = FIXP22_FLOAT_TO(g0 / z0);
	invzB0 = FIXP22_FLOAT_TO(b0 / z0);

	invzR1 = FIXP22_FLOAT_TO(r1 / z1);
	invzG1 = FIXP22_FLOAT_TO(g1 / z1);
	invzB1 = FIXP22_FLOAT_TO(b1 / z1);

	invzR2 = FIXP22_FLOAT_TO(r2 / z2);
	invzG2 = FIXP22_FLOAT_TO(g2 / z2);
	invzB2 = FIXP22_FLOAT_TO(b2 / z2);





	//drawTraingleGouraud16 has ignored  y0 == y1 == y2 so just div
	dy = y1 - y0;
	//dxdyl = FIXP16_FLOAT_TO((x1 - x0) / dy);
	dxdyl = FIXP19_FLOAT_TO((x1 - x0) / dy);
	//dxdyl = (x1 - x0) / dy;
	dzdyl = (invz1 - invz0) / dy;
	drdyl = (invzR1 - invzR0) / dy;
	dgdyl = (invzG1 - invzG0) / dy;
	dbdyl = (invzB1 - invzB0) / dy;




	dy = y2 - y0;
	/*dxdyr = FIXP16_FLOAT_TO((x2 - x0) / dy);*/
	dxdyr = FIXP19_FLOAT_TO((x2 - x0) / dy);
	//dxdyr = (x2 - x0) / dy;
	dzdyr = (invz2 - invz0) / dy;
	drdyr = (invzR2 - invzR0) / dy;
	dgdyr = (invzG2 - invzG0) / dy;
	dbdyr = (invzB2 - invzB0) / dy;



	//clip bottom
	if (y1>cliper.bottom)
	{
		iye = cliper.bottom;
	}
	else
	{
		iye = ceil(y1) - 1;
	}

	//clip top
	if (y0 < cliper.top) {
		iys = cliper.top;

	}
	else
	{
		iys = ceil(y0);
	}

	dy = iys - y0;
	/*xl = FIXP16_FLOAT_TO(x0) + dy*dxdyl;*/
	xl = FIXP19_FLOAT_TO(x0) + dy*dxdyl;
	//xl = x0 + dy*dxdyl;
	zl = invz0 + dy*dzdyl;
	rl = invzR0 + dy*drdyl;
	gl = invzG0 + dy*dgdyl;
	bl = invzB0 + dy*dbdyl;


	/*xr = FIXP16_FLOAT_TO(x0) + dy*dxdyr;*/
	xr = FIXP19_FLOAT_TO(x0) + dy*dxdyr;
	//xr = x0 + dy*dxdyr;
	zr = invz0 + dy*dzdyr;
	rr = invzR0 + dy*drdyr;
	gr = invzG0 + dy*dgdyr;
	br = invzB0 + dy*dbdyr;



	//draw

	//pre round up
	/*xl = FIXP16_ROUNDUP(xl);
	xr = FIXP16_ROUNDUP(xr);*/
	xl = FIXP19_ROUNDUP(xl);
	xr = FIXP19_ROUNDUP(xr);
	/*xl += 0.5f;
	xr += 0.5f;*/

	destBuffer = (uint16 *)drawBuffer + iys*destBufferLingLength;
	destZBuffer = zBuffer->zBuffer + iys*zBufferLineLength;

	if ((x0 < cliper.left) || (x0 > cliper.right) ||
		(x1 <cliper.left) || (x1 > cliper.right) ||
		(x2 < cliper.left) || (x2 >cliper.right)) {
		//should clip


		for (int yidx = iys;yidx <= iye;yidx++) {
			/*ixs = FIXP16_TO_INT(xl);
			ixe = FIXP16_TO_INT(xr);*/
			ixs = FIXP19_TO_INT(xl);
			ixe = FIXP19_TO_INT(xr);
			/*ixs = xl;
			ixe = xr;*/

			if (ixe < cliper.left || ixs>cliper.right) {

				xl += dxdyl;
				xr += dxdyr;

				zl += dzdyl;
				zr += dzdyr;

				rl += drdyl;
				rr += drdyr;

				gl += dgdyl;
				gr += dgdyr;

				bl += dbdyl;
				br += dbdyr;



				destBuffer += destBufferLingLength;
				destZBuffer += zBufferLineLength;
				continue;
			}

			idx = ixe - ixs;
			if (idx>0) {
				drdx = (rr - rl) / idx;
				dgdx = (gr - gl) / idx;
				dbdx = (br - bl) / idx;
				dzdx = (zr - zl) / idx;
			}
			else
			{
				drdx = (rr - rl);
				dgdx = (gr - gl);
				dbdx = (br - gl);
				dzdx = (zr - zl);
			}


			invzR = rl;
			invzG = gl;
			invzB = bl;
			invz = zl;


			if (ixs< cliper.left) {
				idx = cliper.left - ixs;
				invzR += idx*drdx;
				invzG += idx*dgdx;
				invzB += idx*dbdx;
				invz += idx*dzdx;

				ixs = cliper.left;
			}
			if (ixe > cliper.right) {
				ixe = cliper.right;
			}


			for (int xidx = ixs;xidx < ixe;xidx++) {
				if (invz > destZBuffer[xidx]) {
					invzFixp22 = FIXP28_TO_FIXP22(invz);
					ir = invzR / invzFixp22;
					ig = invzG / invzFixp22;
					ib = invzB / invzFixp22;
					drawColor = R_G_B_TO_RGB565(ir, ig, ib);

					destBuffer[xidx] = drawColor;

					destZBuffer[xidx] = invz;
				}
				invzR += drdx;
				invzG += dgdx;
				invzB += dbdx;
				invz += dzdx;
			}

			xl += dxdyl;
			xr += dxdyr;

			zl += dzdyl;
			zr += dzdyr;

			rl += drdyl;
			rr += drdyr;

			gl += dgdyl;
			gr += dgdyr;

			bl += dbdyl;
			br += dbdyr;



			destBuffer += destBufferLingLength;
			destZBuffer += zBufferLineLength;
		}
	}
	else
	{
		//no clip
		for (int yidx = iys;yidx <= iye;yidx++) {
		/*	ixs = FIXP16_TO_INT(xl);
			ixe = FIXP16_TO_INT(xr);*/
			ixs = FIXP19_TO_INT(xl);
			ixe = FIXP19_TO_INT(xr);
		/*	ixs = xl;
			ixe = xr;*/

			idx = ixe - ixs;
			if (idx > 0) {
				drdx = (rr - rl) / idx;
				dgdx = (gr - gl) / idx;
				dbdx = (br - bl) / idx;
				dzdx = (zr - zl) / idx;
			}
			else
			{
				drdx = (rr - rl);
				dgdx = (gr - gl);
				dbdx = (br - gl);
				dzdx = (zr - zl);
			}


			invzR = rl;
			invzG = gl;
			invzB = bl;
			invz = zl;


			for (int xidx = ixs;xidx < ixe;xidx++) {
				if (invz > destZBuffer[xidx]) {
					invzFixp22 = FIXP28_TO_FIXP22(invz);
					ir = invzR / invzFixp22;
					ig = invzG / invzFixp22;
					ib = invzB / invzFixp22;
					drawColor = R_G_B_TO_RGB565(ir, ig, ib);

					destBuffer[xidx] = drawColor;
					destZBuffer[xidx] = invz;
				}
				invzR += drdx;
				invzG += dgdx;
				invzB += dbdx;
				invz += dzdx;
			}

			xl += dxdyl;
			xr += dxdyr;

			zl += dzdyl;
			zr += dzdyr;

			rl += drdyl;
			rr += drdyr;

			gl += dgdyl;
			gr += dgdyr;

			bl += dbdyl;
			br += dbdyr;



			destBuffer += destBufferLingLength;
			destZBuffer += zBufferLineLength;
		}
	}


	return 1;
}

int Canvas::drawPolyTriangleSolidInvZ16(PolySelfTriangle4D * poly)
{
	return drawTriangleSolidInvZ16(poly->transVList[0].x, poly->transVList[0].y, poly->transVList[0].z,
								   poly->transVList[1].x, poly->transVList[1].y, poly->transVList[1].z, 
								   poly->transVList[2].x, poly->transVList[2].y, poly->transVList[2].z, 
								   poly->lightColor[0]);
}

int Canvas::drawTriangleSolidInvZ16(float x0, float y0, float z0, float x1, float y1, float z1, float x2, float y2, float z2, Color c)
{
	float  tmp  ,interX,interZ,fraction;
	//line so return 0 as trivially reject
	if ((FEQUALE1(y0, y1) && FEQUALE1(y1, y2)) || (FEQUALE1(x0, x1) && FEQUALE1(x1, x2))) {
		return 0;
	}

	//select sort for y small -> big
	if (y1 < y0) {
		SWAP(y1, y0, tmp);
		SWAP(x1, x0, tmp);
		SWAP(z1, z0, tmp);
	}
	if (y2 < y0) {
		SWAP(y2, y0, tmp);
		SWAP(x2, x0, tmp);
		SWAP(z2, z0, tmp);
	}
	if (y2 < y1) {
		SWAP(y1, y2, tmp);
		SWAP(x1, x2, tmp);
		SWAP(z1, z2, tmp);
	}

	//trivially reject
	if (y0 > cliper.bottom || y2<cliper.top || (x0<cliper.left&&x1<cliper.left&&x2<cliper.left) || (x0>cliper.right&&x1>cliper.right&&x2>cliper.right)) {
		return 0;
	}

	if (FEQUALE1(y0, y1)) {
		drawTopTriangleSolidInvZ16(x2, y2,z2, x0, y0,z0, x1, y1,z1, c);
	}
	else if (FEQUALE1(y1, y2))
	{
		drawBottomTriangleSolidInvZ(x0, y0,z0, x1, y1,z1, x2, y2,z2,c);
	}
	else
	{
		fraction = (y1 - y0) / (y2 - y0);
		interX = x0 +(x2-x0)*fraction ;
		interZ = z0 +(z2 - z0)*fraction;
		drawBottomTriangleSolidInvZ(x0, y0,z0 ,x1, y1,z1, interX, y1,interZ, c);
		drawTopTriangleSolidInvZ16(x2, y2,z2, x1, y1,z1, interX, y1,interZ, c);
	}

	return 1;
}

int Canvas::drawTopTriangleSolidInvZ16(float x0, float y0, float z0, float x1, float y1, float z1, float x2, float y2, float z2, Color c)
{	
	fixp28 zl, zr,invz0,invz1,invz2, dzdyl, dzdyr,dzdx,fxpz;


	/*fixp16 xl, xr,dxdyl, dxdyr;*/

	fixp19 xl, xr, dxdyl, dxdyr;

	int  iys,iye,
		 ixs,ixe,
		 idx;
	uint16 *destBuffer,drawColor;
	uint32 *destZBuffer;
	int destLineLength,zBufferLineLength;
	float dy,tmpf;

	destLineLength = pitch >> 1;
	zBufferLineLength = zBuffer->pitch >> 2;

	drawColor = R_G_B_TO_RGB565(c.r, c.g, c.b);

	if (x2 < x1) {
		SWAP(x2, x1, tmpf);
		SWAP(y2, y1, tmpf);
		SWAP(z2, z1, tmpf);
	}
	
	invz0 = FIXP28_FLOAT_TO(1/ z0);
	invz1 = FIXP28_FLOAT_TO(1 / z1);
	invz2 = FIXP28_FLOAT_TO(1 / z2);

	dy = y0 - y1;
	/*dxdyl = FIXP16_FLOAT_TO((x0 - x1) / dy);*/
	dxdyl = FIXP19_FLOAT_TO((x0 - x1) / dy);
	dzdyl = (invz0 - invz1) / dy;

	dy = y0 - y2;
	dxdyr = FIXP19_FLOAT_TO((x0 - x2) / dy);
	dzdyr = (invz0 - invz2) / dy;

	//clip bottom
	if (y0>cliper.bottom)
	{
		iye = cliper.bottom;
	}
	else
	{
		iye = ceil(y0) - 1;
	}

	//clip top
	if (y1 < cliper.top) {
		iys = cliper.top;

	}
	else
	{
		iys = ceil(y1);
	}

	dy = iys - y1;
	/*xl = FIXP16_FLOAT_TO(x1) + dy*dxdyl;*/
	xl = FIXP19_FLOAT_TO(x1) + dy*dxdyl;
	zl = invz1 + dy*dzdyl;

	dy = iys - y2;
	/*xr = FIXP16_FLOAT_TO(x2) + dy*dxdyr;*/
	xr = FIXP19_FLOAT_TO(x2) + dy*dxdyr;
	zr = invz2 + dy*dzdyr;

	//draw


	//xl = FIXP16_ROUNDUP(xl);
	//xr = FIXP16_ROUNDUP(xr);

	xl = FIXP19_ROUNDUP(xl);
	xr = FIXP19_ROUNDUP(xr);

	destBuffer =(uint16* )drawBuffer + iys*destLineLength;
	destZBuffer =  this->zBuffer->zBuffer+ iys*zBufferLineLength;

	if ((x0 < cliper.left) || (x0 > cliper.right) ||
		(x1 <cliper.left) || (x1 > cliper.right) ||
		(x2 < cliper.left) || (x2 >cliper.right))
	{	

		for (int yidx = iys;yidx <= iye;yidx++) {
			/*ixs = FIXP16_TO_INT(xl);
			ixe = FIXP16_TO_INT(xr);*/

			ixs = FIXP19_TO_INT(xl);
			ixe = FIXP19_TO_INT(xr);

			if (ixs > cliper.right || ixe < cliper.left) {
				xl += dxdyl;
				xr += dxdyr;

				zl += dzdyl;
				zr += dzdyr;

				destBuffer += destLineLength;
				destZBuffer += zBufferLineLength;
				continue;
			}

			idx = ixe - ixs;
			if (idx > 0) {
				dzdx = (zr - zl) / idx;
			}
			else
			{
				dzdx = (zr - zl);
			}

			fxpz = xl;

			//clip
			if (ixs < cliper.left) {
				idx = cliper.left - ixs;

				fxpz += idx*dzdx;

				ixs = cliper.left;
			}

			if (ixe > cliper.right) {
				ixe = cliper.right;
			}

			for (int xidx = ixs;xidx <= ixe;xidx++) {
				if (fxpz > destZBuffer[xidx]) {
					destBuffer[xidx] = drawColor;
				}

			}
			xl += dxdyl;
			xr += dxdyr;

			zl += dzdyl;
			zr += dzdyr;

			destBuffer += destLineLength;
			destZBuffer += zBufferLineLength;
		}
	}else
	{
		for (int yidx = iys;yidx <= iye;yidx++) {
		/*	ixs = FIXP16_TO_INT(xl);
			ixe = FIXP16_TO_INT(xr);*/

			ixs = FIXP19_TO_INT(xl);
			ixe = FIXP19_TO_INT(xr);

			idx = ixe - ixs;
			if (idx > 0) {
				dzdx = (zr - zl) / idx;
			}
			else
			{
				dzdx = (zr - zl);
			}

			fxpz = xl;

			for (int xidx = ixs;xidx <= ixe;xidx++) {
				if (fxpz > destZBuffer[xidx]) {
					destBuffer[xidx] = drawColor;
				}

			}
			xl += dxdyl;
			xr += dxdyr;

			zl += dzdyl;
			zr += dzdyr;

			destBuffer += destLineLength;
			destZBuffer += zBufferLineLength;
		}
	}


	return 1;
}

int Canvas::drawBottomTriangleSolidInvZ(float x0, float y0, float z0, float x1, float y1, float z1, float x2, float y2, float z2, Color c)
{
	fixp28 zl, zr, invz0, invz1, invz2,dzdyl, dzdyr, dzdx, fxpz;

	/*fixp16 xl, xr, dxdyl, dxdyr;*/

	fixp19 xl, xr, dxdyl, dxdyr;

	int  iys, iye,
		 ixs, ixe,
		 idx;
	uint16 *destBuffer, drawColor;
	uint32 *destZBuffer;
	int destLineLength, zBufferLineLength;
	float dy, tmpf;

	destLineLength = pitch >> 1;
	zBufferLineLength =zBuffer->pitch >> 2;

	drawColor = R_G_B_TO_RGB565(c.r, c.g, c.b);

	if (x2 < x1) {
		SWAP(x2, x1, tmpf);
		SWAP(y2, y1, tmpf);
		SWAP(z2, z1, tmpf);
	}

	invz0 = FIXP28_FLOAT_TO(1 / z0);
	invz1 = FIXP28_FLOAT_TO(1 / z1);
	invz2 = FIXP28_FLOAT_TO(1 / z2);

	dy = y1 - y0;
	/*dxdyl = FIXP16_FLOAT_TO((x1 - x0) / dy);*/
	dxdyl = FIXP19_FLOAT_TO((x1 - x0) / dy);
	dzdyl = (invz1 - invz0) / dy;

	dy = y2 - y0;
	/*dxdyr = FIXP16_FLOAT_TO((x2 - x0) / dy);*/
	dxdyr = FIXP19_FLOAT_TO((x2 - x0) / dy);

	dzdyr = (invz2 - invz0) / dy;

	//clip bottom
	if (y1>cliper.bottom)
	{
		iye = cliper.bottom;
	}
	else
	{
		iye = ceil(y0) - 1;
	}

	//clip top
	if (y0 < cliper.top) {
		iys = cliper.top;

	}
	else
	{
		iys = ceil(y1);
	}

	dy = iys - y0;

	/*xl = FIXP16_FLOAT_TO(x0) + dy*dxdyl;*/
	xl = FIXP19_FLOAT_TO(x0) + dy*dxdyl;
	zl = invz0 + dy*dzdyl;

	
	xr = FIXP19_FLOAT_TO(x0) + dy*dxdyr;
	zr = invz0 + dy*dzdyr;

	//draw


	/*xl = FIXP16_ROUNDUP(xl);
	xr = FIXP16_ROUNDUP(xr);*/

	xl = FIXP19_ROUNDUP(xl);
	xr = FIXP19_ROUNDUP(xr);

	destBuffer = (uint16*)drawBuffer + iys*destLineLength;
	destZBuffer = this->zBuffer->zBuffer + iys*zBufferLineLength;

	if ((x0 < cliper.left) || (x0 > cliper.right) ||
		(x1 <cliper.left) || (x1 > cliper.right) ||
		(x2 < cliper.left) || (x2 >cliper.right))
	{

		for (int yidx = iys;yidx <= iye;yidx++) {
			/*ixs = FIXP16_TO_INT(xl);
			ixe = FIXP16_TO_INT(xr);*/

			ixs = FIXP19_TO_INT(xl);
			ixe = FIXP19_TO_INT(xr);

			if (ixs > cliper.right || ixe < cliper.left) {
				xl += dxdyl;
				xr += dxdyr;

				zl += dzdyl;
				zr += dzdyr;

				destBuffer += destLineLength;
				destZBuffer += zBufferLineLength;
				continue;
			}

			idx = ixe - ixs;
			if (idx > 0) {
				dzdx = (zr - zl) / idx;
			}
			else
			{
				dzdx = (zr - zl);
			}

			fxpz = xl;

			//clip
			if (ixs < cliper.left) {
				idx = cliper.left - ixs;

				fxpz += idx*dzdx;

				ixs = cliper.left;
			}

			if (ixe > cliper.right) {
				ixe = cliper.right;
			}

			for (int xidx = ixs;xidx <= ixe;xidx++) {
				if (fxpz > destZBuffer[xidx]) {
					destBuffer[xidx] = drawColor;
					destZBuffer[xidx] = fxpz;
				}
				fxpz += dzdx;

			}
			xl += dxdyl;
			xr += dxdyr;

			zl += dzdyl;
			zr += dzdyr;

			destBuffer += destLineLength;
			destZBuffer += zBufferLineLength;
		}
	}
	else
	{
		for (int yidx = iys;yidx <= iye;yidx++) {
			/*ixs = FIXP16_TO_INT(xl);
			ixe = FIXP16_TO_INT(xr);*/

			ixs = FIXP19_TO_INT(xl);
			ixe = FIXP19_TO_INT(xr);

			idx = ixe - ixs;
			if (idx > 0) {
				dzdx = (zr - zl) / idx;
			}
			else
			{
				dzdx = (zr - zl);
			}

			fxpz = xl;

			for (int xidx = ixs;xidx <= ixe;xidx++) {
				if (fxpz > destZBuffer[xidx]) {
					destBuffer[xidx] = drawColor;
				}

			}
			xl += dxdyl;
			xr += dxdyr;

			zl += dzdyl;
			zr += dzdyr;

			destBuffer += destLineLength;
			destZBuffer += zBufferLineLength;
		}
	}


	return 1;
}

int Canvas::drawPolyTriangleGouraudDirect16(PolySelfTriangle4D * poly)
{	
	return drawTriangleGouraudDirect16(poly->transVList[0].x, poly->transVList[0].y, //
		poly->transVList[1].x, poly->transVList[1].y,//
		poly->transVList[2].x, poly->transVList[2].y,//
		poly->lightColor[0],poly->lightColor[1],poly->lightColor[2]);
}

int Canvas::drawTriangleGouraudDirect16(float x0, float y0, float x1, float y1, float x2, float y2, Color c0, Color c1, Color c2)
{
	fixp16 drdyl, dgdyl, dbdyl, drdyr, dgdyr, dbdyr, rl, rr, gl, gr, bl, br, drdx, dgdx, dbdx, fxpr, fxpg, fxpb,tmpfxp;
	float xl, xr,dy,dxdyl,dxdyr,tmpf;
	int interY,iys,iye,ixs,ixe,idx,ix,r0,g0,b0,r1,b1,g1,r2,g2,b2,rgbi,tmpi;
	int traType;
	int isLeftInter;
	int lineLength;
	uint16 *buffer;
	lineLength = pitch >> 1;

	if ((FEQUALE1(y0, y1) && FEQUALE1(y1, y2)) || (FEQUALE1(x0, x1) && FEQUALE1(x1, x2))) {
		return 0;
	}


	r0 = c0.r;
	g0 = c0.g;
	b0 = c0.b;


	r1 = c1.r;
	g1 = c1.g;
	b1 = c1.b;


	r2 = c2.r;
	g2 = c2.g;
	b2 = c2.b;

	if (y1 < y0) {
		SWAP(y1, y0, tmpf);
		SWAP(x1, x0, tmpf);
		SWAP(r1, r0, tmpi);
		SWAP(g1, g0, tmpi);
		SWAP(b1, b0, tmpi);
	}
	if (y2 < y0) {
		SWAP(y2, y0, tmpf);
		SWAP(x2, x0, tmpf);
		SWAP(r2, r0, tmpi);
		SWAP(g2, g0, tmpi);
		SWAP(b2, b0, tmpi);
	}
	if (y2 < y1) {
		SWAP(y1, y2, tmpf);
		SWAP(x1, x2, tmpf);
		SWAP(r1, r2, tmpi);
		SWAP(g1, g2, tmpi);
		SWAP(b1, b2, tmpi);
	}


	//trivially reject
	if (y0 > cliper.bottom || y2<cliper.top || (x0<cliper.left&&x1<cliper.left&&x2<cliper.left) || (x0>cliper.right&&x1>cliper.right&&x2>cliper.right)) {
		return 0;
	}





	if (FEQUALE1(y0, y1)) {
		traType = CANVAS_TRA_TYPE_TOP;
		if (x0 > x1) {
			SWAP(y1, y0, tmpf);
			SWAP(x1, x0, tmpf);
			SWAP(r1, r0, tmpi);
			SWAP(g1, g0, tmpi);
			SWAP(b1, b0, tmpi);
		}
	}
	else if (FEQUALE1(y1, y2))
	{
		traType = CANVAS_TRA_TYPE_BOTTOM;
		if (x1 > x2) {
			SWAP(y1, y2, tmpf);
			SWAP(x1, x2, tmpf);
			SWAP(r1, r2, tmpi);
			SWAP(g1, g2, tmpi);
			SWAP(b1, b2, tmpi);
		}
	}
	else
	{
		traType = CANVAS_TRA_TYPE_GENERAL;
	}

	if (traType == CANVAS_TRA_TYPE_BOTTOM || traType == CANVAS_TRA_TYPE_TOP) {
		//clip bottom
		if (y2 > cliper.bottom) {
			iye = cliper.bottom;
		}
		else
		{
			iye = ceil(y2) - 1;
		}

		//clip top
		if (traType == CANVAS_TRA_TYPE_BOTTOM) {
			//flat bottom
			dy = y1 - y0;
			dxdyl = (x1 - x0) / dy;
			drdyl = FIXP16_FLOAT_TO((r1 - r0)/dy);
			dgdyl = FIXP16_FLOAT_TO((g1 - g0) / dy);
			dbdyl = FIXP16_FLOAT_TO((b1 - b0) / dy);

			dxdyr = (x2 - x0) / dy;
			drdyr = FIXP16_FLOAT_TO((r2 - r0) / dy);
			dgdyr = FIXP16_FLOAT_TO((g2 - g0) / dy);
			dbdyr = FIXP16_FLOAT_TO((b2 - b0) / dy);

			//clip top
			if (y0 < cliper.top) {
				iys = cliper.top;
				dy = (iys - y0);
				xl = x0 + dy*dxdyl;
				xr=x0+ dy*dxdyr;
				rl = FIXP16_INT_TO(r0) + dy*drdyl;
				rr=  FIXP16_INT_TO(r0) + dy*drdyr;
				gl = FIXP16_INT_TO(g0) + dy*dgdyl;
				gr = FIXP16_INT_TO(g0) + dy*dgdyr;
				bl = FIXP16_INT_TO(b0) + dy*dbdyl;
				br = FIXP16_INT_TO(b0) + dy*dbdyr;
				
			}
			else
			{	
				iys = ceil(y0);
				dy = (iys - y0);
				xl = x0 + dy*dxdyl;
				xr = x0 + dy*dxdyr;
				rl = FIXP16_INT_TO(r0) + dy*drdyl;
				rr = FIXP16_INT_TO(r0) + dy*drdyr;
				gl = FIXP16_INT_TO(g0) + dy*dgdyl;
				gr = FIXP16_INT_TO(g0) + dy*dgdyr;
				bl = FIXP16_INT_TO(b0) + dy*dbdyl;
				br = FIXP16_INT_TO(b0) + dy*dbdyr;

			}

		}
		else {
			//flat top
			dy = y2 - y0;

			dxdyl = (x2 - x0) / dy;
			dxdyr = (x2 - x1) / dy;

			drdyl = FIXP16_FLOAT_TO((r2 - r0) / dy);
			drdyr= FIXP16_FLOAT_TO((r2 - r1) / dy);

			dgdyl = FIXP16_FLOAT_TO((g2 - g0) / dy);
			dgdyr = FIXP16_FLOAT_TO((g2 - g1) / dy);

			dbdyl = FIXP16_FLOAT_TO((b2 - b0) / dy);
			dbdyr = FIXP16_FLOAT_TO((b2 - b1) / dy);
			
			//clip top
			if (y0 < cliper.top) {
				iys = cliper.top;
				dy = iys - y0;

				xl = x0 + dy*dxdyl;
				xr = x1 + dy*dxdyr;

				rl = r0 + dy*drdyl;
				rr = r1 + dy*drdyr;

				gl = g0 + dy*dgdyl;
				gr = g1 + dy*dgdyr;

				bl = b0 + dy*dbdyl;
				br = b1 + dy*dbdyr;
			}
			else
			{
				iys = ceil(y0);
				dy = iys - y0;

				xl = x0 + dy*dxdyl;
				xr = x1 + dy*dxdyr;

				rl = FIXP16_INT_TO(r0) + dy*drdyl;
				rr = FIXP16_INT_TO(r1) + dy*drdyr;

				gl = FIXP16_INT_TO(g0) + dy*dgdyl;
				gr = FIXP16_INT_TO(g1) + dy*dgdyr;

				bl = FIXP16_INT_TO(b0) + dy*dbdyl;
				br = FIXP16_INT_TO(b1) + dy*dbdyr;
			}
		}

	




		//now  draw
		buffer = (uint16 *)drawBuffer + lineLength*iys;
	


		if ((x0 < cliper.left) || (x0 > cliper.right) ||
			(x1 <cliper.left) || (x1 > cliper.right) ||
			(x2 < cliper.left) || (x2 >cliper.right))
		{
			//so should clip x
			for (int yidx = iys;yidx <= iye;yidx++, buffer += lineLength) {
				idx = xr - xl + 0.5;

				if (idx > 0) {
					drdx = (rr - rl) / idx;
					dgdx = (gr - gl) / idx;
					dbdx = (br - bl) / idx;
				}
				else
				{
					drdx = (rr - rl);
					dgdx = (gr - rl);
					dbdx = (br - bl);
				}

				if (xl < cliper.left)
				{
					
					idx = cliper.left - xl;

			
					ixs = cliper.left;
					rl += idx*drdx;
					gl += idx*dgdx;
					bl += idx*dbdx;

				}
				else
				{
					ixs = xl;
				}

				if (xr > cliper.right) {
					ixe = cliper.right;
				}
				else
				{
					ixe = xr;
				}

				// pre roundup
				fxpr = FIXP16_ROUNDUP(rl);
				fxpg = FIXP16_ROUNDUP(gl);
				fxpb = FIXP16_ROUNDUP(bl);

				for (ix = ixs;ix <= ixe;ix++) {

					buffer[ix] = R_G_B_TO_RGB565_FIXP16(fxpr, fxpg, fxpb);

					fxpr += drdx;
					fxpg += dgdx;
					fxpb += dbdx;
				}

				xl += dxdyl;
				xr += dxdyr;
				rl += drdyl;
				rr += drdyr;

				gl += dgdyl;
				gr += dgdyr;
				bl += dbdyl;
				br += dbdyr;

			}
		}
		else
		{

			for (int yidx = iys;yidx <= iye;yidx++, buffer += lineLength) {

				idx = xr - xl + 0.5;

				if (idx > 0) {
					drdx = (rr - rl) / idx;
					dgdx = (gr - gl) / idx;
					dbdx = (br - bl) / idx;
				}
				else
				{
					drdx = (rr - rl);
					dgdx = (gr - rl);
					dbdx = (br - bl);
				}

				ixs = xl;
				ixe = xr;
				//pre roundup
				fxpr = FIXP16_ROUNDUP(rl);
				fxpg = FIXP16_ROUNDUP(gl);
				fxpb = FIXP16_ROUNDUP(bl);

				for (ix = ixs;ix <= ixe;ix++) {

					buffer[ix] = R_G_B_TO_RGB565_FIXP16(fxpr, fxpg, fxpb);

					fxpr += drdx;
					fxpg += dgdx;
					fxpb += dbdx;
				}

				xl += dxdyl;
				xr += dxdyr;
				rl += drdyl;
				rr += drdyr;

				gl += dgdyl;
				gr += dgdyr;
				bl += dbdyl;
				br += dbdyr;
			}

		}
	
	}
	else if(traType==CANVAS_TRA_TYPE_GENERAL){
		//clip bottom
		if (y2 > cliper.bottom) {
			iye = cliper.bottom;
		}
		else
		{
			iye = ceil(y2) - 1;
		}


		//clip top
		if (y1 < cliper.top) {
			iys = cliper.top;

			//left
			//assume y1 left inter
			dy = y2 - y1;
			dxdyl = (x2 - x1) / dy;
			drdyl = FIXP16_FLOAT_TO((r2 - r1) / dy);
			dgdyl = FIXP16_FLOAT_TO((g2 - g1) / dy);
			dbdyl = FIXP16_FLOAT_TO((b2 - b1) / dy);
			


			//right
			dy = y2 - y0;
			dxdyr = (x2 - x0) / dy;
			drdyr = FIXP16_FLOAT_TO((r2 - r0) / dy);
			dgdyr = FIXP16_FLOAT_TO((g2 - r0) / dy);
			dbdyr = FIXP16_FLOAT_TO((b2 - r0) / dy);

			//left start
			dy = iys - y1;
			xl = x1+dy*dxdyl;
			rl = FIXP16_INT_TO(r1) + dy*drdyl;
			gl = FIXP16_INT_TO(g1) + dy*dgdyl;
			bl = FIXP16_INT_TO(b1) + dy*dbdyl;

			//right start
			dy = iys - y0;
			xr = x0+ dy*dxdyr;
			rr = FIXP16_INT_TO(r0) + dy*drdyr;
			gr = FIXP16_INT_TO(g0) + dy*dgdyr;
			br = FIXP16_INT_TO(b0) + dy*dbdyr;

			if (dxdyr > dxdyl) {
				//assume y1=left but  y1=right so reverse
				SWAP(dxdyl, dxdyr, tmpf);
				SWAP(drdyl, drdyr, tmpfxp);
				SWAP(dgdyl, dgdyr, tmpfxp);
				SWAP(dbdyl, dbdyr, tmpfxp);
				SWAP(xl, xr, tmpf);
				SWAP(rl, rr, tmpfxp);
				SWAP(gl, gr, tmpfxp);
				SWAP(bl, br, tmpfxp);
				isLeftInter = 0;
			}
			else
			{
				isLeftInter = 1;
			}
		}
		else if(y0<cliper.top)
		{	
			iys = cliper.top;
			
			//left
			//assume x1==left,
			dy = y1 - y0;
			dxdyl = (x1 - x0) / dy;
			drdyl = FIXP16_FLOAT_TO((r1 - r0) / dy);
			dgdyl = FIXP16_FLOAT_TO((g1 - g0) / dy);
			dbdyl = FIXP16_FLOAT_TO((b1 - b0) / dy);

			//right
			dy = y2 - y0;
			dxdyr = (x2 - x0) / dy;
			drdyr = FIXP16_FLOAT_TO((r2 - r0) / dy);
			dgdyr = FIXP16_FLOAT_TO((g2 - g0) / dy);
			dbdyr = FIXP16_FLOAT_TO((b2 - b0) / dy);

			dy = iys - y0;
			//left start
			xl = x0 + dy*dxdyl;
			rl = FIXP16_INT_TO(r0) + dy*drdyl;
			gl = FIXP16_INT_TO(g0) + dy*dgdyl;
			bl = FIXP16_INT_TO(b0) + dy*dbdyl;

			//right start
			xr = x0 + dy*dxdyr;
			rr = FIXP16_INT_TO(r0) + dy*drdyr;
			gr = FIXP16_INT_TO(g0) + dy*dgdyr;
			br = FIXP16_INT_TO(b0) + dy*dbdyr;

			if (dxdyl > dxdyr) {
				//assume x1==left  but xi==right ,so reverse
				SWAP(dxdyl, dxdyr, tmpf);
				SWAP(drdyl, drdyr, tmpfxp);
				SWAP(dgdyl, dgdyr, tmpfxp);
				SWAP(dbdyl, dbdyr, tmpfxp);
				SWAP(xl, xr, tmpf);
				SWAP(rl, rr, tmpfxp);
				SWAP(gl, gr, tmpfxp);
				SWAP(bl, br, tmpfxp);
				isLeftInter = 0;
			}
			else
			{
				isLeftInter = 1;
			}
		}
		else
		{
			iys =ceil(y0);

			//left
			//assume x1==left,
			dy = y1 - y0;
			dxdyl = (x1 - x0) / dy;
		

			drdyl = FIXP16_FLOAT_TO((r1 - r0) / dy);
			dgdyl = FIXP16_FLOAT_TO((g1 - g0) / dy);
			dbdyl = FIXP16_FLOAT_TO((b1 - b0) / dy);

			//right
			dy = y2 - y0;
			dxdyr = (x2 - x0) / dy;

		
			

			drdyr = FIXP16_FLOAT_TO((r2 - r0) / dy);
			dgdyr = FIXP16_FLOAT_TO((g2 - g0) / dy);
			dbdyr = FIXP16_FLOAT_TO((b2 - b0) / dy);

			dy = iys - y0;
			//left start
			xl = x0 + dy*dxdyl;
		

			rl = FIXP16_INT_TO(r0) + dy*drdyl;
			gl = FIXP16_INT_TO(g0) + dy*dgdyl;
			bl = FIXP16_INT_TO(b0) + dy*dbdyl;

			//right start
			xr = x0 + dy*dxdyr;



			rr = FIXP16_INT_TO(r0) + dy*drdyr;
			gr = FIXP16_INT_TO(g0) + dy*dgdyr;
			br = FIXP16_INT_TO(b0) + dy*dbdyr;

			if (dxdyl > dxdyr) {
				//assume x1==left  but xi==right ,so reverse
				SWAP(dxdyl, dxdyr, tmpf);
				SWAP(drdyl, drdyr, tmpfxp);
				SWAP(dgdyl, dgdyr, tmpfxp);
				SWAP(dbdyl, dbdyr, tmpfxp);
				SWAP(xl, xr, tmpf);
				SWAP(rl, rr, tmpfxp);
				SWAP(gl, gr, tmpfxp);
				SWAP(bl, br, tmpfxp);
				isLeftInter = 0;
			}
			else
			{
				isLeftInter = 1;
			}
		}


		//now  draw
		buffer = (uint16 *)drawBuffer + lineLength*iys;

		interY = ceil(y1)-1;

		if ((x0 < cliper.left) || (x0 > cliper.right) ||
			(x1 <cliper.left) || (x1 > cliper.right) ||
			(x2 < cliper.left) || (x2 >cliper.right))
		{
			//so should clip 	
			for (int yidx = iys;yidx <= iye;yidx++, buffer += lineLength) {
				idx = xr - xl + 0.5;

				if (idx > 0) {
					drdx = (rr - rl) / idx;
					dgdx = (gr - gl) / idx;
					dbdx = (br - bl) / idx;
				}
				else
				{
					drdx = (rr - rl);
					dgdx = (gr - rl);
					dbdx = (br - bl);
				}
				if (xl < cliper.left)
				{
					idx = cliper.left - xl;

					ixs = cliper.left;
					rl += idx*drdx;
					gl += idx*dgdx;
					bl += idx*dbdx;

				}
				else
				{
					ixs = xl;
				}

			
				if (xr > cliper.right) {
					ixe = cliper.right;
				}
				else
				{
					ixe = xr;
				}

				// pre roundup
				fxpr = FIXP16_ROUNDUP(rl);
				fxpg = FIXP16_ROUNDUP(gl);
				fxpb = FIXP16_ROUNDUP(bl);

				for (ix = ixs;ix <= ixe;ix++) {

					buffer[ix] = R_G_B_TO_RGB565_FIXP16(fxpr, fxpg, fxpb);

					fxpr += drdx;
					fxpg += dgdx;
					fxpb += dbdx;
				}

				xl += dxdyl;
				xr += dxdyr;

				rl += drdyl;
				rr += drdyr;

				gl += dgdyl;
				gr += dgdyr;

				bl += dbdyl;
				br += dbdyr;

				if (yidx == interY) {

					if (isLeftInter) {
						//left inter
						dy = y2 - y1;
						dxdyl = (x2 - x1) / dy;
						
					    drdyl = FIXP16_FLOAT_TO((r2 - r1) / dy);
						dgdyl = FIXP16_FLOAT_TO((g2 - g1) / dy);
						dbdyl = FIXP16_FLOAT_TO((b2 - b1) / dy);

						dy = interY+1-y1;//ceil(y1)-y1;
						xl = x1 + dy*dxdyl;
						rl = FIXP16_INT_TO(r1) + dy*drdyl;
						gl = FIXP16_INT_TO(g1) + dy*dgdyl;
						bl = FIXP16_INT_TO(b1) + dy*dbdyl;

					}
					else
					{	
						//right inter
						dy = y2 - y1;
						dxdyr = (x2 - x1) / dy;

						drdyr = FIXP16_FLOAT_TO((r2 - r1) / dy);
						dgdyr = FIXP16_FLOAT_TO((g2 - g1) / dy);
						dbdyr = FIXP16_FLOAT_TO((b2 - b1) / dy);

						dy = interY + 1 - y1;//ceil(y1)-y1;
						xr = x1 + dy*dxdyl;
						rr = FIXP16_INT_TO(r1) + dy*drdyl;
						gr = FIXP16_INT_TO(g1) + dy*dgdyl;
						br = FIXP16_INT_TO(b1) + dy*dbdyl;
					}
					
				}

			}
		}
		else
		{

		

			for (int yidx = iys;yidx <= iye;yidx++, buffer += lineLength) {
		

				idx = xr - xl + 0.5;

				if (idx > 0) {
					drdx = (rr - rl) / idx;
					dgdx = (gr - gl) / idx;
					dbdx = (br - bl) / idx;
				}
				else
				{
					drdx = (rr - rl);
					dgdx = (gr - rl);
					dbdx = (br - bl);
				}

				
				ixs = xl;
				ixe = xr;

		
				

				// pre roundup
				fxpr = FIXP16_ROUNDUP(rl);
				fxpg = FIXP16_ROUNDUP(gl);
				fxpb = FIXP16_ROUNDUP(bl);

				

				for (ix = ixs;ix <= ixe;ix++) {
					

					buffer[ix] = R_G_B_TO_RGB565_FIXP16(fxpr, fxpg, fxpb);

					fxpr += drdx;
					fxpg += dgdx;
					fxpb += dbdx;
				}

				

				xl += dxdyl;
				xr += dxdyr;

				rl += drdyl;
				rr += drdyr;

				gl += dgdyl;
				gr += dgdyr;

				bl += dbdyl;
				br += dbdyr;

				if (yidx == interY) {

					if (isLeftInter) {
						//left inter
						dy = y2 - y1;
						dxdyl = (x2 - x1) / dy;

						drdyl = FIXP16_FLOAT_TO((r2 - r1) / dy);
						dgdyl = FIXP16_FLOAT_TO((g2 - g1) / dy);
						dbdyl = FIXP16_FLOAT_TO((b2 - b1) / dy);

						dy = interY + 1 - y1;//ceil(y1)-y1;
						xl = x1 + dy*dxdyl;
						rl = FIXP16_INT_TO(r1) + dy*drdyl;
						gl = FIXP16_INT_TO(g1) + dy*dgdyl;
						bl = FIXP16_INT_TO(b1) + dy*dbdyl;

					}
					else
					{
						//right inter
						dy = y2 - y1;
						dxdyr = (x2 - x1) / dy;

						drdyr = FIXP16_FLOAT_TO((r2 - r1) / dy);
						dgdyr = FIXP16_FLOAT_TO((g2 - g1) / dy);
						dbdyr = FIXP16_FLOAT_TO((b2 - b1) / dy);

						dy = interY + 1 - y1;//ceil(y1)-y1;
						xr = x1 + dy*dxdyl;
						rr = FIXP16_INT_TO(r1) + dy*drdyl;
						gr = FIXP16_INT_TO(g1) + dy*dgdyl;
						br = FIXP16_INT_TO(b1) + dy*dbdyl;
					}

				}

			}
		}
	}


	return 1;
}

int Canvas::drawPolyTriangleTextureGouraud16(PolySelfTriangle4D * poly)
{	
	Bitmap *texture = poly->texure;
	return drawTriangleTextureGouraud16(poly->transVList[0].x, poly->transVList[0].y, poly->transVList[0].u, poly->transVList[0].v,poly->lightColor[0],
		poly->transVList[1].x, poly->transVList[1].y, poly->transVList[1].u, poly->transVList[1].v, poly->lightColor[1],
		poly->transVList[2].x, poly->transVList[2].y, poly->transVList[2].u, poly->transVList[2].v, poly->lightColor[2],
		(uint16 *)texture->getBuffer(), texture->getWidth(), texture->getHeight(), texture->getPitch());
}

int Canvas::drawTriangleTextureGouraud16(float x0, float y0, float u0, float v0, Color c0, float x1, float y1, float u1, float v1, Color c1, float x2, float y2, float u2, float v2, Color c2, uint16 * textBuffer, int width, int height, int pitch)
{
	float  tmp, interX, interU, interV, fraction;
	Color tmpC, interC;
	//line so return 0 as trivially reject
	if ((FEQUALE1(y0, y1) && FEQUALE1(y1, y2)) || (FEQUALE1(x0, x1) && FEQUALE1(x1, x2))) {
		return 0;
	}

	//select sort for y small -> big
	if (y1 < y0) {
		SWAP(y1, y0, tmp);
		SWAP(x1, x0, tmp);
		SWAP(u1, u0,tmp);
		SWAP(v1, v0, tmp);

		SWAP(c1, c0, tmpC);
	}
	if (y2 < y0) {
		SWAP(y2, y0, tmp);
		SWAP(x2, x0, tmp);
		SWAP(u2, u0, tmp);
		SWAP(v2, v0, tmp);

		SWAP(c2, c0, tmpC);
	}
	if (y2 < y1) {
		SWAP(y1, y2, tmp);
		SWAP(x1, x2, tmp);
		SWAP(u1, u2, tmp);
		SWAP(v1, v2, tmp);

		SWAP(c1, c2, tmpC);
	}

	//trivially reject
	if (y0 > cliper.bottom || y2<cliper.top || (x0<cliper.left&&x1<cliper.left&&x2<cliper.left) || (x0>cliper.right&&x1>cliper.right&&x2>cliper.right)) {
		return 0;
	}

	if (FEQUALE1(y0, y1)) {
		drawTopTriangleTextureGouraud16(x2, y2,u2,v2,c2, x0, y0,u0,v0,c0, x1, y1,u1,v1,c1,textBuffer,width,height,pitch);
	}
	else if (FEQUALE1(y1, y2))
	{
		drawBottomTriangleTextureGouraud16(x0, y0,u0,v0,c0, x1, y1,u1,v1,c1, x2, y2, u2, v2, c2, textBuffer, width, height, pitch);
	}
	else
	{
		fraction = (y1 - y0) / (y2 - y0);
		interX = x0 + (x2 - x0)*fraction;
		interU = u0 + (u2 - u0)*fraction;
		interV = v0 + (v2 - v0)*fraction;

		interC.r = c0.r + (c2.r - c0.r)*fraction + 0.5;
		interC.g = c0.g + (c2.g - c0.g)*fraction + 0.5;
		interC.b = c0.b + (c2.b - c0.b)*fraction + 0.5;
		drawBottomTriangleTextureGouraud16(x0, y0,u0,v0,c0,
										   x1, y1,u1,v1,c1,
										   interX, y1,interU,interV,interC,
										    textBuffer, width, height, pitch);
		drawTopTriangleTextureGouraud16(x2, y2,u2,v2,c2,
										x1, y1,u1,v1,c1,
										interX, y1, interU,interV,interC,
										textBuffer, width, height, pitch);
	}

	return 1;
}

int Canvas::drawTopTriangleTextureGouraud16(float x0, float y0, float u0, float v0, Color c0, float x1, float y1, float u1, float v1, Color c1, float x2, float y2, float u2, float v2, Color c2, uint16 * textBuffer, int width, int height, int pitch)
{
	fixp16 xl, xr,
		ul, ur,
		vl, vr,
		dxdyl, dxdyr,
		dudyl, dudyr,
		dvdyl, dvdyr,
		dudx, dvdx,
		fxpu, fxpv,
		rl, rr,
		gl, gr,
		bl, br,
		drdyl, drdyr,
		dgdyl, dgdyr,
		dbdyl, dbdyr,
		drdx,  dgdx ,dbdx,
		fxpr,fxpg,fxpb;
	int iu, iv,
		idx, idy,
		iys, iye,
		ixs, ixe,
		r0,b0,g0,
		r1,g1,b1,
		r2,g2,b2,
		ir,ig,ib,
		rDraw,gDraw,bDraw,
		tmpi,
		destLineLength;

	float dx, dy, tmpf;
	uint16 *srcBuffer, *destBuffer,drawColor;
	int  lineShift, textShift;



	r0 = c0.r;
	g0 = c0.g;
	b0 = c0.b;

	r1 = c1.r;
	g1 = c1.g;
	b1 = c1.b;

	r2 = c2.r;
	g2 = c2.g;
	b2 = c2.b;


	destLineLength = this->pitch >> 1;
	textShift = Log_Base2(width);
	lineShift = Log_Base2(width);

	if (x2 < x1) {
		SWAP(x2, x1, tmpf);
		SWAP(y2, y1, tmpf);
		SWAP(u2, u1, tmpf);
		SWAP(v2, v1, tmpf);
		SWAP(r2, r1, tmpi);
		SWAP(g2, g1, tmpi);
		SWAP(b2, b1, tmpi);
	}


	//pre decor has ignored  y0 == y1 == y2 so just div
	dy = y0 - y1;
	dxdyl = FIXP16_FLOAT_TO((x0 - x1) / dy);
	dxdyr = FIXP16_FLOAT_TO((x0 - x2) / dy);

	dudyl = FIXP16_FLOAT_TO((u0 - u1) / dy);
	dudyr = FIXP16_FLOAT_TO((u0 - u2) / dy);

	dvdyl = FIXP16_FLOAT_TO((v0 - v1) / dy);
	dvdyr = FIXP16_FLOAT_TO((v0 - v2) / dy);

	drdyl = FIXP16_FLOAT_TO((r0 - r1) / dy);
	drdyr = FIXP16_FLOAT_TO((r0 - r2) / dy);

	dgdyl = FIXP16_FLOAT_TO((g0 - g1) / dy);
	dgdyr = FIXP16_FLOAT_TO((g0 - g2) / dy);

	dbdyl = FIXP16_FLOAT_TO((b0 - b1) / dy);
	dbdyr = FIXP16_FLOAT_TO((b0 - b2) / dy);


	//clip bottom
	if (y0>cliper.bottom)
	{
		iye = cliper.bottom;
	}
	else
	{
		iye = ceil(y0) - 1;
	}

	//clip top
	if (y1 < cliper.top) {
		iys = cliper.top;

	}
	else
	{
		iys = ceil(y1);
	}

	dy = iys - y1;

	xl = FIXP16_FLOAT_TO(x1) + dy*dxdyl;
	ul = FIXP16_FLOAT_TO(u1) + dy*dudyl;
	vl = FIXP16_FLOAT_TO(v1) + dy*dvdyl;

	rl = FIXP16_INT_TO(r1) + dy*drdyl;
	gl = FIXP16_INT_TO(g1) + dy*dgdyl;
	bl = FIXP16_INT_TO(b1) + dy*dbdyl;

	xr = FIXP16_FLOAT_TO(x2) + dy*dxdyr;
	ur = FIXP16_FLOAT_TO(u2) + dy*dudyr;
	vr = FIXP16_FLOAT_TO(v2) + dy*dvdyr;

	rr = FIXP16_INT_TO(r2) + dy*drdyr;
	gr = FIXP16_INT_TO(g2) + dy*dgdyr;
	br = FIXP16_INT_TO(b2) + dy*dbdyr;


	//draw

	//pre round up  uv must later ,because need mul
	xl = FIXP16_ROUNDUP(xl);
	xr = FIXP16_ROUNDUP(xr);
	rl = FIXP16_ROUNDUP(rl);
	rr = FIXP16_ROUNDUP(rr);
	gl = FIXP16_ROUNDUP(gl);
	gr = FIXP16_ROUNDUP(gr);
	bl = FIXP16_ROUNDUP(bl);
	br = FIXP16_ROUNDUP(br);

	//ready for buffer
	destBuffer = (uint16 *)drawBuffer + iys*destLineLength;
	srcBuffer = textBuffer;
	if ((x0 < cliper.left) || (x0 > cliper.right) ||
		(x1 <cliper.left) || (x1 > cliper.right) ||
		(x2 < cliper.left) || (x2 >cliper.right)) {
		//should clip


		for (int yidx = iys;yidx <= iye;yidx++) {
			ixs = FIXP16_TO_INT(xl);
			ixe = FIXP16_TO_INT(xr);


			if (ixe < cliper.left || ixs>cliper.right) {


				xl += dxdyl;
				xr += dxdyr;

				ul += dudyl;
				ur += dudyr;

				vl += dvdyl;
				vr += dvdyr;


				rl += drdyl;
				rr += drdyr;

				gl += dgdyl;
				gr += dgdyr;

				bl += dbdyl;
				br += dbdyr;

				destBuffer += destLineLength;

				continue;
			}

			fxpu = ul;
			fxpv = vl;
			fxpr = rl;
			fxpg = gl;
			fxpb = bl;


			idx = ixe - ixs;
			if (idx>0) {
				dudx = (ur - ul) / idx;
				dvdx = (vr - vl) / idx;
				drdx = (rr - rl) / idx;
				dgdx = (gr - gl) / idx;
				dbdx = (br - bl) / idx;
			}
			else
			{
				dudx = (ur - ul);
				dvdx = (vr - vl);
				drdx = (rr - rl);
				dgdx = (gr - gl);
				dbdx = (br - bl);
			}



			if (ixs< cliper.left) {
				idx = cliper.left - ixs;
				fxpu += idx*dudx;
				fxpv += idx*dvdx;
				fxpr += idx*drdx;
				fxpg += idx*dgdx;
				fxpb += idx*dbdx;
				ixs = cliper.left;
			}
			if (ixe > cliper.right) {
				ixe = cliper.right;
			}


			for (int xidx = ixs;xidx <= ixe;xidx++) {

				iu = FIXP16_TO_INT(FIXP16_ROUNDUP(fxpu << textShift));
				iv = FIXP16_TO_INT(FIXP16_ROUNDUP(fxpv << textShift)) << lineShift;
				ir = FIXP16_TO_INT(fxpr);
				ig = FIXP16_TO_INT(fxpg);
				ib = FIXP16_TO_INT(fxpb);

				//modulate
				drawColor = srcBuffer[iu + iv];

				RGB565_TO_R5_G6_B5(drawColor, rDraw, gDraw, bDraw);
				rDraw = (rDraw*ir) >> 8;
				gDraw = (gDraw*ig) >> 8;
				bDraw = (bDraw*ib) >> 8;
				drawColor = R5_G6_B5_TO_RGB565(rDraw, gDraw, bDraw);

				destBuffer[xidx] = drawColor;

				fxpu += dudx;
				fxpv += dvdx;
				fxpr += drdx;
				fxpg += dgdx;
				fxpb += dbdx;
			}

			xl += dxdyl;
			xr += dxdyr;

			ul += dudyl;
			ur += dudyr;

			vl += dvdyl;
			vr += dvdyr;
			
			rl += drdyl;
			rr += drdyr;

			gl += dgdyl;
			gr += dgdyr;

			bl += dbdyl;
			br += dbdyr;


			destBuffer += destLineLength;
		}
	}
	else
	{
		//no clip
		for (int yidx = iys;yidx <= iye;yidx++) {
			ixs = FIXP16_TO_INT(xl);
			ixe = FIXP16_TO_INT(xr);

			fxpu = ul;
			fxpv = vl;

			fxpr = rl;
			fxpg = gl;
			fxpb = bl;

			idx = ixe - ixs;
			if (idx>0) {
				dudx = (ur - ul) / idx;
				dvdx = (vr - vl) / idx;

				drdx = (rr - rl) / idx;
				dgdx = (gr - gl) / idx;
				dbdx = (br - bl) / idx;
			}
			else
			{
				dudx = (ur - ul);
				dvdx = (vr - vl);

				drdx = (rr - rl);
				dgdx = (gr - gl);
				dbdx = (br - bl);
			}

			for (int xidx = ixs;xidx <= ixe;xidx++) {

				iu = FIXP16_TO_INT(FIXP16_ROUNDUP(fxpu << textShift));
				iv = FIXP16_TO_INT(FIXP16_ROUNDUP(fxpv << textShift)) << lineShift;
				ir = FIXP16_TO_INT(fxpr);
				ig = FIXP16_TO_INT(fxpg);
				ib = FIXP16_TO_INT(fxpb);

				//modulate
				drawColor = srcBuffer[iu + iv];

				RGB565_TO_R5_G6_B5(drawColor, rDraw, gDraw, bDraw);
				rDraw = (rDraw*ir) >> 8;
				gDraw = (gDraw*ig) >> 8;
				bDraw = (bDraw*ib) >> 8;
				drawColor = R5_G6_B5_TO_RGB565(rDraw, gDraw, bDraw);

				destBuffer[xidx] = drawColor;

				fxpu += dudx;
				fxpv += dvdx;
				fxpr += drdx;
				fxpg += dgdx;
				fxpb += dbdx;

			}
			xl += dxdyl;
			xr += dxdyr;

			ul += dudyl;
			ur += dudyr;

			vl += dvdyl;
			vr += dvdyr;


			rl += drdyl;
			rr += drdyr;

			gl += dgdyl;
			gr += dgdyr;

			bl += dbdyl;
			br += dbdyr;

			destBuffer += destLineLength;
		}

	}


	return 1;
}

int Canvas::drawBottomTriangleTextureGouraud16(float x0, float y0, float u0, float v0, Color c0, float x1, float y1, float u1, float v1, Color c1, float x2, float y2, float u2, float v2, Color c2, uint16 * textBuffer, int width, int height, int pitch)
{
	fixp16 xl, xr,
		ul, ur,
		vl, vr,
		dxdyl, dxdyr,
		dudyl, dudyr,
		dvdyl, dvdyr,
		dudx, dvdx,
		fxpu, fxpv,
		rl, rr,
		gl, gr,
		bl, br,
		drdyl, drdyr,
		dgdyl, dgdyr,
		dbdyl, dbdyr,
		drdx, dgdx, dbdx,
		fxpr, fxpg, fxpb;
	int iu, iv,
		idx, idy,
		iys, iye,
		ixs, ixe,
		r0, b0, g0,
		r1, g1, b1,
		r2, g2, b2,
		ir, ig, ib,
		rDraw, gDraw, bDraw,
		tmpi,
		destLineLength;

	float dx, dy, tmpf;
	uint16 *srcBuffer, *destBuffer,drawColor;
	int  lineShift, textShift;


	r0 = c0.r;
	g0 = c0.g;
	b0 = c0.b;

	r1 = c1.r;
	g1 = c1.g;
	b1 = c1.b;

	r2 = c2.r;
	g2 = c2.g;
	b2 = c2.b;


	destLineLength = this->pitch >> 1;
	textShift = Log_Base2(width);
	lineShift = Log_Base2(width);

	if (x2 < x1) {
		SWAP(x2, x1, tmpf);
		SWAP(y2, y1, tmpf);
		SWAP(u2, u1, tmpf);
		SWAP(v2, v1, tmpf);

		SWAP(r2, r1, tmpi);
		SWAP(g2, g1, tmpi);
		SWAP(b2, b1, tmpi);
	}


	//pre decor has ignored  y0 == y1 == y2 so just div
	dy = y1 - y0;
	dxdyl = FIXP16_FLOAT_TO((x1 - x0) / dy);
	dxdyr = FIXP16_FLOAT_TO((x2 - x0) / dy);

	dudyl = FIXP16_FLOAT_TO((u1 - u0) / dy);
	dudyr = FIXP16_FLOAT_TO((u2 - v0) / dy);

	dvdyl = FIXP16_FLOAT_TO((v1 - v0) / dy);
	dvdyr = FIXP16_FLOAT_TO((v2 - v0) / dy);


	drdyl = FIXP16_FLOAT_TO((r1 - r0) / dy);
	drdyr = FIXP16_FLOAT_TO((r2 - r0) / dy);

	dgdyl = FIXP16_FLOAT_TO((g1 - g0) / dy);
	dgdyr = FIXP16_FLOAT_TO((g2 - g0) / dy);

	dbdyl = FIXP16_FLOAT_TO((b1 - b0) / dy);
	dbdyr = FIXP16_FLOAT_TO((b2 - b0) / dy);



	//clip bottom
	if (y1>cliper.bottom)
	{
		iye = cliper.bottom;
	}
	else
	{
		iye = ceil(y1) - 1;
	}

	//clip top
	if (y0 < cliper.top) {
		iys = cliper.top;

	}
	else
	{
		iys = ceil(y0);
	}

	dy = iys - y0;

	xl = FIXP16_FLOAT_TO(x0) + dy*dxdyl;
	ul = FIXP16_FLOAT_TO(u0) + dy*dudyl;
	vl = FIXP16_FLOAT_TO(v0) + dy*dvdyl;


	rl = FIXP16_INT_TO(r0) + dy*drdyl;
	gl = FIXP16_INT_TO(g0) + dy*dgdyl;
	bl = FIXP16_INT_TO(b0) + dy*dbdyl;

	xr = FIXP16_FLOAT_TO(x0) + dy*dxdyr;
	ur = FIXP16_FLOAT_TO(u0) + dy*dudyr;
	vr = FIXP16_FLOAT_TO(v0) + dy*dvdyr;


	rr = FIXP16_INT_TO(r0) + dy*drdyr;
	gr = FIXP16_INT_TO(g0) + dy*dgdyr;
	br = FIXP16_INT_TO(b0) + dy*dbdyr;

	//draw

	//pre round up ,uv should roudup later because uv need mul ;
	xl = FIXP16_ROUNDUP(xl);
	xr = FIXP16_ROUNDUP(xr);
	rl = FIXP16_ROUNDUP(rl);
	rr = FIXP16_ROUNDUP(rr);
	gl = FIXP16_ROUNDUP(gl);
	gr = FIXP16_ROUNDUP(gr);
	bl = FIXP16_ROUNDUP(bl);
	br = FIXP16_ROUNDUP(br);

	//ready for buffer
	destBuffer = (uint16 *)drawBuffer + iys*destLineLength;
	srcBuffer = textBuffer;
	if ((x0 < cliper.left) || (x0 > cliper.right) ||
		(x1 <cliper.left) || (x1 > cliper.right) ||
		(x2 < cliper.left) || (x2 >cliper.right)) {
		//should clip


		for (int yidx = iys;yidx <= iye;yidx++) {
			ixs = FIXP16_TO_INT(xl);
			ixe = FIXP16_TO_INT(xr);

			if (ixe < cliper.left || ixs>cliper.right) {


				xl += dxdyl;
				xr += dxdyr;

				ul += dudyl;
				ur += dudyr;

				vl += dvdyl;
				vr += dvdyr;


				rl += drdyl;
				rr += drdyr;

				gl += dgdyl;
				gr += dgdyr;

				bl += dbdyl;
				br += dbdyr;

				destBuffer += destLineLength;

				continue;
			}

			fxpu = ul;
			fxpv = vl;
			fxpr = rl;
			fxpg = gl;
			fxpb = bl;


			idx = ixe - ixs;
			if (idx>0) {
				dudx = (ur - ul) / idx;
				dvdx = (vr - vl) / idx;
				drdx = (rr - rl) / idx;
				dgdx = (gr - gl) / idx;
				dbdx = (br - bl) / idx;
			}
			else
			{
				dudx = (ur - ul);
				dvdx = (vr - vl);
				drdx = (rr - rl);
				dgdx = (gr - gl);
				dbdx = (br - bl);
			}



			if (ixs< cliper.left) {
				idx = cliper.left - ixs;

				fxpu += idx*dudx;
				fxpv += idx*dvdx;
				fxpr += idx*drdx;
				fxpg += idx*dgdx;
				fxpb += idx*dbdx;

				ixs = cliper.left;
			}
			if (ixe > cliper.right) {
				ixe = cliper.right;
			}


			for (int xidx = ixs;xidx <= ixe;xidx++) {

				iu = FIXP16_TO_INT(FIXP16_ROUNDUP(fxpu << textShift));
				iv = FIXP16_TO_INT(FIXP16_ROUNDUP(fxpv << textShift)) << lineShift;
				ir = FIXP16_TO_INT(fxpr);
				ig = FIXP16_TO_INT(fxpg);
				ib = FIXP16_TO_INT(fxpb);

				//modulate
				drawColor = srcBuffer[iu + iv];

				RGB565_TO_R5_G6_B5(drawColor, rDraw, gDraw, bDraw);
				rDraw = (rDraw*ir) >> 8;
				gDraw = (gDraw*ig) >> 8;
				bDraw = (bDraw*ib) >> 8;
				drawColor = R5_G6_B5_TO_RGB565(rDraw, gDraw, bDraw);

				destBuffer[xidx] = drawColor;

				fxpu += dudx;
				fxpv += dvdx;
				fxpr += drdx;
				fxpg += dgdx;
				fxpb += dbdx;
			}

			xl += dxdyl;
			xr += dxdyr;

			ul += dudyl;
			ur += dudyr;

			vl += dvdyl;
			vr += dvdyr;

			rl += drdyl;
			rr += drdyr;

			gl += dgdyl;
			gr += dgdyr;

			bl += dbdyl;
			br += dbdyr;


			destBuffer += destLineLength;
		}
	}
	else
	{
		//no clip
		for (int yidx = iys;yidx <= iye;yidx++) {
			ixs = FIXP16_TO_INT(xl);
			ixe = FIXP16_TO_INT(xr);

			fxpu = ul;
			fxpv = vl;

			fxpr = rl;
			fxpg = gl;
			fxpb = bl;

			idx = ixe - ixs;
			if (idx>0) {
				dudx = (ur - ul) / idx;
				dvdx = (vr - vl) / idx;

				drdx = (rr - rl) / idx;
				dgdx = (gr - gl) / idx;
				dbdx = (br - bl) / idx;
			}
			else
			{
				dudx = (ur - ul);
				dvdx = (vr - vl);

				drdx = (rr - rl);
				dgdx = (gr - gl);
				dbdx = (br - bl);
			}

			for (int xidx = ixs;xidx <= ixe;xidx++) {

				iu = FIXP16_TO_INT(FIXP16_ROUNDUP(fxpu << textShift));
				iv = FIXP16_TO_INT(FIXP16_ROUNDUP(fxpv << textShift)) << lineShift;
				ir = FIXP16_TO_INT(fxpr);
				ig = FIXP16_TO_INT(fxpg);
				ib = FIXP16_TO_INT(fxpb);

				//modulate
				drawColor = srcBuffer[iu + iv];

				RGB565_TO_R5_G6_B5(drawColor, rDraw, gDraw, bDraw);
				rDraw = (rDraw*ir) >> 8;
				gDraw = (gDraw*ig) >> 8;
				bDraw = (bDraw*ib) >> 8;
				drawColor = R5_G6_B5_TO_RGB565(rDraw, gDraw, bDraw);

				destBuffer[xidx] = drawColor;

				fxpu += dudx;
				fxpv += dvdx;
				fxpr += drdx;
				fxpg += dgdx;
				fxpb += dbdx;

			}
			xl += dxdyl;
			xr += dxdyr;

			ul += dudyl;
			ur += dudyr;

			vl += dvdyl;
			vr += dvdyr;


			rl += drdyl;
			rr += drdyr;

			gl += dgdyl;
			gr += dgdyr;

			bl += dbdyl;
			br += dbdyr;

			destBuffer += destLineLength;
		}

	}


	return 1;
}

int Canvas::drawPolyTriangleTextureFlat16(PolySelfTriangle4D * poly)
{
	Bitmap *texture = poly->texure;
	return drawTriangleTextureFlat16(poly->transVList[0].x, poly->transVList[0].y, poly->transVList[0].u, poly->transVList[0].v,
		poly->transVList[1].x, poly->transVList[1].y, poly->transVList[1].u, poly->transVList[1].v,
		poly->transVList[2].x, poly->transVList[2].y, poly->transVList[2].u, poly->transVList[2].v,
		poly->lightColor[0],
		(uint16 *)texture->getBuffer(), texture->getWidth(), texture->getHeight(), texture->getPitch());
}

int Canvas::drawTriangleTextureFlat16(float x0, float y0, float u0, float v0, float x1, float y1, float u1, float v1, float x2, float y2, float u2, float v2, Color c, uint16 * textBuffer, int width, int height, int pitch)
{

	float  tmp, interX, interU, interV, fraction;
	//line so return 0 as trivially reject
	if ((FEQUALE1(y0, y1) && FEQUALE1(y1, y2)) || (FEQUALE1(x0, x1) && FEQUALE1(x1, x2))) {
		return 0;
	}

	//select sort for y small -> big
	if (y1 < y0) {
		SWAP(y1, y0, tmp);
		SWAP(x1, x0, tmp);
		SWAP(u1, u0, tmp);
		SWAP(v1, v0, tmp);
	}
	if (y2 < y0) {
		SWAP(y2, y0, tmp);
		SWAP(x2, x0, tmp);
		SWAP(u2, u0, tmp);
		SWAP(v2, v0, tmp);
	}
	if (y2 < y1) {
		SWAP(y1, y2, tmp);
		SWAP(x1, x2, tmp);
		SWAP(u1, u2, tmp);
		SWAP(v1, v2, tmp);
	}

	//trivially reject
	if (y0 > cliper.bottom || y2<cliper.top || (x0<cliper.left&&x1<cliper.left&&x2<cliper.left) || (x0>cliper.right&&x1>cliper.right&&x2>cliper.right)) {
		return 0;
	}

	if (FEQUALE1(y0, y1)) {
		drawTopTriangleTextureFlat16(x2, y2, u2, v2,
			x0, y0, u0, v0,
			x1, y1, u1, v1,
			c,
			textBuffer, width, height, pitch);
	}
	else if (FEQUALE1(y1, y2))
	{
		drawBottomTriangleTextureFlat16(x0, y0, u0, v0,
			x1, y1, u1, v1,
			x2, y2, u2, v2,
			c,
			textBuffer, width, height, pitch);
	}
	else
	{
		fraction = (y1 - y0) / (y2 - y0);
		interX = x0 + (x2 - x0)*fraction;
		interU = u0 + (u2 - u0)*fraction;
		interV = v0 + (v2 - v0)*fraction;
		drawBottomTriangleTextureFlat16(x0, y0, u0, v0,
			x1, y1, u1, v1,
			interX, y1, interU, interV,
			c,
			textBuffer, width, height, pitch);
		drawTopTriangleTextureFlat16(x2, y2, u2, v2,
			x1, y1, u1, v1,
			interX, y1, interU, interV,
			c,
			textBuffer, width, height, pitch);
	}

	return 1;
}

int Canvas::drawTopTriangleTextureFlat16(float x0, float y0, float u0, float v0, float x1, float y1, float u1, float v1, float x2, float y2, float u2, float v2, Color c, uint16 * textBuffer, int width, int height, int pitch)
{
	fixp16 xl, xr,
		ul, ur,
		vl, vr,
		dxdyl, dxdyr,
		dudyl, dudyr,
		dvdyl, dvdyr,
		dudx, dvdx,
		fxpu, fxpv;
	int rBase, gBase, bBase,
		rDraw, gDraw,bDraw,
		iu, iv,
		idx, idy,
		iys, iye,
		ixs, ixe,
		destLineLength;

	float dx, dy, tmpf;
	uint16 *srcBuffer, *destBuffer,drawColor;
	int  lineShift, textShift;


	rBase = c.r;
	gBase = c.g;
	bBase = c.b;

	destLineLength = this->pitch >> 1;
	textShift = Log_Base2(width);
	lineShift = Log_Base2(width);

	if (x2 < x1) {
		SWAP(x2, x1, tmpf);
		SWAP(y2, y1, tmpf);
		SWAP(u2, u1, tmpf);
		SWAP(v2, v1, tmpf);

	}


	//pre decor has ignored  y0 == y1 == y2 so just div
	dy = y0 - y1;
	dxdyl = FIXP16_FLOAT_TO((x0 - x1) / dy);
	dxdyr = FIXP16_FLOAT_TO((x0 - x2) / dy);

	dudyl = FIXP16_FLOAT_TO((u0 - u1) / dy);
	dudyr = FIXP16_FLOAT_TO((u0 - u2) / dy);

	dvdyl = FIXP16_FLOAT_TO((v0 - v1) / dy);
	dvdyr = FIXP16_FLOAT_TO((v0 - v2) / dy);



	//clip bottom
	if (y0>cliper.bottom)
	{
		iye = cliper.bottom;
	}
	else
	{
		iye = ceil(y0) - 1;
	}

	//clip top
	if (y1 < cliper.top) {
		iys = cliper.top;

	}
	else
	{
		iys = ceil(y1);
	}

	dy = iys - y1;

	xl = FIXP16_FLOAT_TO(x1) + dy*dxdyl;
	ul = FIXP16_FLOAT_TO(u1) + dy*dudyl;
	vl = FIXP16_FLOAT_TO(v1) + dy*dvdyl;


	xr = FIXP16_FLOAT_TO(x2) + dy*dxdyr;
	ur = FIXP16_FLOAT_TO(u2) + dy*dudyr;
	vr = FIXP16_FLOAT_TO(v2) + dy*dvdyr;

	//draw

	//pre round up
	xl = FIXP16_ROUNDUP(xl);
	xr = FIXP16_ROUNDUP(xr);

	//ready for buffer
	destBuffer = (uint16 *)drawBuffer + iys*destLineLength;
	srcBuffer = textBuffer;
	if ((x0 < cliper.left) || (x0 > cliper.right) ||
		(x1 <cliper.left) || (x1 > cliper.right) ||
		(x2 < cliper.left) || (x2 >cliper.right)) {
		//should clip


		for (int yidx = iys;yidx <= iye;yidx++) {
			ixs = FIXP16_TO_INT(xl);
			ixe = FIXP16_TO_INT(xr);
			fxpu = ul;
			fxpv = vl;

			if (ixe < cliper.left || ixs>cliper.right) {


				xl += dxdyl;
				xr += dxdyr;

				ul += dudyl;
				ur += dudyr;

				vl += dvdyl;
				vr += dvdyr;

				destBuffer += destLineLength;
				continue;
			}

			idx = ixe - ixs;
			if (idx>0) {
				dudx = (ur - ul) / idx;
				dvdx = (vr - vl) / idx;
			}
			else
			{
				dudx = (ur - ul);
				dvdx = (vr - vl);
			}



			if (ixs< cliper.left) {
				idx = cliper.left - ixs;
				fxpu += idx*dudx;
				fxpv += idx*dvdx;
				ixs = cliper.left;
			}
			if (ixe > cliper.right) {
				ixe = cliper.right;
			}


			for (int xidx = ixs;xidx <= ixe;xidx++) {
				iu = FIXP16_TO_INT(FIXP16_ROUNDUP(fxpu << textShift));
				iv = FIXP16_TO_INT(FIXP16_ROUNDUP(fxpv << textShift)) << lineShift;

				drawColor = srcBuffer[iu + iv];
				RGB565_TO_R5_G6_B5(drawColor, rDraw, gDraw, bDraw);
				rDraw = (rDraw*rBase) >> 8;
				gDraw = (gDraw*gBase) >> 8;
				bDraw = (bDraw*bBase) >> 8;
				drawColor = R5_G6_B5_TO_RGB565(rDraw, gDraw, bDraw);


				destBuffer[xidx] = drawColor;

				fxpu += dudx;
				fxpv += dvdx;
			}

			xl += dxdyl;
			xr += dxdyr;

			ul += dudyl;
			ur += dudyr;

			vl += dvdyl;
			vr += dvdyr;

			destBuffer += destLineLength;
		}
	}
	else
	{
		//no clip
		for (int yidx = iys;yidx <= iye;yidx++) {
			ixs = FIXP16_TO_INT(xl);
			ixe = FIXP16_TO_INT(xr);
			fxpu = ul;
			fxpv = vl;

			idx = ixe - ixs;
			if (idx>0) {
				dudx = (ur - ul) / idx;
				dvdx = (vr - vl) / idx;
			}
			else
			{
				dudx = (ur - ul);
				dvdx = (vr - vl);
			}

			for (int xidx = ixs;xidx <= ixe;xidx++) {

				iu = FIXP16_TO_INT(FIXP16_ROUNDUP(fxpu << textShift));
				iv = FIXP16_TO_INT(FIXP16_ROUNDUP(fxpv << textShift)) << lineShift;


				drawColor = srcBuffer[iu + iv];
				RGB565_TO_R5_G6_B5(drawColor, rDraw, gDraw, bDraw);
				rDraw = (rDraw*rBase) >> 8;
				gDraw = (gDraw*gBase) >> 8;
				bDraw = (bDraw*bBase) >> 8;
				drawColor = R5_G6_B5_TO_RGB565(rDraw, gDraw, bDraw);

				destBuffer[xidx] = drawColor;

				fxpu += dudx;
				fxpv += dvdx;

			}
			xl += dxdyl;
			xr += dxdyr;

			ul += dudyl;
			ur += dudyr;

			vl += dvdyl;
			vr += dvdyr;

			destBuffer += destLineLength;
		}

	}


	return 1;
}

int Canvas::drawBottomTriangleTextureFlat16(float x0, float y0, float u0, float v0, float x1, float y1, float u1, float v1, float x2, float y2, float u2, float v2, Color c, uint16 * textBuffer, int width, int height, int pitch)
{

	fixp16 xl, xr,
		ul, ur,
		vl, vr,
		dxdyl, dxdyr,
		dudyl, dudyr,
		dvdyl, dvdyr,
		dudx, dvdx,
		fxpu, fxpv;
	int rBase, gBase, bBase,
		rDraw, gDraw, bDraw,
		iu, iv,
		idx, idy,
		iys, iye,
		ixs, ixe,
		destLineLength;

	float dx, dy, tmpf;
	uint16 *srcBuffer, *destBuffer,drawColor;
	int  lineShift, textShift;


	rBase = c.r;
	gBase = c.g;
	bBase = c.b;

	destLineLength = this->pitch >> 1;
	textShift = Log_Base2(width);
	lineShift = Log_Base2(width);

	if (x2 < x1) {
		SWAP(x2, x1, tmpf);
		SWAP(y2, y1, tmpf);
		SWAP(u2, u1, tmpf);
		SWAP(v2, v1, tmpf);

	}


	//pre decor has ignored  y0 == y1 == y2 so just div
	dy = y1 - y0;
	dxdyl = FIXP16_FLOAT_TO((x1 - x0) / dy);
	dxdyr = FIXP16_FLOAT_TO((x2 - x0) / dy);

	dudyl = FIXP16_FLOAT_TO((u1 - u0) / dy);
	dudyr = FIXP16_FLOAT_TO((u2 - v0) / dy);

	dvdyl = FIXP16_FLOAT_TO((v1 - v0) / dy);
	dvdyr = FIXP16_FLOAT_TO((v2 - v0) / dy);



	//clip bottom
	if (y1>cliper.bottom)
	{
		iye = cliper.bottom;
	}
	else
	{
		iye = ceil(y1) - 1;
	}

	//clip top
	if (y0 < cliper.top) {
		iys = cliper.top;

	}
	else
	{
		iys = ceil(y0);
	}

	dy = iys - y0;

	xl = FIXP16_FLOAT_TO(x0) + dy*dxdyl;
	ul = FIXP16_FLOAT_TO(u0) + dy*dudyl;
	vl = FIXP16_FLOAT_TO(v0) + dy*dvdyl;


	xr = FIXP16_FLOAT_TO(x0) + dy*dxdyr;
	ur = FIXP16_FLOAT_TO(u0) + dy*dudyr;
	vr = FIXP16_FLOAT_TO(v0) + dy*dvdyr;

	//draw

	//pre round up
	xl = FIXP16_ROUNDUP(xl);
	xr = FIXP16_ROUNDUP(xr);

	//ready for buffer
	destBuffer = (uint16 *)drawBuffer + iys*destLineLength;
	srcBuffer = textBuffer;
	if ((x0 < cliper.left) || (x0 > cliper.right) ||
		(x1 <cliper.left) || (x1 > cliper.right) ||
		(x2 < cliper.left) || (x2 >cliper.right)) {
		//should clip


		for (int yidx = iys;yidx <= iye;yidx++) {
			ixs = FIXP16_TO_INT(xl);
			ixe = FIXP16_TO_INT(xr);
			fxpu = ul;
			fxpv = vl;

			if (ixe < cliper.left || ixs>cliper.right) {


				xl += dxdyl;
				xr += dxdyr;

				ul += dudyl;
				ur += dudyr;

				vl += dvdyl;
				vr += dvdyr;

				destBuffer += destLineLength;
				continue;
			}

			idx = ixe - ixs;
			if (idx>0) {
				dudx = (ur - ul) / idx;
				dvdx = (vr - vl) / idx;
			}
			else
			{
				dudx = (ur - ul);
				dvdx = (vr - vl);
			}



			if (ixs< cliper.left) {
				idx = cliper.left - ixs;
				fxpu += idx*dudx;
				fxpv += idx*dvdx;
				ixs = cliper.left;
			}
			if (ixe > cliper.right) {
				ixe = cliper.right;
			}


			for (int xidx = ixs;xidx <= ixe;xidx++) {
				iu = FIXP16_TO_INT(FIXP16_ROUNDUP(fxpu << textShift));
				iv = FIXP16_TO_INT(FIXP16_ROUNDUP(fxpv << textShift)) << lineShift;
			
				drawColor = srcBuffer[iu + iv];
				RGB565_TO_R5_G6_B5(drawColor, rDraw, gDraw, bDraw);
				rDraw = (rDraw*rBase) >> 8;
				gDraw = (gDraw*gBase) >> 8;
				bDraw = (bDraw*bBase) >> 8;
				drawColor = R5_G6_B5_TO_RGB565(rDraw, gDraw, bDraw);

				destBuffer[xidx] = drawColor;
				fxpu += dudx;
				fxpv += dvdx;
			}

			xl += dxdyl;
			xr += dxdyr;

			ul += dudyl;
			ur += dudyr;

			vl += dvdyl;
			vr += dvdyr;

			destBuffer += destLineLength;
		}
	}
	else
	{
		//no clip
		for (int yidx = iys;yidx <= iye;yidx++) {
			ixs = FIXP16_TO_INT(xl);
			ixe = FIXP16_TO_INT(xr);
			fxpu = ul;
			fxpv = vl;

			idx = ixe - ixs;
			if (idx>0) {
				dudx = (ur - ul) / idx;
				dvdx = (vr - vl) / idx;
			}
			else
			{
				dudx = (ur - ul);
				dvdx = (vr - vl);
			}

			for (int xidx = ixs;xidx <= ixe;xidx++) {

				iu = FIXP16_TO_INT(FIXP16_ROUNDUP(fxpu << textShift));
				iv = FIXP16_TO_INT(FIXP16_ROUNDUP(fxpv << textShift)) << lineShift;

				drawColor = srcBuffer[iu + iv];
				RGB565_TO_R5_G6_B5(drawColor, rDraw, gDraw, bDraw);
				rDraw = (rDraw*rBase) >> 8;
				gDraw = (gDraw*gBase) >> 8;
				bDraw = (bDraw*bBase) >> 8;
				drawColor = R5_G6_B5_TO_RGB565(rDraw, gDraw, bDraw);


				destBuffer[xidx] = drawColor;

				fxpu += dudx;
				fxpv += dvdx;

			}
			xl += dxdyl;
			xr += dxdyr;

			ul += dudyl;
			ur += dudyr;

			vl += dvdyl;
			vr += dvdyr;

			destBuffer += destLineLength;
		}

	}


	return 1;
}

int Canvas::drawPolyTriangleTexture16(PolySelfTriangle4D * poly)
{	
	Bitmap *texture = poly->texure;
	return drawTriangleTexture16(poly->transVList[0].x,poly->transVList[0].y,poly->transVList[0].u,poly->transVList[0].v,
								 poly->transVList[1].x, poly->transVList[1].y, poly->transVList[1].u, poly->transVList[1].v,
								 poly->transVList[2].x, poly->transVList[2].y, poly->transVList[2].u, poly->transVList[2].v,
								 (uint16 *)texture->getBuffer(),texture->getWidth(),texture->getHeight(),texture->getPitch());
}

int Canvas::drawTriangleTexture16(float x0, float y0, float u0, float v0, float x1, float y1, float u1, float v1, float x2, float y2, float u2, float v2, uint16* textBuffer, int width, int height, int pitch)
{
	float  tmp, interX,interU,interV,fraction;
	//line so return 0 as trivially reject
	if ((FEQUALE1(y0, y1) && FEQUALE1(y1, y2)) || (FEQUALE1(x0, x1) && FEQUALE1(x1, x2))) {
		return 0;
	}

	//select sort for y small -> big
	if (y1 < y0) {
		SWAP(y1, y0, tmp);
		SWAP(x1, x0, tmp);
		SWAP(u1, u0, tmp);
		SWAP(v1, v0, tmp);
	}
	if (y2 < y0) {
		SWAP(y2, y0, tmp);
		SWAP(x2, x0, tmp);
		SWAP(u2, u0, tmp);
		SWAP(v2, v0, tmp);
	}
	if (y2 < y1) {
		SWAP(y1, y2, tmp);
		SWAP(x1, x2, tmp);
		SWAP(u1, u2, tmp);
		SWAP(v1, v2, tmp);
	}

	//trivially reject
	if (y0 > cliper.bottom || y2<cliper.top || (x0<cliper.left&&x1<cliper.left&&x2<cliper.left) || (x0>cliper.right&&x1>cliper.right&&x2>cliper.right)) {
		return 0;
	}

	if (FEQUALE1(y0, y1)) {
		drawTopTriangleTexture16(x2, y2,u2,v2, 
								 x0, y0,u0,v0,
								 x1, y1, u1,v1,
								 textBuffer,width,height,pitch);
	}
	else if (FEQUALE1(y1, y2))
	{
		drawBottomTriangleTexture16(x0, y0,u0,v0, 
									x1, y1,u1,v1 ,
									x2, y2,u2,v2,
									textBuffer, width, height, pitch);
	}
	else
	{
		fraction = (y1 - y0) / (y2 - y0);
		interX = x0 + (x2 - x0)*fraction;
		interU = u0 + (u2 - u0)*fraction;
		interV = v0 + (v2 - v0)*fraction;
		drawBottomTriangleTexture16(x0, y0,u0,v0,
									x1, y1,u1,v1,
									interX, y1, interU,interV,
									textBuffer, width, height, pitch);
		drawTopTriangleTexture16(x2, y2,u2,v2,
								 x1, y1,u1,v1, 
								interX, y1,interU,interV, 
								textBuffer, width, height, pitch);
	}

	return 1;
}

int Canvas::drawTopTriangleTexture16(float x0, float y0, float u0, float v0, float x1, float y1, float u1, float v1, float x2, float y2, float u2, float v2, uint16* textBuffer, int width, int height, int pitch)
{	
	fixp16 xl, xr,  
		ul, ur,
		vl, vr,
		dxdyl, dxdyr,
		dudyl, dudyr,
		dvdyl, dvdyr, 
		dudx,dvdx,
		fxpu,fxpv;
	int iu,iv,
		idx, idy, 
		iys, iye, 
		ixs, ixe,
		destLineLength;

	float dx, dy, tmpf;
	uint16 *srcBuffer,*destBuffer;
	int  lineShift, textShift;


	destLineLength = this->pitch >> 1;
	textShift = Log_Base2(width);
	lineShift = Log_Base2(width);

	if (x2 < x1) {
		SWAP(x2, x1, tmpf);
		SWAP(y2, y1, tmpf);
		SWAP(u2, u1, tmpf);
		SWAP(v2, v1, tmpf);

	}


	//pre decor has ignored  y0 == y1 == y2 so just div
	dy = y0 - y1;
	dxdyl = FIXP16_FLOAT_TO((x0 - x1) / dy);
	dxdyr = FIXP16_FLOAT_TO((x0 - x2) / dy);

	dudyl = FIXP16_FLOAT_TO((u0 - u1) / dy);
	dudyr = FIXP16_FLOAT_TO((u0 - u2) / dy);

	dvdyl = FIXP16_FLOAT_TO((v0 - v1) / dy);
	dvdyr = FIXP16_FLOAT_TO((v0 - v2) / dy);



	//clip bottom
	if (y0>cliper.bottom)
	{
		iye = cliper.bottom;
	}
	else
	{
		iye = ceil(y0) - 1;
	}

	//clip top
	if (y1 < cliper.top) {
		iys = cliper.top;

	}
	else
	{
		iys = ceil(y1);
	}

	dy = iys - y1;

	xl = FIXP16_FLOAT_TO(x1) + dy*dxdyl;
	ul = FIXP16_FLOAT_TO(u1) + dy*dudyl;
	vl = FIXP16_FLOAT_TO(v1) + dy*dvdyl;


	xr = FIXP16_FLOAT_TO(x2) + dy*dxdyr;
	ur = FIXP16_FLOAT_TO(u2) + dy*dudyr;
	vr = FIXP16_FLOAT_TO(v2) + dy*dvdyr;
	
	//draw

	//pre round up
	xl = FIXP16_ROUNDUP(xl);
	xr = FIXP16_ROUNDUP(xr);

	//ready for buffer
	destBuffer = (uint16 *)drawBuffer + iys*destLineLength;
	srcBuffer = textBuffer;
	if ((x0 < cliper.left) || (x0 > cliper.right) ||
		(x1 <cliper.left) || (x1 > cliper.right) ||
		(x2 < cliper.left) || (x2 >cliper.right)) {
		//should clip


		for (int yidx = iys;yidx <= iye;yidx++) {
			ixs = FIXP16_TO_INT(xl);
			ixe = FIXP16_TO_INT(xr);
			fxpu = ul;
			fxpv = vl;

			if (ixe < cliper.left || ixs>cliper.right) {

			
				xl += dxdyl;
				xr += dxdyr;

				ul += dudyl;
				ur += dudyr;

				vl += dvdyl;
				vr += dvdyr;

				destBuffer += destLineLength;
				continue;
			}

			idx = ixe - ixs;
			if (idx>0) {
				dudx = (ur - ul) / idx;
				dvdx = (vr - vl) / idx;
			}
			else
			{
				dudx = (ur - ul);
				dvdx = (vr - vl);
			}



			if (ixs< cliper.left) {
				idx = cliper.left - ixs;
				fxpu += idx*dudx;
				fxpv += idx*dvdx;
				ixs = cliper.left;
			}
			if (ixe > cliper.right) {
				ixe = cliper.right;
			}


			for (int xidx = ixs;xidx <= ixe;xidx++) {

				iu = FIXP16_TO_INT(FIXP16_ROUNDUP(fxpu << textShift));
				iv = FIXP16_TO_INT(FIXP16_ROUNDUP(fxpv << textShift)) << lineShift;
				destBuffer[xidx] = srcBuffer[iu+iv];
				fxpu += dudx;
				fxpv += dvdx;
			}

			xl += dxdyl;
			xr += dxdyr;
			
			ul += dudyl;
			ur += dudyr;

			vl += dvdyl;
			vr += dvdyr;

			destBuffer += destLineLength;
		}
	}
	else
	{
		//no clip
		for (int yidx = iys;yidx <= iye;yidx++) {
			ixs = FIXP16_TO_INT(xl);
			ixe = FIXP16_TO_INT(xr);
			fxpu = ul;
			fxpv = vl;
			
			idx = ixe - ixs;
			if (idx>0) {
				dudx = (ur - ul) / idx;
				dvdx = (vr - vl) / idx;
			}
			else
			{
				dudx = (ur - ul);
				dvdx = (vr - vl);
			}

			for (int xidx = ixs;xidx <= ixe;xidx++) {

				iu = FIXP16_TO_INT(FIXP16_ROUNDUP(fxpu << textShift));
				iv = FIXP16_TO_INT(FIXP16_ROUNDUP(fxpv << textShift))<<lineShift; 
				destBuffer[xidx] = srcBuffer[iu + iv];
				fxpu += dudx;
				fxpv += dvdx;

			}
			xl += dxdyl;
			xr += dxdyr;

			ul += dudyl;
			ur += dudyr;

			vl += dvdyl;
			vr += dvdyr;

			destBuffer += destLineLength;
		}

	}


	return 1;
}

int Canvas::drawBottomTriangleTexture16(float x0, float y0, float u0, float v0, float x1, float y1, float u1, float v1, float x2, float y2, float u2, float v2, uint16* textBuffer, int width, int height, int pitch)
{
	fixp16 xl, xr,
		ul, ur,
		vl, vr,
		dxdyl, dxdyr,
		dudyl, dudyr,
		dvdyl, dvdyr,
		dudx, dvdx,
		fxpu, fxpv;
	int iu, iv,
		idx, idy,
		iys, iye,
		ixs, ixe,
		destLineLength;

	float dx, dy, tmpf;
	uint16 *srcBuffer, *destBuffer;
	int  lineShift,textShift;


	destLineLength = this->pitch >> 1;
	textShift = Log_Base2(width);
	lineShift = Log_Base2(width);

	if (x2 < x1) {
		SWAP(x2, x1, tmpf);
		SWAP(y2, y1, tmpf);
		SWAP(u2, u1, tmpf);
		SWAP(v2, v1, tmpf);

	}


	//pre decor has ignored  y0 == y1 == y2 so just div
	dy = y1 - y0;
	dxdyl = FIXP16_FLOAT_TO((x1 - x0) / dy);
	dxdyr = FIXP16_FLOAT_TO((x2 - x0) / dy);

	dudyl = FIXP16_FLOAT_TO((u1 - u0) / dy);
	dudyr = FIXP16_FLOAT_TO((u2 - v0) / dy);

	dvdyl = FIXP16_FLOAT_TO((v1 - v0) / dy);
	dvdyr = FIXP16_FLOAT_TO((v2 - v0) / dy);



	//clip bottom
	if (y1>cliper.bottom)
	{
		iye = cliper.bottom;
	}
	else
	{
		iye = ceil(y1) - 1;
	}

	//clip top
	if (y0 < cliper.top) {
		iys = cliper.top;

	}
	else
	{
		iys = ceil(y0);
	}

	dy = iys - y0;

	xl = FIXP16_FLOAT_TO(x0) + dy*dxdyl;
	ul = FIXP16_FLOAT_TO(u0) + dy*dudyl;
	vl = FIXP16_FLOAT_TO(v0) + dy*dvdyl;


	xr = FIXP16_FLOAT_TO(x0) + dy*dxdyr;
	ur = FIXP16_FLOAT_TO(u0) + dy*dudyr;
	vr = FIXP16_FLOAT_TO(v0) + dy*dvdyr;

	//draw

	//pre round up
	xl = FIXP16_ROUNDUP(xl);
	xr = FIXP16_ROUNDUP(xr);

	//ready for buffer
	destBuffer = (uint16 *)drawBuffer + iys*destLineLength;
	srcBuffer = textBuffer;
	if ((x0 < cliper.left) || (x0 > cliper.right) ||
		(x1 <cliper.left) || (x1 > cliper.right) ||
		(x2 < cliper.left) || (x2 >cliper.right)) {
		//should clip


		for (int yidx = iys;yidx <= iye;yidx++) {
			ixs = FIXP16_TO_INT(xl);
			ixe = FIXP16_TO_INT(xr);
			fxpu = ul;
			fxpv = vl;

			if (ixe < cliper.left || ixs>cliper.right) {


				xl += dxdyl;
				xr += dxdyr;

				ul += dudyl;
				ur += dudyr;

				vl += dvdyl;
				vr += dvdyr;

				destBuffer += destLineLength;
				continue;
			}

			idx = ixe - ixs;
			if (idx>0) {
				dudx = (ur - ul) / idx;
				dvdx = (vr - vl) / idx;
			}
			else
			{
				dudx = (ur - ul);
				dvdx = (vr - vl);
			}



			if (ixs< cliper.left) {
				idx = cliper.left - ixs;
				fxpu += idx*dudx;
				fxpv += idx*dvdx;
				ixs = cliper.left;
			}
			if (ixe > cliper.right) {
				ixe = cliper.right;
			}


			for (int xidx = ixs;xidx <= ixe;xidx++) {

				iu = FIXP16_TO_INT(FIXP16_ROUNDUP(fxpu << textShift));
				iv = FIXP16_TO_INT(FIXP16_ROUNDUP(fxpv << textShift)) << lineShift;
				destBuffer[xidx] = srcBuffer[iu + iv];
				fxpu += dudx;
				fxpv += dvdx;
			}

			xl += dxdyl;
			xr += dxdyr;

			ul += dudyl;
			ur += dudyr;

			vl += dvdyl;
			vr += dvdyr;

			destBuffer += destLineLength;
		}
	}
	else
	{
		//no clip
		for (int yidx = iys;yidx <= iye;yidx++) {
			ixs = FIXP16_TO_INT(xl);
			ixe = FIXP16_TO_INT(xr);
			fxpu = ul;
			fxpv = vl;

			idx = ixe - ixs;
			if (idx>0) {
				dudx = (ur - ul) / idx;
				dvdx = (vr - vl) / idx;
			}
			else
			{
				dudx = (ur - ul);
				dvdx = (vr - vl);
			}

			for (int xidx = ixs;xidx <= ixe;xidx++) {

				iu = FIXP16_TO_INT(FIXP16_ROUNDUP(fxpu << textShift));
				iv = FIXP16_TO_INT(FIXP16_ROUNDUP(fxpv << textShift)) << lineShift;
				destBuffer[xidx] = srcBuffer[iu + iv];
				fxpu += dudx;
				fxpv += dvdx;

			}
			xl += dxdyl;
			xr += dxdyr;

			ul += dudyl;
			ur += dudyr;

			vl += dvdyl;
			vr += dvdyr;

			destBuffer += destLineLength;
		}

	}


	return 1;
}

int Canvas::drawPolyTriangleGouraud16(PolySelfTriangle4D * poly)
{
	return drawTriangleGouraud16(poly->transVList[0].x, poly->transVList[0].y, //
		poly->transVList[1].x, poly->transVList[1].y,//
		poly->transVList[2].x, poly->transVList[2].y,//
		poly->lightColor[0], poly->lightColor[1], poly->lightColor[2]);
}

int Canvas::drawTriangleGouraud16(float x0, float y0, float x1, float y1, float x2, float y2, Color c0, Color c1, Color c2)
{
	float  tmp,interX,fraction;
	Color tmpC, interC;
	//line so return 0 as trivially reject
	if ((FEQUALE1(y0, y1) && FEQUALE1(y1, y2)) || (FEQUALE1(x0, x1) && FEQUALE1(x1, x2))) {
		return 0;
	}

	//select sort for y small -> big
	if (y1 < y0) {
		SWAP(y1, y0, tmp);
		SWAP(x1, x0, tmp);
		SWAP(c1, c0, tmpC);
	}
	if (y2 < y0) {
		SWAP(y2, y0, tmp);
		SWAP(x2, x0, tmp);
		SWAP(c2, c0, tmpC);
	}
	if (y2 < y1) {
		SWAP(y1, y2, tmp);
		SWAP(x1, x2, tmp);
		SWAP(c1, c2, tmpC);
	}

	//trivially reject
	if (y0 > cliper.bottom || y2<cliper.top || (x0<cliper.left&&x1<cliper.left&&x2<cliper.left) || (x0>cliper.right&&x1>cliper.right&&x2>cliper.right)) {
		return 0;
	}

	if (FEQUALE1(y0, y1)) {
		drawTopTriangleGouraud16(x2, y2, x0, y0, x1, y1, c2,c0,c1);
	}
	else if (FEQUALE1(y1, y2))
	{
		drawBottomTriangleGouraud16(x0, y0, x1, y1, x2, y2, c0,c1,c2);
	}
	else
	{
		fraction = (y1 - y0) / (y2 - y0);
		interX = x0 + (x2 - x0)*fraction;
		interC.r = c0.r + (c2.r - c0.r)*fraction+0.5;
		interC.g = c0.g + (c2.g - c0.g)*fraction+0.5;
		interC.b = c0.b + (c2.b - c0.b)*fraction +0.5;
		drawBottomTriangleGouraud16(x0, y0, x1, y1, interX, y1,c0,c1 ,interC);
		drawTopTriangleGouraud16(x2, y2, x1, y1, interX, y1,c2,c1,interC);
	}

	return 1;
}

int Canvas::drawTopTriangleGouraud16(float x0, float y0, float x1, float y1, float x2, float y2, Color c0, Color c1, Color c2)
{
	fixp16 xl, xr, dxdyl, dxdyr, rl,rr,gl,gr,bl,br,drdyl,drdyr,dgdyl,dgdyr,dbdyl,dbdyr,drdx,dgdx,dbdx,fxpr, fxpg, fxpb;
	int idx,idy,iys,iye,ixs,ixe,r0, g0, b0, r1, g1, b1, r2, g2, b2,tmpi,lingLength;
	float dx,dy,tmpf;
	uint16 *buffer;

	lingLength = pitch >> 1;

	r0 = c0.r;
	r1 = c1.r;
	r2 = c2.r;

	g0 = c0.g;
	g1 = c1.g;
	g2 = c2.g;

	b0 = c0.b;
	b1 = c1.b;
	b2 = c2.b;

	if (x2 < x1) {
		SWAP(x2, x1, tmpf);
		SWAP(y2, y1, tmpf);
		SWAP(r2, r1, tmpi);
		SWAP(g2, g1, tmpi);
		SWAP(b2, b1, tmpi);
	}




	//drawTraingleGouraud16 has ignored  y0 == y1 == y2 so just div
	dy = y0 - y1;
	dxdyl = FIXP16_FLOAT_TO((x0-x1) / dy);
	dxdyr = FIXP16_FLOAT_TO((x0-x2) / dy);
	
	drdyl = FIXP16_FLOAT_TO((r0 - r1) / dy);
	drdyr = FIXP16_FLOAT_TO((r0 - r2) / dy);

	dgdyl = FIXP16_FLOAT_TO((g0 - g1) / dy);
	dgdyr = FIXP16_FLOAT_TO((g0 - g2) / dy);

	dbdyl = FIXP16_FLOAT_TO((b0 - b1) / dy);
	dbdyr = FIXP16_FLOAT_TO((b0 - b2) / dy);

	//clip bottom
	if(y0>cliper.bottom)
	{
		iye = cliper.bottom;
	}else
	{
		iye = ceil(y0) - 1;
	}

	//clip top
	if (y1 < cliper.top) {
		iys = cliper.top;
		
	}
	else
	{
		iys = ceil(y1);
	}

	dy = iys - y1;

	xl = FIXP16_FLOAT_TO(x1) + dy*dxdyl;
	rl = FIXP16_INT_TO(r1) + dy*drdyl;
	gl = FIXP16_INT_TO(g1) + dy*dgdyl;
	bl = FIXP16_INT_TO(b1) + dy*dbdyl;

	xr = FIXP16_FLOAT_TO(x2) + dy*dxdyr;
	rr = FIXP16_INT_TO(r2) + dy*drdyr;
	gr = FIXP16_INT_TO(g2) + dy*dgdyr;
	br = FIXP16_INT_TO(b2) + dy*dbdyr;

	//draw

	//pre round up
	xl = FIXP16_ROUNDUP(xl);
	xr = FIXP16_ROUNDUP(xr);
	rl = FIXP16_ROUNDUP(rl);
	rr = FIXP16_ROUNDUP(rr);
	gl = FIXP16_ROUNDUP(gl);
	gr = FIXP16_ROUNDUP(gr);
	bl = FIXP16_ROUNDUP(bl);
	br = FIXP16_ROUNDUP(br);

	buffer = (uint16 *)drawBuffer+iys*lingLength;
	if ((x0 < cliper.left) || (x0 > cliper.right) ||
		(x1 <cliper.left) || (x1 > cliper.right) ||
		(x2 < cliper.left) || (x2 >cliper.right)) {
		//should clip


		for (int yidx = iys;yidx <= iye;yidx++, buffer += lingLength) {
			ixs = FIXP16_TO_INT(xl);
			ixe = FIXP16_TO_INT(xr);
			fxpr = rl;
			fxpg = gl;
			fxpb = bl;

			if (ixe < cliper.left || ixs>cliper.right) {

				xl += dxdyl;
				xr += dxdyr;

				rl += drdyl;
				rr += drdyr;

				gl += dgdyl;
				gr += dgdyr;

				bl += dbdyl;
				br += dbdyr;


				continue;
			}

			idx = ixe- ixs;
			if (idx>0) {
				drdx = (rr - rl) / idx;
				dgdx = (gr - gl) / idx;
				dbdx = (br - bl) / idx;
			}
			else
			{
				drdx = rr - rl;
				dgdx = gr - gl;
				dbdx = br - gl;
			}





			if (ixs< cliper.left) {
				idx = cliper.left-ixs;
				fxpr += idx*drdx;
				fxpg += idx*dgdx;
				fxpb += idx*dbdx;

				ixs = cliper.left;
			}
			if (ixe > cliper.right) {
				ixe = cliper.right;
			}

			
			for (int xidx = ixs;xidx <= ixe;xidx++) {
				uint16 color = R_G_B_TO_RGB565_FIXP16(fxpr, fxpg, fxpb);
				buffer[xidx] = color;
				fxpr += drdx;
				fxpg += dgdx;
				fxpb += dbdx;
			}

			xl += dxdyl;
			xr += dxdyr;

			rl += drdyl;
			rr += drdyr;

			gl += dgdyl;
			gr += dgdyr;

			bl += dbdyl;
			br += dbdyr;
		}
	}
	else
	{
		//no clip
		for (int yidx = iys;yidx <= iye;yidx++, buffer += lingLength) {
			ixs = FIXP16_TO_INT(xl);
			ixe = FIXP16_TO_INT(xr);
			fxpr = rl;
			fxpg = gl;
			fxpb = bl;

			idx = ixe - ixs;
			if (idx>0) {
				drdx = (rr - rl) / idx;
				dgdx = (gr - gl) / idx;
				dbdx = (br - bl) / idx;
			}
			else
			{
				drdx = rr - rl;
				dgdx = gr - gl;
				dbdx = br - gl;
			}
	
			for (int xidx = ixs;xidx <= ixe;xidx++) {
				uint16 color = R_G_B_TO_RGB565_FIXP16(fxpr, fxpg, fxpb);
				if (color > 0xdddd) {
					int d = 0;
					d = 220;
				}
				buffer[xidx] = color;
				fxpr += drdx;
				fxpg += dgdx;
				fxpb += dbdx;
			}

			xl += dxdyl;
			xr += dxdyr;

			rl += drdyl;
			rr += drdyr;

			gl += dgdyl;
			gr += dgdyr;

			bl += dbdyl;
			br += dbdyr;
		}

	}


	return 1;
}

int Canvas::drawBottomTriangleGouraud16(float x0, float y0, float x1, float y1, float x2, float y2, Color c0, Color c1, Color c2)
{
	fixp16 xl, xr, rl, rr, gl, gr, bl, br, dxdyl, dxdyr, drdyl, drdyr, dgdyl, dgdyr, dbdyl, dbdyr, drdx, dgdx, dbdx, fxpr, fxpg, fxpb;
	int idx, iys, iye, ixs, ixe,ix0,ix1,ix2,iy0,iy1,iy2,r0, g0, b0, r1, g1, b1, r2, g2, b2, tmpi, lingLength;
	float dx, dy,tmpf;
	uint16 *buffer;

	lingLength = pitch >> 1;



	r0 = c0.r;
	r1 = c1.r;
	r2 = c2.r;

	g0 = c0.g;
	g1 = c1.g;
	g2 = c2.g;

	b0 = c0.b;
	b1 = c1.b;
	b2 = c2.b;

	if (x2 < x1) {
		SWAP(x2, x1, tmpf);
		SWAP(y2, y1, tmpf);
		SWAP(r2, r1, tmpi);
		SWAP(g2, g1, tmpi);
		SWAP(b2, b1, tmpi);
	}
	//drawTraingleGouraud16 has ignored  y0 == y1 == y2 so just div
	dy = y1 - y0;
	dxdyl = FIXP16_FLOAT_TO((x1 - x0) / dy);
	dxdyr = FIXP16_FLOAT_TO((x2 - x0) / dy);

	drdyl = FIXP16_FLOAT_TO((r1 - r0) / dy);
	drdyr = FIXP16_FLOAT_TO((r2 - r0) / dy);

	dgdyl = FIXP16_FLOAT_TO((g1 - g0) / dy);
	dgdyr = FIXP16_FLOAT_TO((g2 - g0) / dy);

	dbdyl = FIXP16_FLOAT_TO((b1 - b0) / dy);
	dbdyr = FIXP16_FLOAT_TO((b2 - b0) / dy);

	//clip bottom
	if (y1>cliper.bottom)
	{
		iye = cliper.bottom;
	}
	else
	{
		iye = ceil(y1) - 1;
	}

	//clip top
	if (y0 < cliper.top) {
		iys = cliper.top;

	}
	else
	{
		iys = ceil(y0);
	}

	dy = iys - y0;

	xl = FIXP16_FLOAT_TO(x0) + dy*dxdyl;
	rl = FIXP16_INT_TO(r0) + dy*drdyl;
	gl = FIXP16_INT_TO(g0) + dy*dgdyl;
	bl = FIXP16_INT_TO(b0) + dy*dbdyl;

	xr = FIXP16_FLOAT_TO(x0) + dy*dxdyr;
	rr = FIXP16_INT_TO(r0) + dy*drdyr;
	gr = FIXP16_INT_TO(g0) + dy*dgdyr;
	br = FIXP16_INT_TO(b0) + dy*dbdyr;

	//draw

	//pre round up
	xl = FIXP16_ROUNDUP(xl);
	xr = FIXP16_ROUNDUP(xr);
	rl = FIXP16_ROUNDUP(rl);
	rr = FIXP16_ROUNDUP(rr);
	gl = FIXP16_ROUNDUP(gl);
	gr = FIXP16_ROUNDUP(gr);
	bl = FIXP16_ROUNDUP(bl);
	br = FIXP16_ROUNDUP(br);

	buffer = (uint16 *)drawBuffer + iys*lingLength;

	if ((x0 < cliper.left) || (x0 > cliper.right) ||
		(x1 <cliper.left) || (x1 > cliper.right) ||
		(x2 < cliper.left) || (x2 >cliper.right)) {
		//should clip

		for (int yidx = iys;yidx <= iye;yidx++, buffer += lingLength) {
			ixs = FIXP16_TO_INT(xl);
			ixe = FIXP16_TO_INT(xr);
			fxpr = rl;
			fxpg = gl;
			fxpb = bl;

			if (ixe < cliper.left || ixs>cliper.right) {

				xl += dxdyl;
				xr += dxdyr;

				rl += drdyl;
				rr += drdyr;

				gl += dgdyl;
				gr += dgdyr;

				bl += dbdyl;
				br += dbdyr;

				continue;
			}

			idx = ixe - ixs;
			if (idx>0) {
				drdx = (rr - rl) / idx;
				dgdx = (gr - gl) / idx;
				dbdx = (br - bl) / idx;
			}
			else
			{
				drdx = rr - rl;
				dgdx = gr - gl;
				dbdx = br - gl;
			}
			if (ixs < cliper.left) {
				idx = cliper.left - ixs;
				fxpr += idx*drdx;
				fxpg += idx*dgdx;
				fxpb += idx*dbdx;
				ixs = cliper.left;
			}
			if (ixe  > cliper.right) {
				ixe = cliper.right;
			}

			for (int xidx = ixs;xidx <= ixe;xidx++) {
				buffer[xidx] = R_G_B_TO_RGB565_FIXP16(fxpr, fxpg, fxpb);
				fxpr += drdx;
				fxpg += dgdx;
				fxpb += dbdx;
			}

			xl += dxdyl;
			xr += dxdyr;

			rl += drdyl;
			rr += drdyr;

			gl += dgdyl;
			gr += dgdyr;

			bl += dbdyl;
			br += dbdyr;
		}
	}
	else
	{
		//no clip
		for (int yidx = iys;yidx <= iye;yidx++, buffer += lingLength) {
			ixs = FIXP16_TO_INT(xl);
			ixe = FIXP16_TO_INT(xr);
			fxpr = rl;
			fxpg = gl;
			fxpb = bl;

			idx = ixe - ixs;
			if (idx>0) {
				drdx = (rr - rl) / idx;
				dgdx = (gr - gl) / idx;
				dbdx = (br - bl) / idx;
			}
			else
			{
				drdx = rr - rl;
				dgdx = gr - gl;
				dbdx = br - gl;
			}

			
			for (int xidx = ixs;xidx <= ixe;xidx++) {
				buffer[xidx] = R_G_B_TO_RGB565_FIXP16(fxpr, fxpg, fxpb);
				fxpr += drdx;
				fxpg += dgdx;
				fxpb += dbdx;
			}

			xl += dxdyl;
			xr += dxdyr;

			rl += drdyl;
			rr += drdyr;

			gl += dgdyl;
			gr += dgdyr;

			bl += dbdyl;
			br += dbdyr;
		}

	}


	return 1;
}

int Canvas::drawPolyTriangleSolidF16(PolySelfTriangle4D * poly)
{
	return drawTriangleSolidF16(poly->transVList[0].x, poly->transVList[0].y, //
		poly->transVList[1].x, poly->transVList[1].y,//
		poly->transVList[2].x, poly->transVList[2].y,//
		poly->lightColor[0]);
}


int Canvas::drawTriangleSolidF16(float x0, float y0, float x1, float y1, float x2, float y2, Color color)
{
	float  tmp;
	//line so return 0 as trivially reject
	if ((FEQUALE1(y0,y1)&& FEQUALE1(y1, y2)) ||(FEQUALE1(x0, x1) && FEQUALE1(x1, x2))){
		return 0;
	}

	//select sort for y small -> big
	if (y1 < y0) {
		SWAP(y1, y0,tmp);
		SWAP(x1, x0,tmp);
	}
	if (y2 < y0) {
		SWAP(y2, y0, tmp);
		SWAP(x2, x0, tmp);
	}
	if (y2 < y1) {
		SWAP(y1, y2, tmp);
		SWAP(x1, x2, tmp);
	}

	//trivially reject
	if (y0 > cliper.bottom || y2<cliper.top || (x0<cliper.left&&x1<cliper.left&&x2<cliper.left) || (x0>cliper.right&&x1>cliper.right&&x2>cliper.right)) {
		return 0;
	}

	if (FEQUALE1(y0, y1)) {
		drawTopTriangleF16(x2, y2, x0, y0, x1, y1, color);
	}
	else if (FEQUALE1(y1, y2))
	{
		drawBottomTriangleF16(x0, y0, x1, y1, x2, y2, color);
	}
	else
	{
		float newX;
		newX = x0 + (y1 - y0)*(x2 - x0) / (y2 - y0);
		drawBottomTriangleF16(x0, y0, x1, y1, newX, y1, color);
		drawTopTriangleF16(x2, y2, x1, y1, newX, y1, color);
	}

	return 1;
}

int Canvas::drawTopTriangleF16(float x0, float y0, float x1, float y1, float x2, float y2, Color c)
{
	float dy,leftXStep,rightXStep,xs,xe;
	int iys, iye,ixs,ixe,lineLength;
	uint16 color = c.getRGB566();

	float tmp;
	if (x1 > x2) {
		SWAP(x1, x2, tmp);
		SWAP(y1, y2, tmp);
	}

	lineLength = pitch >> 1;
	dy = y0 - y1;
	leftXStep = (x0 - x1) / dy; //by ystep 1
	rightXStep = (x0 - x2) / dy;
	
	xs = x1+0.5;//pre roundup
	xe = x2+0.5;
	if (y1 < cliper.top) {
		iys = cliper.top;
		xs = xs + (iys - y1)*leftXStep;
		xe = xe + (iys - y1)*rightXStep;
	}
	else
	{
		iys = ceil(y1);
		xs = xs + (iys - y1)*leftXStep;
		xe = xe + (iys - y1)*rightXStep;
	}

	if (y0 > cliper.bottom) {
		iye = cliper.bottom; //if clip don't left-top
	}
	else
	{
		iye = ceil(y0) - 1;
	}


	uint16 *buffer = (uint16*)drawBuffer + iys*lineLength;
	//test clip needed
	if (x0 >= cliper.left&&x0 <= cliper.right&&xs >= cliper.left&&xs <= cliper.right&&xe >= cliper.left&&xe <= cliper.right) {
		//not need clip
		for (int i = iys;i <= iye;i++, buffer += lineLength) {
			ixs = xs;
			ixe = xe;
			memSet16(buffer + ixs, color, (ixe - ixs + 1));
			xs += leftXStep;
			xe += rightXStep;
		}
	}
	else {
		for (int i = iys; i <= iye; i++, buffer += lineLength)//step
		{
			ixs = xs;
			ixe = xe;
			//pre step for continue
			xs += leftXStep;
			xe += rightXStep;
			//trivially reject
			if (ixe < cliper.left || ixs>cliper.right) {
				continue;
			}

			//clip
			ixs = ixs < cliper.left ? cliper.left : ixs;
			ixe = ixe> cliper.right ? cliper.right : ixe;
			memSet16(buffer + ixs, color, (ixe - ixs + 1));

		}
	}


	return 1;
}

int Canvas::drawBottomTriangleF16(float x0, float y0, float x1, float y1, float x2, float y2, Color c)
{	
	float dy, leftXStep, rightXStep, xs, xe;
	int iys, iye, ixs, ixe, lineLength;
	float tmp;
	uint16 color = c.getRGB566();
	if (x1 > x2) {
		SWAP(x1, x2, tmp);
		SWAP(y1, y2, tmp);
	}

	lineLength = pitch >> 1;
	dy = y1 - y0;
	leftXStep = (x1 - x0) / dy;
	rightXStep = (x2 - x0) / dy;

	xs = x0+0.5f;//pre round up
	xe = x0+0.5f;
	if (y0 < cliper.top) {
		iys = cliper.top;
		xs = xs + (iys - y0)*leftXStep;
		xe = xe + (iys - y0)*rightXStep;
	}
	else
	{
		iys = ceil(y0);
		xs = xs + (iys - y0)*leftXStep;
		xe = xe + (iys - y0)*rightXStep;
	}

	if (y1 > cliper.bottom) {
		iye = cliper.bottom;
	}
	else
	{
		iye = ceil(y1) - 1;
	}

	
	uint16 *buffer = (uint16*)drawBuffer + iys*lineLength;
	//test clip needed
	if ((x0 >= cliper.left&&x0 <= cliper.right)&&(x1 >= cliper.left&&x1 <= cliper.right)&&(x2 >= cliper.left&&x2 <= cliper.right)) {
		for (int i = iys;i <= iye;i++, buffer += lineLength) {
			ixs = xs;
			ixe = xe;
			memSet16(buffer + ixs, color, (ixe - ixs + 1));
			xs += leftXStep;
			xe += rightXStep;
		}
	}
	else {
		for (int i = iys; i <= iye; i++, buffer += lineLength)//step
		{
			ixs = xs;
			ixe = xe;
			//pre step for continue
			xs += leftXStep;
			xe += rightXStep;
			//trivially reject
			if (ixe < cliper.left || ixs>cliper.right) {
				continue;
			}

			//clip
			ixs = ixs < cliper.left ? cliper.left : ixs;
			ixe = ixe> cliper.right ? cliper.right : ixe;
			memSet16(buffer + ixs, color, (ixe - ixs + 1));

		}
	}

	return 1;
}






int Canvas::drawHLine16(int xl, int xr, int y, Color c)
{
	uint16 *buffer = (uint16 *)drawBuffer;
	int lineLength = pitch >> 1;
	int tmp;
	uint16 color = c.getRGB566();

	if (OUTBOUND(cliper.top, cliper.bottom, y)) {
		return -1;
	}

	if (xl > xr) {
		SWAP(xl, xr, tmp);
	}

	xl = MAX(cliper.left, xl);
	xr = MIN(xr, cliper.right);

	memSet16(buffer + xl + y*lineLength, color, xr - xl + 1);

	return 0;
}

int Canvas::drawVLine16(int y1, int y2, int x, Color c)
{	
	int tmp;
	uint16 color = c.getRGB566();
	//clip
	if (OUTBOUND(cliper.left, cliper.right, x)) {
		return -1;
	}
	if (y1 > y2) {
		SWAP(y1, y2, tmp);
	}
	y1 = MAX(cliper.top, y1);
	y2 = MIN(cliper.bottom, y2);


	//draw
	int lineLength = pitch >> 1;
	uint16 *buffer = (uint16 *)drawBuffer + x + y1*lineLength;

	for (int i = 0;i <= y2 - y1;i++) {
		*buffer = color;
		buffer += lineLength;
	}
	return 0;
}

int Canvas::drawLine16(int x1, int y1, int x2, int y2, Color c)
{
	uint16 color = c.getRGB566();
	if (!clipLine(x1, y1, x2, y2)) {
		return 0;
	}

	int dx, dy, dxX2, dyX2, xStep, yStep, stepIndicator;
	int lineLength = pitch >> 1;
	
	dx = x2 - x1;
	dy = y2 - y1;
	if (dx > 0) {
		xStep = 1;
	}
	else
	{
		xStep = -1;
		dx = -dx;
	}
	if (dy > 0) {
		yStep = lineLength;

	}
	else
	{
		yStep = -lineLength;
		dy = -dy;
	}

	dxX2 = dx << 1;
	dyX2 = dy << 1;

	uint16 *bufferPos = (uint16*)drawBuffer + x1 + y1*lineLength;


	if (dx > dy) {
		stepIndicator = dyX2 - dx;
		for (int i = 0;i <= dx;i++) {
			*bufferPos = color;
			if (stepIndicator > 0) {
				stepIndicator -= dxX2;
				bufferPos += yStep;
			}
			bufferPos += xStep;
			stepIndicator += dyX2;
		}
	}
	else
	{
		stepIndicator = dxX2 - dy;
		for (int i = 0;i <= dy;i++) {
			*bufferPos = color;
			if (stepIndicator > 0) {
				stepIndicator -= dyX2;
				bufferPos += xStep;
			}
			bufferPos += yStep;
			stepIndicator += dxX2;
		}
	}


	return 1;
}

int Canvas::drawTriangleWireI16(int x1, int y1, int x2, int y2, int x3, int y3, Color c)
{
	drawLine16(x1, y1, x2, y2, c);
	drawLine16(x2, y2, x3, y3, c);
	drawLine16(x3, y3, x1, y1, c);
	return 0;
}

int Canvas::drawTriangleSolidI16(int x1, int y1, int x2, int y2, int x3, int y3, Color c)
{	
	int tmp;
	//line so return 0 as trivially reject
	if ((y1 == y2&&y2 == y3) || (x1 == x2&&x2 == x3)) {
		return 0;
	}

	//select sort for y small -> big
	if (y2 < y1) {
		SWAP(y2, y1, tmp);
		SWAP(x2, x1, tmp);
	}
	if (y3 < y1) {
		SWAP(y3, y1, tmp);
		SWAP(x3, x1, tmp);
	}
	if (y3 < y2) {
		SWAP(y2, y3, tmp);
		SWAP(x2, x3, tmp);
	}

	//trivially reject
	if (y1 > cliper.bottom || y3<cliper.top || (x1<cliper.left&&x2<cliper.left&&x3<cliper.left) || (x1>cliper.right&&x2>cliper.right&&x3>cliper.right)) {
		return 0;
	}

	if (y1 == y2) {
		drawTopTriangleI16(x3, y3, x1, y1, x2, y2, c);
	}
	else if (y2 == y3)
	{
		drawBottomTriangleI16(x1, y1, x2, y2, x3, y3, c);
	}
	else
	{
		int newX;
		newX = x1 + 0.5f + (y2 - y1)*(x3 - x1) / (float)(y3 - y1);
		drawBottomTriangleI16(x1, y1, x2, y2, newX, y2, c);
		drawTopTriangleI16(x3, y3, x2, y2, newX, y2, c);
	}

	return 1;
}

int Canvas::drawTopTriangleI16(int x1, int y1, int x2, int y2, int x3, int y3, Color c)
{
	float leftXStep, rightXStep, xs, xe, dy;
	uint16 color = c.getRGB566();
	int tmp;

	int lineLength = pitch >> 1;

	if (x2 > x3) {
		SWAP(x2, x3, tmp);
	}
	dy = y1 - y2;
	leftXStep = (x1 - x2) / dy;
	rightXStep = (x1 - x3) / dy;
	xs = x2;
	xe = x3 + 0.5f;//for after rounding

	if (y2 < cliper.top) {
		xs = xs + leftXStep*(cliper.top - y2);
		xe = xe + rightXStep*(cliper.top - y2);
		y2 = y3 = cliper.top;
	}

	if (y1 > cliper.bottom) {
		y1 = cliper.bottom;
	}

	uint16 *buffer = (uint16*)drawBuffer + y2*lineLength;
	//test clip needed
	if (x1 >= cliper.left&&x1 <= cliper.right&&xs >= cliper.left&&xs <= cliper.right&&xe >= cliper.left&&xe <= cliper.right) {
		for (int i = y2;i <= y1;i++) {
			memSet16(buffer + (unsigned int)xs, color, (xe - xs + 1));
			xs += leftXStep;
			xe += rightXStep;
			buffer += lineLength;
		}
	}
	else {
		int left, right;

		for (int i = y2; i <= y1; i++, buffer += lineLength)//step
		{
			left = xs;
			right = xe;
			//pre step for continue
			xs += leftXStep;
			xe += rightXStep;
			//trivially reject
			if (right < cliper.left || left>cliper.right) {
				continue;
			}

			//clip
			left = left < cliper.left ? cliper.left : left;
			right = right > cliper.right ? cliper.right : right;
			memSet16(buffer + left, color, (right - left + 1));

		}
	}


	return 1;
}

int Canvas::drawBottomTriangleI16(int x1, int y1, int x2, int y2, int x3, int y3, Color c)
{
	float leftXStep, rightXStep, xs, xe, dy;
	uint16 color = c.getRGB566();
	int tmp;

	int lineLength = pitch >> 1;

	if (x2 > x3) {
		SWAP(x2, x3, tmp);
	}
	dy = y1 - y2;
	leftXStep = (x2 - x1) / dy;
	rightXStep = (x3 - x1) / dy;
	xs = x2;
	xe = x3 + 0.5f;//for after rounding

	if (y2 > cliper.bottom) {
		xs = xs + leftXStep*(y2 - cliper.bottom);
		xe = xe + rightXStep*(y2 - cliper.bottom);
		y2 = y3 = cliper.bottom;
	}

	if (y1 < cliper.top) {
		y1 = cliper.top;
	}

	uint16 *buffer = (uint16*)drawBuffer + y2*lineLength;
	//test clip needed
	if (x1 >= cliper.left&&x1 <= cliper.right&&xs >= cliper.left&&xs <= cliper.right&&xe >= cliper.left&&xe <= cliper.right) {
		for (int i = y2;i >= y1;i--) {
			memSet16(buffer + (unsigned int)xs, color, (xe - xs + 1));
			xs += leftXStep;
			xe += rightXStep;
			buffer -= lineLength;
		}
	}
	else {
		int left, right;

		for (int i = y2; i >= y1; i--,buffer -= lineLength)//step
		{
			left = xs;
			right = xe;
			//pre step for continue
			xs += leftXStep;
			xe += rightXStep;
			//trivially reject
			if (right < cliper.left || left>cliper.right) {
				
				continue;
			}

			//clip
			left = left < cliper.left ? cliper.left : left;
			right = right > cliper.right ? cliper.right : right;
			memSet16(buffer + left, color, (right - left + 1));

			
		}
	}


	return 1;
}

int Canvas::drawBitmap(float left, float top, float right, float bottom, Bitmap * bitmap, int sampleMode)
{	//trivially reject
	if (FEQUALE1(left, right) || FEQUALE1(top, bottom)) {
		return 0;
	}
	if (right<cliper.left || left>cliper.right || top > cliper.bottom || bottom < cliper.top) {
		return 0;
	}

	//accept
	if (bitmap->getBitCount() == 16) {
		switch (sampleMode)
		{
		case CANVAS_SAMPLE_POINT:
		default:
			drawBitmapPoint16(left, top, right, bottom, (uint16 *)bitmap->getBuffer(), bitmap->getWidth(), bitmap->getHeight(), bitmap->getPitch());
			break;
		}
	}
	else
	{
		return 0;
	}

	return 1;
}

int Canvas::drawBitmapPoint16(float left, float top, float right, float bottom, uint16 * buffer, int width, int height, int pitch)
{	

	fixp16 us,vs,fpxu,fpxv,dvdy,dudx;
	int ixs, ixe, iys, iye;
	uint16 *srcBuffer, *destBuffer;
	int srcLineLength,destLineLength;
	int u, v;
	int bitmapR, bitmapB;

	bitmapR = width - 1;
	bitmapB = height - 1;

	srcLineLength = pitch >> 1;
	destLineLength = this->pitch >> 1;

	//init 

	iys = top + 0.5f;
	iye = bottom+0.5f;
	ixs = left + 0.5f;
	ixe = right + 0.5f;

	us = 0;
	vs = 0;
	
	
	dudx = FIXP16_FLOAT_TO((1.0f) / (ixe - ixs));
	dvdy = FIXP16_FLOAT_TO((1.0f) / (iye - iys));

	//clip

	if (iys < cliper.top) {
		vs = (cliper.top - iys)*dvdy;
		iys = cliper.top;
	}
	
	if (iye > cliper.bottom) {
		iye = cliper.bottom;
	}
	
	if (ixs < cliper.left) {
		us = (cliper.left - ixs)*dudx;
		ixs = cliper.left;
	}
	if (ixe > cliper.right) {
		ixe = cliper.right;
	}
	
	//draw
	destBuffer = (uint16 *)this->drawBuffer + destLineLength*iys;
	srcBuffer = buffer;
    fpxv = vs;
	for (int yidx = iys;yidx <= iye;yidx++ ){
		fpxu = us;
		v = (FIXP16_TO_INT(FIXP16_ROUNDUP(fpxv*bitmapR)))*srcLineLength;
		for (int xidx = ixs;xidx <= ixe;xidx++) {
			u = FIXP16_TO_INT(FIXP16_ROUNDUP(fpxu*bitmapB));
			destBuffer[xidx] = srcBuffer[v + u];
			fpxu += dudx;
			
		}
		destBuffer += destLineLength;
		fpxv +=dvdy;
	}

	return 1;
}






int Canvas::setClipper(int left, int top, int right, int bottom)
{
	left = MAX(0, left);
	top = MAX(0, top);
	right = MIN(width - 1, right);
	bottom = MIN(height - 1, bottom);
	cliper = { left,top, right,bottom };
	return 0;
}

int Canvas::clearClipper()
{
	cliper = { 0,0,width - 1,height - 1 };
	return 0;
}

//private 



int Canvas::clipLine(int & x1, int & y1, int & x2, int & y2)
{

	//define point location
#define POINT_LOC_C 0x0 //0000
#define POINT_LOC_N 0x1 //0001
#define POINT_LOC_S 0x2 //0010
#define POINT_LOC_W 0x4 //0100
#define POINT_LOC_E 0x8 //1000
#define POINT_LOC_NW 0x5 //0101
#define POINT_LOC_NE 0x9  //1001
#define POINT_LOC_SW 0x6 //0110
#define POINT_LOC_SE 0xa //1010


	int xt1 = x1, yt1 = y1, xt2 = x2, yt2 = y2, dx, dy;
	int pointLoc1 = 0;
	int pointLoc2 = 0;

	//check 

	if (x1 < cliper.left) {
		pointLoc1 |= POINT_LOC_W;
	}
	else if (x1>cliper.right)
	{
		pointLoc1 |= POINT_LOC_E;
	}
	if (y1<cliper.top) {
		pointLoc1 |= POINT_LOC_N;
	}
	else if (y1>cliper.bottom)
	{
		pointLoc1 |= POINT_LOC_S;
	}

	if (x2 < cliper.left) {
		pointLoc2 |= POINT_LOC_W;
	}
	else if (x2>cliper.right)
	{
		pointLoc2 |= POINT_LOC_E;
	}
	if (y2<cliper.top) {
		pointLoc2 |= POINT_LOC_N;
	}
	else if (y2>cliper.bottom)
	{
		pointLoc2 |= POINT_LOC_S;
	}

	//trivially reject
	if (pointLoc1&pointLoc2) {
		return 0;
	}
	//trivially accept
	if (!(pointLoc1 | pointLoc2)) {
		return 1;
	}


	//do clip

	dx = x2 - x1;
	dy = y2 - y1;


	switch (pointLoc1)
	{
	case POINT_LOC_N:
		yt1 = cliper.top;
		xt1 = x1 + 0.5f + ((cliper.top - y1)*dx) / dy;
		break;
	case POINT_LOC_S:
		yt1 = cliper.bottom;
		xt1 = x1 + 0.5f + ((cliper.bottom - y1)*dx) / dy;
		break;
	case POINT_LOC_W:
		xt1 = cliper.left;
		yt1 = y1 + 0.5f + ((cliper.left - x1)*dy) / dx;
		break;
	case POINT_LOC_E:
		xt1 = cliper.right;
		yt1 = y1 + 0.5f + ((cliper.right - x1)*dy) / dx;
		break;
	case POINT_LOC_NW:
		yt1 = cliper.top;
		xt1 = x1 + 0.5f + ((cliper.top - y1)*dx) / dy;
		if (OUTBOUND(cliper.left, cliper.right, xt1)) {
			xt1 = cliper.left;
			yt1 = y1 + 0.5f + ((cliper.left - x1)*dy) / dx;
		}
		break;
	case POINT_LOC_NE:
		yt1 = cliper.top;
		xt1 = x1 + 0.5f + ((cliper.top - y1)*dx) / dy;
		if (OUTBOUND(cliper.left, cliper.right, xt1)) {
			xt1 = cliper.right;
			yt1 = y1 + 0.5f + ((cliper.right - x1)*dy) / dx;
		}
	case POINT_LOC_SW:
		yt1 = cliper.bottom;
		xt1 = x1 + 0.5f + ((cliper.bottom - y1)*dx) / dy;
		if (OUTBOUND(cliper.left, cliper.right, xt1)) {
			xt1 = cliper.left;
			yt1 = y1 + 0.5f + ((cliper.left - x1)*dy) / dx;
		}
		break;
	case POINT_LOC_SE:
		yt1 = cliper.bottom;
		xt1 = x1 + 0.5f + ((cliper.bottom - y1)*dx) / dy;
		if (OUTBOUND(cliper.left, cliper.right, xt1)) {
			xt1 = cliper.right;
			yt1 = y1 + 0.5f + ((cliper.right - x1)*dy) / dx;
		}
		break;

	default:
		break;
	}

	dx = x1 - x2;
	dy = y1 - y2;


	switch (pointLoc2)
	{
	case POINT_LOC_N:
		yt2 = cliper.top;
		xt2 = x2 + 0.5f + ((cliper.top - y2)*dx) / dy;
		break;
	case POINT_LOC_S:
		yt2 = cliper.bottom;
		xt2 = x2 + 0.5f + ((cliper.bottom - y2)*dx) / dy;
		break;
	case POINT_LOC_W:
		xt2 = cliper.left;
		yt2 = y2 + 0.5f + ((cliper.left - x2)*dy) / dx;
		break;
	case POINT_LOC_E:
		xt2 = cliper.right;
		yt2 = y2 + 0.5f + ((cliper.right - x2)*dy) / dx;
		break;
	case POINT_LOC_NW:
		yt2 = cliper.top;
		xt2 = x2 + 0.5f + ((cliper.top - y2)*dx) / dy;
		if (OUTBOUND(cliper.left, cliper.right, xt2)) {
			xt2 = cliper.left;
			yt2 = y2 + 0.5f + ((cliper.left - x2)*dy) / dx;
		}
		break;
	case POINT_LOC_NE:
		yt2 = cliper.top;
		xt2 = x2 + 0.5f + ((cliper.top - y2)*dx) / dy;
		if (OUTBOUND(cliper.left, cliper.right, xt2)) {
			xt2 = cliper.right;
			yt2 = y2 + 0.5f + ((cliper.right - x2)*dy) / dx;
		}
	case POINT_LOC_SW:
		yt2 = cliper.bottom;
		xt2 = x2 + 0.5f + ((cliper.bottom - y2)*dx) / dy;
		if (OUTBOUND(cliper.left, cliper.right, xt2)) {
			xt2 = cliper.left;
			yt2 = y2 + 0.5f + ((cliper.left - x2)*dy) / dx;
		}
		break;
	case POINT_LOC_SE:
		yt2 = cliper.bottom;
		xt2 = x2 + 0.5f + ((cliper.bottom - y2)*dx) / dy;
		if (OUTBOUND(cliper.left, cliper.right, xt2)) {
			xt2 = cliper.right;
			yt2 = y2 + 0.5f + ((cliper.right - x2)*dy) / dx;
		}
		break;

	default:
		break;
	}

	if (OUTBOUND(cliper.left, cliper.right, xt1) || OUTBOUND(cliper.left, cliper.right, xt2) || OUTBOUND(cliper.top, cliper.bottom, yt1) || OUTBOUND(cliper.top, cliper.bottom, yt2)) {
		return 0;
	}
	x1 = xt1;
	y1 = yt1;
	x2 = xt2;
	y2 = yt2;


	return 1;
}