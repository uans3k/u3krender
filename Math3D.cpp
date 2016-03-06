#include "Math3D.h"
#include <math.h>


int logbase2[513] =
{
	0,0,1,1,2,2,2,2,3,3,3,3,3,3,3,3, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5, 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
	6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6, 6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
	6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6, 6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
	7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7, 7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
	7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7, 7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
	7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7, 7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
	7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7, 7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,

	8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8, 8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
	8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8, 8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
	8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8, 8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
	8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8, 8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
	8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8, 8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
	8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8, 8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
	8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8, 8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
	8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8, 8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,

};
float Cos_Look[361];
float Sin_Look[361];


int Log_Base2(int x)
{
	if (x > 512) {
		return 0;
	}

	return logbase2[x];
}

void SinCosTables_Build()
{

	for (int ang = 0; ang <= 360; ang++)
	{
		float theta = (float)ang*PI / (float)180;
		Cos_Look[ang] = cos(theta);
		Sin_Look[ang] = sin(theta);

	}
}
float Sin_Fast(float theta)
{
	int floor;
	theta = fmodf(theta, 360);
	if (theta < 0) {
		theta += 360;
	}
	floor = theta;

	return Sin_Look[floor] + (theta - floor)*(Sin_Look[floor + 1] - Sin_Look[floor]);
}

float Cos_Fast(float theta)
{
	int floor;
	theta = fmodf(theta, 360);
	if (theta < 0) {
		theta += 360;
	}
	floor = theta;

	return Cos_Look[floor] + (theta - floor)*(Cos_Look[floor + 1] - Cos_Look[floor]);
}

float Distance2D_Fast(int x, int y)
{
	x = abs(x);
	y = abs(y);

	// compute the minimum of x,y
	int mn = MIN(x, y);

	// return the distance
	return(x + y - (mn >> 1) - (mn >> 2) + (mn >> 4));
}

float Distance3D_Fast(float fx, float fy, float fz)
{
	int temp;  // used for swaping
	int x, y, z; // used for algorithm

				 // make sure values are all positive
	x = fabs(fx) * 1024;
	y = fabs(fy) * 1024;
	z = fabs(fz) * 1024;

	// sort values
	if (y < x) SWAP(x, y, temp)

		if (z < y) SWAP(y, z, temp)

			if (y < x) SWAP(x, y, temp)

				int dist = (z + 11 * (y >> 5) + (x >> 2));

	// compute distance with 8% error
	return((float)(dist >> 10));
}

void Vector2D_Build(Point2D * p0, Point2D * p1, Vector2D * vResult)
{
	vResult->x = p1->x - p0->x;
	vResult->y = p1->y - p0->y;
}

void Vector2D_Add(Vector2D * v1, Vector2D * v2, Vector2D * vSum)
{
	vSum->x = v1->x + v2->x;
	vSum->y = v1->y + v2->y;
}

void Vector2D_Sub(Vector2D * v1, Vector2D * v2, Vector2D * vDiff)
{
	vDiff->x = v1->x - v2->x;
	vDiff->y = v1->y - v2->y;
}

void Vector2D_Scalar(float scalar, Vector2D * v, Vector2D * vResult)
{
	vResult->x = v->x*scalar;
	vResult->y = v->y*scalar;
}

void Vector2D_Scalar(float scalar, Vector2D * vSelf)
{
	vSelf->x *= scalar;
	vSelf->y *= scalar;
}

float Vector2D_Dot(Vector2D * v1, Vector2D * v2)
{
	return v1->x*v2->x + v1->y*v2->y;
}

float Vector2D_Length(Vector2D * v1)
{
	return sqrtf(v1->x*v1->x + v1->y*v1->y);
}

float Vector2D_Length_Fast(Vector2D * v1)
{
	return Distance2D_Fast(v1->x, v1->y);
}

void Vector2D_Normalize(Vector2D * vSelf)
{
	float length = sqrtf(vSelf->x*vSelf->x + vSelf->y*vSelf->y);

	if (length < EPSILON_E5)
	{
		return;
	}
	float divLength = 1 / length;

	vSelf->x *= divLength;
	vSelf->y *= divLength;
}

void Vector2D_Normalize(Vector2D * v, Vector2D * vResult)
{
	float length = sqrtf(v->x*v->x + v->y*v->y);

	if (length < EPSILON_E5)
	{
		Vector2D_Zero(v);
		return;
	}
	float divLength = 1 / length;

	vResult->x = v->x * divLength;
	vResult->y = v->y * divLength;
}

float Vector2D_CosTh(Vector2D * v1, Vector2D * v2)
{
	return Vector2D_Dot(v1, v2) / (Vector2D_Length(v1)*Vector2D_Length(v2));
}

void Vector3D_Build(Point3D * p0, Point3D * p1, Vector3D * vResult)
{
	vResult->x = p1->x - p0->x;
	vResult->y = p1->y - p0->y;
	vResult->z = p1->z - p0->z;
}

void Vector3D_Add(Vector3D * v1, Vector3D * v2, Vector3D * vSum)
{

	vSum->x = v1->x + v2->x;
	vSum->y = v1->y + v2->y;
	vSum->z = v1->z + v2->z;
}

void Vector3D_Sub(Vector3D * v1, Vector3D * v2, Vector3D * vDiff)
{
	vDiff->x = v1->x - v2->x;
	vDiff->y = v1->y - v2->y;
	vDiff->z = v1->z - v2->z;
}

void Vector3D_Mul_Scalar(float scalar, Vector3D * v, Vector3D * vResult)
{
	vResult->x = v->x*scalar;
	vResult->y = v->y*scalar;
	vResult->z = v->z*scalar;
}

void Vector3D_Mul_Scalar(float scalar, Vector3D * vSelf)
{
	vSelf->x *= scalar;
	vSelf->y *= scalar;
	vSelf->z *= scalar;
}

float Vector3D_Dot(Vector3D * v1, Vector3D * v2)
{
	return v1->x*v2->x + v1->y*v2->y + v1->z*v2->z;
}

void Vector3D_Cross(Vector3D * v1, Vector3D * v2, Vector3D * vResult)
{
	vResult->x = v1->y*v2->z - v1->z*v2->y;
	vResult->y = v1->z*v2->x - v1->x*v2->z;
	vResult->z = v1->x*v2->y - v1->y*v2->x;
	
}

float Vector3D_Length(Vector3D * v1)
{
	return sqrtf(v1->x*v1->x + v1->y*v1->y + v1->z*v1->z);
}

float Vector3D_Length_Fast(Vector3D * v1)
{
	return Distance3D_Fast(v1->x, v1->y, v1->z);
}

void Vector3D_Normalize(Vector3D * vSelf)
{
	float length = sqrtf(vSelf->x*vSelf->x + vSelf->y*vSelf->y + vSelf->z*vSelf->z);

	if (length < EPSILON_E5)
	{
		return;
	}
	float divLength = 1 / length;

	vSelf->x *= divLength;
	vSelf->y *= divLength;
	vSelf->z *= divLength;
}

void Vector3D_Normalize(Vector3D * v, Vector3D * vResult)
{
	float length = sqrtf(v->x*v->x + v->y*v->y + v->z*v->z);

	if (length < EPSILON_E5)
	{
		Vector3D_Zero(vResult);
		return;
	}
	float divLength = 1 / length;

	vResult->x = v->x * divLength;
	vResult->y = v->y * divLength;
	vResult->z = v->z*divLength;
}

float Vector3D_CosTh(Vector3D * v1, Vector3D * v2)
{
	return Vector3D_Dot(v1, v2) / (Vector3D_Length(v1)*Vector3D_Length(v2));
}

void Vector4D_Build(Point4D * p0, Point4D * p1, Vector4D * vResult)
{

	vResult->x = p1->x - p0->x;
	vResult->y = p1->y - p0->y;
	vResult->z = p1->z - p0->z;
	vResult->w = 1;
}

void Vector4D_Add(Vector4D * pSelf, Vector4D * p1)
{
	pSelf->x += p1->x;
	pSelf->y += p1->y;
	pSelf->z += p1->z;
}


void Vector4D_Add(Vector4D * p0, Vector4D * p1, Vector4D * pResult)
{
	pResult->x = p0->x + p1->x;
	pResult->y = p0->x + p1->y;
	pResult->z = p0->x + p1->z;
	pResult->w = 0.0f;
}

void Vector4D_Sub(Vector4D * v0, Vector4D * v1, Vector4D * vResult)
{
	vResult->x = v0->x - v1->x;
	vResult->y = v0->y - v1->y;
	vResult->z = v0->z - v1->z;
	vResult->w = 0.0f;
}

void Vector4D_Sub(Vector4D * vSelf, Vector4D * v)
{
	vSelf->x -= v->x;
	vSelf->y -= v->y;
	vSelf->z -= v->z;
	vSelf->w = 0.0f;
}

void Vector4D_Mul_Scalar(float scalar, Vector4D * vSelf)
{
	vSelf->x *= scalar;
	vSelf->y *= scalar;
	vSelf->z *= scalar;
	vSelf->w= 0.0f;
}

void Vector4D_Mul_Scalar(float scalar, Vector4D * v, Vector4D * vResult)
{
	vResult->x = v->x*scalar;
	vResult->y = v->y*scalar;
	vResult->z = v->z*scalar;
	vResult->w = 0.0f;
}

float Vector4D_Length(Vector4D * v0)
{
	return sqrtf(v0->x*v0->x + v0->y*v0->y + v0->z*v0->z);
}

float Vector4D_Length_Fast(Vector4D * v0)
{
	return Distance3D_Fast(v0->x, v0->y, v0->z);
}

void Vector4D_Normalize(Vector4D * v, Vector4D * vResult)
{
	float length = sqrtf(vResult->x*vResult->x + vResult->y*vResult->y + vResult->z*vResult->z);

	if (length < EPSILON_E5)
	{
		return;
	}
	float divLength = 1 / length;

	vResult->x *= divLength;
	vResult->y *= divLength;
	vResult->z *= divLength;
	vResult->w = 0.0f;
}

void Vector4D_Normalize(Vector4D * vSelf)
{
	float length = sqrtf(vSelf->x*vSelf->x + vSelf->y*vSelf->y + vSelf->z*vSelf->z);

	if (length < EPSILON_E5)
	{
		return;
	}
	float divLength = 1 / length;

	vSelf->x *= divLength;
	vSelf->y *= divLength;
	vSelf->z *= divLength;
	vSelf->w = 0.0f;
}

float Vector4D_Dot(Vector4D * v0, Vector4D * v1)
{
	return v0->x*v1->x + v0->y*v1->y + v0->z*v1->z;
}

void Vector4D_Cross(Vector4D * v0, Vector4D * v1, Vector4D * vResult)
{
	vResult->x = v0->y*v1->z - v0->z*v1->y;
	vResult->y = v0->z*v1->x - v0->x*v1->z;
	vResult->z = v0->x*v1->y - v0->y*v1->x;
	vResult->w = 0.0f;
}

void Point4D_Move(Point4D * p, Vector4D * v)
{
	p->x += v->x;
	p->y += v->y;
	p->z += v->z;
	p->w = 1.0f;
}

void Point4D_Move(Point4D * src, Point4D * dest, Vector4D * v)
{
	dest->x =src->x+ v->x;
	dest->y = src->y+ v->y;
	dest->z = src->z+ v->z;
	dest->w = 1.0f;
}

void Point4D_Move(Point4D * p, Vector3D * v)
{
	p->x += v->x;
	p->y += v->y;
	p->z += v->z;
	p->w = 1.0f;
}

void Matrix_Mul_2x2(Matrix2x2 *m1, Matrix2x2 *m2, Matrix2x2 *mResult)
{
	mResult->M00 = m1->M00*m2->M00 + m1->M01*m2->M10;
	mResult->M01 = m1->M00*m2->M01 + m1->M01*m2->M11;

	mResult->M10 = m1->M10*m2->M00 + m1->M11*m2->M10;
	mResult->M11 = m1->M10*m2->M01 + m1->M11*m2->M11;

}

void Matrix_Add_4x4(Matrix4x4 * m1, Matrix4x4 * m2, Matrix4x4 * mResult)
{

	for (int row = 0; row < 4; row++)
	{
		for (int col = 0; col < 4; col++)
		{

			mResult->M[row][col] = m1->M[row][col] + m2->M[row][col];
		}

	}
}

/*
  M3_ij= M1_ik*M2_kj  (einstein)
*/
void Matrix_Mul_4x4(Matrix4x4 * m1, Matrix4x4 * m2, Matrix4x4 * mResult)
{
	for (int row = 0; row < 4; row++)
	{
		for (int col = 0; col < 4; col++)
		{

			float sum = 0;

			for (int index = 0; index < 4; index++)
			{

				sum += (m1->M[row][index] * m2->M[index][col]);
			}


			mResult->M[row][col] = sum;

		}

	}
}

void Matrix_Mul_1x4_4x4(Matrix1x4 * m1, Matrix4x4 * m2, Matrix1x4 * mResult)
{
	for (int col = 0; col < 4; col++)
	{
		float sum = 0;

		for (int row = 0; row < 4; row++)
		{
			sum += (m1->M[row] * m2->M[row][col]);
		}


		mResult->M[col] = sum;

	}
}

void Matrix_Mul_Vector3D_4x4(Vector3D * v, Matrix4x4 * m, Vector3D * mResult)
{
	float sum;
	int row;
	for (int col = 0;col < 3;col++) {
		sum = 0.0f;
		row = 0;
		for (;row < 3;row++) {
			sum += v->M[row] * m->M[row][col];
		}
		sum += m->M[row][col];
		mResult->M[col] = sum;
	}
}

void Matrix_Mul_Vecotr3D_4x3(Vector3D * v, Matrix4x3 * m, Vector3D * mResult)
{
	float sum;
	int row;
	for (int col = 0;col < 3;col++) {
		sum = 0;
		row = 0;
		for (;row < 3;row++) {
			sum += v->M[row] * m->M[row][col];
		}
		sum += m->M[row][col];//1*m[3][col]
		mResult->M[col] = sum;
	}
}

void Matrix_Mul_Vector4D_4x4(Vector4D * v, Matrix4x4 * m, Vector4D * vResult)
{
	float sum;
	for (int col = 0;col < 4;col++) {
		sum = 0;
		for (int row = 0; row < 4; row++)
		{
			sum += v->M[row] * m->M[row][col];
		}
		vResult->M[col] = sum;
	}
}

void Matrix_Mul_Point4D_4x4(Point4D * v, Matrix4x4 * m, Point4D * vResult)
{	
	float sum;
	for (int col = 0;col < 4;col++) {
		sum = 0;
		for (int row = 0; row < 4; row++)
		{
			sum += v->M[row] * m->M[row][col];
		}
		vResult->M[col] = sum;
	}
}

void Matrix_Build_Rotation_4x4(Matrix4x4 *m, float thX, float thY, float thZ)
{	
	Matrix4x4 mRotate, mTmp[2];
	float sinTh, cosTh;

	int pos=0,posTmp;
	int state=0;
	MATRIX_IDENTITY_4X4(m);
	if (fabs(thX) > EPSILON_E5) {
		sinTh = Sin_Fast(thX);
		cosTh = Cos_Fast(thX);
		Matrix_Init_4x4(&mTmp[pos], 1, 0, 0, 0,//
			0, cosTh, sinTh, 0,//
			0, -sinTh, cosTh, 0,//
			0, 0, 0, 1
			);
		state = 1;
	}
	if (fabs(thY) > EPSILON_E5) {
		sinTh = Sin_Fast(thY);
		cosTh = Cos_Fast(thY);
		if (state) {
			Matrix_Init_4x4(&mRotate, cosTh, 0, -sinTh, 0,//
				0, 1, 0, 0,//
				sinTh, 0, cosTh, 0,//
				0, 0, 0, 1);
			posTmp = pos^1;
			Matrix_Mul_4x4(&mTmp[pos], &mRotate, &mTmp[posTmp]);
			pos = posTmp;
		}
		else
		{
			Matrix_Init_4x4(&mTmp[pos], cosTh, 0, -sinTh, 0,//
				0, 1, 0, 0,//
				sinTh, 0, cosTh, 0,//
				0, 0, 0, 1);
			state = 1;
		}
	}
	if (fabs(thZ) > EPSILON_E5) {
		sinTh = Sin_Fast(thZ);
		cosTh = Cos_Fast(thZ);
		if (state) {
			Matrix_Init_4x4(&mRotate ,cosTh, sinTh, 0, 0,//
				-sinTh, cosTh, 0, 0,//
				0, 0, 1, 0,//
				0, 0, 0, 1);
			posTmp = pos ^ 1;
			Matrix_Mul_4x4(&mTmp[pos], &mRotate, &mTmp[posTmp]);
			pos = posTmp;
		}
		else
		{
			Matrix_Init_4x4(&mTmp[pos],  cosTh,sinTh,0,0,//
				                    -sinTh,cosTh, 0, 0,//
				                    0, 0, 1, 0,//
			                     	0, 0, 0, 1);
			state = 1;
		}
	}
	if (state) {
		MATRIX_COPY_4X4(&mTmp[pos], m);
	}
}


void Matrix_Build_Translation_4x4(Matrix4x4 *m, float x, float y, float z)
{
	Matrix_Init_4x4(m,  1, 0, 0, 0,//
					    0, 1, 0, 0,//
						0, 0, 1, 0,//
						x, y, z, 1);
}

void Matrix_Mul_Vector4D_4x3(Vector4D * v, Matrix4x4 * m, Vector4D * vResult)
{
	float sum;
	for (int col = 0;col < 3;col++) {
		sum = 0;
		for (int row = 0; row < 4; row++)
		{
			sum += v->M[row] * m->M[row][col];
		}
		vResult->M[col] = sum;
	}
	vResult->M[3] = v->M[3];
}


int Matrix_Inverse_4x4(Matrix4x4 * m, Matrix4x4 * mResult)
{
	float det = (m->M00 * (m->M11 * m->M22 - m->M12 * m->M21) -
		m->M01 * (m->M10 * m->M22 - m->M12 * m->M20) +
		m->M02 * (m->M10 * m->M21 - m->M11 * m->M20));


	if (fabs(det) < EPSILON_E5)
		return 0;

	float divDet = 1 / det;


	mResult->M00 = divDet * (m->M11 * m->M22 - m->M12 * m->M21);
	mResult->M01 = divDet * (m->M02 * m->M21 - m->M01 * m->M22);
	mResult->M02 = divDet * (m->M01 * m->M12 - m->M02 * m->M11);
	mResult->M03 = 0.0f;
	mResult->M10 = divDet * (m->M12 * m->M20 - m->M10 * m->M22);
	mResult->M11 = divDet * (m->M00 * m->M22 - m->M02 * m->M20);
	mResult->M12 = divDet * (m->M02 * m->M10 - m->M00 * m->M12);
	mResult->M13 = 0.0f;

	mResult->M20 = divDet * (m->M10 * m->M21 - m->M11 * m->M20);
	mResult->M21 = divDet * (m->M01 * m->M20 - m->M00 * m->M21);
	mResult->M22 = divDet * (m->M00 * m->M11 - m->M01 * m->M10);
	mResult->M23 = 0.0f;

	mResult->M30 = -(m->M30 * mResult->M00 + m->M31 * mResult->M10 + m->M32 * mResult->M20);
	mResult->M31 = -(m->M30 * mResult->M01 + m->M31 * mResult->M11 + m->M32 * mResult->M21);
	mResult->M32 = -(m->M30 * mResult->M02 + m->M31 * mResult->M12 + m->M32 * mResult->M22);
	mResult->M33 = 1.0f;



	return 1;
}

void Matrix_Init_4x4(Matrix4x4 * m, float m00, float m01, float m02, float m03, float m10, float m11, float m12, float m13, float m20, float m21, float m22, float m23, float m30, float m31, float m32, float m33)
{
	m->M00 = m00; m->M01 = m01; m->M02 = m02; m->M03 = m03;
	m->M10 = m10; m->M11 = m11; m->M12 = m12; m->M13 = m13;
	m->M20 = m20; m->M21 = m21; m->M22 = m22; m->M23 = m23;
	m->M30 = m30; m->M31 = m31; m->M32 = m32; m->M33 = m33;

}

void Line2D_Init(Point2D * p0, Point2D * p1, Line2D * line)
{
	Point2D_Init(&line->p0, p0);
	Point2D_Init(&line->p1, p1);
	Vector2D_Build(p0, p1, &line->v);
}

void Line2D_GetFractionPoint(Line2D * line, float fraction, Point2D * point)
{
	point->x = line->p0.x + fraction*line->v.x;
	point->y = line->p0.y + fraction*line->v.y;
}

IntersectResult  Line2D_GetIntersectPoint(Line2D * line1, Line2D * line2, Point2D *pResult)
{	//    t1     t2
	// | l1.p.x  -l2.p.x |   |l2.p0.x-l1.p0.x|
	// | l1.p.y  -l2.p.y |   |l2.p0.y-l1.p0.y|

	float t1, t2, det, divDet;
	det = line2->v.x*line1->v.y - line1->v.x*line2->v.y;
	if (det < EPSILON_E5) {
		return INTERSECT_NONE;
	}
	divDet = 1 / det;
	t1 = divDet*((line1->p0.x - line2->p0.x)*line2->v.y - (line1->p0.y - line2->p0.y)*line2->v.x);
	t2 = divDet*(line1->v.x*(line2->p0.y - line1->p0.y) - line1->v.y*(line2->p0.x - line1->p0.x));

	pResult->x = line1->p0.x + line1->v.x*t1;
	pResult->y = line1->p0.y + line1->v.y*t1;

	if ((t1 >= 0) && (t1 <= 1) && (t2 >= 0) && (t2 <= 1))
		return INTERSECT_IN;
	else
		return INTERSECT_OUT;

}

IntersectResult Line2D_GetIntersectFraction(Line2D * line1, Line2D * line2, float * line1Fraction, float * line2Fraction)
{	// p1+v1t2=p2+v2t2;
	//     
	//     
	//    t1     t2
	// | l1.p.x  -l2.p.x |   |l2.p0.x-l1.p0.x|
	// | l1.p.y  -l2.p.y |   |l2.p0.y-l1.p0.y|

	float det, divDet;
	det = line2->v.x*line1->v.y - line1->v.x*line2->v.y;
	if (det < EPSILON_E5) {
		return INTERSECT_NONE;
	}
	divDet = 1 / det;
	*line1Fraction = divDet*((line1->p0.x - line2->p0.x)*line2->v.y - (line1->p0.y - line2->p0.y)*line2->v.x);
	*line2Fraction = divDet*(line1->v.x*(line2->p0.y - line1->p0.y) - line1->v.y*(line2->p0.x - line1->p0.x));

	if ((*line1Fraction >= 0) && (*line1Fraction <= 1) && (*line1Fraction >= 0) && (*line2Fraction <= 1))
		return INTERSECT_IN;
	else
		return INTERSECT_OUT;
}

void Line3D_GetFractionPoint(Line3D * line, float fraction, Point3D * point)
{
	point->x = line->p0.x + fraction*line->v.x;
	point->y = line->p0.y + fraction*line->v.y;
	point->z = line->p0.z + fraction*line->v.z;

}

void Plane3D_Init(Plane3D * plane, Point3D * p, Vector3D * n, int normalize)
{
	Point3D_Init(&plane->p0, p);
	if (normalize) {
		Vector3D_Normalize(n, &plane->p0);
	}
	else
	{
		Vector3D_Init(&plane->n, n);
	}
}

int Plane3D_IsPointIn(Plane3D * plane, Point3D * p)
{
	float dot = plane->n.x*(p->x - plane->p0.x)//
		+ plane->n.y*(p->y - plane->p0.y)//
		+ plane->n.z*(p->z - plane->p0.z);
	if (dot < EPSILON_E5) {
		return 1;
	}
	return 0;
}

IntersectResult Plane3D_GetIntersectPoint(Plane3D * plane, Line3D * line, float *fraction, Point3D * pResult)
{
	// (p0+vt-s.p0)*s.n=0;
	// vx*nx+vy*ny+vz*nz!=0 ==>
	//t={(s.x0-p.x0)*nx+(s.y0-p.y0)*ny+(s.z0+p.z0)*nz}/(vx*nx+vy*ny+vz*nz);
	//
	// vx*nx+vy*ny+vz*nz=0==>
	// p0 is in plane ? 
	// 
	float dot = Vector3D_Dot(&plane->n, &line->v);
	if (dot < EPSILON_E5) {
		if (Plane3D_IsPointIn(plane, &line->p0)) {
			return INTERSECT_ALLIN;
		}
		else
		{
			return INTERSECT_NONE;
		}
	}
	float divDot = 1 / dot;
	*fraction = divDot*((plane->p0.x - line->p0.x)*plane->n.x + (plane->p0.y - line->p0.y)*plane->n.y + (plane->p0.z - line->p0.z)*plane->n.z);
	pResult->x = line->p0.x + (*fraction)*line->v.x;
	pResult->y = line->p0.y + (*fraction)*line->v.y;
	pResult->z = line->p0.z + (*fraction)*line->v.z;
	if (((*fraction) >= 0) && ((*fraction) <= 1)) {
		return INTERSECT_IN;
	}
	else {
		return INTERSECT_OUT;
	}
}

void Line3D_Init(Point3D * p0, Point3D * p1, Line3D * line)
{
	Point3D_Init(&line->p0, p0);
	Point3D_Init(&line->p1, p1);
	Vector3D_Build(p0, p1, &line->v);
}
