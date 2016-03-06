#pragma once
#include <memory.h>

//use anonymous union and anonymous struct  for struct to let sturct flexible
typedef struct {
	union
	{
		float M[4][4];
		struct
		{
			float M00, M01, M02, M03;
			float M10, M11, M12, M13;
			float M20, M21, M22, M23;
			float M30, M31, M32, M33;
		};
	};
}Matrix4x4;



typedef struct {
	union
	{
		float M[4][3];
		struct
		{
			float M00, M01, M02;
			float M10, M11, M12;
			float M20, M21, M22;
			float M30, M31, M32;

		};
	};
}Matrix4x3;

typedef struct {
	union
	{
		float M[4];
		struct
		{
			float M00;
			float M10;
			float M20;
			float M30;

		};
	};
}Matrix4x1;

typedef struct {
	union
	{
		float M[4];
		struct
		{
			float M00;
			float M01;
			float M02;
			float M03;

		};
	};
}Matrix1x4;

typedef struct {
	union
	{
		float M[3][3];
		struct
		{
			float M00, M01, M02;
			float M10, M11, M12;
			float M20, M21, M22;

		};
	};
}Matrix3x3;

typedef struct {
	union
	{
		float M[2][2];
		struct
		{
			float M00, M01;
			float M10, M11;
		};
	};
}Matrix2x2;



typedef struct {
	union
	{
		float M[2];
		struct
		{
			float x, y;
		};
	};
}Vector2D, Point2D;

typedef struct {
	union
	{
		float M[3];
		struct
		{
			float x, y, z;
		};
	};
}Vector3D, Point3D;

typedef struct {
	union
	{
		float M[4];
		struct
		{
			float x, y, z, w;
		};
	};
}Vector4D, Point4D;


typedef struct {
	Point3D p0;
	Point3D p1;
	Vector3D v;//needed to know p0->p1 or p0<-p1 
}Line3D;

typedef struct {
	Point2D p0;
	Point2D p1;
	Vector2D v;//needed to know p0->p1 or p0<-p1 
}Line2D;

typedef struct {
	Point3D p0;
	Vector3D n;// normal to the plane
}Plane3D;

typedef struct {
	Point2D p0;
	Vector2D n;// normal to the plane
}Plane2D;

const Matrix4x4 E_4x4 = { 1,0,0,0,
						  0,1,0,0,
						  0,0,1,0,
						  0,0,0,1 };
const Matrix4x3 E_4x3 = { 1,0,0,
						  0,1,0,
						  0,0,1,
						  0,0,0, };

const Matrix3x3 E_3x3 = { 1,0,0,
						  0,1,0,
						  0,0,1 };

const Matrix2x2 E_2x2 = { 1,0,
						  0,1 };


/*
	macros
*/


#define PI (3.141596f)
#define DEG_TO_RAD(degree) ((degree)*PI/180.0)
#define RAD_TO_DEG(rad) ((rad)*180.0/PI)

#define EPSILON_E1 (0.1f) 
#define EPSILON_E2 (0.01f) 
#define EPSILON_E3 (0.001f) 
#define EPSILON_E4 (0.0001f) 
#define EPSILON_E5 (0.00001f)
#define EPSILON_E6 (0.000001f)

#define FIXP16_SHIFT     16
#define FIXP16_MAG       65536
#define FIXP16_F_MASK   0x0000ffff
#define FIXP16_I_MASK   0xffff0000
#define FIXP16_ROUNDUP_ADD  0x00008000




#define FIXP30_SHIFT 30   
#define FIXP30_MAG  1073741824 
#define FIXP30_F_MASK 0x3fffffff
#define FIXP30_I_MASK 0xc0000000

#define FIXP29_SHIFT 29
#define FIXP29_MAG 536870912
#define FIXP29_F_MASK 0x1fffffff
#define FIXP29_I_MASK 0xe0000000

#define FIXP28_SHIFT 28
#define FIXP28_MAG 268435456 

#define FIXP22_SHIFT 22
#define FIXP22_MAG  4194304

#define FIXP20_SHIFT 20
#define FIXP20_MAG 1048576

#define FIXP19_SHIFT 19
#define FIXP19_MAG 524288 
#define FIXP19_ROUNDUP_ADD 0x00040000

typedef int fixp16;// 1.15.16
typedef int fixp19;//1.12.19 4k
typedef int fixp20;// 1.11.20
typedef int fixp22;
typedef int fixp28; //1.3.28
typedef int fixp29;//1.2.29  
typedef int fixp30;//1.1.30


#define FIXP20_INT_TO(i) ((i)<<FIXP20_SHIFT)
#define FIXP20_FLOAT_TO(f) ((float)(f) * (float)FIXP20_MAG+0.5)
#define FIXP20_TO_INT(fp) ((fp)>>FIXP20_SHIFT)

#define FIXP19_INT_TO(i) ((i)<<FIXP19_SHIFT)
#define FIXP19_FLOAT_TO(f) ((float)(f) * (float)FIXP19_MAG+0.5)
#define FIXP19_TO_INT(fp) ((fp)>>FIXP19_SHIFT)
#define FIXP19_ROUNDUP(fp) ((fp)+FIXP19_ROUNDUP_ADD)



#define FIXP16_INT_TO(i) ((i)<<FIXP16_SHIFT)
#define FIXP16_FLOAT_TO(f) ((float)(f) * (float)FIXP16_MAG+0.5)
#define FIXP16_TO_FLOAT(fp) ((float)(fp)/(float)FIXP16_MAG)
#define FIXP16_TO_INT(fp) ((fp)>>FIXP16_SHIFT)
#define FIXP16_ROUNDUP(fp) ((fp)+FIXP16_ROUNDUP_ADD)

#define FIXP22_INT_TO(i)  ((i)<<FIXP22_SHIFT)
#define FIXP22_FLOAT_TO(f) ((float)(f) * (float)FIXP22_MAG+0.5)

#define FIXP28_TO_FIXP22(i) ((i)>>6)
#define FIXP28_TO_FIXP20(i) ((i)>>8)

#define FIXP28_INT_TO(i)  ((i)<<FIXP28_SHIFT)
#define FIXP28_FLOAT_TO(f) ((float)(f) * (float)FIXP28_MAG+0.5)

#define FIXP30_INT_T0(i) ((i)<<FIXP30_SHIFT)
#define FIXP30_FLOAT_TO(f) ((float)(f) * (float)FIXP30_MAG+0.5)

#define FIXP29_INT_T0(i) ((i)<<FIXP29_SHIFT)
#define FIXP29_FLOAT_TO(f) ((float)(f) * (float)FIXP29_MAG+0.5)


#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)>(b)?(b):(a))
#define SWAP(a,b,tmp) {tmp=a;a=b;b=tmp;}

#define FEQUALE3(a,b) ( (fabs((a)-(b)) < EPSILON_E3) ? 1 : 0)
#define FEQUALE2(a,b) ( (fabs((a)-(b)) < EPSILON_E2) ? 1 : 0)
#define FEQUALE1(a,b) ( (fabs((a)-(b)) < EPSILON_E1) ? 1 : 0)
#define OUTBOUND(mi,ma,a) (((a)<(mi))||((a)>(ma)))


#define SET_BIT(dest,bit) ((dest)|=(bit))
#define CLEAR_BIT(dest,bit) ((dest)&= ~(bit))
#define IS_BIT(dest,bit) ((dest)&(bit))


#define MATRIX_ZERO_4X4(pm) {memset((void *)pm,0,sizeof(Matrix4x4));}
#define MATRIX_ZERO_4X3(pm) {memset((void *)pm,0,sizeof(Matrix4x3));}
#define MATRIX_ZERO_4X1(pm) {memset((void *)pm,0,sizeof(Matrix4x1));}
#define MATRIX_ZERO_3X3(pm) {memset((void *)pm,0,sizeof(Matrix3x3));}
#define MATRIX_ZERO_2X2(pm) {memset((void *)pm,0,sizeof(Matrix2x2));}

#define MATRIX_IDENTITY_2X2(pm) {memcpy((void *)(pm), (void *)&E_2x2, sizeof(Matrix2x2));}
#define MATRIX_IDENTITY_3X3(pm) {memcpy((void *)(pm), (void *)&E_3x3, sizeof(Matrix3x3));}
#define MATRIX_IDENTITY_4X4(pm) {memcpy((void *)(pm), (void *)&E_4x4, sizeof(Matrix4x4));}
#define MATRIX_IDENTITY_4X3(pm) {memcpy((void *)(pm), (void *)&E_4x3, sizeof(Matrix4x3));}
#define MATRIX_IDENTITY_4X1(pm) {memcpy((void *)(pm), (void *)&E_4x3, sizeof(Matrix4x1));}

#define MATRIX_COPY_2X2(src, dest) {memcpy((void *)(dest), (void *)(src), sizeof(Matrix2x2) ); }
#define MATRIX_COPY_3X3(src, dest) {memcpy((void *)(dest), (void *)(src), sizeof(Matrix3x3) ); }
#define MATRIX_COPY_4X4(src, dest) {memcpy((void *)(dest), (void *)(src), sizeof(Matrix4x4) ); }
#define MATRIX_COPY_4X3(src, dest) {memcpy((void *)(dest), (void *)(src), sizeof(Matrix4x3) ); }
#define MATRIX_COPY_4X1(src, dest) {memcpy((void *)(dest), (void *)(src), sizeof(Matrix4x1) ); }

/*
	inline method
*/



inline void Matrix_Transpose_4x4(Matrix4x4 *src) {
	Matrix4x4 m;
	m.M00 = src->M00; m.M01 = src->M10; m.M02 = src->M20;m.M03 = src->M30;
	m.M10 = src->M01; m.M11 = src->M11; m.M12 = src->M21;m.M13 = src->M31;
	m.M20 = src->M02; m.M21 = src->M12; m.M22 = src->M22;m.M23 = src->M32;
	m.M30 = src->M03; m.M31 = src->M13; m.M32 = src->M32;m.M33 = src->M33;
	MATRIX_COPY_4X4(&m, src);
}
inline void Matrix_Transpose_3x3(Matrix3x3 *src) {
	Matrix3x3 m;
	m.M00 = src->M00; m.M01 = src->M10; m.M02 = src->M20;
	m.M10 = src->M01; m.M11 = src->M11; m.M12 = src->M21;
	m.M20 = src->M02; m.M21 = src->M12; m.M22 = src->M22;
	MATRIX_COPY_3X3(&m, src);
}

inline void Matrix_Transpose_4x4(Matrix4x4 *src, Matrix4x4 *dest) {
	dest->M00 = src->M00; dest->M01 = src->M10; dest->M02 = src->M20;dest->M03 = src->M30;
	dest->M10 = src->M01; dest->M11 = src->M11; dest->M12 = src->M21;dest->M13 = src->M31;
	dest->M20 = src->M02; dest->M21 = src->M12; dest->M22 = src->M22;dest->M23 = src->M32;
	dest->M30 = src->M03; dest->M31 = src->M13; dest->M32 = src->M32;dest->M33 = src->M33;
}

inline void Matrix_Transpose_4x4(Matrix3x3 *src, Matrix3x3 *dest) {
	dest->M00 = src->M00; dest->M01 = src->M10; dest->M02 = src->M20;
	dest->M10 = src->M01; dest->M11 = src->M11; dest->M12 = src->M21;
	dest->M20 = src->M02; dest->M21 = src->M12; dest->M22 = src->M22;
}

inline void Vector2D_Zero(Vector2D *src) {
	src->x = 0;
	src->y = 0;
}

inline void Vector3D_Zero(Vector3D *src) {
	src->x = 0;
	src->y = 0;
	src->z = 0;
}


inline void Vector4D_Zero(Vector4D *src) {
	src->x = 0;
	src->y = 0;
	src->z = 0;
	src->w = 0.0f;
}



inline void Vector2D_Init(Vector2D *src, float x, float y) {
	src->x = x;
	src->y = y;
}
inline void Vector2D_Init(Vector2D *dest, Vector2D *src) {
	dest->x = src->x;
	dest->y = src->y;
}

inline void Vector3D_Init(Vector3D *src, float x, float y, float z) {
	src->x = x;
	src->y = y;
	src->z = z;
}

inline void Vector3D_Init(Vector3D *dest, Vector3D *src) {
	dest->x = src->x;
	dest->y = src->y;
	dest->z = src->z;
}


inline void Vector4D_Init(Vector4D *src, float x, float y, float z) {
	src->x = x;
	src->y = y;
	src->z = z;
	src->w = 0.0f;
}

inline void Vector4D_Init(Vector4D *dest, Vector4D *src) {
	dest->x = src->x;
	dest->y = src->y;
	dest->z = src->z;
	dest->w = 0.0f;
}




inline void Vector2D_Copy(Vector2D *dest,Vector2D *src) {
	dest->x = src->x;
	dest->y = src->y;
}

inline void Vector3D_Copy(Vector3D *dest, Vector3D *src) {
	dest->x = src->x;
	dest->y = src->y;
	dest->z = src->z;
}

inline void Vector4D_Copy(Vector4D *dest, Vector4D *src) {
	dest->x = src->x;
	dest->y = src->y;
	dest->z = src->z;
	dest->w = src->w;
}






inline void Point2D_Zero(Point2D *src) {
	src->x = 0;
	src->y = 0;
}

inline void  Point3D_Zero(Point3D *src) {
	src->x = 0;
	src->y = 0;
	src->z = 0;
}


inline void  Point4D_Zero(Point4D *src) {
	src->x = 0;
	src->y = 0;
	src->z = 0;
	src->w = 1.0f;
}



inline void  Point2D_Init(Point2D *src, float x, float y) {
	src->x = x;
	src->y = y;
}

inline void  Point2D_Init(Point2D *dest, Point2D *src) {
	dest->x = src->x;
	dest->y = src->y;
}

inline void Point3D_Init(Point3D *src, float x, float y, float z) {
	src->x = x;
	src->y = y;
	src->z = z;
}


inline void  Point3D_Init(Point3D *dest, Point3D *src) {
	dest->x = src->x;
	dest->y = src->y;
	dest->z = src->z;
}



inline void Point4D_Init(Point4D *src, float x, float y, float z) {
	src->x = x;
	src->y = y;
	src->z = z;
	src->w = 1.0f;
}

inline void  Point4D_Init(Point4D *dest, Point4D *src) {
	dest->x = src->x;
	dest->y = src->y;
	dest->z = src->z;
	dest->w = 1.0f;
}

inline void Point2D_Copy(Point2D *dest, Point2D *src) {
	dest->x = src->x;
	dest->y = src->y;
}

inline void Point3D_Copy(Point3D *dest, Point3D *src) {
	dest->x = src->x;
	dest->y = src->y;
	dest->z = src->z;
}

inline void Point4D_Copy(Point4D *dest, Point4D *src) {
	dest->x = src->x;
	dest->y = src->y;
	dest->z = src->z;
	dest->w = src->w;
}



inline void Point4D_DivByW(Point4D *src) {
	src->x /= src->w;
	src->y /= src->w;
	src->z /= src->w;

}
inline void Point4D_DivByW(Point4D*src, Point3D *dest) {
	dest->x = src->x /= src->w;
	dest->y = src->y /= src->w;
	dest->z = src->z /= src->w;
}




/*
	function
*/

int Log_Base2(int x);
void SinCosTables_Build();
float Sin_Fast(float theta);
float Cos_Fast(float theta);
float Distance2D_Fast(int x,int y);
float Distance3D_Fast(float x, float y, float z);


//vector point math, vector accuracy ensure E5;
//cross   v0xv1
//add v0+v1
//sub v0-v1
//build p0->p1 // v=p0-p1

void Vector2D_Build(Point2D *p0, Point2D *p1, Vector2D *vResult);
void Vector2D_Add(Vector2D *v1, Vector2D *v2, Vector2D *vSum);
void Vector2D_Sub(Vector2D *v1, Vector2D *v2, Vector2D  *vDiff);
void Vector2D_Scalar(float scalar,Vector2D *v,Vector2D *vResult);
void Vector2D_Scalar(float scalar, Vector2D *vSelf);
float Vector2D_Dot(Vector2D *v1, Vector2D *v2);
float Vector2D_Length(Vector2D *v1);
float Vector2D_Length_Fast(Vector2D *v1);
void Vector2D_Normalize(Vector2D *vSelf);
void Vector2D_Normalize(Vector2D *v,Vector2D *vResult);
float Vector2D_CosTh(Vector2D *v1, Vector2D *v2);

void Vector3D_Build(Point3D *p0, Point3D *p1, Vector3D *vResult);
void Vector3D_Add(Vector3D *v1, Vector3D *v2, Vector3D *vSum);
void Vector3D_Sub(Vector3D *v1, Vector3D *v2, Vector3D  *vDiff);
void Vector3D_Mul_Scalar(float scalar, Vector3D *v, Vector3D *vResult);
void Vector3D_Mul_Scalar(float scalar, Vector3D *vSelf);
float Vector3D_Dot(Vector3D *v1, Vector3D *v2);
void Vector3D_Cross(Vector3D *v1, Vector3D *v2, Vector3D *vResult);
float Vector3D_Length(Vector3D *v1);
float Vector3D_Length_Fast(Vector3D *v1VECTOR3D_CosTh);
void Vector3D_Normalize(Vector3D *vSelf);
void Vector3D_Normalize(Vector3D *v, Vector3D *vResult);
float Vector3D_CosTh(Vector3D *v1, Vector3D *v2);

void Vector4D_Build(Point4D *p0, Point4D *p1, Vector4D *vResult);

void Vector4D_Add(Vector4D * vSelf, Vector4D * p);
void Vector4D_Add(Vector4D *p0, Vector4D *p1, Vector4D *pResult);
void Vector4D_Sub(Vector4D *v0, Vector4D *v1, Vector4D  *vResult);
void Vector4D_Sub(Vector4D *vSelf, Vector4D *v);
void Vector4D_Mul_Scalar(float scalar, Vector4D *vSelf);
void Vector4D_Mul_Scalar(float scalar, Vector4D *v, Vector4D *vResult);
float Vector4D_Length(Vector4D *v0);
float Vector4D_Length_Fast(Vector4D *v0);

void Vector4D_Normalize(Vector4D *v, Vector4D *vResult);
void Vector4D_Normalize(Vector4D *vSelf);
float Vector4D_Dot(Vector4D *v0, Vector4D *v1);
void Vector4D_Cross(Vector4D *v0, Vector4D *v1,Vector4D *vResult);//v0 X v1

void Point4D_Move(Point4D *p, Vector4D * v);
void Point4D_Move(Point4D *src,Point4D *dest, Vector4D * v);
void Point4D_Move(Point4D *p, Vector3D * v);


//matrix math
void Matrix_Add_2x2(Matrix2x2 *m1, Matrix2x2 *m2, Matrix2x2 *mResult);
void Matrix_Sub_2x2(Matrix2x2 *m1, Matrix2x2 *m2, Matrix2x2 *mResult);
void Matrix_Mul_2x2(Matrix2x2 *m1, Matrix2x2 *m2, Matrix2x2 *mResult);
float Matrix_Det_2x2(Matrix2x2 m);
void Matrix_Inverse(Matrix2x2 m, Matrix2x2 mResult);






//true Matrix4x4 math,because of speed,no 4x4 inverse  
void Matrix_Add_4x4(Matrix4x4 * m1, Matrix4x4 * m2, Matrix4x4 *mResult);
void Matrix_Mul_4x4(Matrix4x4 * m1, Matrix4x4 * m2, Matrix4x4 * mResult);
void Matrix_Mul_1x4_4x4(Matrix1x4 * m1, Matrix4x4 * m2, Matrix1x4 * mResult);
void Matrix_Mul_Vector4D_4x4(Vector4D *  v, Matrix4x4 * m, Vector4D *  vResult);
void Matrix_Mul_Point4D_4x4(Point4D *  v, Matrix4x4 * m, Point4D *  vResult);
void Matrix_Build_Rotation_4x4(Matrix4x4 *m,float thX, float thY, float thZ);
void Matrix_Build_Translation_4x4(Matrix4x4 *m,float x,float y,float z);

//not true Matrix4x4 math,but for homogeneous Matrix math,so keep Matrix.col4 [0,0,0,1] and Vector.col4[0]  in mind,
//homogeneous matrix mul add sub inverse (no transpose) always homogeneous matrix ,and so homogeneous matrix not a orthogonal matrix
void Matrix_Mul_Vector3D_4x4(Vector3D *  v, Matrix4x4 * m, Vector3D *  vResult);
void Matrix_Mul_Vecotr3D_4x3(Vector3D *  v, Matrix4x3 *m, Vector3D *  vResult);
void Matrix_Mul_Vector4D_4x3(Vector4D *  v, Matrix4x4 * m, Vector4D *  vResult);
int Matrix_Inverse_4x4(Matrix4x4 * m, Matrix4x4 * mResult);
void Matrix_Init_4x4(Matrix4x4 * m,
	float m00, float m01, float m02, float m03,
	float m10, float m11, float m12, float m13,
	float m20, float m21, float m22, float m23,
	float m30, float m31, float m32, float m33);

//line 

enum IntersectResult
{
	INTERSECT_NONE = 0,
	INTERSECT_IN = 1,
	INTERSECT_OUT = 2,
	INTERSECT_ALLIN = 3
};

void Line2D_Init(Point2D *p0,Point2D *p1,Line2D *line);  
void Line2D_GetFractionPoint(Line2D *line, float fraction, Point2D *point);
IntersectResult Line2D_GetIntersectPoint(Line2D *line1, Line2D *line2, Point2D *pResult);
IntersectResult Line2D_GetIntersectFraction(Line2D *line1, Line2D *line2, float *line1Fraction,float *line2Fraction);


void Line3D_Init(Point3D *p0, Point3D *p1, Line3D *line);
void Line3D_GetFractionPoint(Line3D *line, float fraction, Point3D *point);

//plane

void Plane3D_Init(Plane3D *plane, Point3D *p, Vector3D *n,int normalize);
int Plane3D_IsPointIn(Plane3D *plane, Point3D *p);
IntersectResult Plane3D_GetIntersectPoint(Plane3D *plane, Line3D * line, float *fraction, Point3D *pResult);

