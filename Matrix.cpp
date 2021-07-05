#include <d3dx9.h>
#include <d3d9.h>
#pragma comment(lib, "d3d9.lib")
#pragma comment(lib, "d3dx9.lib")
#include <dwmapi.h>
#pragma comment(lib, "dwmapi.lib")
#include<windows.h>
#include<cassert>
#include<vector>
#include <iostream>
#include <cmath>

#include "./Blackbone/src/BlackBone/Process/Process.h"


using namespace std;

#define PI ((float)3.141592)
#define DEG2RAD( x ) ( ( float )( x ) * ( float )( ( float )( PI ) / 180.0f ) )
#define RAD2DEG( x ) ( ( float )( x ) * ( float )( 180.0f / ( float )( PI ) ) )


IDirect3D9Ex* p_Object = nullptr;
D3DPRESENT_PARAMETERS p_Params;
IDirect3DDevice9Ex* p_Device = nullptr;

int DirectXInit(HWND hWnd);

int Width = GetSystemMetrics(SM_CXSCREEN);
int Height = GetSystemMetrics(SM_CYSCREEN);

typedef float vec_t;
class QAngle
{
public:
	// Members
	vec_t x, y, z;
};

class Vector
{
public:
	// Members
	vec_t x, y, z;
	Vector operator*(float v);
	float operator[] (int i) const ;
	
};

 void VectorMultiply(const Vector& a, vec_t b, Vector& c)
{
	c.x = a.x * b;
	c.y = a.y * b;
	c.z = a.z * b;
}

void AngleVectors(const Vector& angles, Vector* forward);

void MatrixMultiply(const D3DMATRIX& s1, const D3DMATRIX& s2, D3DMATRIX& dst);
void MatrixRotate(D3DMATRIX& dst, const Vector& vAxisOfRot, float angleDegrees);
void MatrixBuildRotationAboutAxis(const Vector& vAxisOfRot, float angleDegrees, D3DMATRIX& dst);
void MatrixBuildRotationAboutAxis(D3DMATRIX& dst, const Vector& vAxisOfRot, float angleDegrees);
void MatrixSetIdentity(D3DMATRIX& dst);
void MatrixBuildTranslation(D3DMATRIX& dst, const Vector& translation);
void MatrixTranslate(D3DMATRIX& dst, const Vector& translation);

static Vector g_CamerOrigin;
static Vector g_ViewAngles;
static float g_Fovx;
static float g_AsppectRatio;
static float g_Znear;
static float g_Zfar;
static float g_Fovy;

static D3DXMATRIX ModelView;

static D3DXMATRIX MatrixProjection;

static D3DXMATRIX Tview;

static D3DXMATRIX Rview;

static D3DXMATRIX martix_view;

static D3DXMATRIX ModelViewProtection;

static D3DXMATRIX Rotation;

static D3DMATRIX baseRotation;

int main()
{
    
	blackbone::Process csgo;
	csgo.Attach(L"csgo.exe");

	blackbone::ProcessMemory& mm = csgo.memory();
	blackbone::ProcessModules& modules = csgo.modules();

	auto client = modules.GetModule(L"client.dll");

	

	
		//Camera_t	camera.h (public\mathlib)
		mm.Read<Vector>(client.get()->baseAddress + 0x4D937F0, g_CamerOrigin);
		mm.Read<Vector>(client.get()->baseAddress + 0x4D937F0 + 12, g_ViewAngles);
		mm.Read<float>(client.get()->baseAddress + 0x4D937F0 + 24, g_Fovx);
		mm.Read<float>(client.get()->baseAddress + 0x4D937F0 + 28, g_AsppectRatio);
		mm.Read<float>(client.get()->baseAddress + 0x4D937F0 + 32, g_Znear);
		mm.Read<float>(client.get()->baseAddress + 0x4D937F0 + 36, g_Zfar);

		ModelView.m[0][0] = 1;
		ModelView.m[1][1] = 1;
		ModelView.m[2][2] = 1;
		ModelView.m[3][0] = 1;

		ModelView.m[0][3] = -g_CamerOrigin.x;
		ModelView.m[1][3] = -g_CamerOrigin.y;
		ModelView.m[2][3] = -g_CamerOrigin.z;

		Tview.m[0][0] = 1;
		Tview.m[1][1] = 1;
		Tview.m[2][2] = 1;
		Tview.m[3][3] = 1;
		Tview.m[0][3] = -g_CamerOrigin.x;
		Tview.m[1][3] = -g_CamerOrigin.y;
		Tview.m[2][3] = -g_CamerOrigin.z;


		static D3DMATRIX baseRotation;

		MatrixBuildRotationAboutAxis(baseRotation, Vector(1, 0, 0), -90);
		MatrixRotate(baseRotation, Vector(0, 0, 1), 90);
		martix_view = baseRotation;

		//roll绕x轴转
		MatrixRotate(martix_view, Vector(1, 0, 0), -g_ViewAngles[2]);
		//Pitch绕y轴转
		MatrixRotate(martix_view, Vector(0, 1, 0), -g_ViewAngles[0]);
		//Yaw是绕着z轴转
		MatrixRotate(martix_view, Vector(0, 0, 1), -g_ViewAngles[1]);

		//将摄像机移到原点
		g_CamerOrigin.x = -g_CamerOrigin.x;
		g_CamerOrigin.y = -g_CamerOrigin.y;
		g_CamerOrigin.z = -g_CamerOrigin.z;
		MatrixTranslate(martix_view, g_CamerOrigin);

		/*Vector forward;
		AngleVectors(g_ViewAngles, &forward);
		auto target_point = forward * 100;*/

		


		//D3DXVECTOR3 vEye(g_CamerOrigin.x, g_CamerOrigin.y, g_CamerOrigin.z);
		//D3DXVECTOR3 vAt(target_point.x + g_CamerOrigin.x, 
		//	target_point.y+g_CamerOrigin.y,
		//	target_point.z+g_CamerOrigin.z);
		////坐标系是可以旋转的
		////vUp可以是(1,0,0)(-1,0,0)(0,1,0)(0,-1,0)
		//D3DXVECTOR3 vUp(1.0f, 0.0f, 0.0f);
		//D3DXMatrixLookAtLH(&martix_view, &vEye, &vAt, &vUp);




		float halfWidth = tan(g_Fovx * PI / 360.0);
		float halfHeight = halfWidth / g_AsppectRatio;
		MatrixProjection.m[0][0] = 1.0f / halfWidth;
		MatrixProjection.m[1][1] = 1.0f / halfHeight;
		MatrixProjection.m[2][2] = g_Zfar / (g_Znear - g_Zfar);
		MatrixProjection.m[3][2] = -1.0f;
		MatrixProjection.m[2][3] = g_Znear * g_Zfar / (g_Znear - g_Zfar);

		ModelViewProtection = MatrixProjection * martix_view ;

	printf("%8.3f %8.3f %8.3f %8.3f \n%8.3f %8.3f %8.3f %8.3f\n %8.3f %8.3f %8.3f %8.3f\n %8.3f %8.3f %8.3f %8.3f\n",
		ModelViewProtection.m[0][0], ModelViewProtection.m[0][1], ModelViewProtection.m[0][2], ModelViewProtection.m[0][3], ModelViewProtection.m[1][0], ModelViewProtection.m[1][1], ModelViewProtection.m[1][2], ModelViewProtection.m[1][3], ModelViewProtection.m[2][0], ModelViewProtection.m[2][1], ModelViewProtection.m[2][2], ModelViewProtection.m[2][3], ModelViewProtection.m[3][0], ModelViewProtection.m[3][1], ModelViewProtection.m[3][2], ModelViewProtection.m[3][3] );



	return 0;
}

int DirectXInit(HWND hWnd)
{
	if (FAILED(Direct3DCreate9Ex(D3D_SDK_VERSION, &p_Object)))
		exit(1);

	ZeroMemory(&p_Params, sizeof(p_Params));
	p_Params.Windowed = TRUE;
	p_Params.SwapEffect = D3DSWAPEFFECT_DISCARD;
	p_Params.hDeviceWindow = hWnd;
	p_Params.MultiSampleQuality = D3DMULTISAMPLE_NONE;
	p_Params.BackBufferFormat = D3DFMT_A8R8G8B8;
	p_Params.BackBufferWidth = Width;
	p_Params.BackBufferHeight = Height;
	p_Params.EnableAutoDepthStencil = TRUE;
	p_Params.AutoDepthStencilFormat = D3DFMT_D16;
	if (FAILED(p_Object->CreateDeviceEx(D3DADAPTER_DEFAULT,
		D3DDEVTYPE_HAL,
		hWnd,
		D3DCREATE_SOFTWARE_VERTEXPROCESSING,
		&p_Params,
		0,
		&p_Device)))
		exit(1);
	return 0;

}


void AngleVectors(const Vector& angles, Vector* forward) // forward.lenth(3D) is 1
{

	float	sp, sy, cp, cy;

	//yaw
	sy = sin(DEG2RAD(angles.y));//sin(0)=0
	cy = cos(DEG2RAD(angles.y));//cos(0)=1

	//pitch
	sp = sin(DEG2RAD(angles.x));//sin(pi/4) = √2/2
	cp = cos(DEG2RAD(angles.x));//cos(pi/4) = √2/2

	forward->x = cp * cy;//√2/2
	forward->y = cp * sy;//0
	forward->z = -sp;//-√2/2
}

Vector Vector::operator*(float v)
{
	Vector res;
	VectorMultiply(*this, v, res);
	return res;
}

float Vector::operator[](int i) const 
{
	return ((vec_t*)this)[i];
}


void MatrixBuildRotationAboutAxis(D3DMATRIX& dst, const Vector& vAxisOfRot, float angleDegrees)
{
	MatrixBuildRotationAboutAxis(vAxisOfRot, angleDegrees, dst);
	dst.m[3][0] = 0;
	dst.m[3][1] = 0;
	dst.m[3][2] = 0;
	dst.m[3][3] = 1;
}

void MatrixBuildRotationAboutAxis(const Vector& vAxisOfRot, float angleDegrees, D3DMATRIX& dst)
{
	float radians;
	float axisXSquared;
	float axisYSquared;
	float axisZSquared;
	float fSin;
	float fCos;

	radians = angleDegrees * (PI / 180.0);
	fSin = sin(radians);
	fCos = cos(radians);

	axisXSquared = vAxisOfRot[0] * vAxisOfRot[0];
	axisYSquared = vAxisOfRot[1] * vAxisOfRot[1];
	axisZSquared = vAxisOfRot[2] * vAxisOfRot[2];

	// Column 0:
	dst.m[0][0] = axisXSquared + (1 - axisXSquared) * fCos;
	dst.m[1][0] = vAxisOfRot[0] * vAxisOfRot[1] * (1 - fCos) + vAxisOfRot[2] * fSin;
	dst.m[2][0] = vAxisOfRot[2] * vAxisOfRot[0] * (1 - fCos) - vAxisOfRot[1] * fSin;

	// Column 1:
	dst.m[0][1] = vAxisOfRot[0] * vAxisOfRot[1] * (1 - fCos) - vAxisOfRot[2] * fSin;
	dst.m[1][1] = axisYSquared + (1 - axisYSquared) * fCos;
	dst.m[2][1] = vAxisOfRot[1] * vAxisOfRot[2] * (1 - fCos) + vAxisOfRot[0] * fSin;

	// Column 2:
	dst.m[0][2] = vAxisOfRot[2] * vAxisOfRot[0] * (1 - fCos) + vAxisOfRot[1] * fSin;
	dst.m[1][2] = vAxisOfRot[1] * vAxisOfRot[2] * (1 - fCos) - vAxisOfRot[0] * fSin;
	dst.m[2][2] = axisZSquared + (1 - axisZSquared) * fCos;

	// Column 3:
	dst.m[0][3] = 0;
	dst.m[1][3] = 0;
	dst.m[2][3] = 0;
}

void MatrixRotate(D3DMATRIX& dst, const Vector& vAxisOfRot, float angleDegrees)
{
	D3DMATRIX rotation, temp;
	MatrixBuildRotationAboutAxis(rotation, vAxisOfRot, angleDegrees);
	MatrixMultiply(dst, rotation, temp);
	dst = temp;
}

void MatrixMultiply(const D3DMATRIX& s1, const D3DMATRIX& s2, D3DMATRIX& dst)
{
	// Make sure it works if src1 == dst.m.m or src2 == dst.m

	dst.m[0][0] = s1.m[0][0] * s2.m[0][0] + s1.m[0][1] * s2.m[1][0] + s1.m[0][2] * s2.m[2][0] + s1.m[0][3] * s2.m[3][0];
	dst.m[0][1] = s1.m[0][0] * s2.m[0][1] + s1.m[0][1] * s2.m[1][1] + s1.m[0][2] * s2.m[2][1] + s1.m[0][3] * s2.m[3][1];
	dst.m[0][2] = s1.m[0][0] * s2.m[0][2] + s1.m[0][1] * s2.m[1][2] + s1.m[0][2] * s2.m[2][2] + s1.m[0][3] * s2.m[3][2];
	dst.m[0][3] = s1.m[0][0] * s2.m[0][3] + s1.m[0][1] * s2.m[1][3] + s1.m[0][2] * s2.m[2][3] + s1.m[0][3] * s2.m[3][3];

	dst.m[1][0] = s1.m[1][0] * s2.m[0][0] + s1.m[1][1] * s2.m[1][0] + s1.m[1][2] * s2.m[2][0] + s1.m[1][3] * s2.m[3][0];
	dst.m[1][1] = s1.m[1][0] * s2.m[0][1] + s1.m[1][1] * s2.m[1][1] + s1.m[1][2] * s2.m[2][1] + s1.m[1][3] * s2.m[3][1];
	dst.m[1][2] = s1.m[1][0] * s2.m[0][2] + s1.m[1][1] * s2.m[1][2] + s1.m[1][2] * s2.m[2][2] + s1.m[1][3] * s2.m[3][2];
	dst.m[1][3] = s1.m[1][0] * s2.m[0][3] + s1.m[1][1] * s2.m[1][3] + s1.m[1][2] * s2.m[2][3] + s1.m[1][3] * s2.m[3][3];

	dst.m[2][0] = s1.m[2][0] * s2.m[0][0] + s1.m[2][1] * s2.m[1][0] + s1.m[2][2] * s2.m[2][0] + s1.m[2][3] * s2.m[3][0];
	dst.m[2][1] = s1.m[2][0] * s2.m[0][1] + s1.m[2][1] * s2.m[1][1] + s1.m[2][2] * s2.m[2][1] + s1.m[2][3] * s2.m[3][1];
	dst.m[2][2] = s1.m[2][0] * s2.m[0][2] + s1.m[2][1] * s2.m[1][2] + s1.m[2][2] * s2.m[2][2] + s1.m[2][3] * s2.m[3][2];
	dst.m[2][3] = s1.m[2][0] * s2.m[0][3] + s1.m[2][1] * s2.m[1][3] + s1.m[2][2] * s2.m[2][3] + s1.m[2][3] * s2.m[3][3];

	dst.m[3][0] = s1.m[3][0] * s2.m[0][0] + s1.m[3][1] * s2.m[1][0] + s1.m[3][2] * s2.m[2][0] + s1.m[3][3] * s2.m[3][0];
	dst.m[3][1] = s1.m[3][0] * s2.m[0][1] + s1.m[3][1] * s2.m[1][1] + s1.m[3][2] * s2.m[2][1] + s1.m[3][3] * s2.m[3][1];
	dst.m[3][2] = s1.m[3][0] * s2.m[0][2] + s1.m[3][1] * s2.m[1][2] + s1.m[3][2] * s2.m[2][2] + s1.m[3][3] * s2.m[3][2];
	dst.m[3][3] = s1.m[3][0] * s2.m[0][3] + s1.m[3][1] * s2.m[1][3] + s1.m[3][2] * s2.m[2][3] + s1.m[3][3] * s2.m[3][3];
}

void MatrixSetIdentity(D3DMATRIX& dst)
{
	dst.m[0][0] = 1.0f; dst.m[0][1] = 0.0f; dst.m[0][2] = 0.0f; dst.m[0][3] = 0.0f;
	dst.m[1][0] = 0.0f; dst.m[1][1] = 1.0f; dst.m[1][2] = 0.0f; dst.m[1][3] = 0.0f;
	dst.m[2][0] = 0.0f; dst.m[2][1] = 0.0f; dst.m[2][2] = 1.0f; dst.m[2][3] = 0.0f;
	dst.m[3][0] = 0.0f; dst.m[3][1] = 0.0f; dst.m[3][2] = 0.0f; dst.m[3][3] = 1.0f;
}

void MatrixBuildTranslation(D3DMATRIX& dst, const Vector& translation)
{
	MatrixSetIdentity(dst);
	dst.m[0][3] = translation[0];
	dst.m[1][3] = translation[1];
	dst.m[2][3] = translation[2];
}

void MatrixTranslate(D3DMATRIX& dst, const Vector& translation)
{
	D3DMATRIX matTranslation, temp;
	MatrixBuildTranslation(matTranslation, translation);
	MatrixMultiply(dst, matTranslation, temp);
	dst = temp;
}
