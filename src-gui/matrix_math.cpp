/*
Copyright 2010 Etay Meiri

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifdef WIN32

#define _USE_MATH_DEFINES
#endif // WIN32

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <matrix_math.h>

using namespace std;

Quaternion::Quaternion(float _x, float _y, float _z, float _w){
  x = _x;
  y = _y;
  z = _z;
  w = _w;
}

Quaternion::Quaternion(Vector3f axis, float angle){
  float s = sinf(angle/2.0f);
  x = axis.x * s;
  y = axis.y * s;
  z = axis.z * s;
  w = cosf(angle/2.0f);
}

void Quaternion::Normalize(){
  float Length = sqrtf(x*x + y*y + z*z + w*w);
  x /= Length;
  y /= Length;
  z /= Length;
  w /= Length;
}

Quaternion Quaternion::Conjugate(){
  Quaternion ret(-x, -y, -z, w);
  return ret;
}

Quaternion operator*(const Quaternion& l, const Quaternion& r){
  const float w = (l.w * r.w) - (l.x * r.x) - (l.y * r.y) - (l.z * r.z);
  const float x = (l.x * r.w) + (l.w * r.x) + (l.y * r.z) - (l.z * r.y);
  const float y = (l.y * r.w) + (l.w * r.y) + (l.z * r.x) - (l.x * r.z);
  const float z = (l.z * r.w) + (l.w * r.z) + (l.x * r.y) - (l.y * r.x);
  Quaternion ret(x, y, z, w);
  return ret;
}

Quaternion operator*(const Quaternion& q, const float * v){
  const float w = - (q.x * v[0]) - (q.y * v[1]) - (q.z * v[2]);
  const float x =   (q.w * v[0]) - (q.y * v[2]) - (q.z * v[1]);
  const float y =   (q.w * v[1]) - (q.z * v[0]) - (q.x * v[2]);
  const float z =   (q.w * v[2]) - (q.x * v[1]) - (q.y * v[0]);
  Quaternion ret(x, y, z, w);
  return ret;
}

void Quaternion::ToDegrees(float * f){
  f[0] = atan2(x*z + y*w, x*w - y*z);
  f[1] = acos(-x*x - y*y -z*z -w*w);
  f[2] = atan2(x*z - y*w, x*w + y*z);
}

Vector3f Vector3f::Cross(const Vector3f& v) const {
  const float _x = y * v.z - z * v.y;
  const float _y = z * v.x - x * v.z;
  const float _z = x * v.y - y * v.x;
  return Vector3f(_x, _y, _z);
}

float Vector3f::Dot(const Vector3f& v) const {
  const float _x = x * v.x;
  const float _y = y * v.y;
  const float _z = z * v.z;
  return _x + _y + _z;
}

float Vector3f::Length() {
  return sqrt(x*x + y*y + z*z);
}

Vector3f& Vector3f::Normalize(){
  const float Length = sqrtf(x*x + y*y + z*z);
  if (Length == 0){
    x = 0; y = 0; z = 0;
  } else {
    x /= Length;
    y /= Length;
    z /= Length;
  }
  return *this;
}

void Vector3f::Rotate(float Angle, const Vector3f& Axe){
  const float SinHalfAngle = sinf(ToRadian(Angle/2));
  const float CosHalfAngle = cosf(ToRadian(Angle/2));

  const float Rx = Axe.x * SinHalfAngle;
  const float Ry = Axe.y * SinHalfAngle;
  const float Rz = Axe.z * SinHalfAngle;
  const float Rw = CosHalfAngle;
  Quaternion RotationQ(Rx, Ry, Rz, Rw);
  Quaternion ConjugateQ = RotationQ.Conjugate();
  Quaternion W = RotationQ * (*this) * ConjugateQ;

  x = W.x;
  y = W.y;
  z = W.z;
}

void Matrix4f::InitScaleTransform(float ScaleX, float ScaleY, float ScaleZ)
{
    m[0][0] = ScaleX; m[0][1] = 0.0f;   m[0][2] = 0.0f;   m[0][3] = 0.0f;
    m[1][0] = 0.0f;   m[1][1] = ScaleY; m[1][2] = 0.0f;   m[1][3] = 0.0f;
    m[2][0] = 0.0f;   m[2][1] = 0.0f;   m[2][2] = ScaleZ; m[2][3] = 0.0f;
    m[3][0] = 0.0f;   m[3][1] = 0.0f;   m[3][2] = 0.0f;   m[3][3] = 1.0f;
}

void Matrix4f::InitRotateTransform(float RotateX, float RotateY, float RotateZ)
{
    Matrix4f rx, ry, rz;

    const float x = ToRadian(RotateX);
    const float y = ToRadian(RotateY);
    const float z = ToRadian(RotateZ);

    rx.m[0][0] = 1.0f; rx.m[0][1] = 0.0f   ; rx.m[0][2] = 0.0f    ; rx.m[0][3] = 0.0f;
    rx.m[1][0] = 0.0f; rx.m[1][1] = cosf(x); rx.m[1][2] = -sinf(x); rx.m[1][3] = 0.0f;
    rx.m[2][0] = 0.0f; rx.m[2][1] = sinf(x); rx.m[2][2] = cosf(x) ; rx.m[2][3] = 0.0f;
    rx.m[3][0] = 0.0f; rx.m[3][1] = 0.0f   ; rx.m[3][2] = 0.0f    ; rx.m[3][3] = 1.0f;

    ry.m[0][0] = cosf(y); ry.m[0][1] = 0.0f; ry.m[0][2] = -sinf(y); ry.m[0][3] = 0.0f;
    ry.m[1][0] = 0.0f   ; ry.m[1][1] = 1.0f; ry.m[1][2] = 0.0f    ; ry.m[1][3] = 0.0f;
    ry.m[2][0] = sinf(y); ry.m[2][1] = 0.0f; ry.m[2][2] = cosf(y) ; ry.m[2][3] = 0.0f;
    ry.m[3][0] = 0.0f   ; ry.m[3][1] = 0.0f; ry.m[3][2] = 0.0f    ; ry.m[3][3] = 1.0f;

    rz.m[0][0] = cosf(z); rz.m[0][1] = -sinf(z); rz.m[0][2] = 0.0f; rz.m[0][3] = 0.0f;
    rz.m[1][0] = sinf(z); rz.m[1][1] = cosf(z) ; rz.m[1][2] = 0.0f; rz.m[1][3] = 0.0f;
    rz.m[2][0] = 0.0f   ; rz.m[2][1] = 0.0f    ; rz.m[2][2] = 1.0f; rz.m[2][3] = 0.0f;
    rz.m[3][0] = 0.0f   ; rz.m[3][1] = 0.0f    ; rz.m[3][2] = 0.0f; rz.m[3][3] = 1.0f;

    *this = rz * ry * rx;
}

void Matrix4f::InitRotateAxisTransform(Vector3f u, float th)
{
  Matrix4f R;

  R.m[0][0] = cosf(th) + u.x*u.x*(1.0f - cosf(th));
  R.m[0][1] = u.x*u.y*(1 - cosf(th)) - u.z*sinf(th);
  R.m[0][2] = u.x*u.z*(1 - cosf(th)) + u.y*sinf(th);
  R.m[0][3] = 0;

  R.m[1][0] = u.x*u.y*(1 - cosf(th)) + u.z*sinf(th);
  R.m[1][1] = cosf(th) + u.y*u.y*(1.0f - cosf(th));
  R.m[1][2] = u.y*u.z*(1 - cosf(th)) - u.x*sinf(th);
  R.m[1][3] = 0;

  R.m[2][0] = u.x*u.z*(1 - cosf(th)) - u.y*sinf(th);
  R.m[2][1] = u.z*u.y*(1 - cosf(th)) + u.x*sinf(th);
  R.m[2][2] = cosf(th) + u.z*u.z*(1.0f - cosf(th));
  R.m[2][3] = 0;

  R.m[3][0] = 0;
  R.m[3][1] = 0;
  R.m[3][2] = 0;
  R.m[3][3] = 1;
  
    *this = R;
}


void Matrix4f::InitRotateTransform(const Quaternion& quat)
{
  float yy2 = 2.0f * quat.y * quat.y;
  float xy2 = 2.0f * quat.x * quat.y;
  float xz2 = 2.0f * quat.x * quat.z;
  float yz2 = 2.0f * quat.y * quat.z;
  float zz2 = 2.0f * quat.z * quat.z;
  float wz2 = 2.0f * quat.w * quat.z;
  float wy2 = 2.0f * quat.w * quat.y;
  float wx2 = 2.0f * quat.w * quat.x;
  float xx2 = 2.0f * quat.x * quat.x;
  m[0][0] = 1 - yy2 - zz2;
  m[0][1] = xy2 - wz2;
  m[0][2] = xz2 + wy2;
  m[0][3] = 0;
  m[1][0] = xy2 + wz2;
  m[1][1] = 1 - xx2 - zz2;
  m[1][2] = yz2 + wx2;
  m[1][3] = 0;
  m[2][0] = xz2 - wy2;
  m[2][1] = yz2 - wx2;
  m[2][2] = 1 - xx2 - yy2;
  m[2][3] = 0.0f;
  m[3][0] = m[3][1] = m[3][2] = 0;
  m[3][3] = 1.0f;
}


void Matrix4f::InitTranslateTransform(float x, float y, float z)
{
    m[0][0] = 1.0f; m[0][1] = 0.0f; m[0][2] = 0.0f; m[0][3] = x;
    m[1][0] = 0.0f; m[1][1] = 1.0f; m[1][2] = 0.0f; m[1][3] = y;
    m[2][0] = 0.0f; m[2][1] = 0.0f; m[2][2] = 1.0f; m[2][3] = z;
    m[3][0] = 0.0f; m[3][1] = 0.0f; m[3][2] = 0.0f; m[3][3] = 1.0f;
}

void Matrix4f::InitCameraTransform(const float Target[3], const float Up[3])
{
    Normalize((float *)Target);
    float * N = (float *)Target;
    float U[3];
    Cross(Up, N, U);
    float V[3];
    Cross(N, U, V);

    m[0][0] = U[0];   m[0][1] = U[1];   m[0][2] = U[2];   m[0][3] = 0.0f;
    m[1][0] = V[0];   m[1][1] = V[1];   m[1][2] = V[2];   m[1][3] = 0.0f;
    m[2][0] = N[0];   m[2][1] = N[1];   m[2][2] = N[2];   m[2][3] = 0.0f;
    m[3][0] = 0.0f;  m[3][1] = 0.0f;  m[3][2] = 0.0f;  m[3][3] = 1.0f;
}

void Matrix4f::InitPersProjTransform(const PersProjInfo& p)
{
    const float ar         = p.Width / p.Height;
    const float zRange     = p.zNear - p.zFar;
    const float tanHalfFOV = tanf(ToRadian(p.FOV / 2.0f));

    m[0][0] = 1.0f/(tanHalfFOV * ar); m[0][1] = 0.0f;            m[0][2] = 0.0f;            m[0][3] = 0.0;
    m[1][0] = 0.0f;                   m[1][1] = 1.0f/tanHalfFOV; m[1][2] = 0.0f;            m[1][3] = 0.0;
    m[2][0] = 0.0f;                   m[2][1] = 0.0f;            m[2][2] = (-p.zNear - p.zFar)/zRange ; m[2][3] = 2.0f*p.zFar*p.zNear/zRange;
    m[3][0] = 0.0f;                   m[3][1] = 0.0f;            m[3][2] = 1.0f;            m[3][3] = 0.0;
}

void Matrix4f::InitOrthoProjTransform(const OrthoProjInfo& p)
{
  float l = p.Left;
  float r = p.Right;
  float b = p.Bottom;
  float t = p.Top;
  float n = p.zNear;
  float f = p.zFar;

  m[0][0] = 2.0f/(r-l); m[0][1] = 0.f; m[0][2] = 0.f; m[0][3] = -(r+l)/(r-l);
  m[1][0] = 0.f; m[1][1] = 2.0f/(t-b); m[1][2] = 0.f; m[1][3] = -(t+b)/(t-b);
  m[2][0] = 0.f; m[2][1] = 0.f; m[2][2] = 2.0f/(f-n); m[2][3] = -(f+n)/(f-n);
  m[3][0] = 0.f; m[3][1] = 0.f; m[3][2] = 0.f; m[3][3] = 1.f;
}

#ifdef WIN32
std::istream& safeGetline(std::istream& is, std::string& t){
	t.clear();

	// The characters in the stream are read one-by-one using a std::streambuf.
	// That is faster than reading them one-by-one using the std::istream.
	// Code that uses streambuf this way must be guarded by a sentry object.
	// The sentry object performs various tasks,
	// such as thread synchronization and updating the stream state.

	std::istream::sentry se(is, true);
	std::streambuf* sb = is.rdbuf();

	for (;;) {
		int c = sb->sbumpc();
		//cout << c << endl; //Debug code to see what is being parsed
		switch (c) {
		case '\n':
			return is;
		case '\r':
			if (sb->sgetc() == '\n')
				sb->sbumpc();
			return is;
		case EOF:
			// Also handle the case when the last line has no line ending
			if (t.empty())
				is.setstate(std::ios::eofbit);
			return is;
		default:
			t += (char)c;
		}
	}
}

void ReadMesh(GLfloat *v, unsigned int* i, const char * v_file, const char * i_file) {
	int v_i = 0;
	int i_i = 0;
	GLfloat x, y, z;
	unsigned int a, b, c;
	//new code


	std::string path = v_file; // path to file (v first)
	path = path.erase(0, 2); // remove unix ./ file path

	std::ifstream vfile(path.c_str());

	//making sure the file stream is open
	if (!vfile.is_open()) { //is_open reaturns true if working
		std::cout << ("Failed to open the file.  " + path) << std::endl;
		return;
	}

	int n = 0;
	std::string t;

	while (!safeGetline(vfile, t).eof()) { // read the .v file
		sscanf(t.c_str(), "%f %f %f\n", &x, &y, &z); //string is converted to constant char
		v[v_i] = x;
		v_i += 1;
		v[v_i] = y;
		v_i += 1;
		v[v_i] = z;
		v_i += 1;
	}
	//this will print garbage data if v[0 to 2] is not set
	std::cout << "printing x,y,z " << v[0] << "," << v[1] << "," << v[2] << "," << std::endl;

	path = i_file; // path to file (v first)
	path = path.erase(0, 2); // remove unix ./ file path

	std::ifstream ifile(path.c_str());

	//making sure the file stream is open
	if (!ifile.is_open()) { //is_open reaturns true if working
		std::cout << ("Failed to open the file.  " + path) << std::endl;
		return;
	}


	n = 0;
	while (!safeGetline(ifile, t).eof()) { // read the .v file
		sscanf(t.c_str(), "%d %d %d\n", &a, &b, &c); //string is converted to constant char
		i[i_i] = a;
		i_i += 1;
		i[i_i] = b;
		i_i += 1;
		i[i_i] = c;
		i_i += 1;
	}
	//std::cout << "printing x,y,z " << i[0] << "," << i[1] << "," << i[2] << "," << std::endl;

	ifile.close();
	vfile.close();

}

#endif // WIN32


// #if defined(LINUX) || defined(__APPLE__)
void ReadMesh(GLfloat *v, unsigned int* i, const char * v_file, const char * i_file) {
	FILE * fpv = NULL;
	FILE * fpi = NULL;
	char * line = NULL;
	size_t len = 0;
	ssize_t read;

	fpv = fopen(v_file, "r");
	if (fpv == NULL)
		exit(EXIT_FAILURE);


	fpi = fopen(i_file, "r");
	if (fpv == NULL)
		exit(EXIT_FAILURE);

	int v_i = 0;
	int i_i = 0;
	GLfloat x, y, z;
	unsigned int a, b, c;
	while ((read = getline(&line, &len, fpv)) != -1) {
		sscanf(line, "%f %f %f\n", &x, &y, &z);
		v[v_i] = x;
		v_i += 1;
		v[v_i] = y;
		v_i += 1;
		v[v_i] = z;
		v_i += 1;
	}

	while ((read = getline(&line, &len, fpi)) != -1) {
		sscanf(line, "%d %d %d\n", &a, &b, &c);
		i[i_i] = a;
		i_i += 1;
		i[i_i] = b;
		i_i += 1;
		i[i_i] = c;
		i_i += 1;
	}

	fclose(fpv);
	fclose(fpi);
	if (line) free(line);
}

// #endif // LINUX

void Cross(const float left[3], const float right[3], float * result)
{
 result[0] = left[1]*right[2] - left[2]*right[1];
 result[1] = left[2]*right[0] - left[0]*right[2];
 result[2] = left[0]*right[1] - left[1]*right[0];
}

void Normalize(float * v)
{
  float d = sqrt(pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2));
  v[0] = v[0]/d;
  v[1] = v[1]/d;
  v[2] = v[2]/d;
}

void Pipeline::SetRotationMatrix(const Matrix4f _m)
{
  m_rotate_trans.m[0][0] = _m.m[0][0];
  m_rotate_trans.m[0][1] = _m.m[0][1];
  m_rotate_trans.m[0][2] = _m.m[0][2];
  m_rotate_trans.m[0][3] = _m.m[0][3];
  m_rotate_trans.m[1][0] = _m.m[1][0];
  m_rotate_trans.m[1][1] = _m.m[1][1];
  m_rotate_trans.m[1][2] = _m.m[1][2];
  m_rotate_trans.m[1][3] = _m.m[1][3];
  m_rotate_trans.m[2][0] = _m.m[2][0];
  m_rotate_trans.m[2][1] = _m.m[2][1];
  m_rotate_trans.m[2][2] = _m.m[2][2];
  m_rotate_trans.m[2][3] = _m.m[2][3];
  m_rotate_trans.m[3][0] = _m.m[3][0];
  m_rotate_trans.m[3][1] = _m.m[3][1];
  m_rotate_trans.m[3][2] = _m.m[3][2];
  m_rotate_trans.m[3][3] = _m.m[3][3];
}

void Pipeline::SetPostRotationMatrix(const Matrix4f _m)
{
  m_post_rotate_trans.m[0][0] = _m.m[0][0];
  m_post_rotate_trans.m[0][1] = _m.m[0][1];
  m_post_rotate_trans.m[0][2] = _m.m[0][2];
  m_post_rotate_trans.m[0][3] = _m.m[0][3];
  m_post_rotate_trans.m[1][0] = _m.m[1][0];
  m_post_rotate_trans.m[1][1] = _m.m[1][1];
  m_post_rotate_trans.m[1][2] = _m.m[1][2];
  m_post_rotate_trans.m[1][3] = _m.m[1][3];
  m_post_rotate_trans.m[2][0] = _m.m[2][0];
  m_post_rotate_trans.m[2][1] = _m.m[2][1];
  m_post_rotate_trans.m[2][2] = _m.m[2][2];
  m_post_rotate_trans.m[2][3] = _m.m[2][3];
  m_post_rotate_trans.m[3][0] = _m.m[3][0];
  m_post_rotate_trans.m[3][1] = _m.m[3][1];
  m_post_rotate_trans.m[3][2] = _m.m[3][2];
  m_post_rotate_trans.m[3][3] = _m.m[3][3];
}

const Matrix4f * Pipeline::GetProjTrans(){
  m_ProjTransformation.InitPersProjTransform(m_projInfo);
//  m_ProjTransformation.InitOrthoProjTransform(m_orthoInfo);
  return &m_ProjTransformation;
}

const Matrix4f * Pipeline::GetViewTrans(){
  Matrix4f CamTranslateTrans, CamRotateTrans;
  CamTranslateTrans.InitTranslateTransform(-m_camera.Pos[0], -m_camera.Pos[1], -m_camera.Pos[2]);
  CamRotateTrans.InitCameraTransform(m_camera.Target, m_camera.Up);
  m_Vtransformation = CamRotateTrans * CamTranslateTrans;
  return &m_Vtransformation;
}

const Matrix4f * Pipeline::GetWorldTrans(){
  Matrix4f ScaleTrans, RotateTrans, TranslateTrans, PostRotateTrans;
  ScaleTrans.InitScaleTransform(m_scale[0], m_scale[1], m_scale[2]);
  TranslateTrans.InitTranslateTransform(m_pos[0], m_pos[1], m_pos[2]);
  m_Wtransformation = m_post_rotate_trans * TranslateTrans * m_rotate_trans * ScaleTrans;
  return &m_Wtransformation;
}

const Matrix4f * Pipeline::GetVPTrans(){
  GetViewTrans();
  GetProjTrans();
  m_VPtransformation = m_ProjTransformation * m_Vtransformation;
  return &m_VPtransformation;
}

const Matrix4f * Pipeline::GetWVPTrans(){
  GetWorldTrans();
  GetVPTrans();
  m_WVPtransformation = m_VPtransformation * m_Wtransformation;
  return &m_WVPtransformation;
}


class Camera
{
public:
  Camera();
  Camera(const float Pos[3], const float Target[3], const float Up[3]);
  bool OnKeyboard(int key);
  const float * GetPos();
  const float * GetTarget();
  const float * GetUp();

private:
  float m_pos[3];
  float m_target[3];
  float m_up[3];

};
