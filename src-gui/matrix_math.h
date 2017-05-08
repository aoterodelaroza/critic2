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
#include <math.h>
#include <stdio.h>
#include <string.h>

#define ToRadian(x) ((x) * M_PI / 180.0f)
#define ToDegree(x) ((x) * 180.0f / M_PI)

void Normalize(float * v);
void Cross(const float left[3], const float right[3], float * result);

struct CameraInfo{
  float Pos[3];
  float Target[3];
  float Up[3];
};

struct PersProjInfo {
  float FOV;
  float Width;
  float Height;
  float zNear;
  float zFar;
};

struct OrthoProjInfo
{
  float Right;
  float Left;
  float Bottom;
  float Top;
  float zNear;
  float zFar;
};

struct Vector3f
{
  float x;
  float y;
  float z;
  Vector3f(){}
  Vector3f(float _x, float _y, float _z){
    x = _x;
    y = _y;
    z = _z;
  }

  Vector3f(const float * pFloat){
    x = pFloat[0];
    y = pFloat[1];
    z = pFloat[2];
  }

  Vector3f (float f){
    x = y = z= f;
  }

  inline Vector3f operator+(const Vector3f& Right) const
  {
    Vector3f Ret;
    Ret.x = x + Right.x;
    Ret.y = y + Right.y;
    Ret.z = z + Right.z;
    return Ret;
  }

  inline Vector3f operator-(const Vector3f& Right) const
  {
    Vector3f Ret;
    Ret.x = x - Right.x;
    Ret.y = y - Right.y;
    Ret.z = z - Right.z;
    return Ret;
  }



  inline Vector3f operator*(float Right) const
  {
    Vector3f Ret;
    Ret.x = x * Right;
    Ret.y = y * Right;
    Ret.z = z * Right;
    return Ret;
  }

 inline Vector3f operator/(float Right) const
  {
    Vector3f Ret;
    Ret.x = x / Right;
    Ret.y = y / Right;
    Ret.z = z / Right;
    return Ret;
  }

  Vector3f& operator += (const Vector3f& r){
    x += r.x;
    y += r.y;
    z += r.z;
    return *this;
  }

  Vector3f& operator -= (const Vector3f& r){
    x -= r.x;
    y -= r.y;
    z -= r.z;
    return *this;
  }

  Vector3f& operator *= (const Vector3f& r){
    x *= r.x;
    y *= r.y;
    z *= r.z;
    return *this;
  }

  operator const float*() const{
    return &(x);
  }

  Vector3f Cross(const Vector3f &v) const;
  Vector3f& Normalize();
  float Dot(const Vector3f &v) const;
  float Length();
  void Rotate(float Angle, const Vector3f& Axis);
  void Print() const {
    printf("(%.02f, %.02f, %.02f)\n", x, y ,z);
  }
};

struct Quaternion
{
  float x, y, z, w;
  Quaternion(float _x, float _y, float _z, float _q);
  Quaternion(Vector3f axis, float angle);
  void Normalize();
  Quaternion Conjugate();
  void ToDegrees(float * result);
};

class Matrix4f
{
public:
  float m[4][4];
  Matrix4f()
  {
    m[0][0] = 1.f; m[0][1] = 0.f; m[0][2] = 0.f; m[0][3] = 0.f;
    m[1][0] = 0.f; m[1][1] = 1.f; m[1][2] = 0.f; m[1][3] = 0.f;
    m[2][0] = 0.f; m[2][1] = 0.f; m[2][2] = 1.f; m[2][3] = 0.f;
    m[3][0] = 0.f; m[3][1] = 0.f; m[3][2] = 0.f; m[3][3] = 1.f;
  }
  inline void InitIdentity()
  {
    m[0][0] = 1.f; m[0][1] = 0.f; m[0][2] = 0.f; m[0][3] = 0.f;
    m[1][0] = 0.f; m[1][1] = 1.f; m[1][2] = 0.f; m[1][3] = 0.f;
    m[2][0] = 0.f; m[2][1] = 0.f; m[2][2] = 1.f; m[2][3] = 0.f;
    m[3][0] = 0.f; m[3][1] = 0.f; m[3][2] = 0.f; m[3][3] = 1.f;
  }
  inline Matrix4f operator*(const Matrix4f& Right) const
  {
    Matrix4f Ret;
    for (unsigned int i=0; i<4; i++){
      for (unsigned int j=0; j<4; j++){
        Ret.m[i][j] = m[i][0] * Right.m[0][j] +
                      m[i][1] * Right.m[1][j] +
                      m[i][2] * Right.m[2][j] +
                      m[i][3] * Right.m[3][j];
      }
    }
    return Ret;
  }

  inline Vector3f operator*(const Vector3f& Right) const
  {
    Vector3f Ret;
    Ret.x = m[0][0]*Right.x + m[0][1]*Right.y + m[0][2]*Right.z;
    Ret.y = m[1][0]*Right.x + m[1][1]*Right.y + m[1][2]*Right.z;
    Ret.z = m[2][0]*Right.x + m[2][1]*Right.y + m[2][2]*Right.z;
    return Ret;
  }

  inline Matrix4f operator*(const float s) const
  {
    Matrix4f Ret;
    for (unsigned int i=0; i<4; i++){
      for (unsigned int j=0; j<4; j++){
        Ret.m[i][j] = m[i][j] * s;
      }
    }
    return Ret;
  }

  inline Matrix4f operator+(const Matrix4f& Right) const
  {
    Matrix4f Ret;
    for (unsigned int i=0; i<4; i++){
      for (unsigned int j=0; j<4; j++){
        Ret.m[i][j] = m[i][j] + Right.m[i][j];
      }
    }
    return Ret;
  }


  void InitScaleTransform(float sx, float sy, float sz);
  void InitRotateTransform(float rx, float ry, float rz);
  void InitRotateAxisTransform(Vector3f u, float th);
  void InitRotateTransform(const Quaternion& quat);
  void InitTranslateTransform(float x, float y, float z);
  void InitCameraTransform(const float Target[3], const float Up[3]);
  void InitPersProjTransform(const PersProjInfo& p);
  void InitOrthoProjTransform(const OrthoProjInfo& p);
};

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

class Pipeline
{
public:
  Pipeline(){
    m_scale[0] = 1.f; m_scale[1] = 1.f, m_scale[2] = 1.f;
    m_pos[0] = 0.f; m_pos[1] = 0.f, m_pos[2] = 0.f;
    m_rotate[0] = 0.f; m_rotate[1] = 0.f, m_rotate[2] = 0.f;
    m_post_rotate[0] = 0.f; m_post_rotate[1] = 0.f, m_post_rotate[2] = 0.f;
  }

  void Scale(float x, float y, float z){
    m_scale[0] = x; m_scale[1] = y; m_scale[2] = z;
  }

  void Translate(float x, float y, float z){
    m_pos[0] = x; m_pos[1] = y; m_pos[2] = z;
  }

  void Rotate(float x, float y, float z){
    m_rotate_trans.InitRotateTransform(x, y, z);
    m_rotate[0] = x; m_rotate[1] = y; m_rotate[2] = z;
  }

  void RotateAxis(Vector3f u, float th){
    m_rotate_trans.InitRotateAxisTransform(u, th);
  }

  void PostRotate(float x, float y, float z){
    m_post_rotate_trans.InitRotateTransform(x, y, z);
    m_post_rotate[0] = x; m_post_rotate[1] = y; m_post_rotate[2] = z;
  }

  void SetCamera(float Pos[3], float Target[3], float Up[3]){
    memcpy(&m_camera.Pos, Pos, sizeof(float)*3);
    memcpy(&m_camera.Target, Target, sizeof(float)*3);
    memcpy(&m_camera.Up, Up, sizeof(float)*3);
  }

  void SetCamera(CameraInfo cam){
    memcpy(&m_camera.Pos, cam.Pos, sizeof(float)*3);
    memcpy(&m_camera.Target, cam.Target, sizeof(float)*3);
    memcpy(&m_camera.Up, cam.Up, sizeof(float)*3);
  }

  void SetPersProjInfo(float FOV, float Width, float Height, float zNear, float zFar){
    m_projInfo.FOV = FOV;
    m_projInfo.Width = Width;
    m_projInfo.Height = Height;
    m_projInfo.zNear = zNear;
    m_projInfo.zFar = zFar;
  }

  void SetOrthoProjInfo(const OrthoProjInfo& p){
    m_orthoInfo.Left = p.Left;
    m_orthoInfo.Right = p.Right;
    m_orthoInfo.Bottom = p.Bottom;
    m_orthoInfo.Top = p.Top;
    m_orthoInfo.zNear = p.zNear;
    m_orthoInfo.zFar = p.zFar;
  }

  void SetOrthoProjInfo(float l, float r, float b, float t, float n, float f){
    m_orthoInfo.Left = l;
    m_orthoInfo.Right = r;
    m_orthoInfo.Bottom = b;
    m_orthoInfo.Top = t;
    m_orthoInfo.zNear = n;
    m_orthoInfo.zFar = f;
  }

  const Matrix4f * GetProjTrans();
  const Matrix4f * GetViewTrans();
  const Matrix4f * GetWorldTrans();
  const Matrix4f * GetVPTrans();
  const Matrix4f * GetWVPTrans();

  void SetRotationMatrix(const Matrix4f _m);
  void SetRotationMatrix(const float _m[4][4]);
  void SetPostRotationMatrix(const Matrix4f _m);

private:
  float m_scale[3];
  float m_pos[3];
  float m_rotate[3];
  Matrix4f m_rotate_trans;
  Matrix4f m_transform;
  float m_post_rotate[3];
  Matrix4f m_post_rotate_trans;

  PersProjInfo m_projInfo;
  OrthoProjInfo m_orthoInfo;

  CameraInfo m_camera;

  Matrix4f m_WVPtransformation;
  Matrix4f m_VPtransformation;
  Matrix4f m_ProjTransformation;
  Matrix4f m_Vtransformation;
  Matrix4f m_Wtransformation;
  Matrix4f m_WVtransformation;
  Matrix4f m_WPtransformation;
};

