/*
Copyright (c) 2017 Alberto Otero de la Roza
<aoterodelaroza@gmail.com>, Robin Myhr <x@example.com>, Isaac
Visintainer <x@example.com>, Richard Greaves <x@example.com>, Ángel
Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
<victor@fluor.quimica.uniovi.es>.

critic2 is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

critic2 is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "imgui.h"
#include "imgui_impl_glfw_gl3.h"
#include "draw.h"
#include "critic2.h"
#include "geometry.h"
#include "settings.h"

// GUI global variables //
// Bond and atom resolutions (0 = coarse -> 3 = smooth)
static void drawstick(Pipeline *p, GLuint shad, const c_stick *s);
static void drawball(Pipeline *p, GLuint shad, const c_ball *b, float scal);

// draw all the scene elements
void draw_all_elements(Pipeline *p, GLuint shad, GLFWwindow* window){

  glEnableVertexAttribArray(0);

  if (settings.show_bonds){
    for (int i=0; i<nbond; i++){
      drawstick(p, shad, &(bond[i].s));
    }
  }
  if (settings.show_atoms){
    for (int i=0; i<nat; i++){
      drawball(p, shad, &(at[i].b), settings.atomsize);
    }
  }
  if (settings.show_cps){
    for (int i=0; i<ncritp; i++) {
      drawball(p, shad, &(critp[i].b), settings.cpsize);
    }
  }
  if (settings.show_cell){
    for (int i=0; i<cell_nstick; i++) {
      drawstick(p, shad, &(cell_s[i]));
    }
  }

  glDisableVertexAttribArray(0);
}

// draw a bond between 2 atoms defined in the bond struct
static void drawstick(Pipeline *p, GLuint shad, const c_stick *s){
  p->Scale(s->thick, s->thick, s->length);
  p->Translate(s->r2[0], s->r2[1], s->r2[2]);
  p->SetRotationMatrix(s->rot);

  GLuint gWorldLocation = glGetUniformLocation(shad, "gWorld");
  GLuint gWVPLocation = glGetUniformLocation(shad, "gWVP");
  GLuint vColorLocation = glGetUniformLocation(shad, "vColor");

  glUniformMatrix4fv(gWVPLocation, 1, GL_TRUE, (const GLfloat *)p->GetWVPTrans());
  glUniformMatrix4fv(gWorldLocation, 1, GL_TRUE, (const GLfloat *)p->GetWorldTrans());
  glUniform4fv(vColorLocation, 1, (const GLfloat *)&(s->rgb));

  glBindBuffer(GL_ARRAY_BUFFER, bufcylv[settings.bondresolution]);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufcyli[settings.bondresolution]);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
  glDrawElements(GL_TRIANGLES, 3*ncyli[settings.bondresolution], GL_UNSIGNED_INT, 0);
}

// Draw a ball, with optional scaling
static void drawball(Pipeline *p, GLuint shad, const c_ball *b, float scal = 1.0){
  p->Scale(b->rad * scal,b->rad * scal,b->rad * scal);
  p->Translate(b->r[0], b->r[1], b->r[2]);

  GLuint gWorldLocation = glGetUniformLocation(shad, "gWorld");
  GLuint gWVPLocation = glGetUniformLocation(shad, "gWVP");
  GLuint vColorLocation = glGetUniformLocation(shad, "vColor");

  glUniformMatrix4fv(gWVPLocation, 1, GL_TRUE, (const GLfloat *)p->GetWVPTrans());
  glUniformMatrix4fv(gWorldLocation, 1, GL_TRUE, (const GLfloat *)p->GetWorldTrans());
  glUniform4fv(vColorLocation, 1, b->rgb);

  glBindBuffer(GL_ARRAY_BUFFER, bufsphv[settings.atomresolution]);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufsphi[settings.atomresolution]);
  glDrawElements(GL_TRIANGLES, 3*nsphi[settings.atomresolution], GL_UNSIGNED_INT, 0);
}

// Mouse scroll callback
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
  if (!ImGui::GetIO().WantCaptureMouse)
    settings.cam_pos[2] += yoffset * box_xmaxlen * 0.2f;
  else
    ImGui_ImplGlfwGL3_ScrollCallback(window, xoffset, yoffset);
}

// process mouse input
void process_mouse_input(GLFWwindow* window, Matrix4f *rot){
  static bool firstpass = true;

  // c means for current loop, l means last loop, p means last pressed
  static int cLMB, cRMB, lLMB, lRMB;
  static double cMPosX, cMPosY, lMPosX, lMPosY, pMPosX, pMPosY;
  static Matrix4f lastRot;

  float camPanFactor = fabs(0.00115f * settings.cam_pos[2]);
  float camRotateFactor = 0.015f;

  if (!firstpass){
    lLMB = cLMB;
    lRMB = cRMB;
    lMPosX = cMPosX;
    lMPosY = cMPosY;
  }

  cLMB = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
  cRMB = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT);
  glfwGetCursorPos(window, &cMPosX, &cMPosY);

  if (firstpass){
    firstpass = false;
    return;
  }

  // Process mouse input
  if (!ImGui::GetIO().WantCaptureMouse) {
    if (cRMB == GLFW_PRESS){
      settings.cam_pos[0] -= camPanFactor * (cMPosX - lMPosX);
      settings.cam_pos[1] += camPanFactor * (cMPosY - lMPosY);
    }
    if (cLMB == GLFW_PRESS){
      if (lLMB != GLFW_PRESS){
	pMPosX = cMPosX;
	pMPosY = cMPosY;
	lastRot = *rot;
      } else {
	Vector3f curRotAxis = Vector3f((float)(cMPosX-pMPosX), (float)(pMPosY-cMPosY), 0);
	curRotAxis = curRotAxis.Cross(Vector3f(0, 0, 1));
	float curRotAng = curRotAxis.Length() * camRotateFactor;
	curRotAxis.Normalize();
 
	Matrix4f curRot;
	curRot.InitRotateAxisTransform(curRotAxis, curRotAng);
	*rot = curRot * lastRot;
      }
    }
  }
}
