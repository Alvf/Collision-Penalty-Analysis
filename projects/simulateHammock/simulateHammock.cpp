#include <iostream>
#include <random>
#include "util/FFMPEG_MOVIE.h"

#include "Scenes/Shell/HAMMOCK_SCENE.h"
#include "util/DRAW_GL.h"

#include <snapshot.h>

using namespace HOBAK;
using namespace std;

GLVU glvu;
FFMPEG_MOVIE movie;

bool animate = false;
bool singleStep = false;
int timestep = -1;

// animation should pause on this frame;
int pauseFrame = -1;

// scene being simulated
SIMULATION_SCENE* scene = NULL;

///////////////////////////////////////////////////////////////////////
// GL and GLUT callbacks
///////////////////////////////////////////////////////////////////////
void glutDisplay()
{
  glvu.BeginFrame();
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    scene->drawScene();
    if (animate)
      movie.addFrameGL();
  glvu.EndFrame();
}

///////////////////////////////////////////////////////////////////////
// animate and display new result
///////////////////////////////////////////////////////////////////////
void glutIdle()
{
  if (animate) 
  {
    scene->stepSimulation();
    // recordShellFrameToJSON(*scene);
    if (singleStep) 
    {
      animate = false;
      singleStep = false;
    }
    timestep++;
  }
  
  if (timestep == scene->pauseFrame() && animate)
  {
    cout << " Hit pause frame specific in scene file " << endl;
    animate = false;
  }

  glutPostRedisplay();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutKeyboard(unsigned char key, int x, int y)
{
  switch (key) {
  case 27:  // escape key
  case 'Q':
    TIMER::printTimings();
    exit(0);
    break;
  case 'q':
    TIMER::printTimings();
    movie.streamWriteMovie(scene->movieName().c_str());
    exit(0);
    break;
  case 'a':
    animate = !animate;
    break;
  case 'd':
    scene->drawFeature() = !scene->drawFeature();
    break;
  case ' ':
    animate = true;
    singleStep = true;
    TIMER::printTimings();
    break;
  case 'v': {
      Camera *camera = glvu.GetCurrentCam();
      glvuVec3f eye;
      glvuVec3f ref;
      glvuVec3f up;
      camera->GetLookAtParams(&eye, &ref, &up);
      cout << "    _eye    = VECTOR3(" << eye[0] << ", " << eye[1] << ", " << eye[2] << ");" << endl;
      cout << "    _lookAt = VECTOR3(" << ref[0] << ", " << ref[1] << ", " << ref[2] << ");" << endl;
      cout << "    _up     = VECTOR3(" << up[0] << ", " << up[1] << ", " << up[2] << ");" << endl;
    } 
    break;
  default:
    break;
  }

  glvu.Keyboard(key, x, y);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutSpecial(int key, int x, int y)
{
  switch (key) {
  case GLUT_KEY_LEFT:
    scene->leftArrow() = true;
    break;
  case GLUT_KEY_RIGHT:
    scene->rightArrow() = true;
    break;
  case GLUT_KEY_UP:
    scene->arrowCounter()++;
    cout << "Arrow counter: " << scene->arrowCounter() << endl;
    break;
  case GLUT_KEY_DOWN:
    scene->arrowCounter()--;
    cout << "Arrow counter: " << scene->arrowCounter() << endl;
    break;
  case GLUT_KEY_HOME:
    break;
  case GLUT_KEY_END:
    break;
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutMouseClick(int button, int state, int x, int y) 
{ 
  glvu.Mouse(button, state, x, y);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutMouseMotion(int x, int y) 
{ 
  glvu.Motion(x, y);
}

//////////////////////////////////////////////////////////////////////////////
// open the GLVU window
//////////////////////////////////////////////////////////////////////////////
int glvuWindow()
{
  const string title = string("Simulating scene: ") + scene->sceneName();
  char buffer[title.length() + 1];
  strcpy(buffer, title.c_str());
  glvu.Init(buffer, GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_MULTISAMPLE, 0, 0, 800, 800);

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  GLfloat lightZeroPosition[] = { 10.0, 4.0, 10.0, 1.0 };
  GLfloat lightZeroColor[] = { 0.2, 0.2, 0.1, 1.0 };
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, 1);
  glLightfv(GL_LIGHT0, GL_POSITION, lightZeroPosition);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor);

  GLfloat lightPosition1[] = { -4.0, 1.0, 1.0, 1.0 };
  GLfloat lightColor1[] = { 0.2, 0.0, 0.0, 1.0 };
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, 1);
  glLightfv(GL_LIGHT1, GL_POSITION, lightPosition1);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, lightColor1);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  glEnable(GL_COLOR_MATERIAL);
  glShadeModel(GL_SMOOTH);
  glClearColor(1, 1, 1, 0);

  glutDisplayFunc(&glutDisplay);
  glutIdleFunc(&glutIdle);
  glutKeyboardFunc(&glutKeyboard);
  glutSpecialFunc(&glutSpecial);
  glutMouseFunc(&glutMouseClick);
  glutMotionFunc(&glutMouseMotion);

  glEnable(GL_MULTISAMPLE);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glvuVec3f ModelMin(-10, -10, -10), ModelMax(10, 10, 10);
  glvuVec3f Eye, LookAtCntr, Up, WorldCenter;

  const VECTOR3& eye    = scene->eye();
  const VECTOR3& lookAt = scene->lookAt();
  const VECTOR3& up     = scene->up();
  const VECTOR3& worldCenter = scene->worldCenter();

  for (unsigned int x = 0; x < 3; x++)
  {
    Eye[x] = eye[x];
    LookAtCntr[x] = lookAt[x];
    Up[x] = up[x];
    WorldCenter[x] = worldCenter[x];
  }

  float Yfov = 45;
  float Aspect = 1;
  float Near = 0.001f;
  float Far = 10.0f;
  glvu.SetAllCams(ModelMin, ModelMax, Eye, LookAtCntr, Up, Yfov, Aspect, Near, Far);
  glvu.SetWorldCenter(WorldCenter);

  glutMainLoop();

  // Control flow will never reach here
  return EXIT_SUCCESS;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  scene = new HAMMOCK_SCENE();

  // make sure to print out what we should expect from this scene
  scene->printSceneDescription();

  bool success = scene->buildScene();

  // prepare scene for JSON output
  // initializeShellSceneJSON(*scene);

  if (!success)
  {
    cout << " Failed to build scene-> Exiting ..." << endl;
    return 1;
  }


  glutInit(&argc, argv);
  glvuWindow();

  return 0;
}
