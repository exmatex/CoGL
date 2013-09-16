/*
Copyright (c) 2013, Los Alamos National Security, LLC
All rights reserved.
Copyright 2013. Los Alamos National Security, LLC. This software was produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL),
which is operated by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S. Government has rights to use, reproduce, and distribute this software.

NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.

If software is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
·         Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
·         Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other
          materials provided with the distribution.
·         Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, LANL, the U.S. Government, nor the names of its contributors may be used
          to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "CoGLRender.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#define DEG_TO_RAD 3.14159/180.0


// Global variables
int mouse_prev_x, mouse_prev_y;
int mouse_buttons = 0;
CoGLRender* renderer;


//------------------------------------------------------------------------------
// mouse()
//
// param    int a_button    Active button index
// param    int a_state     Active button state
// param    int a_x         X coordinate of mouse location
// param    int a_y         Y coordinate of mouse location
//
// Set the mouse state variables
//------------------------------------------------------------------------------
void mouse(int a_button, int a_state, int a_x, int a_y)
{
    if (a_state == GLUT_DOWN) mouse_buttons |= 1<<a_button;
    else if (a_state == GLUT_UP) mouse_buttons = 0;

    mouse_prev_x = a_x;
    mouse_prev_y = a_y;
    glutPostRedisplay();
}


//------------------------------------------------------------------------------
// motion()
//
// param    int a_x         X coordinate of mouse location
// param    int a_y         Y coordinate of mouse location
//
// Update the view rotation or zoom based on mouse movement and mouse button
//------------------------------------------------------------------------------
void motion(int a_x, int a_y)
{
    float dx = a_x - mouse_prev_x;
    float dy = a_y - mouse_prev_y;

    if (mouse_buttons == 1) renderer->rotate(-0.2*dx*DEG_TO_RAD, -0.2*dy*DEG_TO_RAD);
    else if (mouse_buttons == 4) renderer->zoom(dy/1000.0);

    mouse_prev_x = a_x;
    mouse_prev_y = a_y;
    glutPostRedisplay();
}


//------------------------------------------------------------------------------
// display()
//
// Call the instance of the CoGLRender class to update and render the simulation
//------------------------------------------------------------------------------
void display()
{
    renderer->display();
    glutSwapBuffers();
}


//------------------------------------------------------------------------------
// idle()
//
// Idle callback fuction
//------------------------------------------------------------------------------
void idle()
{
    glutPostRedisplay();
}


//------------------------------------------------------------------------------
// initGL()
//
// param    int argc         Number of command-line arguments
// param    char** argv      Command-line arguments
//
// Set up the GLUT callbacks, initialize the simulation, and start the GLUT loop
//------------------------------------------------------------------------------
void initGL(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(1024, 1024);
    glutCreateWindow("CoGL");

    renderer->initialize();

    glutDisplayFunc(display);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutIdleFunc(idle);
    glutMainLoop();
}


//------------------------------------------------------------------------------
// main()
//
// param    int argc         Number of command-line arguments
// param    char** argv      Command-line arguments
//
// If called with a command-line argument, run and time the simulation the
// specified number of iterations; otherwise, initialize the main GLUT callback
// loop to iteratively update and render the simulation
//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
    renderer = new CoGLRender();
    if (argc > 1) renderer->time_simulation(atoi(argv[1]));
    else initGL(argc, argv);
    return 0;
}

