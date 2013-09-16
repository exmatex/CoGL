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

#include <fstream>
#include <iostream>
#include <sstream>
#include <float.h>

#include "CoGLRender.h"

#define STRINGIZE(x) #x
#define STRINGIZE_VALUE_OF(x) STRINGIZE(x)


//------------------------------------------------------------------------------
// CoGLRender::CoGLRender()
//
// Constructor for CoGLRender class
//------------------------------------------------------------------------------ 
CoGLRender::CoGLRender()
{
}


//---------------------------------------------------------------------------
// CoGLRender::~CoGLRender()
//
// Destructor of CoGLRender class
//---------------------------------------------------------------------------
CoGLRender::~CoGLRender()
{
    #ifdef USE_INTEROP
      if (vbo_buffers[0])
      {
	for (int i=0; i<3; i++) cudaGraphicsUnregisterResource(vbo_resources[i]);
	for (int i=0; i<3; i++)
	{
	  glBindBuffer(1, vbo_buffers[i]);
	  glDeleteBuffers(1, &(vbo_buffers[i]));
	  vbo_buffers[i] = 0;
	}
      }     
    #else
      vertices.clear(); normals.clear(); colors.clear(); 
    #endif
}


//---------------------------------------------------------------------------
// CoGLRender::initialize()
//
// Initialize OpenGL settings, and call initialize_simulation
//---------------------------------------------------------------------------
void CoGLRender::initialize()
{
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);

    float white[] = { 0.5, 0.5, 0.5, 1.0 };
    float black[] = { 0.0, 0.0, 0.0, 1.0 };
    float light_pos[] = { 0.0, 0.0, 10.0, 1.0 };
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, white);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100);
    glLightfv(GL_LIGHT0, GL_AMBIENT, white);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, white);
    glLightfv(GL_LIGHT0, GL_SPECULAR, black);
    glLightfv(GL_LIGHT0, GL_POSITION, light_pos);

    glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, 1);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);

    #ifdef USE_INTEROP
      glewInit();
      cudaGLSetGLDevice(0);

      // initialize contour buffer objects
      glGenBuffers(3, vbo_buffers);
      for (int i=0; i<3; i++)
      {
        unsigned int buffer_size = DIM3*VERTICES_PER_CELL*sizeof(float4);
        glBindBuffer(GL_ARRAY_BUFFER, vbo_buffers[i]);
        glBufferData(GL_ARRAY_BUFFER, buffer_size, 0, GL_DYNAMIC_DRAW);
      }

      glBindBuffer(GL_ARRAY_BUFFER, 0);
      for (int i=0; i<3; i++) cudaGraphicsGLRegisterBuffer(&(vbo_resources[i]), vbo_buffers[i], cudaGraphicsMapFlagsWriteDiscard);
    #endif

    initialize_simulation();

    rotate(0.2, -0.2);
    max_camera_fov = 40.0;  zoom_pct = 0.5;
}


//---------------------------------------------------------------------------
// CoGLRender::display()
//
// Render the structured grid, with colors based on the deviatoric strain
//---------------------------------------------------------------------------
void CoGLRender::display()
{
    time_step();

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(max_camera_fov*zoom_pct, 1.0, 1.0, 1000.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0, 0, 8.0f*DIM, 0, 0, 0, 0, 1, 0);
    glPushMatrix();

    glTranslatef(-DIM/2.0f, -DIM/2.0f, -DIM/2.0f);
    qrot.getRotMat(rotation_matrix);
    glMultMatrixf(rotation_matrix);

    float center_x = DIM/2.0f;  float center_y = DIM/2.0f;  float center_z = DIM/2.0f;
    GLfloat matrix[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
    float offset_x = matrix[0]*center_x + matrix[1]*center_y + matrix[2]*center_z;
    float offset_y = matrix[4]*center_x + matrix[5]*center_y + matrix[6]*center_z;
    float offset_z = matrix[8]*center_x + matrix[9]*center_y + matrix[10]*center_z;
    offset_x = center_x - offset_x; offset_y = center_y - offset_y; offset_z = center_z - offset_z;
    glTranslatef(-offset_x, -offset_y, -offset_z);

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
  
    #ifdef USE_INTEROP
      glBindBuffer(GL_ARRAY_BUFFER, vbo_buffers[0]);     
      glVertexPointer(3, GL_FLOAT, 0, 0);     
      glBindBuffer(GL_ARRAY_BUFFER, vbo_buffers[1]);
      glColorPointer(4, GL_FLOAT, 0, 0);
      glBindBuffer(GL_ARRAY_BUFFER, vbo_buffers[2]);
      glNormalPointer(GL_FLOAT, 0, 0); 
      glDrawArrays(GL_QUADS, 0, simulation->n_vertices);   
      glBindBuffer(GL_ARRAY_BUFFER, 0);
    #else       
      if (cur_timestep < max_timesteps) colors.assign(simulation->colors_begin(), simulation->colors_end());
          
      glNormalPointer(GL_FLOAT, 0, &normals[0]);
      glColorPointer(4, GL_FLOAT, 0, &colors[0]);
      glVertexPointer(3, GL_FLOAT, 0, &vertices[0]);
      glDrawArrays(GL_QUADS, 0, simulation->n_vertices);   
    #endif

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);

    glPopMatrix();
}


//------------------------------------------------------------------------------
// CoGLRender::rotate()
//
// param    float a_angle1   Rotation about first axis
// param    float a_angle2   Rotation about second axis
//
// Rotate view based on input angles
//------------------------------------------------------------------------------
void CoGLRender::rotate(float a_angle1, float a_angle2)
{
    Quaternion new_rot_x;
    new_rot_x.setEulerAngles(a_angle1, 0.0, 0.0);
    qrot.mul(new_rot_x);

    Quaternion new_rot_y;
    new_rot_y.setEulerAngles(0.0, 0.0, a_angle2);
    qrot.mul(new_rot_y);
}


//------------------------------------------------------------------------------
// CoGLRender::zoom()
//
// param    float a_delta    Change to percentage of max FOV to use for field of view angle
//
// Increment the zoom level based on input parameter
//------------------------------------------------------------------------------
void CoGLRender::zoom(float a_delta)
{
    zoom_pct += a_delta;
    if (zoom_pct > 1.0) zoom_pct = 1.0;  if (zoom_pct < 0.0) zoom_pct = 0.0;
}


//---------------------------------------------------------------------------
// CoGLRender::initialize_simulation()
//
// Initialize simulation parameters, read input data, and pass to CoGLSim class
//---------------------------------------------------------------------------
void CoGLRender::initialize_simulation()
{
    simulation = new CoGLSim<double>(DIM);

    #ifdef USE_INTEROP
      for (int i=0; i<3; i++) simulation->vbo_resources[i] = vbo_resources[i];
    #endif

    simulation->create_mesh();

    #ifndef USE_INTEROP
      normals.assign(simulation->normals_begin(), simulation->normals_end());
      vertices.assign(simulation->vertices_begin(), simulation->vertices_end());
      colors.assign(simulation->colors_begin(), simulation->colors_end());
    #endif 
  
    h_ux.resize(DIM3);    h_uy.resize(DIM3);    h_uz.resize(DIM3);
    h_uxdt.resize(DIM3);  h_uydt.resize(DIM3);  h_uzdt.resize(DIM3);
    h_uxx_applied.resize(DIM3);
    
    max_timesteps=200000;  temp=255.0;  d2=1.0;
    eta=10.0;  delta=0.01;  h=1.0;  
    shear_a=28.0*pow(10.0,10);  bulk_a=14.0*pow(10.0,10);
    ao=1.97*pow(10.0,10);
    
    as = shear_a/ao;  ac = bulk_a/ao;
    printf("AS=%lf, AC=%lf\n", as, ac);

    tau=(temp-270.)/(295.-270.);
    eta0_xx=(3.795-3.756)/3.756;
    eta0_yy=eta0_xx;
    eta0_zz=(3.725-3.756)/3.756;
    e1_0=(eta0_xx+eta0_yy+eta0_zz)/sqrt(3.);
    e3_0=(eta0_xx+eta0_yy-2.0*eta0_zz)/sqrt(6.);
    ao=1.97*pow(10.0,10);
    eo=-0.5000*ac*e1_0/e3_0;
   
    char fuxFilename[1024], fuyFilename[1024], fuzFilename[1024];
    sprintf(fuxFilename, "%s/fort.31", STRINGIZE_VALUE_OF(DATA_DIRECTORY));
    sprintf(fuyFilename, "%s/fort.32", STRINGIZE_VALUE_OF(DATA_DIRECTORY));
    sprintf(fuzFilename, "%s/fort.33", STRINGIZE_VALUE_OF(DATA_DIRECTORY));

    FILE* fux = fopen(fuxFilename, "r");
    FILE* fuy = fopen(fuyFilename, "r");
    FILE* fuz = fopen(fuzFilename, "r");
    if ((!fux)) { printf("Initial conditions file not found\n"); exit(-1); }
    for (unsigned int i=0; i<DIM; i++)
      for (unsigned int j=0; j<DIM; j++)
        for (unsigned l=0; l<DIM; l++)
        {
          fscanf(fux, "%lf ", &(h_ux[i*DIM2+j*DIM+l]));
          fscanf(fuy, "%lf ", &(h_uy[i*DIM2+j*DIM+l]));
          fscanf(fuz, "%lf ", &(h_uz[i*DIM2+j*DIM+l]));
        }
    fclose(fux); fclose(fuy); fclose(fuz);

    thrust::fill(h_uxdt.begin(), h_uxdt.end(), 0.0);
    thrust::fill(h_uydt.begin(), h_uydt.end(), 0.0);
    thrust::fill(h_uzdt.begin(), h_uzdt.end(), 0.0);
    thrust::fill(h_uxx_applied.begin(), h_uxx_applied.end(), 0.0);

    simulation->initialize_simulation(h_ux, h_uy, h_uz, h_uxdt, h_uydt, h_uzdt, h_uxx_applied, ac, as, d2, delta, eo, eta, h, tau);

    printf("Initial conditions set\n");
    cur_timestep=0;
}


//---------------------------------------------------------------------------
// CoGLRender::time_step()
//
// Advance the simulation by calling the CoGLSim class, and report timings
//---------------------------------------------------------------------------
void CoGLRender::time_step()
{
    if (cur_timestep == 0) gettimeofday(&begin, 0);
    
    if (cur_timestep < max_timesteps)
    {
      simulation->advance_simulation(); 
      simulation->color_mesh();
      cur_timestep++;
    }

    if (cur_timestep == max_timesteps)
    {
      gettimeofday(&end, 0);
      timersub(&end, &begin, &diff);
      float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
      std::cout << "Time: " << seconds << std::endl;
    }

}


//---------------------------------------------------------------------------
// CoGLRender::time_simulation()
//
// param    int a_iters    Number of iterations to time   
//
// Time the simulation for a given number of iterations, without rendering
//---------------------------------------------------------------------------
void CoGLRender::time_simulation(int a_iters)
{
    #ifdef USE_INTEROP
      printf("You must disable interop in CMake configuration in order to run timing tests without rendering\n");
    #else
      initialize_simulation();
      cur_timestep = 0;
      max_timesteps = a_iters;
      for (unsigned int i=0; i<max_timesteps; i++) time_step();
    #endif
}



