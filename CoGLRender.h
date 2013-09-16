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

#ifndef CoGLRender_H
#define CoGLRender_H

#if defined (__APPLE__) || defined(MACOSX)
#include <mac/glew.h>
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#include <OpenGL/OpenGL.h>
#else
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glut.h>
#endif

#ifdef USE_INTEROP
#include <cuda_gl_interop.h>
#endif

#include <thrust/host_vector.h>
#include <thrust/functional.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/discard_iterator.h>

#include <sys/time.h>
#include <stdio.h>

#include "piston/util/quaternion.h"
#include "CoGLSim.h"
using namespace piston;


//!************************************************************************************************************
// Class CoGLRender
//
// This class provides the OpenGL rendering infrastructure for visualizing the results of the simulation.  It 
// contains a pointer to an instance of the CoGLSim class, which performs the actual simulation.
//!************************************************************************************************************
class CoGLRender
{
    //-----------------------
    // METHODS
    //-----------------------
public:

    //! Constructor and destructor
    CoGLRender();
    ~CoGLRender();

    //! Rendering methods 
    void initialize();
    void display();
    void rotate(float a_angle1, float a_angle2);
    void zoom(float a_delta);

    //! Methods to initialize and update the simulation
    void initialize_simulation();
    void time_step();
    void time_simulation(int a_iters);

    //-----------------------
    // PROPERTIES
    //-----------------------
protected:

    //! Pointer to instance of CoGLSim class, which performs the simulation
    CoGLSim<double>* simulation;

    //! Simulation variables for initializations
    int max_timesteps, cur_timestep;
    double delta,h,as,ac,ao,d2,eo,tau,eta,temp,shear_a,bulk_a;
    double eta0_xx,eta0_yy,eta0_zz,e1_0,e3_0;

    //! Host vectors for simulation initial values
    thrust::host_vector<double> h_ux, h_uy, h_uz;
    thrust::host_vector<double> h_uxdt, h_uydt, h_uzdt;
    thrust::host_vector<double> h_uxx_applied;

    //! Rendering variables
    float rotation_matrix[16];
    Quaternion qrot;
    float max_camera_fov, zoom_pct;

    //! Vectors for rendering vertices, normals, and colors
    thrust::host_vector<float3> vertices;
    thrust::host_vector<float3> normals;
    thrust::host_vector<float4> colors;

    //! Vertex buffers for interop
    #ifdef USE_INTEROP
      GLuint vbo_buffers[3];  struct cudaGraphicsResource* vbo_resources[3];
    #endif

    //! Timing variables
    struct timeval begin, end, diff;
};

#endif
