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

#ifndef CoGLSim_H
#define CoGLSim_H

#include <piston/piston_math.h>
#include <piston/choose_container.h>

using namespace std;

#define DIM 64
#define DIM2 DIM*DIM
#define DIM3 DIM*DIM*DIM
#define VERTICES_PER_CELL 24

//#define PRECOMPUTE_INDICES
//#define OPTIMIZE_KERNEL_DIVISIONS
//#define OPTIMIZE_KERNEL_BRANCHES
#define PRINT_SIM_OUTPUT


#ifdef PRECOMPUTE_INDICES
#define INDICES indices.begin(), indices.begin()+n_cells
#else
#define INDICES CountingIterator(0), CountingIterator(0)+n_cells
#endif


#ifdef PRECOMPUTE_INDICES
#define COMPUTE_INDICES \
        const int i = index >> 16; \
	const int j = (index << 16) >> 24; \
	const int l = (index << 24) >> 24;
#define COMPUTE_INDEX \
        const int i = index >> 16;
#else
#ifdef OPTIMIZE_KERNEL_DIVISIONS
#define COMPUTE_INDICES \
        unsigned divl = ((unsigned long long)index * reciprocal) >> 35; \
	const int l = index - divl*n_dim; \
	int num = index*n_dim_r; \
	unsigned divj = ((unsigned long long)num * reciprocal) >> 35; \
	const int j = num - divj*n_dim; \
	const int i = index*dim_sq_r;
#define COMPUTE_INDEX \
	const int i = index*dim_sq_r;        
#else
#define COMPUTE_INDICES \
	const int l = index % (n_dim); \
	const int j = (index / (n_dim)) % (n_dim); \
	const int i = index / (dim_sq);
#define COMPUTE_INDEX \
	const int i = index / (dim_sq);
#endif
#endif


namespace piston {


//!************************************************************************************************************
// Class CoGLSim
//
// This class uses data-parallel vectors and Thrust/PISTON to compute the Ginzburg-Landau simulation. 
//!************************************************************************************************************
template <typename T>
class CoGLSim
{
public:

    //-----------------------
    // TYPE DEFINITIONS
    //-----------------------

    //! Container types
    typedef typename thrust::device_vector<float3> VerticesContainer;
    typedef typename thrust::device_vector<float3> NormalsContainer;
    typedef typename thrust::device_vector<float4> ColorsContainer;
    typedef typename thrust::device_vector<float> ScalarsContainer;
    typedef typename thrust::device_vector<float> TableContainer;

    //! Iterator types
    typedef typename VerticesContainer::iterator VerticesIterator;
    typedef typename NormalsContainer::iterator NormalsIterator;
    typedef typename ColorsContainer::iterator ColorsIterator;
    typedef typename ScalarsContainer::iterator ScalarsIterator;
    typedef typename TableContainer::iterator TableIterator;
    typedef typename thrust::counting_iterator<int> CountingIterator;

    //-----------------------
    // PROPERTIES
    //-----------------------

    //! Constant static variables
    static const int vertices_per_box = VERTICES_PER_CELL;
    static const GLfloat box_vertices[vertices_per_box*3];
    static const GLfloat box_normals[vertices_per_box*3];

    //! Table data
    TableContainer box_verts;
    TableContainer box_norms;

    //! Vertices, normals, colors, and scalar values
    VerticesContainer vertices;
    NormalsContainer normals;
    ColorsContainer colors;
    ScalarsContainer scalars;

    //! Vertex buffers, including normals and colors
    float3 *vertex_buffer_data;
    float3 *normal_buffer_data;
    float4 *color_buffer_data;

    //! Variables for size dimensions
    int n_vertices, n_cells, n_dim;

    //! Data-parallel vectors for simulation intermediate computations
    thrust::device_vector<T> e1, e2, e3;
    thrust::device_vector<T> e4, e5, e6;
    thrust::device_vector<T> ux, uy, uz;
    thrust::device_vector<T> uxdt, uydt, uzdt;
    thrust::device_vector<T> gxx, gyy, gzz;
    thrust::device_vector<T> gxy, gyz, gxz;
    thrust::device_vector<T> g1, g2, g3;
    thrust::device_vector<T> lpe2, lpe3, lpvx, lpvy, lpvz;
    thrust::device_vector<T> sxxx, sxyy, sxzz;
    thrust::device_vector<T> sxyx, syyy, syzz;
    thrust::device_vector<T> sxzx, syzy, szzz;
    thrust::device_vector<T> etaxx, etayy, etazz;
    thrust::device_vector<T> etaxy, etaxz, etayz;
    thrust::device_vector<T> uxx, uyy, uzz;
    thrust::device_vector<T> uxy, uxz, uzy;
    thrust::device_vector<T> uyx, uyz, uzx;
    thrust::device_vector<T> uxx_applied;
    thrust::device_vector<T> uxxdt, uyydt, uzzdt;
    thrust::device_vector<T> uxydt, uxzdt, uzydt;
    thrust::device_vector<T> uyxdt, uyzdt, uzxdt;
    thrust::device_vector<int> indices;
    thrust::host_vector<int> h_indices;

    //! Simulation parameters
    T ac, as, d2, delta, eo, eta, h, tau;
    int sim_time;

    //! Vertex buffers used for interop
    #ifdef USE_INTEROP
      struct cudaGraphicsResource* vbo_resources[3];
    #endif

    //-----------------------
    // METHODS
    //-----------------------

    //------------------------------------------------------------------------------
    // CoGLSim::CoGLSim()
    //
    // Constructor for CoGLSim class
    //------------------------------------------------------------------------------ 
    CoGLSim(int a_n_dim) : n_dim(a_n_dim), box_verts((GLfloat*) box_vertices, (GLfloat*) box_vertices+vertices_per_box*3), 
                          box_norms((GLfloat*) box_normals, (GLfloat*) box_normals+vertices_per_box*3)
    {
      n_cells = n_dim*n_dim*n_dim;
      n_vertices = n_cells*vertices_per_box;

      printf("Parameters: %d %d %d %d\n", n_dim, n_cells, vertices_per_box, n_vertices);
    };


    //---------------------------------------------------------------------------
    // CoGLSim::initialize_simulation()
    //
    // param    thrust::host_vector<T> a_ux              Initial x positions
    // param    thrust::host_vector<T> a_uy              Initial y positions
    // param    thrust::host_vector<T> a_uz              Initial z positions
    // param    thrust::host_vector<T> a_uxdt            Initial x velocities
    // param    thrust::host_vector<T> a_uydt            Initial y velocities
    // param    thrust::host_vector<T> a_uzdt            Initial z velocities
    // param    thrust::host_vector<T> a_uxx_applied     Initial displacements
    // param    T a_ac                                   Simulation parameter
    // param    T a_as                                   Simulation parameter
    // param    T a_d2                                   Simulation parameter
    // param    T a_delta                                Simulation parameter
    // param    T a_eo                                   Simulation parameter
    // param    T a_eta                                  Simulation parameter
    // param    T a_h                                    Simulation parameter
    // param    T a_tau                                  Simulation parameter
    //
    // Initialize simulation parameters, read input data, and pass to CoGLSim class
    //---------------------------------------------------------------------------
    void initialize_simulation(thrust::host_vector<T> a_ux, thrust::host_vector<T> a_uy, thrust::host_vector<T> a_uz,
                              thrust::host_vector<T> a_uxdt, thrust::host_vector<T> a_uydt, thrust::host_vector<T> a_uzdt,
                              thrust::host_vector<T> a_uxx_applied, T a_ac, T a_as, T a_d2, T a_delta, T a_eo, 
                              T a_eta, T a_h, T a_tau)
    {
      ux = a_ux;  uy = a_uy;  uz = a_uz;  uxdt = a_uxdt;  uydt = a_uydt;  uzdt = a_uzdt;  uxx_applied = a_uxx_applied; 

      e1.resize(n_cells);     e2.resize(n_cells);     e3.resize(n_cells);
      e4.resize(n_cells);     e5.resize(n_cells);     e6.resize(n_cells);
      gxx.resize(n_cells);    gyy.resize(n_cells);    gzz.resize(n_cells);
      gxy.resize(n_cells);    gyz.resize(n_cells);    gxz.resize(n_cells);
      g1.resize(n_cells);     g2.resize(n_cells);     g3.resize(n_cells);
      lpe2.resize(n_cells);   lpe3.resize(n_cells);
      lpvx.resize(n_cells);   lpvy.resize(n_cells);   lpvz.resize(n_cells);
      sxxx.resize(n_cells);   sxyy.resize(n_cells);   sxzz.resize(n_cells);
      sxyx.resize(n_cells);   syyy.resize(n_cells);   syzz.resize(n_cells);
      sxzx.resize(n_cells);   syzy.resize(n_cells);   szzz.resize(n_cells);
      etaxx.resize(n_cells);  etayy.resize(n_cells);  etazz.resize(n_cells);
      etaxy.resize(n_cells);  etaxz.resize(n_cells);  etayz.resize(n_cells);
      uxx.resize(n_cells);    uyy.resize(n_cells);    uzz.resize(n_cells);
      uxy.resize(n_cells);    uxz.resize(n_cells);    uzy.resize(n_cells);
      uyx.resize(n_cells);    uyz.resize(n_cells);    uzx.resize(n_cells);
      uxxdt.resize(n_cells);  uyydt.resize(n_cells);  uzzdt.resize(n_cells);
      uxydt.resize(n_cells);  uxzdt.resize(n_cells);  uzydt.resize(n_cells);
      uyxdt.resize(n_cells);  uyzdt.resize(n_cells);  uzxdt.resize(n_cells);

      ac = a_ac;  as = a_as;  d2 = a_d2;  delta = a_delta;  eo = a_eo;  eta = a_eta;  h = a_h;  tau = a_tau;  
      scalars.resize(n_cells);  thrust::fill(scalars.begin(), scalars.end(), 1.0);
      sim_time = 0;

      #ifdef PRECOMPUTE_INDICES
        h_indices.resize(n_cells);
        for (unsigned int i=0; i<n_dim; i++)
        {
	  for (unsigned int j=0; j<n_dim; j++)
	  {
	    for (unsigned int l=0; l<n_dim; l++)
	    {
	      int index = i*n_dim*n_dim+j*n_dim+l;
	      int i_code = i << 16;
	      int j_code = j << 8;
	      h_indices[index] = i_code + j_code + l;
	    }
	  }
        }
        indices = h_indices;
        h_indices.clear();
      #endif
    }

    //---------------------------------------------------------------------------
    // CoGLSim::advance_simulation()
    //
    // Advance the simulation one time step using Ginzburg-Landau method
    //---------------------------------------------------------------------------
    void advance_simulation()
    {
      #ifdef PRINT_SIM_OUTPUT
        sim_time++;
        std::cout.precision(17);
        std::cout << "t=" << sim_time << " e2=" << std::scientific << e2[31*n_dim*n_dim+29*n_dim+39] << " " << e2[0*n_dim*n_dim] << " " << e2[1*n_dim*n_dim] << " " 
                                                                  << e2[30*n_dim*n_dim] << " " << e2[63*n_dim*n_dim] << std::endl;
      #endif

      const T sqrt2 = 1.41421356237;
      const T sqrt3 = 1.73205080757;
      const T sqrt6 = 2.44948974278;

      thrust::transform(INDICES, uxx.begin(), gradient_x(thrust::raw_pointer_cast(&*ux.begin()), n_dim, h));
      thrust::transform(INDICES, uxy.begin(), gradient_y(thrust::raw_pointer_cast(&*ux.begin()), n_dim, h));
      thrust::transform(INDICES, uxz.begin(), gradient_z(thrust::raw_pointer_cast(&*ux.begin()), n_dim, h));
      thrust::transform(INDICES, uyx.begin(), gradient_x(thrust::raw_pointer_cast(&*uy.begin()), n_dim, h));
      thrust::transform(INDICES, uyy.begin(), gradient_y(thrust::raw_pointer_cast(&*uy.begin()), n_dim, h));
      thrust::transform(INDICES, uyz.begin(), gradient_z(thrust::raw_pointer_cast(&*uy.begin()), n_dim, h));
      thrust::transform(INDICES, uzx.begin(), gradient_x(thrust::raw_pointer_cast(&*uz.begin()), n_dim, h));
      thrust::transform(INDICES, uzy.begin(), gradient_y(thrust::raw_pointer_cast(&*uz.begin()), n_dim, h));
      thrust::transform(INDICES, uzz.begin(), gradient_z(thrust::raw_pointer_cast(&*uz.begin()), n_dim, h));

      etayy = uyy;  etazz = uzz;
      thrust::transform(uxx.begin(), uxx.end(), uxx_applied.begin(), etaxx.begin(), thrust::plus<T>());
      thrust::transform(CountingIterator(0), CountingIterator(0)+n_cells, etaxy.begin(), etaf(thrust::raw_pointer_cast(&*uxy.begin()), thrust::raw_pointer_cast(&*uyx.begin()), h));
      thrust::transform(CountingIterator(0), CountingIterator(0)+n_cells, etaxz.begin(), etaf(thrust::raw_pointer_cast(&*uxz.begin()), thrust::raw_pointer_cast(&*uzx.begin()), h));
      thrust::transform(CountingIterator(0), CountingIterator(0)+n_cells, etayz.begin(), etaf(thrust::raw_pointer_cast(&*uyz.begin()), thrust::raw_pointer_cast(&*uzy.begin()), h));

      thrust::transform(CountingIterator(0), CountingIterator(0)+n_cells, e1.begin(),
                        sax3(thrust::raw_pointer_cast(&*etaxx.begin()), thrust::raw_pointer_cast(&*etayy.begin()), thrust::raw_pointer_cast(&*etazz.begin()), 1.0, 1.0, 1.0, sqrt3));
      thrust::transform(CountingIterator(0), CountingIterator(0)+n_cells, e2.begin(),
                        sax2(thrust::raw_pointer_cast(&*etaxx.begin()), thrust::raw_pointer_cast(&*etayy.begin()), 1.0, -1.0, sqrt2));
      thrust::transform(CountingIterator(0), CountingIterator(0)+n_cells, e3.begin(),
                        sax3(thrust::raw_pointer_cast(&*etaxx.begin()), thrust::raw_pointer_cast(&*etayy.begin()), thrust::raw_pointer_cast(&*etazz.begin()), 1.0, 1.0, -2.0, sqrt6));
      e6 = etaxy;  e5 = etaxz;  e4 = etayz;

      thrust::transform(INDICES, lpe2.begin(), laplacian(thrust::raw_pointer_cast(&*e2.begin()), n_dim, h));
      thrust::transform(INDICES, lpe3.begin(), laplacian(thrust::raw_pointer_cast(&*e3.begin()), n_dim, h));
      thrust::transform(INDICES, lpvx.begin(), laplacian(thrust::raw_pointer_cast(&*uxdt.begin()), n_dim, h));
      thrust::transform(INDICES, lpvy.begin(), laplacian(thrust::raw_pointer_cast(&*uydt.begin()), n_dim, h));
      thrust::transform(INDICES, lpvz.begin(), laplacian(thrust::raw_pointer_cast(&*uzdt.begin()), n_dim, h));

      thrust::transform(CountingIterator(0), CountingIterator(0)+n_cells, g1.begin(), g1f(thrust::raw_pointer_cast(&*e1.begin()), ac));
      thrust::transform(CountingIterator(0), CountingIterator(0)+n_cells, g2.begin(), g2f(thrust::raw_pointer_cast(&*e2.begin()), 
                        thrust::raw_pointer_cast(&*e3.begin()), thrust::raw_pointer_cast(&*lpe2.begin()), tau, d2));
      thrust::transform(CountingIterator(0), CountingIterator(0)+n_cells, g3.begin(), g3f(thrust::raw_pointer_cast(&*e2.begin()), 
                        thrust::raw_pointer_cast(&*e3.begin()), thrust::raw_pointer_cast(&*lpe3.begin()), tau, d2));
      thrust::transform(CountingIterator(0), CountingIterator(0)+n_cells, gxx.begin(), sax3(thrust::raw_pointer_cast(&*g1.begin()), thrust::raw_pointer_cast(&*g2.begin()),
                        thrust::raw_pointer_cast(&*g3.begin()), 1.0/sqrt3, 1.0/sqrt2, 1.0/sqrt6, 1.0));
      thrust::transform(CountingIterator(0), CountingIterator(0)+n_cells, gyy.begin(), sax3(thrust::raw_pointer_cast(&*g1.begin()), thrust::raw_pointer_cast(&*g2.begin()),
                        thrust::raw_pointer_cast(&*g3.begin()), 1.0/sqrt3, -1.0/sqrt2, 1.0/sqrt6, 1.0));
      thrust::transform(CountingIterator(0), CountingIterator(0)+n_cells, gzz.begin(), sax2(thrust::raw_pointer_cast(&*g1.begin()), thrust::raw_pointer_cast(&*g3.begin()),
                        1.0/sqrt3, -2.0/sqrt6, 1.0));
      thrust::transform(CountingIterator(0), CountingIterator(0)+n_cells, gxy.begin(), sax1(thrust::raw_pointer_cast(&*e6.begin()), as));
      thrust::transform(CountingIterator(0), CountingIterator(0)+n_cells, gxz.begin(), sax1(thrust::raw_pointer_cast(&*e5.begin()), as));
      thrust::transform(CountingIterator(0), CountingIterator(0)+n_cells, gyz.begin(), sax1(thrust::raw_pointer_cast(&*e4.begin()), as));

      thrust::transform(INDICES, sxxx.begin(), gradient_x(thrust::raw_pointer_cast(&*gxx.begin()), n_dim, h));
      thrust::transform(INDICES, syyy.begin(), gradient_y(thrust::raw_pointer_cast(&*gyy.begin()), n_dim, h));
      thrust::transform(INDICES, szzz.begin(), gradient_z(thrust::raw_pointer_cast(&*gzz.begin()), n_dim, h));
      thrust::transform(INDICES, sxyx.begin(), gradient_x(thrust::raw_pointer_cast(&*gxy.begin()), n_dim, h));
      thrust::transform(INDICES, sxyy.begin(), gradient_y(thrust::raw_pointer_cast(&*gxy.begin()), n_dim, h));
      thrust::transform(INDICES, sxzx.begin(), gradient_x(thrust::raw_pointer_cast(&*gxz.begin()), n_dim, h));
      thrust::transform(INDICES, sxzz.begin(), gradient_z(thrust::raw_pointer_cast(&*gxz.begin()), n_dim, h));
      thrust::transform(INDICES, syzz.begin(), gradient_z(thrust::raw_pointer_cast(&*gyz.begin()), n_dim, h));
      thrust::transform(INDICES, syzy.begin(), gradient_y(thrust::raw_pointer_cast(&*gyz.begin()), n_dim, h));   

      thrust::transform(CountingIterator(0), CountingIterator(0)+n_cells, uxdt.begin(), vff(thrust::raw_pointer_cast(&*uxdt.begin()), thrust::raw_pointer_cast(&*sxxx.begin()),
                        thrust::raw_pointer_cast(&*sxyy.begin()), thrust::raw_pointer_cast(&*sxzz.begin()), thrust::raw_pointer_cast(&*lpvx.begin()), delta, eta));
      thrust::transform(CountingIterator(0), CountingIterator(0)+n_cells, uydt.begin(), vff(thrust::raw_pointer_cast(&*uydt.begin()), thrust::raw_pointer_cast(&*sxyx.begin()),
                        thrust::raw_pointer_cast(&*syyy.begin()), thrust::raw_pointer_cast(&*syzz.begin()), thrust::raw_pointer_cast(&*lpvy.begin()), delta, eta));
      thrust::transform(CountingIterator(0), CountingIterator(0)+n_cells, uzdt.begin(), vff(thrust::raw_pointer_cast(&*uzdt.begin()), thrust::raw_pointer_cast(&*sxzx.begin()),
                        thrust::raw_pointer_cast(&*syzy.begin()), thrust::raw_pointer_cast(&*szzz.begin()), thrust::raw_pointer_cast(&*lpvz.begin()), delta, eta));

      thrust::transform(CountingIterator(0), CountingIterator(0)+n_cells, ux.begin(), sax2(thrust::raw_pointer_cast(&*ux.begin()), thrust::raw_pointer_cast(&*uxdt.begin()),
                        1.0, delta, 1.0));
      thrust::transform(CountingIterator(0), CountingIterator(0)+n_cells, uy.begin(), sax2(thrust::raw_pointer_cast(&*uy.begin()), thrust::raw_pointer_cast(&*uydt.begin()),
                        1.0, delta, 1.0));
      thrust::transform(CountingIterator(0), CountingIterator(0)+n_cells, uz.begin(), sax2(thrust::raw_pointer_cast(&*uz.begin()), thrust::raw_pointer_cast(&*uzdt.begin()),
                        1.0, delta, 1.0)); 
    }


    //---------------------------------------------------------------------------
    // CoGLSim::create_mesh()
    //
    // Create the simulation mesh (a simple cube)
    //---------------------------------------------------------------------------
    void create_mesh()
    {
      #ifdef USE_INTEROP       
        size_t num_bytes;
	cudaGraphicsMapResources(1, &vbo_resources[0], 0);
	cudaGraphicsResourceGetMappedPointer((void **)&vertex_buffer_data, &num_bytes, vbo_resources[0]);

	if (vbo_resources[1])
	{
	  cudaGraphicsMapResources(1, &vbo_resources[1], 0);
	  cudaGraphicsResourceGetMappedPointer((void **)&color_buffer_data, &num_bytes, vbo_resources[1]);
	}

	cudaGraphicsMapResources(1, &vbo_resources[2], 0);
        cudaGraphicsResourceGetMappedPointer((void **)&normal_buffer_data, &num_bytes, vbo_resources[2]);

        thrust::for_each(CountingIterator(0), CountingIterator(0)+n_cells, generate_box(n_dim, vertex_buffer_data, normal_buffer_data, color_buffer_data, 
                         thrust::raw_pointer_cast(&*box_verts.begin()), thrust::raw_pointer_cast(&*box_norms.begin())));

        for (int i=0; i<3; i++) cudaGraphicsUnmapResources(1, &vbo_resources[i], 0);    
      #else    
	normals.resize(n_vertices);
	vertices.resize(n_vertices);
	colors.resize(n_vertices);

        thrust::for_each(CountingIterator(0), CountingIterator(0)+n_cells, generate_box(n_dim, thrust::raw_pointer_cast(&*vertices.begin()),
	                 thrust::raw_pointer_cast(&*normals.begin()), thrust::raw_pointer_cast(&*colors.begin()), 
                         thrust::raw_pointer_cast(&*box_verts.begin()), thrust::raw_pointer_cast(&*box_norms.begin())));  
      #endif
    }


    //---------------------------------------------------------------------------
    // CoGLSim::color_mesh()
    //
    // Color map the simulation mesh according to the e2 deviatoric strain
    //---------------------------------------------------------------------------
    void color_mesh()
    {
        float min_value = *(thrust::min_element(e2.begin(), e2.end()));
        float max_value = *(thrust::max_element(e2.begin(), e2.end()));
        
        #ifdef USE_INTEROP     
          size_t num_bytes;
	  if (vbo_resources[1])
	  {
	    cudaGraphicsMapResources(1, &vbo_resources[1], 0);
	    cudaGraphicsResourceGetMappedPointer((void **)&color_buffer_data, &num_bytes, vbo_resources[1]);
	            
            thrust::for_each(CountingIterator(0), CountingIterator(0)+n_cells,
                             color_box(color_buffer_data, thrust::raw_pointer_cast(&*e2.begin()), min_value, max_value));
            cudaGraphicsUnmapResources(1, &vbo_resources[1], 0);
          }        
        #else    
          thrust::for_each(CountingIterator(0), CountingIterator(0)+n_cells,
                           color_box(thrust::raw_pointer_cast(&*colors.begin()), thrust::raw_pointer_cast(&*etaxx.begin()), min_value, max_value));
        #endif
    }


    //---------------------------------------------------------------------------
    // generate_box functor
    //
    // Create vertices and normals for a cell with given index
    //---------------------------------------------------------------------------
    struct generate_box : public thrust::unary_function<int, void>
    {
      int dim;
      float3* vertices;
      float3* normals;
      float4* colors;
      float *box_verts, *box_norms;

      __host__ __device__
      generate_box(int dim, float3* vertices, float3* normals, float4* colors, float* box_verts, float* box_norms) : 
                     dim(dim), vertices(vertices), normals(normals), colors(colors), box_verts(box_verts), box_norms(box_norms) {};

      __host__ __device__
      void operator() (int id) const
      {
	int offset = id*VERTICES_PER_CELL;
	int x = id % (dim);
        int y = (id / (dim)) % (dim);
	int z = id / ((dim)*(dim));

	float4 cresult;
	cresult.x = 0.0f; cresult.y = 0.0f; cresult.z = 0.8f; cresult.w = 1.0f;
	for (unsigned int i=0; i<VERTICES_PER_CELL; i++)
	{
	  float3 vresult, nresult;
	  vresult.x = box_verts[i*3+0] + 1.0f*x;
	  vresult.y = box_verts[i*3+1] + 1.0f*y;
	  vresult.z = box_verts[i*3+2] + 1.0f*z;
	  nresult.x = box_norms[i*3+0];
	  nresult.y = box_norms[i*3+1];
	  nresult.z = box_norms[i*3+2];
	  *(vertices+offset+i) = vresult;
	  *(normals+offset+i) = nresult;
	  *(colors+offset+i) = cresult;
	}
      }
    };


    //---------------------------------------------------------------------------
    // color_box functor
    //
    // Set vertex colors for cell with given index based on given scalar values and range
    //---------------------------------------------------------------------------
    struct color_box : public thrust::unary_function<int, void>
    {
      float4* colors;
      T* values;
      T min_value, max_value;

      __host__ __device__
      color_box(float4* colors, T* values, T min_value, T max_value) : colors(colors), values(values), min_value(min_value), max_value(max_value) {};

      __host__ __device__
      void operator() (int id) const
      {
    	int offset = id*VERTICES_PER_CELL;
    	float4 result;

        const float V = 0.7f, S = 1.0f;
	float H = (1.0f - static_cast<float> (values[id] - min_value) / (max_value - min_value));
	if (H < 0.0f) H = 0.0f;
	else if (H > 1.0f) H = 1.0f;
	H *= 4.0f;

	float i = floor(H);
	float f = H - i;

	float p = V * (1.0 - S);
	float q = V * (1.0 - S * f);
	float t = V * (1.0 - S * (1 - f));

	float R, G, B;
	if (i == 0.0) { R = V; G = t; B = p; } 
        else if (i == 1.0) { R = q; G = V; B = p; } 
        else if (i == 2.0) { R = p; G = V; B = t; } 
        else if (i == 3.0) { R = p; G = q; B = V; } 
        else if (i == 4.0) { R = t; G = p; B = V; } 
        else { R = V; G = p; B = q; }
        result.x = R; result.y = G; result.z = B; result.w = 1.0f;

    	for (unsigned int i=0; i<VERTICES_PER_CELL; i++)
    	  *(colors+offset+i) = result;
      }
    };

    
    //---------------------------------------------------------------------------
    // gradient_x functor
    //
    // Compute gradient in x dimension at given cell index
    //---------------------------------------------------------------------------
    struct gradient_x : public thrust::unary_function<int, T>
    {
      T* x;
      int n_dim, dim_sq;
      unsigned reciprocal;
      T hi, n_dim_r, dim_sq_r;

      gradient_x(T* x, int n_dim, T h) : x(x), n_dim(n_dim), hi(1.0/h), dim_sq(n_dim*n_dim)
      {
        #ifdef OPTIMIZE_KERNEL_DIVISIONS
          reciprocal = ((34359738368.0) / n_dim + 0.9999); n_dim_r = 1.0/n_dim;  dim_sq_r = 1.0/(n_dim*n_dim);
        #endif
      };

      __host__ __device__
      T operator() (int index) const
      {
        COMPUTE_INDICES

        #ifdef OPTIMIZE_KERNEL_BRANCHES
          int offset1 = i+1; offset1 = offset1 % n_dim;
          int offset2 = n_dim+i-1; offset2 = offset2 % n_dim;
        #else
          int offset1 = i+1; if (offset1 >= n_dim) offset1 = 0;
          int offset2 = i-1; if (offset2 < 0) offset2 = n_dim-1;
        #endif
	T grad = 0.5*hi*(x[offset1*dim_sq+j*n_dim+l]-x[offset2*dim_sq+j*n_dim+l]);
	return grad;
      }
    };


    //---------------------------------------------------------------------------
    // gradient_y functor
    //
    // Compute gradient in y dimension at given cell index
    //---------------------------------------------------------------------------
    struct gradient_y : public thrust::unary_function<int, T>
    {
      T* x;
      int n_dim, dim_sq;
      unsigned reciprocal;
      T hi, n_dim_r, dim_sq_r;

      gradient_y(T* x, int n_dim, T h) : x(x), n_dim(n_dim), hi(1.0/h), dim_sq(n_dim*n_dim)
      {
        #ifdef OPTIMIZE_KERNEL_DIVISIONS
          reciprocal = ((34359738368.0) / n_dim + 0.9999); n_dim_r = 1.0/n_dim;  dim_sq_r = 1.0/(n_dim*n_dim);
        #endif
      };

      __host__ __device__
      T operator() (int index) const
      {
	COMPUTE_INDICES

        #ifdef OPTIMIZE_KERNEL_BRANCHES
          int offset1 = j+1; offset1 = offset1 % n_dim;
          int offset2 = n_dim+j-1; offset2 = offset2 % n_dim;
        #else
    	  int offset1 = j+1; if (offset1 >= n_dim) offset1 = 0;
          int offset2 = j-1; if (offset2 < 0) offset2 = n_dim-1;
        #endif
    	T grad = 0.5*hi*(x[i*dim_sq+offset1*n_dim+l]-x[i*dim_sq+offset2*n_dim+l]);
        return grad;
      }
    };


    //---------------------------------------------------------------------------
    // gradient_z functor
    //
    // Compute gradient in z dimension at given cell index
    //---------------------------------------------------------------------------
    struct gradient_z : public thrust::unary_function<int, T>
    {
      T* x;
      int n_dim, dim_sq;
      unsigned reciprocal;
      T hi, n_dim_r, dim_sq_r;

      gradient_z(T* x, int n_dim, T h) : x(x), n_dim(n_dim), hi(1.0/h), dim_sq(n_dim*n_dim)
      {
        #ifdef OPTIMIZE_KERNEL_DIVISIONS
          reciprocal = ((34359738368.0) / n_dim + 0.9999); n_dim_r = 1.0/n_dim;  dim_sq_r = 1.0/(n_dim*n_dim);
        #endif
      };

      __host__ __device__
      T operator() (int index) const
      {
	COMPUTE_INDICES

        #ifdef OPTIMIZE_KERNEL_BRANCHES
          int offset1 = l+1; offset1 = offset1 % n_dim;
          int offset2 = n_dim+l-1; offset2 = offset2 % n_dim;
        #else
	  int offset1 = l+1; if (offset1 >= n_dim) offset1 = 0;
          int offset2 = l-1; if (offset2 < 0) offset2 = n_dim-1;
        #endif
	T grad = 0.5*hi*(x[i*dim_sq+j*n_dim+offset1]-x[i*dim_sq+j*n_dim+offset2]);
	return grad;
      }
    };


    //---------------------------------------------------------------------------
    // laplacian functor
    //
    // Compute laplacian at given cell index
    //---------------------------------------------------------------------------
    struct laplacian : public thrust::unary_function<int, T>
    {
      T* x;
      int n_dim, dim_sq;
      unsigned reciprocal;
      T hi_sq, n_dim_r, dim_sq_r;

      laplacian(T* x, int n_dim, T h) : x(x), n_dim(n_dim), hi_sq(1.0/(h*h)), dim_sq(n_dim*n_dim)
      {
        #ifdef OPTIMIZE_KERNEL_DIVISIONS
          reciprocal = ((34359738368.0) / n_dim + 0.9999); n_dim_r = 1.0/n_dim;  dim_sq_r = 1.0/(n_dim*n_dim);
        #endif
      };

      __host__ __device__
      T operator() (int index) const
      {
	COMPUTE_INDICES

        #ifdef OPTIMIZE_KERNEL_BRANCHES
	  int i1 = i+1; i1 = i1 % n_dim;   int i2 = n_dim+i-1; i2 = i2 % n_dim;
	  int j1 = j+1; j1 = j1 % n_dim;   int j2 = n_dim+j-1; j2 = j2 % n_dim;
	  int l1 = l+1; l1 = l1 % n_dim;   int l2 = n_dim+l-1; l2 = l2 % n_dim;
        #else
	  int i1 = i+1; if (i1 >= n_dim) i1 = 0;  int i2 = i-1; if (i2 <  0) i2 = n_dim-1;
	  int j1 = j+1; if (j1 >= n_dim) j1 = 0;  int j2 = j-1; if (j2 <  0) j2 = n_dim-1;
	  int l1 = l+1; if (l1 >= n_dim) l1 = 0;  int l2 = l-1; if (l2 <  0) l2 = n_dim-1;
        #endif

	T lap = hi_sq*(x[i1*dim_sq+j*n_dim+l]+x[i2*dim_sq+j*n_dim+l]+x[i*dim_sq+j1*n_dim+l]+x[i*dim_sq+j2*n_dim+l]+
                       x[i*dim_sq+j*n_dim+l1]+x[i*dim_sq+j*n_dim+l2]-6.0*x[i*dim_sq+j*n_dim+l]);
        return lap;
      }
    };


    //---------------------------------------------------------------------------
    // etaf functor
    //
    // Calculate intermediate values for Ginzburg-Landau simulation computations
    //---------------------------------------------------------------------------
    struct etaf : public thrust::unary_function<int, T>
    {
      T *x, *y;
      T hi;

      etaf(T* x, T* y, T h) : x(x), y(y), hi(1.0/h) {};

      __host__ __device__
      T operator() (int index) const
      {
    	T value = 0.5*hi*(x[index]+y[index]);
    	return value;
      }
    };


    //---------------------------------------------------------------------------
    // sax3 functor
    //
    // Calculate intermediate values for Ginzburg-Landau simulation computations
    //---------------------------------------------------------------------------
    struct sax3 : public thrust::unary_function<int, T>
    {
      T *x, *y, *z;
      T a, b, c, hi;

      sax3(T* x, T* y, T* z, T a, T b, T c, T h) : x(x), y(y), z(z), a(a), b(b), c(c), hi(1.0/h) {};

      __host__ __device__
      T operator() (int index) const
      {
        T value = hi*(a*x[index]+b*y[index]+c*z[index]);
        return value;
      }
    };


    //---------------------------------------------------------------------------
    // sax2 functor
    //
    // Calculate intermediate values for Ginzburg-Landau simulation computations
    //---------------------------------------------------------------------------
    struct sax2 : public thrust::unary_function<int, T>
    {
      T *x, *y;
      T a, b, hi;

      sax2(T* x, T* y, T a, T b, T h) : x(x), y(y), a(a), b(b), hi(1.0/h) {};

      __host__ __device__
      T operator() (int index) const
      {
        T value = hi*(a*x[index]+b*y[index]);
        return value;
      }
    };


    //---------------------------------------------------------------------------
    // sax1 functor
    //
    // Calculate intermediate values for Ginzburg-Landau simulation computations
    //---------------------------------------------------------------------------
    struct sax1 : public thrust::unary_function<int, T>
    {
      T *x;
      T a;

      sax1(T* x, T a) : x(x), a(a) {};

      __host__ __device__
      T operator() (int index) const
      {
        T value = a*x[index];
        return value;
      }
    };


    //---------------------------------------------------------------------------
    // g1f functor
    //
    // Calculate intermediate values for Ginzburg-Landau simulation computations
    //---------------------------------------------------------------------------
    struct g1f : public thrust::unary_function<int, T>
    {
      T *x;
      T AC;

      g1f(T* x, T AC) : x(x), AC(AC) {};

      __host__ __device__
      T operator() (int i) const
      {
        T value = AC*x[i];
        return value;
      }
    };


    //---------------------------------------------------------------------------
    // g2f functor
    //
    // Calculate intermediate values for Ginzburg-Landau simulation computations
    //---------------------------------------------------------------------------
    struct g2f : public thrust::unary_function<int, T>
    {
      T *x, *y, *z;
      T tau, D2;

      g2f(T* x, T* y, T* z, T tau, T D2) : 
               x(x), y(y), z(z), tau(tau), D2(D2) {};

      __host__ __device__
      T operator() (int i) const
      {
        T value = 2.*tau*x[i]-D2*z[i]-12.*x[i]*y[i]+4.*x[i]*(x[i]*x[i]+y[i]*y[i]);
        return value;
      }
    };


    //---------------------------------------------------------------------------
    // g3f functor
    //
    // Calculate intermediate values for Ginzburg-Landau simulation computations
    //---------------------------------------------------------------------------
    struct g3f : public thrust::unary_function<int, T>
    {
      T *x, *y, *z;
      T tau, D2;

      g3f(T* x, T* y, T* z, T tau, T D2) : 
               x(x), y(y), z(z), tau(tau), D2(D2) {};

      __host__ __device__
      T operator() (int i) const
      {
        T value = 2.*tau*y[i]-D2*z[i]+6.*(y[i]*y[i]-x[i]*x[i])+4.*y[i]*(x[i]*x[i]+y[i]*y[i]);
        return value;
      }
    };


    //---------------------------------------------------------------------------
    // vff functor
    //
    // Calculate intermediate values for Ginzburg-Landau simulation computations
    //---------------------------------------------------------------------------
    struct vff : public thrust::unary_function<int, T>
    {
      T *v, *w, *x, *y, *z;
      T delta, eta;

      vff(T* v, T* w, T* x, T* y, T* z, T delta, T eta) : v(v), w(w), x(x), y(y), z(z), delta(delta), eta(eta) {};

      __host__ __device__
      T operator() (int i) const
      {
        T value = v[i]+delta*(w[i]+x[i]+y[i]+eta*z[i]);
        return value;
      }
    };


    //---------------------------------------------------------------------------
    // Utility functions
    //
    // Return beginning and ending iterators for vertex, normal, and color vectors
    //---------------------------------------------------------------------------
    VerticesIterator vertices_begin()  { return vertices.begin(); }
    VerticesIterator vertices_end()    { return vertices.end();   }
    NormalsIterator normals_begin()    { return normals.begin();  }
    NormalsIterator normals_end()      { return normals.end();    }
    ColorsIterator colors_begin()      { return colors.begin(); }
    ColorsIterator colors_end()        { return colors.end();   }
};


//---------------------------------------------------------------------------
// CoGL::boxVertices
//
// Table of vertices for unit cube
//---------------------------------------------------------------------------
template <typename T>
const GLfloat CoGLSim<T>::box_vertices[12*6] =
                         {0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0,
                          1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0,
                          1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0,
                          0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0,
                          0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
                          1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0};


//---------------------------------------------------------------------------
// CoGL::boxNormals
//
// Table of normals for unit cube
//---------------------------------------------------------------------------
template <typename T>
const GLfloat CoGLSim<T>::box_normals[12*6] =
                        {0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
                         0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0, 0.0, -1.0,
                         1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                         -1.0, 0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0, 0.0,
                         0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
                         0.0, -1.0, 0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0};


}

#endif /* CoGLSim_H */
