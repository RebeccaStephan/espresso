//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \\file ColorGradientCollideSweepSinglePrecisionAVX.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 3247aa7395049ca5bfb69d34d55e45db19fa439c


#include <cmath>

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "ColorGradientCollideSweepSinglePrecisionAVX.h"


#include <immintrin.h>



#define FUNC_PREFIX

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wfloat-equal"
#   pragma GCC diagnostic ignored "-Wshadow"
#   pragma GCC diagnostic ignored "-Wconversion"
#   pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning push
#pragma warning( disable :  1599 )
#endif

using namespace std;

namespace walberla {
namespace pystencils {


namespace internal_196de4e60a60b4c0f3c29ed6234bf8bc {
static FUNC_PREFIX void colorgradientcollidesweepsingleprecisionavx_colorgradientcollidesweepsingleprecisionavx(float * RESTRICT const _data_force_a, float * RESTRICT const _data_force_b, float * RESTRICT  _data_pdfs_a, float * RESTRICT  _data_pdfs_b, float * RESTRICT const _data_phasefield, float * RESTRICT const _data_rho_a, float * RESTRICT const _data_rho_b, float * RESTRICT const _data_velocity, int64_t const _size_force_a_0, int64_t const _size_force_a_1, int64_t const _size_force_a_2, int64_t const _stride_force_a_1, int64_t const _stride_force_a_2, int64_t const _stride_force_a_3, int64_t const _stride_force_b_1, int64_t const _stride_force_b_2, int64_t const _stride_force_b_3, int64_t const _stride_pdfs_a_1, int64_t const _stride_pdfs_a_2, int64_t const _stride_pdfs_a_3, int64_t const _stride_pdfs_b_1, int64_t const _stride_pdfs_b_2, int64_t const _stride_pdfs_b_3, int64_t const _stride_phasefield_1, int64_t const _stride_phasefield_2, int64_t const _stride_rho_a_1, int64_t const _stride_rho_a_2, int64_t const _stride_rho_b_1, int64_t const _stride_rho_b_2, int64_t const _stride_velocity_1, int64_t const _stride_velocity_2, int64_t const _stride_velocity_3, float beta, float omega_shear_a, float omega_shear_b, float sigma)
{
#ifdef _OPENMP
   #pragma omp parallel
#endif
   {
      const float xi_2 = ((1.0f) / (omega_shear_a + omega_shear_b));
      const float xi_3 = omega_shear_a*omega_shear_b;
      const float xi_4 = xi_2*xi_3;
      const float xi_5 = xi_4*2.0f;
      const float xi_6 = omega_shear_a*-4.0f;
      const float xi_7 = xi_4*8.0f;
      const float xi_8 = xi_6 + xi_7;
      const float xi_9 = -xi_8;
      const float xi_12 = omega_shear_b*-4.0f;
      const float xi_13 = xi_12 + xi_7;
      const float xia_1_collide = xi_2;
      const float xi_14 = xi_3*xia_1_collide;
      const float xi_15 = xi_14*2.0f;
      const float xi_16 = xi_14*8.0f;
      const float xi_17 = xi_16 + xi_6;
      const float xi_18 = xi_12 + xi_16;
      const float rr_0_a_collide = 0.0f;
      const float xi_27 = rr_0_a_collide*0.25f;
      const float xib_1_collide = xi_2;
      const float xi_105 = xi_3*xib_1_collide;
      const float xi_106 = xi_105*2.0f;
      const float xi_107 = xi_105*8.0f;
      const float xi_108 = xi_107 + xi_6;
      const float xi_109 = xi_107 + xi_12;
      const float rr_0_b_collide = 0.0f;
      const float xi_118 = rr_0_b_collide*0.25f;
#ifdef _OPENMP
      #pragma omp for schedule(static)
#endif
      for (int64_t ctr_2 = 1; ctr_2 < _size_force_a_2 - 1; ctr_2 += 1)
      {
         for (int64_t ctr_1 = 1; ctr_1 < _size_force_a_1 - 1; ctr_1 += 1)
         {
            {
               for (int64_t ctr_0 = 1; ctr_0 < (int64_t)((_size_force_a_0 - 2) / (8)) * (8) + 1; ctr_0 += 8)
               {
                  const __m256 xi_212 = _mm256_mul_ps(_mm256_set_ps(0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1]));
                  const __m256 xi_213 = _mm256_mul_ps(_mm256_set_ps(0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1]));
                  const __m256 xi_214 = _mm256_mul_ps(_mm256_set_ps(0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1]));
                  const __m256 xi_215 = _mm256_mul_ps(_mm256_set_ps(-0.013888888888888888f,-0.013888888888888888f,-0.013888888888888888f,-0.013888888888888888f,-0.013888888888888888f,-0.013888888888888888f,-0.013888888888888888f,-0.013888888888888888f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1]));
                  const __m256 xi_216 = _mm256_mul_ps(_mm256_set_ps(0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1]));
                  const __m256 xi_217 = _mm256_mul_ps(_mm256_set_ps(0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1]));
                  const __m256 xi_218 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_214,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_215),xi_216),xi_217);
                  const __m256 xi_219 = _mm256_mul_ps(_mm256_set_ps(-0.055555555555555552f,-0.055555555555555552f,-0.055555555555555552f,-0.055555555555555552f,-0.055555555555555552f,-0.055555555555555552f,-0.055555555555555552f,-0.055555555555555552f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1]));
                  const __m256 xi_220 = _mm256_mul_ps(_mm256_set_ps(0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1]));
                  const __m256 xi_221 = _mm256_mul_ps(_mm256_set_ps(0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1]));
                  const __m256 xi_222 = _mm256_mul_ps(_mm256_set_ps(0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f,0.013888888888888888f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1]));
                  const __m256 xi_223 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_220,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_219),xi_221),xi_222);
                  const __m256 xi_224 = _mm256_mul_ps(_mm256_set_ps(0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1]));
                  const __m256 xi_225 = _mm256_mul_ps(_mm256_set_ps(0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1]));
                  const __m256 xi_226 = _mm256_add_ps(_mm256_mul_ps(xi_224,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_225);
                  const __m256 xi_227 = _mm256_mul_ps(_mm256_set_ps(0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0]));
                  const __m256 xi_228 = _mm256_mul_ps(xi_227,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_229 = _mm256_mul_ps(_mm256_set_ps(0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0]));
                  const __m256 xi_230 = _mm256_add_ps(xi_228,xi_229);
                  const __m256 xi_231 = _mm256_mul_ps(_mm256_set_ps(-0.22222222222222221f,-0.22222222222222221f,-0.22222222222222221f,-0.22222222222222221f,-0.22222222222222221f,-0.22222222222222221f,-0.22222222222222221f,-0.22222222222222221f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0]));
                  const __m256 xi_232 = _mm256_mul_ps(_mm256_set_ps(0.22222222222222221f,0.22222222222222221f,0.22222222222222221f,0.22222222222222221f,0.22222222222222221f,0.22222222222222221f,0.22222222222222221f,0.22222222222222221f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0]));
                  const __m256 xi_233 = _mm256_mul_ps(_mm256_set_ps(0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0]));
                  const __m256 xi_234 = _mm256_mul_ps(_mm256_set_ps(0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0]));
                  const __m256 xi_235 = _mm256_add_ps(_mm256_mul_ps(xi_233,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_234);
                  const __m256 xi_236 = _mm256_add_ps(_mm256_add_ps(xi_231,xi_232),xi_235);
                  const __m256 xi_237 = _mm256_add_ps(xi_230,xi_236);
                  const __m256 xi_238 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_212,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_213),xi_218),xi_223),xi_226),xi_237);
                  const __m256 xi_239 = _mm256_mul_ps(xi_238,xi_238);
                  const __m256 xi_240 = _mm256_mul_ps(_mm256_set_ps(0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1]));
                  const __m256 xi_241 = _mm256_mul_ps(xi_240,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_242 = _mm256_mul_ps(_mm256_set_ps(0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1]));
                  const __m256 xi_243 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_228,xi_229),xi_241),xi_242);
                  const __m256 xi_244 = _mm256_mul_ps(_mm256_set_ps(0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1]));
                  const __m256 xi_245 = _mm256_mul_ps(_mm256_set_ps(0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1]));
                  const __m256 xi_246 = _mm256_add_ps(_mm256_mul_ps(xi_244,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_245);
                  const __m256 xi_247 = _mm256_add_ps(_mm256_mul_ps(xi_213,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_212);
                  const __m256 xi_248 = _mm256_mul_ps(_mm256_set_ps(0.22222222222222221f,0.22222222222222221f,0.22222222222222221f,0.22222222222222221f,0.22222222222222221f,0.22222222222222221f,0.22222222222222221f,0.22222222222222221f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0]));
                  const __m256 xi_249 = _mm256_mul_ps(xi_248,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_250 = _mm256_mul_ps(_mm256_set_ps(0.22222222222222221f,0.22222222222222221f,0.22222222222222221f,0.22222222222222221f,0.22222222222222221f,0.22222222222222221f,0.22222222222222221f,0.22222222222222221f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0]));
                  const __m256 xi_251 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_234,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_233),xi_249),xi_250);
                  const __m256 xi_252 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_222,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_218),xi_220),xi_243),xi_246),xi_247),xi_251);
                  const __m256 xi_253 = _mm256_mul_ps(xi_252,xi_252);
                  const __m256 xi_254 = _mm256_add_ps(_mm256_mul_ps(xi_245,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_244);
                  const __m256 xi_255 = _mm256_add_ps(xi_241,xi_242);
                  const __m256 xi_256 = _mm256_mul_ps(_mm256_set_ps(0.22222222222222221f,0.22222222222222221f,0.22222222222222221f,0.22222222222222221f,0.22222222222222221f,0.22222222222222221f,0.22222222222222221f,0.22222222222222221f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1]));
                  const __m256 xi_257 = _mm256_mul_ps(xi_256,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_258 = _mm256_mul_ps(_mm256_set_ps(0.22222222222222221f,0.22222222222222221f,0.22222222222222221f,0.22222222222222221f,0.22222222222222221f,0.22222222222222221f,0.22222222222222221f,0.22222222222222221f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1]));
                  const __m256 xi_259 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_225,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_224),xi_257),xi_258);
                  const __m256 xi_260 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_217,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_214),xi_215),xi_216),xi_223),xi_247),xi_254),xi_255),xi_259);
                  const __m256 xi_261 = _mm256_mul_ps(xi_260,xi_260);
                  const __m256 xi_262 = _mm256_add_ps(_mm256_add_ps(xi_239,xi_253),xi_261);
                  const __m256 xi_263 = _mm256_sqrt_ps(xi_262);
                  const __m256 xi_264 = _mm256_mul_ps(_mm256_load_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]),_mm256_load_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]));
                  const __m256 xi_265 = _mm256_mul_ps(_mm256_mul_ps(xi_263,_mm256_set_ps(sigma,sigma,sigma,sigma,sigma,sigma,sigma,sigma)),_mm256_blendv_ps(_mm256_blendv_ps(_mm256_blendv_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_264,_mm256_set_ps(xi_13,xi_13,xi_13,xi_13,xi_13,xi_13,xi_13,xi_13)),_mm256_mul_ps(_mm256_set_ps(xi_13,xi_13,xi_13,xi_13,xi_13,xi_13,xi_13,xi_13),_mm256_load_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]))),_mm256_set_ps(xi_5,xi_5,xi_5,xi_5,xi_5,xi_5,xi_5,xi_5)),_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_264,_mm256_set_ps(xi_8,xi_8,xi_8,xi_8,xi_8,xi_8,xi_8,xi_8)),_mm256_mul_ps(_mm256_set_ps(xi_9,xi_9,xi_9,xi_9,xi_9,xi_9,xi_9,xi_9),_mm256_load_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]))),_mm256_set_ps(xi_5,xi_5,xi_5,xi_5,xi_5,xi_5,xi_5,xi_5)),_mm256_cmp_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_load_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]),_CMP_NGE_UQ)),_mm256_set_ps(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b),_mm256_cmp_ps(_mm256_set_ps(-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f),_mm256_load_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]),_CMP_NLE_UQ)),_mm256_set_ps(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a),_mm256_cmp_ps(_mm256_set_ps(0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f),_mm256_load_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]),_CMP_NGE_UQ)));
                  const __m256 xi_266 = _mm256_cmp_ps(xi_263,_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_CMP_NLE_UQ);
                  const __m256 xi_267 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_265,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f)),xi_266);
                  const __m256 xi_279 = _mm256_div_ps(_mm256_set_ps(1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f),xi_262);
                  const __m256 xi_280 = _mm256_mul_ps(xi_265,_mm256_set_ps(1.125f,1.125f,1.125f,1.125f,1.125f,1.125f,1.125f,1.125f));
                  const __m256 xi_281 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_280,_mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(xi_253,xi_279),_mm256_set_ps(0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f)),_mm256_set_ps(-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f))),xi_266);
                  const __m256 xi_285 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_280,_mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(xi_279,_mm256_set_ps(0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f)),_mm256_mul_ps(xi_252,xi_252)),_mm256_set_ps(-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f))),xi_266);
                  const __m256 xi_293 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_280,_mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(xi_279,_mm256_set_ps(0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f)),_mm256_mul_ps(xi_260,xi_260)),_mm256_set_ps(-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f))),xi_266);
                  const __m256 xi_296 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_280,_mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(xi_261,xi_279),_mm256_set_ps(0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f)),_mm256_set_ps(-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f))),xi_266);
                  const __m256 xi_302 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_280,_mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(xi_279,_mm256_set_ps(0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f)),_mm256_mul_ps(xi_238,xi_238)),_mm256_set_ps(-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f))),xi_266);
                  const __m256 xi_304 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_280,_mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(xi_239,xi_279),_mm256_set_ps(0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f)),_mm256_set_ps(-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f))),xi_266);
                  const __m256 xi_312 = _mm256_mul_ps(_mm256_set_ps(0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1]));
                  const __m256 xi_313 = _mm256_mul_ps(_mm256_set_ps(0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1]));
                  const __m256 xi_314 = _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1])),_mm256_mul_ps(_mm256_set_ps(-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1])));
                  const __m256 xi_315 = _mm256_add_ps(xi_219,xi_221);
                  const __m256 xi_316 = _mm256_add_ps(xi_259,xi_315);
                  const __m256 xi_317 = _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1])),_mm256_mul_ps(_mm256_set_ps(-0.027777777777777776f,-0.027777777777777776f,-0.027777777777777776f,-0.027777777777777776f,-0.027777777777777776f,-0.027777777777777776f,-0.027777777777777776f,-0.027777777777777776f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1])));
                  const __m256 xi_318 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_250,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_248),xi_317);
                  const __m256 xi_319 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_229,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_312,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_227),xi_235),xi_313),xi_314),xi_316),xi_318);
                  const __m256 xi_320 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_280,_mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(xi_279,_mm256_set_ps(0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f)),_mm256_mul_ps(xi_319,xi_319)),_mm256_set_ps(-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f))),xi_266);
                  const __m256 xi_324 = _mm256_mul_ps(_mm256_set_ps(0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1]));
                  const __m256 xi_325 = _mm256_mul_ps(_mm256_set_ps(0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1]));
                  const __m256 xi_326 = _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1])),_mm256_mul_ps(_mm256_set_ps(-0.027777777777777776f,-0.027777777777777776f,-0.027777777777777776f,-0.027777777777777776f,-0.027777777777777776f,-0.027777777777777776f,-0.027777777777777776f,-0.027777777777777776f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1])));
                  const __m256 xi_327 = _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1])),_mm256_mul_ps(_mm256_set_ps(-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1])));
                  const __m256 xi_328 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_324,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_230),xi_251),xi_316),xi_325),xi_326),xi_327);
                  const __m256 xi_329 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_280,_mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(xi_279,_mm256_set_ps(0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f)),_mm256_mul_ps(xi_328,xi_328)),_mm256_set_ps(-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f))),xi_266);
                  const __m256 xi_330 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_280,_mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(xi_279,_mm256_set_ps(0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f)),_mm256_mul_ps(xi_328,xi_328)),_mm256_set_ps(-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f))),xi_266);
                  const __m256 xi_332 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_280,_mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(xi_279,_mm256_set_ps(0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f)),_mm256_mul_ps(xi_319,xi_319)),_mm256_set_ps(-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f))),xi_266);
                  const __m256 xi_335 = _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0])),_mm256_mul_ps(_mm256_set_ps(-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0])));
                  const __m256 xi_336 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_226,xi_231),xi_232),xi_315);
                  const __m256 xi_337 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_242,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_325,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_240),xi_324);
                  const __m256 xi_338 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_254,xi_318),xi_335),xi_336),xi_337);
                  const __m256 xi_339 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_280,_mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(xi_279,_mm256_set_ps(0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f)),_mm256_mul_ps(xi_338,xi_338)),_mm256_set_ps(-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f))),xi_266);
                  const __m256 xi_345 = _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0])),_mm256_mul_ps(_mm256_set_ps(-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0])));
                  const __m256 xi_346 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_313,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_246),xi_312);
                  const __m256 xi_347 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_249,xi_250),xi_255),xi_326),xi_336),xi_345),xi_346);
                  const __m256 xi_348 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_280,_mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(xi_279,_mm256_set_ps(0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f)),_mm256_mul_ps(xi_347,xi_347)),_mm256_set_ps(-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f))),xi_266);
                  const __m256 xi_353 = _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1])),_mm256_mul_ps(_mm256_set_ps(-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1])));
                  const __m256 xi_354 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_236,xi_243),xi_254),xi_257),xi_258),xi_317),xi_326),xi_353);
                  const __m256 xi_355 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_280,_mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(xi_279,_mm256_set_ps(0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f)),_mm256_mul_ps(xi_354,xi_354)),_mm256_set_ps(-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f))),xi_266);
                  const __m256 xi_359 = _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f,0.1111111111111111f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1])),_mm256_mul_ps(_mm256_set_ps(-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f),_mm256_loadu_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1])));
                  const __m256 xi_360 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_258,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_237),xi_256),xi_337),xi_346),xi_359);
                  const __m256 xi_361 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_280,_mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(xi_279,_mm256_set_ps(0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f)),_mm256_mul_ps(xi_360,xi_360)),_mm256_set_ps(-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f))),xi_266);
                  const __m256 xi_362 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_280,_mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(xi_279,_mm256_set_ps(0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f)),_mm256_mul_ps(xi_347,xi_347)),_mm256_set_ps(-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f))),xi_266);
                  const __m256 xi_364 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_280,_mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(xi_279,_mm256_set_ps(0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f)),_mm256_mul_ps(xi_338,xi_338)),_mm256_set_ps(-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f))),xi_266);
                  const __m256 xi_365 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_280,_mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(xi_279,_mm256_set_ps(0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f)),_mm256_mul_ps(xi_360,xi_360)),_mm256_set_ps(-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f))),xi_266);
                  const __m256 xi_367 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_280,_mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(xi_279,_mm256_set_ps(0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f)),_mm256_mul_ps(xi_354,xi_354)),_mm256_set_ps(-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f))),xi_266);
                  const __m256 xia_3_collide = _mm256_loadu_ps(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3 + ctr_0]);
                  const __m256 xia_4_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xia_5_collide = _mm256_load_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]);
                  const __m256 xi_0 = _mm256_cmp_ps(xia_5_collide,_mm256_set_ps(0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f),_CMP_NLE_UQ);
                  const __m256 xi_1 = _mm256_cmp_ps(xia_5_collide,_mm256_set_ps(-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f),_CMP_NGE_UQ);
                  const __m256 xi_10 = _mm256_mul_ps(xia_5_collide,xia_5_collide);
                  const __m256 xi_11 = _mm256_cmp_ps(xia_5_collide,_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_CMP_NLE_UQ);
                  const __m256 xia_6_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xi_270 = _mm256_mul_ps(xia_6_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xia_7_collide = _mm256_load_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + ctr_0]);
                  const __m256 xia_8_collide = _mm256_load_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xi_275 = _mm256_mul_ps(xia_8_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xia_9_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xia_10_collide = _mm256_loadu_ps(& _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + _stride_force_a_3 + ctr_0]);
                  const __m256 xi_26 = _mm256_mul_ps(xia_10_collide,_mm256_set_ps(0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f));
                  const __m256 xi_53 = _mm256_mul_ps(xia_10_collide,_mm256_set_ps(0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f));
                  const __m256 xi_57 = _mm256_mul_ps(xia_10_collide,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xia_11_collide = _mm256_load_ps(& _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + ctr_0]);
                  const __m256 xi_41 = _mm256_mul_ps(xia_11_collide,_mm256_set_ps(0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f));
                  const __m256 xi_66 = _mm256_mul_ps(xia_11_collide,_mm256_set_ps(0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f));
                  const __m256 xia_12_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xi_289 = _mm256_mul_ps(xia_12_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xia_13_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3 + ctr_0]);
                  const __m256 xia_14_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xia_15_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xia_16_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xia_17_collide = _mm256_load_ps(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0]);
                  const __m256 xia_18_collide = _mm256_load_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xi_271 = _mm256_add_ps(xi_270,xia_18_collide);
                  const __m256 xia_19_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xia_20_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xi_197 = _mm256_add_ps(xia_16_collide,xia_20_collide);
                  const __m256 xi_198 = _mm256_add_ps(_mm256_add_ps(xi_197,xia_18_collide),xia_6_collide);
                  const __m256 xia_21_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xi_276 = _mm256_add_ps(xi_275,xia_21_collide);
                  const __m256 xia_22_collide = _mm256_loadu_ps(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3 + ctr_0]);
                  const __m256 xia_23_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xi_201 = _mm256_add_ps(xia_23_collide,xia_9_collide);
                  const __m256 xia_24_collide = _mm256_load_ps(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + ctr_0]);
                  const __m256 xia_25_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xi_204 = _mm256_add_ps(xia_14_collide,xia_25_collide);
                  const __m256 xi_205 = _mm256_add_ps(_mm256_add_ps(xi_204,xia_21_collide),xia_8_collide);
                  const __m256 xia_26_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xia_27_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xi_202 = _mm256_add_ps(_mm256_add_ps(xi_201,xia_12_collide),xia_27_collide);
                  const __m256 xi_290 = _mm256_add_ps(xi_289,xia_27_collide);
                  const __m256 xia_28_collide = _mm256_loadu_ps(& _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + 2*_stride_force_a_3 + ctr_0]);
                  const __m256 xi_50 = _mm256_mul_ps(xia_28_collide,_mm256_set_ps(0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f));
                  const __m256 xi_76 = _mm256_mul_ps(xia_28_collide,_mm256_set_ps(0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f));
                  const __m256 xia_29_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xia_0_collide = _mm256_div_ps(_mm256_set_ps(1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f),_mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-0.25f,-0.25f,-0.25f,-0.25f,-0.25f,-0.25f,-0.25f,-0.25f),_mm256_blendv_ps(_mm256_blendv_ps(_mm256_blendv_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_10,_mm256_set_ps(xi_13,xi_13,xi_13,xi_13,xi_13,xi_13,xi_13,xi_13)),_mm256_mul_ps(xia_5_collide,_mm256_set_ps(xi_13,xi_13,xi_13,xi_13,xi_13,xi_13,xi_13,xi_13))),_mm256_set_ps(xi_5,xi_5,xi_5,xi_5,xi_5,xi_5,xi_5,xi_5)),_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_10,_mm256_set_ps(xi_8,xi_8,xi_8,xi_8,xi_8,xi_8,xi_8,xi_8)),_mm256_mul_ps(xia_5_collide,_mm256_set_ps(xi_9,xi_9,xi_9,xi_9,xi_9,xi_9,xi_9,xi_9))),_mm256_set_ps(xi_5,xi_5,xi_5,xi_5,xi_5,xi_5,xi_5,xi_5)),xi_11),_mm256_set_ps(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b),xi_1),_mm256_set_ps(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a),xi_0)),_mm256_set_ps(2.0f,2.0f,2.0f,2.0f,2.0f,2.0f,2.0f,2.0f)));
                  const __m256 xi_19 = _mm256_blendv_ps(_mm256_blendv_ps(_mm256_blendv_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_10,_mm256_set_ps(xi_18,xi_18,xi_18,xi_18,xi_18,xi_18,xi_18,xi_18)),_mm256_mul_ps(xia_5_collide,_mm256_set_ps(xi_18,xi_18,xi_18,xi_18,xi_18,xi_18,xi_18,xi_18))),_mm256_set_ps(xi_15,xi_15,xi_15,xi_15,xi_15,xi_15,xi_15,xi_15)),_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_10,_mm256_set_ps(xi_17,xi_17,xi_17,xi_17,xi_17,xi_17,xi_17,xi_17)),_mm256_mul_ps(_mm256_mul_ps(xia_5_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_set_ps(xi_17,xi_17,xi_17,xi_17,xi_17,xi_17,xi_17,xi_17))),_mm256_set_ps(xi_15,xi_15,xi_15,xi_15,xi_15,xi_15,xi_15,xi_15)),xi_11),_mm256_set_ps(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b),xi_1),_mm256_set_ps(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a),xi_0);
                  const __m256 rr_1_a_collide = xi_19;
                  const __m256 xi_23 = _mm256_mul_ps(rr_1_a_collide,_mm256_set_ps(0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f));
                  const __m256 xi_61 = _mm256_mul_ps(rr_1_a_collide,_mm256_set_ps(0.125f,0.125f,0.125f,0.125f,0.125f,0.125f,0.125f,0.125f));
                  const __m256 xi_62 = _mm256_mul_ps(xi_61,xia_10_collide);
                  const __m256 xi_70 = _mm256_mul_ps(rr_1_a_collide,_mm256_set_ps(0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f));
                  const __m256 xi_306 = _mm256_mul_ps(rr_1_a_collide,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 rr_2_a_collide = _mm256_mul_ps(xia_0_collide,_mm256_add_ps(_mm256_mul_ps(xi_19,_mm256_set_ps(-2.0f,-2.0f,-2.0f,-2.0f,-2.0f,-2.0f,-2.0f,-2.0f)),_mm256_set_ps(4.0f,4.0f,4.0f,4.0f,4.0f,4.0f,4.0f,4.0f)));
                  const __m256 xi_28 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(rr_2_a_collide,xi_26),_mm256_mul_ps(_mm256_mul_ps(xia_10_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_set_ps(xi_27,xi_27,xi_27,xi_27,xi_27,xi_27,xi_27,xi_27))),xi_26);
                  const __m256 xi_42 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(rr_2_a_collide,xi_41),_mm256_mul_ps(_mm256_mul_ps(xia_11_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_set_ps(xi_27,xi_27,xi_27,xi_27,xi_27,xi_27,xi_27,xi_27))),xi_41);
                  const __m256 xi_51 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(rr_2_a_collide,xi_50),_mm256_mul_ps(_mm256_mul_ps(xia_28_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_set_ps(xi_27,xi_27,xi_27,xi_27,xi_27,xi_27,xi_27,xi_27))),xi_50);
                  const __m256 xi_54 = _mm256_mul_ps(rr_2_a_collide,_mm256_set_ps(0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f));
                  const __m256 xi_55 = _mm256_mul_ps(xi_54,xia_10_collide);
                  const __m256 xi_56 = _mm256_add_ps(_mm256_mul_ps(xi_55,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_53);
                  const __m256 xi_67 = _mm256_mul_ps(xi_54,xia_11_collide);
                  const __m256 xi_68 = _mm256_add_ps(_mm256_mul_ps(xi_66,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_67);
                  const __m256 xi_72 = _mm256_add_ps(_mm256_mul_ps(xi_53,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_55);
                  const __m256 xi_74 = _mm256_add_ps(_mm256_mul_ps(xi_67,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_66);
                  const __m256 xi_77 = _mm256_mul_ps(xi_54,xia_28_collide);
                  const __m256 xi_78 = _mm256_add_ps(_mm256_mul_ps(xi_77,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_76);
                  const __m256 xi_89 = _mm256_add_ps(_mm256_mul_ps(xi_76,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_77);
                  const __m256 xi_273 = _mm256_mul_ps(rr_2_a_collide,_mm256_set_ps(0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f));
                  const __m256 xi_309 = _mm256_mul_ps(rr_2_a_collide,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 rho_a_collide = xia_17_collide;
                  const __m256 xi_191 = _mm256_mul_ps(rho_a_collide,_mm256_set_ps(-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f));
                  const __m256 xi_207 = _mm256_mul_ps(rho_a_collide,_mm256_set_ps(-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f));
                  const __m256 xi_208 = _mm256_add_ps(xi_202,xi_207);
                  const __m256 xi_268 = _mm256_mul_ps(rho_a_collide,_mm256_set_ps(0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f));
                  const __m256 u_0_a_collide = xia_24_collide;
                  const __m256 xi_20 = _mm256_mul_ps(u_0_a_collide,xia_11_collide);
                  const __m256 xi_33 = _mm256_mul_ps(xi_20,_mm256_set_ps(0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f));
                  const __m256 xi_34 = _mm256_mul_ps(xi_20,_mm256_set_ps(0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f));
                  const __m256 xi_35 = _mm256_mul_ps(rr_1_a_collide,xi_34);
                  const __m256 xi_36 = _mm256_add_ps(_mm256_mul_ps(xi_33,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_35);
                  const __m256 xi_39 = _mm256_add_ps(_mm256_mul_ps(xi_35,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_33);
                  const __m256 xi_40 = _mm256_mul_ps(rr_1_a_collide,xi_33);
                  const __m256 xi_58 = _mm256_mul_ps(_mm256_mul_ps(u_0_a_collide,xi_57),_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_63 = _mm256_mul_ps(u_0_a_collide,xi_62);
                  const __m256 xi_79 = _mm256_mul_ps(xi_34,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_80 = _mm256_mul_ps(xi_20,xi_70);
                  const __m256 xi_92 = _mm256_mul_ps(u_0_a_collide,xia_28_collide);
                  const __m256 xi_93 = _mm256_mul_ps(xi_92,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_96 = _mm256_mul_ps(xi_61,xi_92);
                  const __m256 xi_196 = _mm256_mul_ps(u_0_a_collide,u_0_a_collide);
                  const __m256 xi_200 = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(rho_a_collide,xi_196),_mm256_set_ps(-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f)),xi_191);
                  const __m256 xi_209 = _mm256_mul_ps(rr_1_a_collide,_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_205,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_208,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xia_15_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xia_26_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(rho_a_collide,xi_196)));
                  const __m256 xi_286 = _mm256_mul_ps(u_0_a_collide,xi_268);
                  const __m256 xi_287 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xia_25_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_276),xi_286),xia_14_collide);
                  const __m256 xi_288 = _mm256_mul_ps(xi_273,xi_287);
                  const __m256 xi_291 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xia_23_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_286),xi_290),xia_9_collide);
                  const __m256 xi_292 = _mm256_mul_ps(xi_273,xi_291);
                  const __m256 xi_310 = _mm256_mul_ps(xi_287,xi_309);
                  const __m256 xi_311 = _mm256_mul_ps(xi_310,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_351 = _mm256_mul_ps(xi_291,xi_309);
                  const __m256 xi_352 = _mm256_mul_ps(xi_351,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 u_1_a_collide = xia_3_collide;
                  const __m256 xi_21 = _mm256_mul_ps(u_1_a_collide,xia_10_collide);
                  const __m256 xi_24 = _mm256_mul_ps(xi_21,_mm256_set_ps(0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f));
                  const __m256 xi_25 = _mm256_mul_ps(rr_1_a_collide,xi_24);
                  const __m256 xi_43 = _mm256_mul_ps(xi_21,_mm256_set_ps(0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f));
                  const __m256 xi_44 = _mm256_mul_ps(rr_1_a_collide,xi_43);
                  const __m256 xi_45 = _mm256_add_ps(_mm256_mul_ps(xi_44,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_24);
                  const __m256 xi_47 = _mm256_add_ps(_mm256_mul_ps(xi_24,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_44);
                  const __m256 xi_49 = _mm256_add_ps(xi_36,xi_47);
                  const __m256 xi_52 = _mm256_add_ps(xi_39,xi_45);
                  const __m256 xi_59 = _mm256_mul_ps(u_1_a_collide,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_60 = _mm256_mul_ps(_mm256_mul_ps(xi_59,xia_11_collide),_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_64 = _mm256_mul_ps(u_1_a_collide,xi_61);
                  const __m256 xi_65 = _mm256_mul_ps(xi_64,xia_11_collide);
                  const __m256 xi_69 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_58,xi_60),xi_63),xi_65),xi_68);
                  const __m256 xi_75 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_58,xi_60),xi_63),xi_65),xi_74);
                  const __m256 xi_82 = _mm256_mul_ps(xi_59,xia_28_collide);
                  const __m256 xi_84 = _mm256_mul_ps(xi_64,xia_28_collide);
                  const __m256 xi_90 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(rr_1_a_collide,u_1_a_collide),xia_10_collide),_mm256_set_ps(-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f));
                  const __m256 xi_194 = _mm256_mul_ps(u_1_a_collide,u_1_a_collide);
                  const __m256 xi_195 = _mm256_mul_ps(_mm256_mul_ps(rho_a_collide,xi_194),_mm256_set_ps(-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f));
                  const __m256 xi_210 = _mm256_mul_ps(rr_1_a_collide,_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_198,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_205,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_207,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xia_13_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xia_4_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(rho_a_collide,xi_194)));
                  const __m256 xi_269 = _mm256_mul_ps(u_1_a_collide,xi_268);
                  const __m256 xi_272 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xia_20_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_269),xi_271),xia_16_collide);
                  const __m256 xi_274 = _mm256_mul_ps(xi_272,xi_273);
                  const __m256 xi_277 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xia_14_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_269),xi_276),xia_25_collide);
                  const __m256 xi_278 = _mm256_mul_ps(xi_273,xi_277);
                  const __m256 xi_305 = _mm256_mul_ps(rho_a_collide,u_1_a_collide);
                  const __m256 xi_307 = _mm256_mul_ps(xi_306,_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xia_21_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(u_0_a_collide,xi_305)),xi_204),xi_275));
                  const __m256 xi_308 = _mm256_mul_ps(xi_307,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_322 = _mm256_mul_ps(xi_277,xi_309);
                  const __m256 xi_334 = _mm256_mul_ps(xi_272,xi_309);
                  const __m256 xi_344 = _mm256_mul_ps(xi_334,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 u_2_a_collide = xia_22_collide;
                  const __m256 xi_22 = _mm256_mul_ps(u_2_a_collide,xia_28_collide);
                  const __m256 xi_29 = _mm256_mul_ps(xi_22,_mm256_set_ps(0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f));
                  const __m256 xi_30 = _mm256_mul_ps(xi_22,_mm256_set_ps(0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f));
                  const __m256 xi_31 = _mm256_mul_ps(rr_1_a_collide,xi_30);
                  const __m256 xi_32 = _mm256_add_ps(_mm256_mul_ps(xi_29,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_31);
                  const __m256 xi_37 = _mm256_add_ps(xi_32,xi_36);
                  const __m256 xi_38 = _mm256_add_ps(_mm256_mul_ps(xi_31,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_29);
                  const __m256 xi_46 = _mm256_add_ps(xi_38,xi_45);
                  const __m256 xi_48 = _mm256_mul_ps(rr_1_a_collide,xi_29);
                  const __m256 xi_71 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_30,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_22,xi_70)),xi_52);
                  const __m256 xi_73 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(rr_1_a_collide,u_2_a_collide),xia_28_collide),_mm256_set_ps(-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f)),xi_30),xi_49);
                  const __m256 xi_81 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_46,xi_56),xi_79),xi_80);
                  const __m256 xi_83 = _mm256_mul_ps(u_2_a_collide,xi_57);
                  const __m256 xi_85 = _mm256_mul_ps(u_2_a_collide,xi_62);
                  const __m256 xi_86 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_84,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_85,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_82),xi_83);
                  const __m256 xi_87 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_82,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_83,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_84),xi_85);
                  const __m256 xi_88 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_46,xi_72),xi_79),xi_80);
                  const __m256 xi_91 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_37,xi_43),xi_74),xi_90);
                  const __m256 xi_94 = _mm256_mul_ps(u_2_a_collide,xia_11_collide);
                  const __m256 xi_95 = _mm256_mul_ps(xi_94,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_97 = _mm256_mul_ps(xi_61,xi_94);
                  const __m256 xi_98 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_96,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_97,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_93),xi_95);
                  const __m256 xi_99 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_93,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_95,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_96),xi_97);
                  const __m256 xi_100 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_37,xi_43),xi_68),xi_90);
                  const __m256 xi_192 = _mm256_mul_ps(u_2_a_collide,u_2_a_collide);
                  const __m256 xi_193 = _mm256_mul_ps(_mm256_mul_ps(rho_a_collide,xi_192),_mm256_set_ps(-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f));
                  const __m256 xi_199 = _mm256_mul_ps(rr_1_a_collide,_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_191,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_193,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_195,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_198,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(_mm256_mul_ps(rho_a_collide,xi_196),_mm256_set_ps(-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f))));
                  const __m256 xi_203 = _mm256_mul_ps(rr_1_a_collide,_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_193,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_200,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_202,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(_mm256_mul_ps(rho_a_collide,xi_194),_mm256_set_ps(-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f))));
                  const __m256 xi_206 = _mm256_mul_ps(rr_1_a_collide,_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_195,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_200,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_205,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(_mm256_mul_ps(rho_a_collide,xi_192),_mm256_set_ps(-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f))));
                  const __m256 xi_211 = _mm256_mul_ps(rr_1_a_collide,_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_198,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_208,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xia_19_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xia_29_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(rho_a_collide,xi_192)));
                  const __m256 xi_282 = _mm256_mul_ps(xi_206,_mm256_set_ps(-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f));
                  const __m256 xi_283 = _mm256_mul_ps(xi_199,_mm256_set_ps(-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f));
                  const __m256 xi_284 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_210,_mm256_set_ps(0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f)),xi_282),xi_283);
                  const __m256 xi_294 = _mm256_mul_ps(xi_203,_mm256_set_ps(-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f));
                  const __m256 xi_295 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_209,_mm256_set_ps(0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f)),xi_282),xi_294);
                  const __m256 xi_297 = _mm256_mul_ps(u_2_a_collide,xi_268);
                  const __m256 xi_298 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xia_16_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_271),xi_297),xia_20_collide);
                  const __m256 xi_299 = _mm256_mul_ps(xi_273,xi_298);
                  const __m256 xi_300 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xia_9_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_290),xi_297),xia_23_collide);
                  const __m256 xi_301 = _mm256_mul_ps(xi_273,xi_300);
                  const __m256 xi_303 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_211,_mm256_set_ps(0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f)),xi_283),xi_294);
                  const __m256 xi_321 = _mm256_mul_ps(xi_206,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_323 = _mm256_add_ps(xi_321,xi_322);
                  const __m256 xi_331 = _mm256_add_ps(_mm256_mul_ps(xi_322,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_321);
                  const __m256 xi_333 = _mm256_mul_ps(xi_306,_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xia_18_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(u_2_a_collide,xi_305)),xi_197),xi_270));
                  const __m256 xi_340 = _mm256_mul_ps(xi_199,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_341 = _mm256_mul_ps(xi_298,xi_309);
                  const __m256 xi_342 = _mm256_add_ps(xi_340,xi_341);
                  const __m256 xi_343 = _mm256_mul_ps(xi_333,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_349 = _mm256_mul_ps(xi_306,_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xia_27_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(_mm256_mul_ps(rho_a_collide,u_0_a_collide),u_2_a_collide)),xi_201),xi_289));
                  const __m256 xi_350 = _mm256_mul_ps(xi_349,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_356 = _mm256_mul_ps(xi_203,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_357 = _mm256_mul_ps(xi_300,xi_309);
                  const __m256 xi_358 = _mm256_add_ps(xi_356,xi_357);
                  const __m256 xi_363 = _mm256_add_ps(_mm256_mul_ps(xi_341,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_340);
                  const __m256 xi_366 = _mm256_add_ps(_mm256_mul_ps(xi_357,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_356);
                  const __m256 forceTerm_0_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_20,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_21,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_22,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_20,xi_23)),_mm256_mul_ps(xi_21,xi_23)),_mm256_mul_ps(xi_22,xi_23));
                  const __m256 forceTerm_1_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_21,_mm256_set_ps(0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f)),_mm256_mul_ps(xi_25,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_28),xi_37);
                  const __m256 forceTerm_2_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_25,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_28,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_38,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_39,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(_mm256_mul_ps(u_1_a_collide,xia_10_collide),_mm256_set_ps(0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f)));
                  const __m256 forceTerm_3_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_40,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_42,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_46,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(_mm256_mul_ps(u_0_a_collide,xia_11_collide),_mm256_set_ps(0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f)));
                  const __m256 forceTerm_4_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_20,_mm256_set_ps(0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f)),_mm256_mul_ps(xi_40,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_32),xi_42),xi_47);
                  const __m256 forceTerm_5_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_22,_mm256_set_ps(0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f)),_mm256_mul_ps(xi_48,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_49),xi_51);
                  const __m256 forceTerm_6_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_48,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_51,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_52,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(_mm256_mul_ps(u_2_a_collide,xia_28_collide),_mm256_set_ps(0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f)));
                  const __m256 forceTerm_7_a_collide = _mm256_add_ps(_mm256_add_ps(xi_56,xi_69),xi_71);
                  const __m256 forceTerm_8_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_69,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_72,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_73,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_9_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_56,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_73,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_75,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_10_a_collide = _mm256_add_ps(_mm256_add_ps(xi_71,xi_72),xi_75);
                  const __m256 forceTerm_11_a_collide = _mm256_add_ps(_mm256_add_ps(xi_78,xi_81),xi_86);
                  const __m256 forceTerm_12_a_collide = _mm256_add_ps(_mm256_add_ps(xi_78,xi_87),xi_88);
                  const __m256 forceTerm_13_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_89,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_91,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_98,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_14_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_100,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_89,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_99,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_15_a_collide = _mm256_add_ps(_mm256_add_ps(xi_81,xi_87),xi_89);
                  const __m256 forceTerm_16_a_collide = _mm256_add_ps(_mm256_add_ps(xi_86,xi_88),xi_89);
                  const __m256 forceTerm_17_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_78,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_91,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_99,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_18_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_100,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_78,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_98,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 xib_3_collide = _mm256_loadu_ps(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3 + ctr_0]);
                  const __m256 xib_4_collide = _mm256_load_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xi_396 = _mm256_mul_ps(xib_4_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xib_5_collide = _mm256_load_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]);
                  const __m256 xi_101 = _mm256_cmp_ps(xib_5_collide,_mm256_set_ps(0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f),_CMP_NLE_UQ);
                  const __m256 xi_102 = _mm256_cmp_ps(xib_5_collide,_mm256_set_ps(-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f),_CMP_NGE_UQ);
                  const __m256 xi_103 = _mm256_mul_ps(xib_5_collide,xib_5_collide);
                  const __m256 xi_104 = _mm256_cmp_ps(xib_5_collide,_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_CMP_NLE_UQ);
                  const __m256 xib_6_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3 + ctr_0]);
                  const __m256 xib_7_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xib_8_collide = _mm256_loadu_ps(& _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + 2*_stride_force_b_3 + ctr_0]);
                  const __m256 xi_141 = _mm256_mul_ps(xib_8_collide,_mm256_set_ps(0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f));
                  const __m256 xi_166 = _mm256_mul_ps(xib_8_collide,_mm256_set_ps(0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f));
                  const __m256 xib_9_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xib_10_collide = _mm256_load_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xib_11_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xi_397 = _mm256_add_ps(xi_396,xib_11_collide);
                  const __m256 xib_12_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xi_391 = _mm256_mul_ps(xib_12_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_392 = _mm256_add_ps(xi_391,xib_10_collide);
                  const __m256 xib_13_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xib_14_collide = _mm256_load_ps(& _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + ctr_0]);
                  const __m256 xi_132 = _mm256_mul_ps(xib_14_collide,_mm256_set_ps(0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f));
                  const __m256 xi_143 = _mm256_mul_ps(xib_14_collide,_mm256_set_ps(0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f));
                  const __m256 xib_15_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xib_16_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xib_17_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xib_18_collide = _mm256_load_ps(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]);
                  const __m256 xib_19_collide = _mm256_loadu_ps(& _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + _stride_force_b_3 + ctr_0]);
                  const __m256 xi_117 = _mm256_mul_ps(xib_19_collide,_mm256_set_ps(0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f));
                  const __m256 xi_148 = _mm256_mul_ps(xib_19_collide,_mm256_set_ps(0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f));
                  const __m256 xi_152 = _mm256_mul_ps(xib_19_collide,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xib_20_collide = _mm256_loadu_ps(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3 + ctr_0]);
                  const __m256 xib_21_collide = _mm256_load_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + ctr_0]);
                  const __m256 xib_22_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xi_378 = _mm256_add_ps(xib_13_collide,xib_22_collide);
                  const __m256 xib_23_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xi_406 = _mm256_mul_ps(xib_23_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xib_24_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xib_25_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xib_26_collide = _mm256_load_ps(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + ctr_0]);
                  const __m256 xib_27_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xi_381 = _mm256_add_ps(xib_24_collide,xib_27_collide);
                  const __m256 xi_382 = _mm256_add_ps(_mm256_add_ps(xi_381,xib_11_collide),xib_4_collide);
                  const __m256 xib_28_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xi_374 = _mm256_add_ps(xib_17_collide,xib_28_collide);
                  const __m256 xi_375 = _mm256_add_ps(_mm256_add_ps(xi_374,xib_10_collide),xib_12_collide);
                  const __m256 xib_29_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xi_379 = _mm256_add_ps(_mm256_add_ps(xi_378,xib_23_collide),xib_29_collide);
                  const __m256 xi_407 = _mm256_add_ps(xi_406,xib_29_collide);
                  const __m256 xib_0_collide = _mm256_div_ps(_mm256_set_ps(1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f),_mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-0.25f,-0.25f,-0.25f,-0.25f,-0.25f,-0.25f,-0.25f,-0.25f),_mm256_blendv_ps(_mm256_blendv_ps(_mm256_blendv_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_103,_mm256_set_ps(xi_13,xi_13,xi_13,xi_13,xi_13,xi_13,xi_13,xi_13)),_mm256_mul_ps(xib_5_collide,_mm256_set_ps(xi_13,xi_13,xi_13,xi_13,xi_13,xi_13,xi_13,xi_13))),_mm256_set_ps(xi_5,xi_5,xi_5,xi_5,xi_5,xi_5,xi_5,xi_5)),_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_103,_mm256_set_ps(xi_8,xi_8,xi_8,xi_8,xi_8,xi_8,xi_8,xi_8)),_mm256_mul_ps(xib_5_collide,_mm256_set_ps(xi_9,xi_9,xi_9,xi_9,xi_9,xi_9,xi_9,xi_9))),_mm256_set_ps(xi_5,xi_5,xi_5,xi_5,xi_5,xi_5,xi_5,xi_5)),xi_104),_mm256_set_ps(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b),xi_102),_mm256_set_ps(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a),xi_101)),_mm256_set_ps(2.0f,2.0f,2.0f,2.0f,2.0f,2.0f,2.0f,2.0f)));
                  const __m256 xi_110 = _mm256_blendv_ps(_mm256_blendv_ps(_mm256_blendv_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_103,_mm256_set_ps(xi_109,xi_109,xi_109,xi_109,xi_109,xi_109,xi_109,xi_109)),_mm256_mul_ps(xib_5_collide,_mm256_set_ps(xi_109,xi_109,xi_109,xi_109,xi_109,xi_109,xi_109,xi_109))),_mm256_set_ps(xi_106,xi_106,xi_106,xi_106,xi_106,xi_106,xi_106,xi_106)),_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_103,_mm256_set_ps(xi_108,xi_108,xi_108,xi_108,xi_108,xi_108,xi_108,xi_108)),_mm256_mul_ps(_mm256_mul_ps(xib_5_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_set_ps(xi_108,xi_108,xi_108,xi_108,xi_108,xi_108,xi_108,xi_108))),_mm256_set_ps(xi_106,xi_106,xi_106,xi_106,xi_106,xi_106,xi_106,xi_106)),xi_104),_mm256_set_ps(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b),xi_102),_mm256_set_ps(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a),xi_101);
                  const __m256 rr_1_b_collide = xi_110;
                  const __m256 xi_114 = _mm256_mul_ps(rr_1_b_collide,_mm256_set_ps(0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f));
                  const __m256 xi_156 = _mm256_mul_ps(rr_1_b_collide,_mm256_set_ps(0.125f,0.125f,0.125f,0.125f,0.125f,0.125f,0.125f,0.125f));
                  const __m256 xi_157 = _mm256_mul_ps(xi_156,xib_19_collide);
                  const __m256 xi_419 = _mm256_mul_ps(rr_1_b_collide,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 rr_2_b_collide = _mm256_mul_ps(xib_0_collide,_mm256_add_ps(_mm256_mul_ps(xi_110,_mm256_set_ps(-2.0f,-2.0f,-2.0f,-2.0f,-2.0f,-2.0f,-2.0f,-2.0f)),_mm256_set_ps(4.0f,4.0f,4.0f,4.0f,4.0f,4.0f,4.0f,4.0f)));
                  const __m256 xi_119 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(rr_2_b_collide,xi_117),_mm256_mul_ps(_mm256_mul_ps(xib_19_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_set_ps(xi_118,xi_118,xi_118,xi_118,xi_118,xi_118,xi_118,xi_118))),xi_117);
                  const __m256 xi_133 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(rr_2_b_collide,xi_132),_mm256_mul_ps(_mm256_mul_ps(xib_14_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_set_ps(xi_118,xi_118,xi_118,xi_118,xi_118,xi_118,xi_118,xi_118))),xi_132);
                  const __m256 xi_142 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(rr_2_b_collide,xi_141),_mm256_mul_ps(_mm256_mul_ps(xib_8_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_set_ps(xi_118,xi_118,xi_118,xi_118,xi_118,xi_118,xi_118,xi_118))),xi_141);
                  const __m256 xi_144 = _mm256_mul_ps(rr_2_b_collide,_mm256_set_ps(0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f));
                  const __m256 xi_145 = _mm256_mul_ps(xi_144,xib_14_collide);
                  const __m256 xi_146 = _mm256_add_ps(_mm256_mul_ps(xi_145,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_143);
                  const __m256 xi_149 = _mm256_mul_ps(xi_144,xib_19_collide);
                  const __m256 xi_150 = _mm256_add_ps(_mm256_mul_ps(xi_148,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_149);
                  const __m256 xi_162 = _mm256_add_ps(_mm256_mul_ps(xi_143,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_145);
                  const __m256 xi_164 = _mm256_add_ps(_mm256_mul_ps(xi_149,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_148);
                  const __m256 xi_167 = _mm256_mul_ps(xi_144,xib_8_collide);
                  const __m256 xi_168 = _mm256_add_ps(_mm256_mul_ps(xi_167,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_166);
                  const __m256 xi_179 = _mm256_add_ps(_mm256_mul_ps(xi_166,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_167);
                  const __m256 xi_394 = _mm256_mul_ps(rr_2_b_collide,_mm256_set_ps(0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f));
                  const __m256 xi_422 = _mm256_mul_ps(rr_2_b_collide,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 rho_b_collide = xib_18_collide;
                  const __m256 xi_371 = _mm256_mul_ps(rho_b_collide,_mm256_set_ps(-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f));
                  const __m256 xi_384 = _mm256_mul_ps(rho_b_collide,_mm256_set_ps(-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f));
                  const __m256 xi_386 = _mm256_add_ps(xi_375,xi_384);
                  const __m256 xi_389 = _mm256_mul_ps(rho_b_collide,_mm256_set_ps(0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f));
                  const __m256 u_0_b_collide = xib_26_collide;
                  const __m256 xi_111 = _mm256_mul_ps(u_0_b_collide,xib_14_collide);
                  const __m256 xi_124 = _mm256_mul_ps(xi_111,_mm256_set_ps(0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f));
                  const __m256 xi_125 = _mm256_mul_ps(xi_111,_mm256_set_ps(0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f));
                  const __m256 xi_126 = _mm256_mul_ps(rr_1_b_collide,xi_125);
                  const __m256 xi_127 = _mm256_add_ps(_mm256_mul_ps(xi_124,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_126);
                  const __m256 xi_130 = _mm256_add_ps(_mm256_mul_ps(xi_126,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_124);
                  const __m256 xi_131 = _mm256_mul_ps(rr_1_b_collide,xi_124);
                  const __m256 xi_153 = _mm256_mul_ps(u_0_b_collide,xi_152);
                  const __m256 xi_158 = _mm256_mul_ps(u_0_b_collide,xi_157);
                  const __m256 xi_169 = _mm256_mul_ps(xi_125,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_170 = _mm256_mul_ps(_mm256_mul_ps(rr_1_b_collide,xi_111),_mm256_set_ps(0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f));
                  const __m256 xi_182 = _mm256_mul_ps(u_0_b_collide,xib_8_collide);
                  const __m256 xi_183 = _mm256_mul_ps(xi_182,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_186 = _mm256_mul_ps(xi_156,xi_182);
                  const __m256 xi_370 = _mm256_mul_ps(u_0_b_collide,u_0_b_collide);
                  const __m256 xi_377 = _mm256_mul_ps(_mm256_mul_ps(rho_b_collide,xi_370),_mm256_set_ps(-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f));
                  const __m256 xi_385 = _mm256_mul_ps(rr_1_b_collide,_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_379,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_382,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_384,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xib_15_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xib_7_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(rho_b_collide,xi_370)));
                  const __m256 xi_403 = _mm256_mul_ps(u_0_b_collide,xi_389);
                  const __m256 xi_404 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xib_24_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_397),xi_403),xib_27_collide);
                  const __m256 xi_405 = _mm256_mul_ps(xi_394,xi_404);
                  const __m256 xi_408 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xib_13_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_403),xi_407),xib_22_collide);
                  const __m256 xi_409 = _mm256_mul_ps(xi_394,xi_408);
                  const __m256 xi_423 = _mm256_mul_ps(xi_404,xi_422);
                  const __m256 xi_424 = _mm256_mul_ps(xi_423,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_438 = _mm256_mul_ps(xi_408,xi_422);
                  const __m256 xi_439 = _mm256_mul_ps(xi_438,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 u_1_b_collide = xib_3_collide;
                  const __m256 xi_112 = _mm256_mul_ps(u_1_b_collide,xib_19_collide);
                  const __m256 xi_115 = _mm256_mul_ps(xi_112,_mm256_set_ps(0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f));
                  const __m256 xi_116 = _mm256_mul_ps(rr_1_b_collide,xi_115);
                  const __m256 xi_134 = _mm256_mul_ps(xi_112,_mm256_set_ps(0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f));
                  const __m256 xi_135 = _mm256_mul_ps(rr_1_b_collide,xi_134);
                  const __m256 xi_136 = _mm256_add_ps(_mm256_mul_ps(xi_135,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_115);
                  const __m256 xi_138 = _mm256_add_ps(_mm256_mul_ps(xi_115,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_135);
                  const __m256 xi_140 = _mm256_add_ps(xi_127,xi_138);
                  const __m256 xi_154 = _mm256_mul_ps(u_1_b_collide,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_155 = _mm256_mul_ps(xi_154,xib_14_collide);
                  const __m256 xi_159 = _mm256_mul_ps(u_1_b_collide,xi_156);
                  const __m256 xi_160 = _mm256_mul_ps(xi_159,xib_14_collide);
                  const __m256 xi_161 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_158,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_160,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_153),xi_155);
                  const __m256 xi_163 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_153,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_155,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_158),xi_160);
                  const __m256 xi_172 = _mm256_mul_ps(xi_154,xib_8_collide);
                  const __m256 xi_174 = _mm256_mul_ps(xi_159,xib_8_collide);
                  const __m256 xi_180 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(rr_1_b_collide,u_1_b_collide),xib_19_collide),_mm256_set_ps(-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f));
                  const __m256 xi_372 = _mm256_mul_ps(u_1_b_collide,u_1_b_collide);
                  const __m256 xi_373 = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(rho_b_collide,xi_372),_mm256_set_ps(-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f)),xi_371);
                  const __m256 xi_387 = _mm256_mul_ps(rr_1_b_collide,_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_382,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_386,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xib_25_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xib_6_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(rho_b_collide,xi_372)));
                  const __m256 xi_390 = _mm256_mul_ps(u_1_b_collide,xi_389);
                  const __m256 xi_393 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xib_28_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_390),xi_392),xib_17_collide);
                  const __m256 xi_395 = _mm256_mul_ps(xi_393,xi_394);
                  const __m256 xi_398 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xib_27_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_390),xi_397),xib_24_collide);
                  const __m256 xi_399 = _mm256_mul_ps(xi_394,xi_398);
                  const __m256 xi_418 = _mm256_mul_ps(rho_b_collide,u_1_b_collide);
                  const __m256 xi_420 = _mm256_mul_ps(xi_419,_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xib_11_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(u_0_b_collide,xi_418)),xi_381),xi_396));
                  const __m256 xi_421 = _mm256_mul_ps(xi_420,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_426 = _mm256_mul_ps(xi_398,xi_422);
                  const __m256 xi_430 = _mm256_mul_ps(xi_393,xi_422);
                  const __m256 xi_435 = _mm256_mul_ps(xi_430,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 u_2_b_collide = xib_20_collide;
                  const __m256 xi_113 = _mm256_mul_ps(u_2_b_collide,xib_8_collide);
                  const __m256 xi_120 = _mm256_mul_ps(xi_113,_mm256_set_ps(0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f));
                  const __m256 xi_121 = _mm256_mul_ps(xi_113,_mm256_set_ps(0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f));
                  const __m256 xi_122 = _mm256_mul_ps(rr_1_b_collide,xi_121);
                  const __m256 xi_123 = _mm256_add_ps(_mm256_mul_ps(xi_120,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_122);
                  const __m256 xi_128 = _mm256_add_ps(xi_123,xi_127);
                  const __m256 xi_129 = _mm256_add_ps(_mm256_mul_ps(xi_122,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_120);
                  const __m256 xi_137 = _mm256_add_ps(xi_129,xi_136);
                  const __m256 xi_139 = _mm256_mul_ps(rr_1_b_collide,xi_120);
                  const __m256 xi_147 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(rr_1_b_collide,u_2_b_collide),xib_8_collide),_mm256_set_ps(-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f));
                  const __m256 xi_151 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_121,xi_140),xi_147),xi_150);
                  const __m256 xi_165 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_121,xi_140),xi_147),xi_164);
                  const __m256 xi_171 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_137,xi_164),xi_169),xi_170);
                  const __m256 xi_173 = _mm256_mul_ps(u_2_b_collide,xi_152);
                  const __m256 xi_175 = _mm256_mul_ps(u_2_b_collide,xi_157);
                  const __m256 xi_176 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_174,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_175,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_172),xi_173);
                  const __m256 xi_177 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_172,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_173,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_174),xi_175);
                  const __m256 xi_178 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_137,xi_150),xi_169),xi_170);
                  const __m256 xi_181 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_128,xi_134),xi_146),xi_180);
                  const __m256 xi_184 = _mm256_mul_ps(u_2_b_collide,xib_14_collide);
                  const __m256 xi_185 = _mm256_mul_ps(xi_184,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_187 = _mm256_mul_ps(xi_156,xi_184);
                  const __m256 xi_188 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_186,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_187,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_183),xi_185);
                  const __m256 xi_189 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_183,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_185,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_186),xi_187);
                  const __m256 xi_190 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_128,xi_134),xi_162),xi_180);
                  const __m256 xi_368 = _mm256_mul_ps(u_2_b_collide,u_2_b_collide);
                  const __m256 xi_369 = _mm256_mul_ps(_mm256_mul_ps(rho_b_collide,xi_368),_mm256_set_ps(-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f));
                  const __m256 xi_376 = _mm256_mul_ps(rr_1_b_collide,_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_369,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_373,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_375,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(_mm256_mul_ps(rho_b_collide,xi_370),_mm256_set_ps(-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f))));
                  const __m256 xi_380 = _mm256_mul_ps(rr_1_b_collide,_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_369,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_371,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_377,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_379,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(_mm256_mul_ps(rho_b_collide,xi_372),_mm256_set_ps(-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f))));
                  const __m256 xi_383 = _mm256_mul_ps(rr_1_b_collide,_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_373,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_377,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_382,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(_mm256_mul_ps(rho_b_collide,xi_368),_mm256_set_ps(-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f))));
                  const __m256 xi_388 = _mm256_mul_ps(rr_1_b_collide,_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_379,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_386,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xib_16_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xib_9_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(rho_b_collide,xi_368)));
                  const __m256 xi_400 = _mm256_mul_ps(xi_376,_mm256_set_ps(-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f));
                  const __m256 xi_401 = _mm256_mul_ps(xi_383,_mm256_set_ps(-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f));
                  const __m256 xi_402 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_387,_mm256_set_ps(0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f)),xi_400),xi_401);
                  const __m256 xi_410 = _mm256_mul_ps(xi_380,_mm256_set_ps(-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f));
                  const __m256 xi_411 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_385,_mm256_set_ps(0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f)),xi_401),xi_410);
                  const __m256 xi_412 = _mm256_mul_ps(u_2_b_collide,xi_389);
                  const __m256 xi_413 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xib_17_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_392),xi_412),xib_28_collide);
                  const __m256 xi_414 = _mm256_mul_ps(xi_394,xi_413);
                  const __m256 xi_415 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xib_22_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_407),xi_412),xib_13_collide);
                  const __m256 xi_416 = _mm256_mul_ps(xi_394,xi_415);
                  const __m256 xi_417 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_388,_mm256_set_ps(0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f)),xi_400),xi_410);
                  const __m256 xi_425 = _mm256_mul_ps(xi_383,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_427 = _mm256_add_ps(xi_425,xi_426);
                  const __m256 xi_428 = _mm256_add_ps(_mm256_mul_ps(xi_426,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_425);
                  const __m256 xi_429 = _mm256_mul_ps(xi_419,_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xib_10_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(u_2_b_collide,xi_418)),xi_374),xi_391));
                  const __m256 xi_431 = _mm256_mul_ps(xi_376,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_432 = _mm256_mul_ps(xi_413,xi_422);
                  const __m256 xi_433 = _mm256_add_ps(xi_431,xi_432);
                  const __m256 xi_434 = _mm256_mul_ps(xi_429,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_436 = _mm256_mul_ps(xi_419,_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xib_29_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(_mm256_mul_ps(rho_b_collide,u_0_b_collide),u_2_b_collide)),xi_378),xi_406));
                  const __m256 xi_437 = _mm256_mul_ps(xi_436,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_440 = _mm256_mul_ps(xi_380,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_441 = _mm256_mul_ps(xi_415,xi_422);
                  const __m256 xi_442 = _mm256_add_ps(xi_440,xi_441);
                  const __m256 xi_443 = _mm256_add_ps(_mm256_mul_ps(xi_432,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_431);
                  const __m256 xi_444 = _mm256_add_ps(_mm256_mul_ps(xi_441,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_440);
                  const __m256 forceTerm_0_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_111,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_112,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_113,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_111,xi_114)),_mm256_mul_ps(xi_112,xi_114)),_mm256_mul_ps(xi_113,xi_114));
                  const __m256 forceTerm_1_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_112,_mm256_set_ps(0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f)),_mm256_mul_ps(xi_116,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_119),xi_128);
                  const __m256 forceTerm_2_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_116,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_119,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_129,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_130,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(_mm256_mul_ps(u_1_b_collide,xib_19_collide),_mm256_set_ps(0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f)));
                  const __m256 forceTerm_3_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_131,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_133,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_137,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(_mm256_mul_ps(u_0_b_collide,xib_14_collide),_mm256_set_ps(0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f)));
                  const __m256 forceTerm_4_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_111,_mm256_set_ps(0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f)),_mm256_mul_ps(xi_131,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_123),xi_133),xi_138);
                  const __m256 forceTerm_5_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_113,_mm256_set_ps(0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f)),_mm256_mul_ps(xi_139,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_140),xi_142);
                  const __m256 forceTerm_6_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_130,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_136,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_139,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_142,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(_mm256_mul_ps(u_2_b_collide,xib_8_collide),_mm256_set_ps(0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f)));
                  const __m256 forceTerm_7_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_146,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_151,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_161,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_8_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_151,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_162,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_163,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_9_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_146,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_163,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_165,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_10_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_161,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_162,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_165,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_11_b_collide = _mm256_add_ps(_mm256_add_ps(xi_168,xi_171),xi_176);
                  const __m256 forceTerm_12_b_collide = _mm256_add_ps(_mm256_add_ps(xi_168,xi_177),xi_178);
                  const __m256 forceTerm_13_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_179,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_181,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_188,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_14_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_179,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_189,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_190,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_15_b_collide = _mm256_add_ps(_mm256_add_ps(xi_171,xi_177),xi_179);
                  const __m256 forceTerm_16_b_collide = _mm256_add_ps(_mm256_add_ps(xi_176,xi_178),xi_179);
                  const __m256 forceTerm_17_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_168,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_181,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_189,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_18_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_168,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_188,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_190,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 tmp_a0 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_209,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_210,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_211,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),forceTerm_0_a_collide),xi_199),xi_203),xi_206),xi_267),xia_7_collide);
                  const __m256 tmp_a1 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_274,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_278,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),forceTerm_1_a_collide),xi_281),xi_284),xia_13_collide);
                  const __m256 tmp_a2 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_2_a_collide,xi_274),xi_278),xi_284),xi_285),xia_4_collide);
                  const __m256 tmp_a3 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_3_a_collide,xi_288),xi_292),xi_293),xi_295),xia_15_collide);
                  const __m256 tmp_a4 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_288,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_292,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),forceTerm_4_a_collide),xi_295),xi_296),xia_26_collide);
                  const __m256 tmp_a5 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_299,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_301,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),forceTerm_5_a_collide),xi_302),xi_303),xia_29_collide);
                  const __m256 tmp_a6 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_6_a_collide,xi_299),xi_301),xi_303),xi_304),xia_19_collide);
                  const __m256 tmp_a7 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_7_a_collide,xi_308),xi_311),xi_320),xi_323),xia_14_collide);
                  const __m256 tmp_a8 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_8_a_collide,xi_307),xi_310),xi_323),xi_329),xia_8_collide);
                  const __m256 tmp_a9 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_9_a_collide,xi_307),xi_311),xi_330),xi_331),xia_21_collide);
                  const __m256 tmp_a10 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_10_a_collide,xi_308),xi_310),xi_331),xi_332),xia_25_collide);
                  const __m256 tmp_a11 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_11_a_collide,xi_333),xi_334),xi_339),xi_342),xia_6_collide);
                  const __m256 tmp_a12 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_12_a_collide,xi_342),xi_343),xi_344),xi_348),xia_16_collide);
                  const __m256 tmp_a13 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_13_a_collide,xi_350),xi_352),xi_355),xi_358),xia_9_collide);
                  const __m256 tmp_a14 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_14_a_collide,xi_349),xi_351),xi_358),xi_361),xia_12_collide);
                  const __m256 tmp_a15 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_15_a_collide,xi_334),xi_343),xi_362),xi_363),xia_20_collide);
                  const __m256 tmp_a16 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_16_a_collide,xi_333),xi_344),xi_363),xi_364),xia_18_collide);
                  const __m256 tmp_a17 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_17_a_collide,xi_349),xi_352),xi_365),xi_366),xia_27_collide);
                  const __m256 tmp_a18 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_18_a_collide,xi_350),xi_351),xi_366),xi_367),xia_23_collide);
                  const __m256 tmp_b0 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_385,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_387,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_388,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),forceTerm_0_b_collide),xi_267),xi_376),xi_380),xi_383),xib_21_collide);
                  const __m256 tmp_b1 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_395,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_399,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),forceTerm_1_b_collide),xi_281),xi_402),xib_6_collide);
                  const __m256 tmp_b2 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_2_b_collide,xi_285),xi_395),xi_399),xi_402),xib_25_collide);
                  const __m256 tmp_b3 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_3_b_collide,xi_293),xi_405),xi_409),xi_411),xib_7_collide);
                  const __m256 tmp_b4 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_405,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_409,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),forceTerm_4_b_collide),xi_296),xi_411),xib_15_collide);
                  const __m256 tmp_b5 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_414,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_416,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),forceTerm_5_b_collide),xi_302),xi_417),xib_16_collide);
                  const __m256 tmp_b6 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_6_b_collide,xi_304),xi_414),xi_416),xi_417),xib_9_collide);
                  const __m256 tmp_b7 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_7_b_collide,xi_320),xi_421),xi_424),xi_427),xib_27_collide);
                  const __m256 tmp_b8 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_8_b_collide,xi_329),xi_420),xi_423),xi_427),xib_4_collide);
                  const __m256 tmp_b9 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_9_b_collide,xi_330),xi_420),xi_424),xi_428),xib_11_collide);
                  const __m256 tmp_b10 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_10_b_collide,xi_332),xi_421),xi_423),xi_428),xib_24_collide);
                  const __m256 tmp_b11 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_11_b_collide,xi_339),xi_429),xi_430),xi_433),xib_12_collide);
                  const __m256 tmp_b12 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_12_b_collide,xi_348),xi_433),xi_434),xi_435),xib_17_collide);
                  const __m256 tmp_b13 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_13_b_collide,xi_355),xi_437),xi_439),xi_442),xib_22_collide);
                  const __m256 tmp_b14 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_14_b_collide,xi_361),xi_436),xi_438),xi_442),xib_23_collide);
                  const __m256 tmp_b15 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_15_b_collide,xi_362),xi_430),xi_434),xi_443),xib_28_collide);
                  const __m256 tmp_b16 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_16_b_collide,xi_364),xi_429),xi_435),xi_443),xib_10_collide);
                  const __m256 tmp_b17 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_17_b_collide,xi_365),xi_436),xi_439),xi_444),xib_29_collide);
                  const __m256 tmp_b18 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_18_b_collide,xi_367),xi_437),xi_438),xi_444),xib_13_collide);
                  const __m256 xirecolor_0 = _mm256_add_ps(tmp_a0,tmp_b0);
                  const __m256 xirecolor_1 = _mm256_add_ps(_mm256_load_ps(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0]),_mm256_load_ps(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]));
                  const __m256 xirecolor_2 = _mm256_div_ps(_mm256_set_ps(1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f),xirecolor_1);
                  const __m256 xi_466 = _mm256_mul_ps(xirecolor_2,_mm256_load_ps(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]));
                  const __m256 xirecolor_3 = _mm256_mul_ps(xirecolor_2,_mm256_load_ps(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0]));
                  const __m256 xirecolor_4 = _mm256_add_ps(tmp_a1,tmp_b1);
                  const __m256 xirecolor_5 = xi_222;
                  const __m256 xirecolor_6 = xi_220;
                  const __m256 xirecolor_7 = xi_240;
                  const __m256 xi_445 = _mm256_mul_ps(xirecolor_7,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xirecolor_8 = xi_445;
                  const __m256 xirecolor_9 = xi_227;
                  const __m256 xirecolor_10 = _mm256_mul_ps(xirecolor_9,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xirecolor_11 = xi_229;
                  const __m256 xi_446 = _mm256_add_ps(xirecolor_10,xirecolor_11);
                  const __m256 xirecolor_12 = xi_242;
                  const __m256 xi_447 = _mm256_add_ps(xirecolor_12,xirecolor_8);
                  const __m256 xirecolor_13 = _mm256_add_ps(xi_446,xi_447);
                  const __m256 xirecolor_14 = xi_244;
                  const __m256 xirecolor_15 = xi_245;
                  const __m256 xi_448 = _mm256_add_ps(_mm256_mul_ps(xirecolor_15,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_14);
                  const __m256 xirecolor_16 = _mm256_mul_ps(xi_448,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xirecolor_17 = xi_214;
                  const __m256 xirecolor_18 = xi_215;
                  const __m256 xirecolor_19 = xi_216;
                  const __m256 xi_449 = _mm256_add_ps(xirecolor_18,xirecolor_19);
                  const __m256 xirecolor_20 = xi_217;
                  const __m256 xirecolor_21 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xirecolor_17,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_449),xirecolor_20);
                  const __m256 xirecolor_22 = xi_213;
                  const __m256 xirecolor_23 = xi_212;
                  const __m256 xi_450 = _mm256_add_ps(_mm256_mul_ps(xirecolor_23,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_22);
                  const __m256 xirecolor_24 = _mm256_mul_ps(xi_450,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xirecolor_25 = xi_248;
                  const __m256 xirecolor_26 = _mm256_mul_ps(xirecolor_25,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xirecolor_27 = xi_234;
                  const __m256 xirecolor_28 = xi_233;
                  const __m256 xirecolor_29 = xi_250;
                  const __m256 xi_451 = _mm256_add_ps(xirecolor_26,xirecolor_29);
                  const __m256 xirecolor_30 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xirecolor_27,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_451),xirecolor_28);
                  const __m256 xirecolor_31 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xirecolor_5,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_13),xirecolor_16),xirecolor_21),xirecolor_24),xirecolor_30),xirecolor_6);
                  const __m256 xirecolor_32 = xi_219;
                  const __m256 xirecolor_33 = xi_221;
                  const __m256 xi_452 = _mm256_add_ps(xirecolor_32,xirecolor_33);
                  const __m256 xirecolor_34 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xirecolor_6,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_452),xirecolor_5);
                  const __m256 xirecolor_35 = xi_224;
                  const __m256 xirecolor_36 = xi_225;
                  const __m256 xi_453 = _mm256_add_ps(_mm256_mul_ps(xirecolor_36,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_35);
                  const __m256 xirecolor_37 = _mm256_mul_ps(xi_453,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xirecolor_38 = xi_446;
                  const __m256 xirecolor_39 = xi_231;
                  const __m256 xirecolor_40 = xi_232;
                  const __m256 xi_454 = _mm256_add_ps(xirecolor_39,xirecolor_40);
                  const __m256 xirecolor_41 = _mm256_add_ps(_mm256_mul_ps(xirecolor_28,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_27);
                  const __m256 xirecolor_42 = _mm256_add_ps(xi_454,xirecolor_41);
                  const __m256 xirecolor_43 = _mm256_add_ps(xirecolor_38,xirecolor_42);
                  const __m256 xirecolor_44 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_450,xirecolor_21),xirecolor_34),xirecolor_37),xirecolor_43);
                  const __m256 xirecolor_45 = xi_448;
                  const __m256 xirecolor_46 = xi_447;
                  const __m256 xirecolor_47 = xi_256;
                  const __m256 xirecolor_48 = _mm256_mul_ps(xirecolor_47,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xirecolor_49 = xi_258;
                  const __m256 xi_455 = _mm256_add_ps(xirecolor_48,xirecolor_49);
                  const __m256 xirecolor_50 = _mm256_add_ps(xi_453,xi_455);
                  const __m256 xirecolor_51 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xirecolor_20,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_449),xirecolor_17),xirecolor_24),xirecolor_34),xirecolor_45),xirecolor_46),xirecolor_50);
                  const __m256 xirecolor_52 = _mm256_sqrt_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xirecolor_31,xirecolor_31),_mm256_mul_ps(xirecolor_44,xirecolor_44)),_mm256_mul_ps(xirecolor_51,xirecolor_51)));
                  const __m256 xirecolor_53 = _mm256_div_ps(_mm256_set_ps(1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f),xirecolor_52);
                  const __m256 xi_456 = _mm256_mul_ps(xirecolor_31,xirecolor_53);
                  const __m256 xi_457 = _mm256_mul_ps(xirecolor_51,xirecolor_53);
                  const __m256 xi_458 = _mm256_mul_ps(xirecolor_44,xirecolor_53);
                  const __m256 xirecolor_54 = _mm256_cmp_ps(xirecolor_52,_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_CMP_NLE_UQ);
                  const __m256 xirecolor_55 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set_ps(beta,beta,beta,beta,beta,beta,beta,beta),_mm256_div_ps(_mm256_set_ps(1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f),_mm256_mul_ps(xirecolor_1,xirecolor_1))),_mm256_load_ps(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0])),_mm256_load_ps(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]));
                  const __m256 xirecolor_56 = _mm256_mul_ps(xirecolor_55,_mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f),_mm256_load_ps(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0])),_mm256_mul_ps(_mm256_set_ps(0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f),_mm256_load_ps(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]))));
                  const __m256 xirecolor_57 = _mm256_mul_ps(xirecolor_56,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),xi_456,xirecolor_54));
                  const __m256 xirecolor_58 = _mm256_add_ps(tmp_a2,tmp_b2);
                  const __m256 xirecolor_59 = _mm256_mul_ps(xirecolor_56,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_456,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_54));
                  const __m256 xirecolor_60 = _mm256_add_ps(tmp_a3,tmp_b3);
                  const __m256 xirecolor_61 = _mm256_mul_ps(xirecolor_56,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_457,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_54));
                  const __m256 xirecolor_62 = _mm256_add_ps(tmp_a4,tmp_b4);
                  const __m256 xirecolor_63 = _mm256_mul_ps(xirecolor_56,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),xi_457,xirecolor_54));
                  const __m256 xirecolor_64 = _mm256_add_ps(tmp_a5,tmp_b5);
                  const __m256 xirecolor_65 = _mm256_mul_ps(xirecolor_56,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_458,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_54));
                  const __m256 xirecolor_66 = _mm256_add_ps(tmp_a6,tmp_b6);
                  const __m256 xirecolor_67 = _mm256_mul_ps(xirecolor_56,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),xi_458,xirecolor_54));
                  const __m256 xirecolor_68 = _mm256_add_ps(tmp_a7,tmp_b7);
                  const __m256 xirecolor_69 = xi_312;
                  const __m256 xirecolor_70 = xi_313;
                  const __m256 xirecolor_71 = xi_452;
                  const __m256 xirecolor_72 = _mm256_add_ps(xirecolor_50,xirecolor_71);
                  const __m256 xirecolor_73 = xi_317;
                  const __m256 xirecolor_74 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xirecolor_29,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_25),xirecolor_73);
                  const __m256 xirecolor_75 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xirecolor_11,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xirecolor_69,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_314),xirecolor_41),xirecolor_70),xirecolor_72),xirecolor_74),xirecolor_9);
                  const __m256 xirecolor_76 = _mm256_mul_ps(xirecolor_53,_mm256_set_ps(0.70710678118654757f,0.70710678118654757f,0.70710678118654757f,0.70710678118654757f,0.70710678118654757f,0.70710678118654757f,0.70710678118654757f,0.70710678118654757f));
                  const __m256 xi_459 = _mm256_mul_ps(xirecolor_75,xirecolor_76);
                  const __m256 xirecolor_77 = _mm256_mul_ps(xirecolor_55,_mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f),_mm256_load_ps(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0])),_mm256_mul_ps(_mm256_set_ps(0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f),_mm256_load_ps(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]))));
                  const __m256 xirecolor_78 = _mm256_mul_ps(xirecolor_77,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_459,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_54));
                  const __m256 xirecolor_79 = _mm256_add_ps(tmp_a8,tmp_b8);
                  const __m256 xirecolor_80 = xi_324;
                  const __m256 xirecolor_81 = xi_325;
                  const __m256 xi_460 = _mm256_add_ps(_mm256_mul_ps(xirecolor_80,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_81);
                  const __m256 xirecolor_82 = xi_326;
                  const __m256 xirecolor_83 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_327,xi_460),xirecolor_30),xirecolor_38),xirecolor_72),xirecolor_82);
                  const __m256 xi_461 = _mm256_mul_ps(xirecolor_76,xirecolor_83);
                  const __m256 xirecolor_84 = _mm256_mul_ps(xirecolor_77,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),xi_461,xirecolor_54));
                  const __m256 xirecolor_85 = _mm256_add_ps(tmp_a9,tmp_b9);
                  const __m256 xirecolor_86 = _mm256_mul_ps(xirecolor_77,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_461,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_54));
                  const __m256 xirecolor_87 = _mm256_add_ps(tmp_a10,tmp_b10);
                  const __m256 xirecolor_88 = _mm256_mul_ps(xirecolor_77,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),xi_459,xirecolor_54));
                  const __m256 xirecolor_89 = _mm256_add_ps(tmp_a11,tmp_b11);
                  const __m256 xirecolor_90 = _mm256_add_ps(_mm256_add_ps(xi_454,xirecolor_37),xirecolor_71);
                  const __m256 xirecolor_91 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_445,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_460,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xirecolor_12,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 xirecolor_92 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_335,xirecolor_45),xirecolor_74),xirecolor_90),xirecolor_91);
                  const __m256 xi_462 = _mm256_mul_ps(xirecolor_76,xirecolor_92);
                  const __m256 xirecolor_93 = _mm256_mul_ps(xirecolor_77,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_462,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_54));
                  const __m256 xirecolor_94 = _mm256_add_ps(tmp_a12,tmp_b12);
                  const __m256 xirecolor_95 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xirecolor_70,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_16),xirecolor_69);
                  const __m256 xirecolor_96 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_345,xi_451),xirecolor_46),xirecolor_82),xirecolor_90),xirecolor_95);
                  const __m256 xi_463 = _mm256_mul_ps(xirecolor_76,xirecolor_96);
                  const __m256 xirecolor_97 = _mm256_mul_ps(xirecolor_77,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_463,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_54));
                  const __m256 xirecolor_98 = _mm256_add_ps(tmp_a13,tmp_b13);
                  const __m256 xirecolor_99 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_353,xi_455),xirecolor_13),xirecolor_42),xirecolor_45),xirecolor_73),xirecolor_82);
                  const __m256 xi_464 = _mm256_mul_ps(xirecolor_76,xirecolor_99);
                  const __m256 xirecolor_100 = _mm256_mul_ps(xirecolor_77,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_464,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_54));
                  const __m256 xirecolor_101 = _mm256_add_ps(tmp_a14,tmp_b14);
                  const __m256 xirecolor_102 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xirecolor_49,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_359),xirecolor_43),xirecolor_47),xirecolor_91),xirecolor_95);
                  const __m256 xi_465 = _mm256_mul_ps(xirecolor_102,xirecolor_76);
                  const __m256 xirecolor_103 = _mm256_mul_ps(xirecolor_77,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_465,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_54));
                  const __m256 xirecolor_104 = _mm256_add_ps(tmp_a15,tmp_b15);
                  const __m256 xirecolor_105 = _mm256_mul_ps(xirecolor_77,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),xi_463,xirecolor_54));
                  const __m256 xirecolor_106 = _mm256_add_ps(tmp_a16,tmp_b16);
                  const __m256 xirecolor_107 = _mm256_mul_ps(xirecolor_77,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),xi_462,xirecolor_54));
                  const __m256 xirecolor_108 = _mm256_add_ps(tmp_a17,tmp_b17);
                  const __m256 xirecolor_109 = _mm256_mul_ps(xirecolor_77,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),xi_465,xirecolor_54));
                  const __m256 xirecolor_110 = _mm256_add_ps(tmp_a18,tmp_b18);
                  const __m256 xirecolor_111 = _mm256_mul_ps(xirecolor_77,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),xi_464,xirecolor_54));
                  _mm256_store_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + ctr_0],_mm256_mul_ps(xirecolor_0,xirecolor_3));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_3,xirecolor_4),xirecolor_57));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_3,xirecolor_58),xirecolor_59));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_3,xirecolor_60),xirecolor_61));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_3,xirecolor_62),xirecolor_63));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_3,xirecolor_64),xirecolor_65));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_3,xirecolor_66),xirecolor_67));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_3,xirecolor_68),xirecolor_78));
                  _mm256_store_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_3,xirecolor_79),xirecolor_84));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_3,xirecolor_85),xirecolor_86));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_3,xirecolor_87),xirecolor_88));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_3,xirecolor_89),xirecolor_93));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_3,xirecolor_94),xirecolor_97));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_3,xirecolor_98),xirecolor_100));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_101,xirecolor_3),xirecolor_103));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_104,xirecolor_3),xirecolor_105));
                  _mm256_store_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_106,xirecolor_3),xirecolor_107));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_108,xirecolor_3),xirecolor_109));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_110,xirecolor_3),xirecolor_111));
                  _mm256_store_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + ctr_0],_mm256_mul_ps(xi_466,xirecolor_0));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_57,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_466,xirecolor_4)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_59,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_466,xirecolor_58)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_61,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_466,xirecolor_60)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_63,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_466,xirecolor_62)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_65,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_466,xirecolor_64)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_67,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_466,xirecolor_66)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_78,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_466,xirecolor_68)));
                  _mm256_store_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_84,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_466,xirecolor_79)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_86,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_466,xirecolor_85)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_88,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_466,xirecolor_87)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_93,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_466,xirecolor_89)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_97,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_466,xirecolor_94)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_100,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_466,xirecolor_98)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_103,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_466,xirecolor_101)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_105,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_466,xirecolor_104)));
                  _mm256_store_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_107,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_466,xirecolor_106)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_109,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_466,xirecolor_108)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_111,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_466,xirecolor_110)));
               }
               for (int64_t ctr_0 = (int64_t)((_size_force_a_0 - 2) / (8)) * (8) + 1; ctr_0 < _size_force_a_0 - 1; ctr_0 += 1)
               {
                  const float xi_212 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1];
                  const float xi_213 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1];
                  const float xi_214 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1];
                  const float xi_215 = -0.013888888888888888f*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1];
                  const float xi_216 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1];
                  const float xi_217 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1];
                  const float xi_218 = -xi_214 + xi_215 + xi_216 + xi_217;
                  const float xi_219 = -0.055555555555555552f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1];
                  const float xi_220 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1];
                  const float xi_221 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1];
                  const float xi_222 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1];
                  const float xi_223 = xi_219 - xi_220 + xi_221 + xi_222;
                  const float xi_224 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1];
                  const float xi_225 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1];
                  const float xi_226 = -xi_224 + xi_225;
                  const float xi_227 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0];
                  const float xi_228 = -xi_227;
                  const float xi_229 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0];
                  const float xi_230 = xi_228 + xi_229;
                  const float xi_231 = -0.22222222222222221f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0];
                  const float xi_232 = 0.22222222222222221f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0];
                  const float xi_233 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0];
                  const float xi_234 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0];
                  const float xi_235 = -xi_233 + xi_234;
                  const float xi_236 = xi_231 + xi_232 + xi_235;
                  const float xi_237 = xi_230 + xi_236;
                  const float xi_238 = -xi_212 + xi_213 + xi_218 + xi_223 + xi_226 + xi_237;
                  const float xi_239 = (xi_238*xi_238);
                  const float xi_240 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1];
                  const float xi_241 = -xi_240;
                  const float xi_242 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1];
                  const float xi_243 = xi_228 + xi_229 + xi_241 + xi_242;
                  const float xi_244 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1];
                  const float xi_245 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1];
                  const float xi_246 = -xi_244 + xi_245;
                  const float xi_247 = xi_212 - xi_213;
                  const float xi_248 = 0.22222222222222221f*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0];
                  const float xi_249 = -xi_248;
                  const float xi_250 = 0.22222222222222221f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0];
                  const float xi_251 = xi_233 - xi_234 + xi_249 + xi_250;
                  const float xi_252 = xi_218 + xi_220 - xi_222 + xi_243 + xi_246 + xi_247 + xi_251;
                  const float xi_253 = (xi_252*xi_252);
                  const float xi_254 = xi_244 - xi_245;
                  const float xi_255 = xi_241 + xi_242;
                  const float xi_256 = 0.22222222222222221f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1];
                  const float xi_257 = -xi_256;
                  const float xi_258 = 0.22222222222222221f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1];
                  const float xi_259 = xi_224 - xi_225 + xi_257 + xi_258;
                  const float xi_260 = xi_214 + xi_215 + xi_216 - xi_217 + xi_223 + xi_247 + xi_254 + xi_255 + xi_259;
                  const float xi_261 = (xi_260*xi_260);
                  const float xi_262 = xi_239 + xi_253 + xi_261;
                  const float xi_263 = powf(xi_262, 0.5f);
                  const float xi_264 = (_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]);
                  const float xi_265 = sigma*xi_263*((0.5f < _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]) ? (omega_shear_a): ((-0.5f > _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]) ? (omega_shear_b): ((0.0f < _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]) ? (xi_264*xi_8 + xi_5 + xi_9*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]): (xi_13*xi_264 + xi_13*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0] + xi_5))));
                  const bool xi_266 = xi_263 > 0.0f;
                  const float xi_267 = ((xi_266) ? (xi_265*0.25f): (0.0f));
                  const float xi_279 = ((1.0f) / (xi_262));
                  const float xi_280 = xi_265*1.125f;
                  const float xi_281 = ((xi_266) ? (xi_280*(xi_253*xi_279*0.055555555555555552f - 0.018518518518518517f)): (0.0f));
                  const float xi_285 = ((xi_266) ? (xi_280*(xi_279*0.055555555555555552f*(xi_252*xi_252) - 0.018518518518518517f)): (0.0f));
                  const float xi_293 = ((xi_266) ? (xi_280*(xi_279*0.055555555555555552f*(xi_260*xi_260) - 0.018518518518518517f)): (0.0f));
                  const float xi_296 = ((xi_266) ? (xi_280*(xi_261*xi_279*0.055555555555555552f - 0.018518518518518517f)): (0.0f));
                  const float xi_302 = ((xi_266) ? (xi_280*(xi_279*0.055555555555555552f*(xi_238*xi_238) - 0.018518518518518517f)): (0.0f));
                  const float xi_304 = ((xi_266) ? (xi_280*(xi_239*xi_279*0.055555555555555552f - 0.018518518518518517f)): (0.0f));
                  const float xi_312 = 0.027777777777777776f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1];
                  const float xi_313 = 0.027777777777777776f*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1];
                  const float xi_314 = -0.1111111111111111f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1] + 0.1111111111111111f*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1];
                  const float xi_315 = xi_219 + xi_221;
                  const float xi_316 = xi_259 + xi_315;
                  const float xi_317 = -0.027777777777777776f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1] + 0.027777777777777776f*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1];
                  const float xi_318 = xi_248 - xi_250 + xi_317;
                  const float xi_319 = xi_227 - xi_229 + xi_235 - xi_312 + xi_313 + xi_314 + xi_316 + xi_318;
                  const float xi_320 = ((xi_266) ? (xi_280*(xi_279*0.027777777777777776f*(xi_319*xi_319) - 0.037037037037037035f)): (0.0f));
                  const float xi_324 = 0.027777777777777776f*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1];
                  const float xi_325 = 0.027777777777777776f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1];
                  const float xi_326 = -0.027777777777777776f*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1] + 0.027777777777777776f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1];
                  const float xi_327 = -0.1111111111111111f*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1] + 0.1111111111111111f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1];
                  const float xi_328 = xi_230 + xi_251 + xi_316 - xi_324 + xi_325 + xi_326 + xi_327;
                  const float xi_329 = ((xi_266) ? (xi_280*(xi_279*0.027777777777777776f*(xi_328*xi_328) - 0.037037037037037035f)): (0.0f));
                  const float xi_330 = ((xi_266) ? (xi_280*(xi_279*0.027777777777777776f*(xi_328*xi_328) - 0.037037037037037035f)): (0.0f));
                  const float xi_332 = ((xi_266) ? (xi_280*(xi_279*0.027777777777777776f*(xi_319*xi_319) - 0.037037037037037035f)): (0.0f));
                  const float xi_335 = -0.1111111111111111f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0] + 0.1111111111111111f*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0];
                  const float xi_336 = xi_226 + xi_231 + xi_232 + xi_315;
                  const float xi_337 = xi_240 - xi_242 + xi_324 - xi_325;
                  const float xi_338 = xi_254 + xi_318 + xi_335 + xi_336 + xi_337;
                  const float xi_339 = ((xi_266) ? (xi_280*(xi_279*0.027777777777777776f*(xi_338*xi_338) - 0.037037037037037035f)): (0.0f));
                  const float xi_345 = -0.1111111111111111f*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0] + 0.1111111111111111f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0];
                  const float xi_346 = xi_246 + xi_312 - xi_313;
                  const float xi_347 = xi_249 + xi_250 + xi_255 + xi_326 + xi_336 + xi_345 + xi_346;
                  const float xi_348 = ((xi_266) ? (xi_280*(xi_279*0.027777777777777776f*(xi_347*xi_347) - 0.037037037037037035f)): (0.0f));
                  const float xi_353 = -0.1111111111111111f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1] + 0.1111111111111111f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1];
                  const float xi_354 = xi_236 + xi_243 + xi_254 + xi_257 + xi_258 + xi_317 + xi_326 + xi_353;
                  const float xi_355 = ((xi_266) ? (xi_280*(xi_279*0.027777777777777776f*(xi_354*xi_354) - 0.037037037037037035f)): (0.0f));
                  const float xi_359 = -0.1111111111111111f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1] + 0.1111111111111111f*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1];
                  const float xi_360 = xi_237 + xi_256 - xi_258 + xi_337 + xi_346 + xi_359;
                  const float xi_361 = ((xi_266) ? (xi_280*(xi_279*0.027777777777777776f*(xi_360*xi_360) - 0.037037037037037035f)): (0.0f));
                  const float xi_362 = ((xi_266) ? (xi_280*(xi_279*0.027777777777777776f*(xi_347*xi_347) - 0.037037037037037035f)): (0.0f));
                  const float xi_364 = ((xi_266) ? (xi_280*(xi_279*0.027777777777777776f*(xi_338*xi_338) - 0.037037037037037035f)): (0.0f));
                  const float xi_365 = ((xi_266) ? (xi_280*(xi_279*0.027777777777777776f*(xi_360*xi_360) - 0.037037037037037035f)): (0.0f));
                  const float xi_367 = ((xi_266) ? (xi_280*(xi_279*0.027777777777777776f*(xi_354*xi_354) - 0.037037037037037035f)): (0.0f));
                  const float xia_3_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3 + ctr_0];
                  const float xia_4_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3 + ctr_0];
                  const float xia_5_collide = _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0];
                  const bool xi_0 = xia_5_collide > 0.5f;
                  const bool xi_1 = xia_5_collide < -0.5f;
                  const float xi_10 = (xia_5_collide*xia_5_collide);
                  const bool xi_11 = xia_5_collide > 0.0f;
                  const float xia_6_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3 + ctr_0];
                  const float xi_270 = -xia_6_collide;
                  const float xia_7_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + ctr_0];
                  const float xia_8_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3 + ctr_0];
                  const float xi_275 = -xia_8_collide;
                  const float xia_9_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3 + ctr_0];
                  const float xia_10_collide = _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + _stride_force_a_3 + ctr_0];
                  const float xi_26 = xia_10_collide*0.16666666666666666f;
                  const float xi_53 = xia_10_collide*0.083333333333333329f;
                  const float xi_57 = xia_10_collide*0.25f;
                  const float xia_11_collide = _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + ctr_0];
                  const float xi_41 = xia_11_collide*0.16666666666666666f;
                  const float xi_66 = xia_11_collide*0.083333333333333329f;
                  const float xia_12_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3 + ctr_0];
                  const float xi_289 = -xia_12_collide;
                  const float xia_13_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3 + ctr_0];
                  const float xia_14_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3 + ctr_0];
                  const float xia_15_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3 + ctr_0];
                  const float xia_16_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3 + ctr_0];
                  const float xia_17_collide = _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0];
                  const float xia_18_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3 + ctr_0];
                  const float xi_271 = xi_270 + xia_18_collide;
                  const float xia_19_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3 + ctr_0];
                  const float xia_20_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3 + ctr_0];
                  const float xi_197 = xia_16_collide + xia_20_collide;
                  const float xi_198 = xi_197 + xia_18_collide + xia_6_collide;
                  const float xia_21_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3 + ctr_0];
                  const float xi_276 = xi_275 + xia_21_collide;
                  const float xia_22_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3 + ctr_0];
                  const float xia_23_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3 + ctr_0];
                  const float xi_201 = xia_23_collide + xia_9_collide;
                  const float xia_24_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + ctr_0];
                  const float xia_25_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3 + ctr_0];
                  const float xi_204 = xia_14_collide + xia_25_collide;
                  const float xi_205 = xi_204 + xia_21_collide + xia_8_collide;
                  const float xia_26_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3 + ctr_0];
                  const float xia_27_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3 + ctr_0];
                  const float xi_202 = xi_201 + xia_12_collide + xia_27_collide;
                  const float xi_290 = xi_289 + xia_27_collide;
                  const float xia_28_collide = _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + 2*_stride_force_a_3 + ctr_0];
                  const float xi_50 = xia_28_collide*0.16666666666666666f;
                  const float xi_76 = xia_28_collide*0.083333333333333329f;
                  const float xia_29_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3 + ctr_0];
                  const float xia_0_collide = ((1.0f) / (-0.25f*((xi_0) ? (omega_shear_a): ((xi_1) ? (omega_shear_b): ((xi_11) ? (xi_10*xi_8 + xi_5 + xi_9*xia_5_collide): (xi_10*xi_13 + xi_13*xia_5_collide + xi_5)))) + 2.0f));
                  const float xi_19 = ((xi_0) ? (omega_shear_a): ((xi_1) ? (omega_shear_b): ((xi_11) ? (xi_10*xi_17 + xi_15 - xi_17*xia_5_collide): (xi_10*xi_18 + xi_15 + xi_18*xia_5_collide))));
                  const float rr_1_a_collide = xi_19;
                  const float xi_23 = rr_1_a_collide*0.5f;
                  const float xi_61 = rr_1_a_collide*0.125f;
                  const float xi_62 = xi_61*xia_10_collide;
                  const float xi_70 = rr_1_a_collide*0.041666666666666664f;
                  const float xi_306 = rr_1_a_collide*0.25f;
                  const float rr_2_a_collide = xia_0_collide*(xi_19*-2.0f + 4.0f);
                  const float xi_28 = rr_2_a_collide*xi_26 + xi_26 - xi_27*xia_10_collide;
                  const float xi_42 = rr_2_a_collide*xi_41 - xi_27*xia_11_collide + xi_41;
                  const float xi_51 = rr_2_a_collide*xi_50 - xi_27*xia_28_collide + xi_50;
                  const float xi_54 = rr_2_a_collide*0.041666666666666664f;
                  const float xi_55 = xi_54*xia_10_collide;
                  const float xi_56 = xi_53 - xi_55;
                  const float xi_67 = xi_54*xia_11_collide;
                  const float xi_68 = -xi_66 + xi_67;
                  const float xi_72 = -xi_53 + xi_55;
                  const float xi_74 = xi_66 - xi_67;
                  const float xi_77 = xi_54*xia_28_collide;
                  const float xi_78 = xi_76 - xi_77;
                  const float xi_89 = -xi_76 + xi_77;
                  const float xi_273 = rr_2_a_collide*0.5f;
                  const float xi_309 = rr_2_a_collide*0.25f;
                  const float rho_a_collide = xia_17_collide;
                  const float xi_191 = rho_a_collide*-0.1111111111111111f;
                  const float xi_207 = rho_a_collide*-0.33333333333333331f;
                  const float xi_208 = xi_202 + xi_207;
                  const float xi_268 = rho_a_collide*0.33333333333333331f;
                  const float u_0_a_collide = xia_24_collide;
                  const float xi_20 = u_0_a_collide*xia_11_collide;
                  const float xi_33 = xi_20*0.16666666666666666f;
                  const float xi_34 = xi_20*0.083333333333333329f;
                  const float xi_35 = rr_1_a_collide*xi_34;
                  const float xi_36 = -xi_33 + xi_35;
                  const float xi_39 = xi_33 - xi_35;
                  const float xi_40 = rr_1_a_collide*xi_33;
                  const float xi_58 = -u_0_a_collide*xi_57;
                  const float xi_63 = u_0_a_collide*xi_62;
                  const float xi_79 = -xi_34;
                  const float xi_80 = xi_20*xi_70;
                  const float xi_92 = u_0_a_collide*xia_28_collide;
                  const float xi_93 = xi_92*0.25f;
                  const float xi_96 = xi_61*xi_92;
                  const float xi_196 = (u_0_a_collide*u_0_a_collide);
                  const float xi_200 = rho_a_collide*xi_196*-0.33333333333333331f + xi_191;
                  const float xi_209 = rr_1_a_collide*(rho_a_collide*xi_196 - xi_205 - xi_208 - xia_15_collide - xia_26_collide);
                  const float xi_286 = u_0_a_collide*xi_268;
                  const float xi_287 = xi_276 + xi_286 + xia_14_collide - xia_25_collide;
                  const float xi_288 = xi_273*xi_287;
                  const float xi_291 = xi_286 + xi_290 - xia_23_collide + xia_9_collide;
                  const float xi_292 = xi_273*xi_291;
                  const float xi_310 = xi_287*xi_309;
                  const float xi_311 = -xi_310;
                  const float xi_351 = xi_291*xi_309;
                  const float xi_352 = -xi_351;
                  const float u_1_a_collide = xia_3_collide;
                  const float xi_21 = u_1_a_collide*xia_10_collide;
                  const float xi_24 = xi_21*0.16666666666666666f;
                  const float xi_25 = rr_1_a_collide*xi_24;
                  const float xi_43 = xi_21*0.083333333333333329f;
                  const float xi_44 = rr_1_a_collide*xi_43;
                  const float xi_45 = xi_24 - xi_44;
                  const float xi_47 = -xi_24 + xi_44;
                  const float xi_49 = xi_36 + xi_47;
                  const float xi_52 = xi_39 + xi_45;
                  const float xi_59 = u_1_a_collide*0.25f;
                  const float xi_60 = -xi_59*xia_11_collide;
                  const float xi_64 = u_1_a_collide*xi_61;
                  const float xi_65 = xi_64*xia_11_collide;
                  const float xi_69 = xi_58 + xi_60 + xi_63 + xi_65 + xi_68;
                  const float xi_75 = xi_58 + xi_60 + xi_63 + xi_65 + xi_74;
                  const float xi_82 = xi_59*xia_28_collide;
                  const float xi_84 = xi_64*xia_28_collide;
                  const float xi_90 = rr_1_a_collide*u_1_a_collide*xia_10_collide*-0.041666666666666664f;
                  const float xi_194 = (u_1_a_collide*u_1_a_collide);
                  const float xi_195 = rho_a_collide*xi_194*-0.33333333333333331f;
                  const float xi_210 = rr_1_a_collide*(rho_a_collide*xi_194 - xi_198 - xi_205 - xi_207 - xia_13_collide - xia_4_collide);
                  const float xi_269 = u_1_a_collide*xi_268;
                  const float xi_272 = xi_269 + xi_271 + xia_16_collide - xia_20_collide;
                  const float xi_274 = xi_272*xi_273;
                  const float xi_277 = xi_269 + xi_276 - xia_14_collide + xia_25_collide;
                  const float xi_278 = xi_273*xi_277;
                  const float xi_305 = rho_a_collide*u_1_a_collide;
                  const float xi_307 = xi_306*(u_0_a_collide*xi_305 + xi_204 + xi_275 - xia_21_collide);
                  const float xi_308 = -xi_307;
                  const float xi_322 = xi_277*xi_309;
                  const float xi_334 = xi_272*xi_309;
                  const float xi_344 = -xi_334;
                  const float u_2_a_collide = xia_22_collide;
                  const float xi_22 = u_2_a_collide*xia_28_collide;
                  const float xi_29 = xi_22*0.16666666666666666f;
                  const float xi_30 = xi_22*0.083333333333333329f;
                  const float xi_31 = rr_1_a_collide*xi_30;
                  const float xi_32 = -xi_29 + xi_31;
                  const float xi_37 = xi_32 + xi_36;
                  const float xi_38 = xi_29 - xi_31;
                  const float xi_46 = xi_38 + xi_45;
                  const float xi_48 = rr_1_a_collide*xi_29;
                  const float xi_71 = xi_22*xi_70 - xi_30 + xi_52;
                  const float xi_73 = rr_1_a_collide*u_2_a_collide*xia_28_collide*-0.041666666666666664f + xi_30 + xi_49;
                  const float xi_81 = xi_46 + xi_56 + xi_79 + xi_80;
                  const float xi_83 = u_2_a_collide*xi_57;
                  const float xi_85 = u_2_a_collide*xi_62;
                  const float xi_86 = xi_82 + xi_83 - xi_84 - xi_85;
                  const float xi_87 = -xi_82 - xi_83 + xi_84 + xi_85;
                  const float xi_88 = xi_46 + xi_72 + xi_79 + xi_80;
                  const float xi_91 = xi_37 + xi_43 + xi_74 + xi_90;
                  const float xi_94 = u_2_a_collide*xia_11_collide;
                  const float xi_95 = xi_94*0.25f;
                  const float xi_97 = xi_61*xi_94;
                  const float xi_98 = xi_93 + xi_95 - xi_96 - xi_97;
                  const float xi_99 = -xi_93 - xi_95 + xi_96 + xi_97;
                  const float xi_100 = xi_37 + xi_43 + xi_68 + xi_90;
                  const float xi_192 = (u_2_a_collide*u_2_a_collide);
                  const float xi_193 = rho_a_collide*xi_192*-0.33333333333333331f;
                  const float xi_199 = rr_1_a_collide*(rho_a_collide*xi_196*-0.16666666666666666f - xi_191 - xi_193 - xi_195 - xi_198);
                  const float xi_203 = rr_1_a_collide*(rho_a_collide*xi_194*-0.16666666666666666f - xi_193 - xi_200 - xi_202);
                  const float xi_206 = rr_1_a_collide*(rho_a_collide*xi_192*-0.16666666666666666f - xi_195 - xi_200 - xi_205);
                  const float xi_211 = rr_1_a_collide*(rho_a_collide*xi_192 - xi_198 - xi_208 - xia_19_collide - xia_29_collide);
                  const float xi_282 = xi_206*-0.5f;
                  const float xi_283 = xi_199*-0.5f;
                  const float xi_284 = xi_210*0.5f + xi_282 + xi_283;
                  const float xi_294 = xi_203*-0.5f;
                  const float xi_295 = xi_209*0.5f + xi_282 + xi_294;
                  const float xi_297 = u_2_a_collide*xi_268;
                  const float xi_298 = xi_271 + xi_297 - xia_16_collide + xia_20_collide;
                  const float xi_299 = xi_273*xi_298;
                  const float xi_300 = xi_290 + xi_297 + xia_23_collide - xia_9_collide;
                  const float xi_301 = xi_273*xi_300;
                  const float xi_303 = xi_211*0.5f + xi_283 + xi_294;
                  const float xi_321 = xi_206*0.25f;
                  const float xi_323 = xi_321 + xi_322;
                  const float xi_331 = xi_321 - xi_322;
                  const float xi_333 = xi_306*(u_2_a_collide*xi_305 + xi_197 + xi_270 - xia_18_collide);
                  const float xi_340 = xi_199*0.25f;
                  const float xi_341 = xi_298*xi_309;
                  const float xi_342 = xi_340 + xi_341;
                  const float xi_343 = -xi_333;
                  const float xi_349 = xi_306*(rho_a_collide*u_0_a_collide*u_2_a_collide + xi_201 + xi_289 - xia_27_collide);
                  const float xi_350 = -xi_349;
                  const float xi_356 = xi_203*0.25f;
                  const float xi_357 = xi_300*xi_309;
                  const float xi_358 = xi_356 + xi_357;
                  const float xi_363 = xi_340 - xi_341;
                  const float xi_366 = xi_356 - xi_357;
                  const float forceTerm_0_a_collide = xi_20*xi_23 - xi_20 + xi_21*xi_23 - xi_21 + xi_22*xi_23 - xi_22;
                  const float forceTerm_1_a_collide = xi_21*0.33333333333333331f - xi_25 + xi_28 + xi_37;
                  const float forceTerm_2_a_collide = u_1_a_collide*xia_10_collide*0.33333333333333331f - xi_25 - xi_28 - xi_38 - xi_39;
                  const float forceTerm_3_a_collide = u_0_a_collide*xia_11_collide*0.33333333333333331f - xi_40 - xi_42 - xi_46;
                  const float forceTerm_4_a_collide = xi_20*0.33333333333333331f + xi_32 - xi_40 + xi_42 + xi_47;
                  const float forceTerm_5_a_collide = xi_22*0.33333333333333331f - xi_48 + xi_49 + xi_51;
                  const float forceTerm_6_a_collide = u_2_a_collide*xia_28_collide*0.33333333333333331f - xi_48 - xi_51 - xi_52;
                  const float forceTerm_7_a_collide = xi_56 + xi_69 + xi_71;
                  const float forceTerm_8_a_collide = -xi_69 - xi_72 - xi_73;
                  const float forceTerm_9_a_collide = -xi_56 - xi_73 - xi_75;
                  const float forceTerm_10_a_collide = xi_71 + xi_72 + xi_75;
                  const float forceTerm_11_a_collide = xi_78 + xi_81 + xi_86;
                  const float forceTerm_12_a_collide = xi_78 + xi_87 + xi_88;
                  const float forceTerm_13_a_collide = -xi_89 - xi_91 - xi_98;
                  const float forceTerm_14_a_collide = -xi_100 - xi_89 - xi_99;
                  const float forceTerm_15_a_collide = xi_81 + xi_87 + xi_89;
                  const float forceTerm_16_a_collide = xi_86 + xi_88 + xi_89;
                  const float forceTerm_17_a_collide = -xi_78 - xi_91 - xi_99;
                  const float forceTerm_18_a_collide = -xi_100 - xi_78 - xi_98;
                  const float xib_3_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3 + ctr_0];
                  const float xib_4_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3 + ctr_0];
                  const float xi_396 = -xib_4_collide;
                  const float xib_5_collide = _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0];
                  const bool xi_101 = xib_5_collide > 0.5f;
                  const bool xi_102 = xib_5_collide < -0.5f;
                  const float xi_103 = (xib_5_collide*xib_5_collide);
                  const bool xi_104 = xib_5_collide > 0.0f;
                  const float xib_6_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3 + ctr_0];
                  const float xib_7_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3 + ctr_0];
                  const float xib_8_collide = _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + 2*_stride_force_b_3 + ctr_0];
                  const float xi_141 = xib_8_collide*0.16666666666666666f;
                  const float xi_166 = xib_8_collide*0.083333333333333329f;
                  const float xib_9_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3 + ctr_0];
                  const float xib_10_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3 + ctr_0];
                  const float xib_11_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3 + ctr_0];
                  const float xi_397 = xi_396 + xib_11_collide;
                  const float xib_12_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3 + ctr_0];
                  const float xi_391 = -xib_12_collide;
                  const float xi_392 = xi_391 + xib_10_collide;
                  const float xib_13_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3 + ctr_0];
                  const float xib_14_collide = _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + ctr_0];
                  const float xi_132 = xib_14_collide*0.16666666666666666f;
                  const float xi_143 = xib_14_collide*0.083333333333333329f;
                  const float xib_15_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3 + ctr_0];
                  const float xib_16_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3 + ctr_0];
                  const float xib_17_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3 + ctr_0];
                  const float xib_18_collide = _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0];
                  const float xib_19_collide = _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + _stride_force_b_3 + ctr_0];
                  const float xi_117 = xib_19_collide*0.16666666666666666f;
                  const float xi_148 = xib_19_collide*0.083333333333333329f;
                  const float xi_152 = xib_19_collide*0.25f;
                  const float xib_20_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3 + ctr_0];
                  const float xib_21_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + ctr_0];
                  const float xib_22_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3 + ctr_0];
                  const float xi_378 = xib_13_collide + xib_22_collide;
                  const float xib_23_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3 + ctr_0];
                  const float xi_406 = -xib_23_collide;
                  const float xib_24_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3 + ctr_0];
                  const float xib_25_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3 + ctr_0];
                  const float xib_26_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + ctr_0];
                  const float xib_27_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3 + ctr_0];
                  const float xi_381 = xib_24_collide + xib_27_collide;
                  const float xi_382 = xi_381 + xib_11_collide + xib_4_collide;
                  const float xib_28_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3 + ctr_0];
                  const float xi_374 = xib_17_collide + xib_28_collide;
                  const float xi_375 = xi_374 + xib_10_collide + xib_12_collide;
                  const float xib_29_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3 + ctr_0];
                  const float xi_379 = xi_378 + xib_23_collide + xib_29_collide;
                  const float xi_407 = xi_406 + xib_29_collide;
                  const float xib_0_collide = ((1.0f) / (-0.25f*((xi_101) ? (omega_shear_a): ((xi_102) ? (omega_shear_b): ((xi_104) ? (xi_103*xi_8 + xi_5 + xi_9*xib_5_collide): (xi_103*xi_13 + xi_13*xib_5_collide + xi_5)))) + 2.0f));
                  const float xi_110 = ((xi_101) ? (omega_shear_a): ((xi_102) ? (omega_shear_b): ((xi_104) ? (xi_103*xi_108 + xi_106 - xi_108*xib_5_collide): (xi_103*xi_109 + xi_106 + xi_109*xib_5_collide))));
                  const float rr_1_b_collide = xi_110;
                  const float xi_114 = rr_1_b_collide*0.5f;
                  const float xi_156 = rr_1_b_collide*0.125f;
                  const float xi_157 = xi_156*xib_19_collide;
                  const float xi_419 = rr_1_b_collide*0.25f;
                  const float rr_2_b_collide = xib_0_collide*(xi_110*-2.0f + 4.0f);
                  const float xi_119 = rr_2_b_collide*xi_117 + xi_117 - xi_118*xib_19_collide;
                  const float xi_133 = rr_2_b_collide*xi_132 - xi_118*xib_14_collide + xi_132;
                  const float xi_142 = rr_2_b_collide*xi_141 - xi_118*xib_8_collide + xi_141;
                  const float xi_144 = rr_2_b_collide*0.041666666666666664f;
                  const float xi_145 = xi_144*xib_14_collide;
                  const float xi_146 = xi_143 - xi_145;
                  const float xi_149 = xi_144*xib_19_collide;
                  const float xi_150 = -xi_148 + xi_149;
                  const float xi_162 = -xi_143 + xi_145;
                  const float xi_164 = xi_148 - xi_149;
                  const float xi_167 = xi_144*xib_8_collide;
                  const float xi_168 = xi_166 - xi_167;
                  const float xi_179 = -xi_166 + xi_167;
                  const float xi_394 = rr_2_b_collide*0.5f;
                  const float xi_422 = rr_2_b_collide*0.25f;
                  const float rho_b_collide = xib_18_collide;
                  const float xi_371 = rho_b_collide*-0.1111111111111111f;
                  const float xi_384 = rho_b_collide*-0.33333333333333331f;
                  const float xi_386 = xi_375 + xi_384;
                  const float xi_389 = rho_b_collide*0.33333333333333331f;
                  const float u_0_b_collide = xib_26_collide;
                  const float xi_111 = u_0_b_collide*xib_14_collide;
                  const float xi_124 = xi_111*0.16666666666666666f;
                  const float xi_125 = xi_111*0.083333333333333329f;
                  const float xi_126 = rr_1_b_collide*xi_125;
                  const float xi_127 = -xi_124 + xi_126;
                  const float xi_130 = xi_124 - xi_126;
                  const float xi_131 = rr_1_b_collide*xi_124;
                  const float xi_153 = u_0_b_collide*xi_152;
                  const float xi_158 = u_0_b_collide*xi_157;
                  const float xi_169 = -xi_125;
                  const float xi_170 = rr_1_b_collide*xi_111*0.041666666666666664f;
                  const float xi_182 = u_0_b_collide*xib_8_collide;
                  const float xi_183 = xi_182*0.25f;
                  const float xi_186 = xi_156*xi_182;
                  const float xi_370 = (u_0_b_collide*u_0_b_collide);
                  const float xi_377 = rho_b_collide*xi_370*-0.33333333333333331f;
                  const float xi_385 = rr_1_b_collide*(rho_b_collide*xi_370 - xi_379 - xi_382 - xi_384 - xib_15_collide - xib_7_collide);
                  const float xi_403 = u_0_b_collide*xi_389;
                  const float xi_404 = xi_397 + xi_403 - xib_24_collide + xib_27_collide;
                  const float xi_405 = xi_394*xi_404;
                  const float xi_408 = xi_403 + xi_407 - xib_13_collide + xib_22_collide;
                  const float xi_409 = xi_394*xi_408;
                  const float xi_423 = xi_404*xi_422;
                  const float xi_424 = -xi_423;
                  const float xi_438 = xi_408*xi_422;
                  const float xi_439 = -xi_438;
                  const float u_1_b_collide = xib_3_collide;
                  const float xi_112 = u_1_b_collide*xib_19_collide;
                  const float xi_115 = xi_112*0.16666666666666666f;
                  const float xi_116 = rr_1_b_collide*xi_115;
                  const float xi_134 = xi_112*0.083333333333333329f;
                  const float xi_135 = rr_1_b_collide*xi_134;
                  const float xi_136 = xi_115 - xi_135;
                  const float xi_138 = -xi_115 + xi_135;
                  const float xi_140 = xi_127 + xi_138;
                  const float xi_154 = u_1_b_collide*0.25f;
                  const float xi_155 = xi_154*xib_14_collide;
                  const float xi_159 = u_1_b_collide*xi_156;
                  const float xi_160 = xi_159*xib_14_collide;
                  const float xi_161 = xi_153 + xi_155 - xi_158 - xi_160;
                  const float xi_163 = -xi_153 - xi_155 + xi_158 + xi_160;
                  const float xi_172 = xi_154*xib_8_collide;
                  const float xi_174 = xi_159*xib_8_collide;
                  const float xi_180 = rr_1_b_collide*u_1_b_collide*xib_19_collide*-0.041666666666666664f;
                  const float xi_372 = (u_1_b_collide*u_1_b_collide);
                  const float xi_373 = rho_b_collide*xi_372*-0.33333333333333331f + xi_371;
                  const float xi_387 = rr_1_b_collide*(rho_b_collide*xi_372 - xi_382 - xi_386 - xib_25_collide - xib_6_collide);
                  const float xi_390 = u_1_b_collide*xi_389;
                  const float xi_393 = xi_390 + xi_392 + xib_17_collide - xib_28_collide;
                  const float xi_395 = xi_393*xi_394;
                  const float xi_398 = xi_390 + xi_397 + xib_24_collide - xib_27_collide;
                  const float xi_399 = xi_394*xi_398;
                  const float xi_418 = rho_b_collide*u_1_b_collide;
                  const float xi_420 = xi_419*(u_0_b_collide*xi_418 + xi_381 + xi_396 - xib_11_collide);
                  const float xi_421 = -xi_420;
                  const float xi_426 = xi_398*xi_422;
                  const float xi_430 = xi_393*xi_422;
                  const float xi_435 = -xi_430;
                  const float u_2_b_collide = xib_20_collide;
                  const float xi_113 = u_2_b_collide*xib_8_collide;
                  const float xi_120 = xi_113*0.16666666666666666f;
                  const float xi_121 = xi_113*0.083333333333333329f;
                  const float xi_122 = rr_1_b_collide*xi_121;
                  const float xi_123 = -xi_120 + xi_122;
                  const float xi_128 = xi_123 + xi_127;
                  const float xi_129 = xi_120 - xi_122;
                  const float xi_137 = xi_129 + xi_136;
                  const float xi_139 = rr_1_b_collide*xi_120;
                  const float xi_147 = rr_1_b_collide*u_2_b_collide*xib_8_collide*-0.041666666666666664f;
                  const float xi_151 = xi_121 + xi_140 + xi_147 + xi_150;
                  const float xi_165 = xi_121 + xi_140 + xi_147 + xi_164;
                  const float xi_171 = xi_137 + xi_164 + xi_169 + xi_170;
                  const float xi_173 = u_2_b_collide*xi_152;
                  const float xi_175 = u_2_b_collide*xi_157;
                  const float xi_176 = xi_172 + xi_173 - xi_174 - xi_175;
                  const float xi_177 = -xi_172 - xi_173 + xi_174 + xi_175;
                  const float xi_178 = xi_137 + xi_150 + xi_169 + xi_170;
                  const float xi_181 = xi_128 + xi_134 + xi_146 + xi_180;
                  const float xi_184 = u_2_b_collide*xib_14_collide;
                  const float xi_185 = xi_184*0.25f;
                  const float xi_187 = xi_156*xi_184;
                  const float xi_188 = xi_183 + xi_185 - xi_186 - xi_187;
                  const float xi_189 = -xi_183 - xi_185 + xi_186 + xi_187;
                  const float xi_190 = xi_128 + xi_134 + xi_162 + xi_180;
                  const float xi_368 = (u_2_b_collide*u_2_b_collide);
                  const float xi_369 = rho_b_collide*xi_368*-0.33333333333333331f;
                  const float xi_376 = rr_1_b_collide*(rho_b_collide*xi_370*-0.16666666666666666f - xi_369 - xi_373 - xi_375);
                  const float xi_380 = rr_1_b_collide*(rho_b_collide*xi_372*-0.16666666666666666f - xi_369 - xi_371 - xi_377 - xi_379);
                  const float xi_383 = rr_1_b_collide*(rho_b_collide*xi_368*-0.16666666666666666f - xi_373 - xi_377 - xi_382);
                  const float xi_388 = rr_1_b_collide*(rho_b_collide*xi_368 - xi_379 - xi_386 - xib_16_collide - xib_9_collide);
                  const float xi_400 = xi_376*-0.5f;
                  const float xi_401 = xi_383*-0.5f;
                  const float xi_402 = xi_387*0.5f + xi_400 + xi_401;
                  const float xi_410 = xi_380*-0.5f;
                  const float xi_411 = xi_385*0.5f + xi_401 + xi_410;
                  const float xi_412 = u_2_b_collide*xi_389;
                  const float xi_413 = xi_392 + xi_412 - xib_17_collide + xib_28_collide;
                  const float xi_414 = xi_394*xi_413;
                  const float xi_415 = xi_407 + xi_412 + xib_13_collide - xib_22_collide;
                  const float xi_416 = xi_394*xi_415;
                  const float xi_417 = xi_388*0.5f + xi_400 + xi_410;
                  const float xi_425 = xi_383*0.25f;
                  const float xi_427 = xi_425 + xi_426;
                  const float xi_428 = xi_425 - xi_426;
                  const float xi_429 = xi_419*(u_2_b_collide*xi_418 + xi_374 + xi_391 - xib_10_collide);
                  const float xi_431 = xi_376*0.25f;
                  const float xi_432 = xi_413*xi_422;
                  const float xi_433 = xi_431 + xi_432;
                  const float xi_434 = -xi_429;
                  const float xi_436 = xi_419*(rho_b_collide*u_0_b_collide*u_2_b_collide + xi_378 + xi_406 - xib_29_collide);
                  const float xi_437 = -xi_436;
                  const float xi_440 = xi_380*0.25f;
                  const float xi_441 = xi_415*xi_422;
                  const float xi_442 = xi_440 + xi_441;
                  const float xi_443 = xi_431 - xi_432;
                  const float xi_444 = xi_440 - xi_441;
                  const float forceTerm_0_b_collide = xi_111*xi_114 - xi_111 + xi_112*xi_114 - xi_112 + xi_113*xi_114 - xi_113;
                  const float forceTerm_1_b_collide = xi_112*0.33333333333333331f - xi_116 + xi_119 + xi_128;
                  const float forceTerm_2_b_collide = u_1_b_collide*xib_19_collide*0.33333333333333331f - xi_116 - xi_119 - xi_129 - xi_130;
                  const float forceTerm_3_b_collide = u_0_b_collide*xib_14_collide*0.33333333333333331f - xi_131 - xi_133 - xi_137;
                  const float forceTerm_4_b_collide = xi_111*0.33333333333333331f + xi_123 - xi_131 + xi_133 + xi_138;
                  const float forceTerm_5_b_collide = xi_113*0.33333333333333331f - xi_139 + xi_140 + xi_142;
                  const float forceTerm_6_b_collide = u_2_b_collide*xib_8_collide*0.33333333333333331f - xi_130 - xi_136 - xi_139 - xi_142;
                  const float forceTerm_7_b_collide = -xi_146 - xi_151 - xi_161;
                  const float forceTerm_8_b_collide = -xi_151 - xi_162 - xi_163;
                  const float forceTerm_9_b_collide = -xi_146 - xi_163 - xi_165;
                  const float forceTerm_10_b_collide = -xi_161 - xi_162 - xi_165;
                  const float forceTerm_11_b_collide = xi_168 + xi_171 + xi_176;
                  const float forceTerm_12_b_collide = xi_168 + xi_177 + xi_178;
                  const float forceTerm_13_b_collide = -xi_179 - xi_181 - xi_188;
                  const float forceTerm_14_b_collide = -xi_179 - xi_189 - xi_190;
                  const float forceTerm_15_b_collide = xi_171 + xi_177 + xi_179;
                  const float forceTerm_16_b_collide = xi_176 + xi_178 + xi_179;
                  const float forceTerm_17_b_collide = -xi_168 - xi_181 - xi_189;
                  const float forceTerm_18_b_collide = -xi_168 - xi_188 - xi_190;
                  const float tmp_a0 = forceTerm_0_a_collide + xi_199 + xi_203 + xi_206 - xi_209 - xi_210 - xi_211 + xi_267 + xia_7_collide;
                  const float tmp_a1 = forceTerm_1_a_collide - xi_274 - xi_278 + xi_281 + xi_284 + xia_13_collide;
                  const float tmp_a2 = forceTerm_2_a_collide + xi_274 + xi_278 + xi_284 + xi_285 + xia_4_collide;
                  const float tmp_a3 = forceTerm_3_a_collide + xi_288 + xi_292 + xi_293 + xi_295 + xia_15_collide;
                  const float tmp_a4 = forceTerm_4_a_collide - xi_288 - xi_292 + xi_295 + xi_296 + xia_26_collide;
                  const float tmp_a5 = forceTerm_5_a_collide - xi_299 - xi_301 + xi_302 + xi_303 + xia_29_collide;
                  const float tmp_a6 = forceTerm_6_a_collide + xi_299 + xi_301 + xi_303 + xi_304 + xia_19_collide;
                  const float tmp_a7 = forceTerm_7_a_collide + xi_308 + xi_311 + xi_320 + xi_323 + xia_14_collide;
                  const float tmp_a8 = forceTerm_8_a_collide + xi_307 + xi_310 + xi_323 + xi_329 + xia_8_collide;
                  const float tmp_a9 = forceTerm_9_a_collide + xi_307 + xi_311 + xi_330 + xi_331 + xia_21_collide;
                  const float tmp_a10 = forceTerm_10_a_collide + xi_308 + xi_310 + xi_331 + xi_332 + xia_25_collide;
                  const float tmp_a11 = forceTerm_11_a_collide + xi_333 + xi_334 + xi_339 + xi_342 + xia_6_collide;
                  const float tmp_a12 = forceTerm_12_a_collide + xi_342 + xi_343 + xi_344 + xi_348 + xia_16_collide;
                  const float tmp_a13 = forceTerm_13_a_collide + xi_350 + xi_352 + xi_355 + xi_358 + xia_9_collide;
                  const float tmp_a14 = forceTerm_14_a_collide + xi_349 + xi_351 + xi_358 + xi_361 + xia_12_collide;
                  const float tmp_a15 = forceTerm_15_a_collide + xi_334 + xi_343 + xi_362 + xi_363 + xia_20_collide;
                  const float tmp_a16 = forceTerm_16_a_collide + xi_333 + xi_344 + xi_363 + xi_364 + xia_18_collide;
                  const float tmp_a17 = forceTerm_17_a_collide + xi_349 + xi_352 + xi_365 + xi_366 + xia_27_collide;
                  const float tmp_a18 = forceTerm_18_a_collide + xi_350 + xi_351 + xi_366 + xi_367 + xia_23_collide;
                  const float tmp_b0 = forceTerm_0_b_collide + xi_267 + xi_376 + xi_380 + xi_383 - xi_385 - xi_387 - xi_388 + xib_21_collide;
                  const float tmp_b1 = forceTerm_1_b_collide + xi_281 - xi_395 - xi_399 + xi_402 + xib_6_collide;
                  const float tmp_b2 = forceTerm_2_b_collide + xi_285 + xi_395 + xi_399 + xi_402 + xib_25_collide;
                  const float tmp_b3 = forceTerm_3_b_collide + xi_293 + xi_405 + xi_409 + xi_411 + xib_7_collide;
                  const float tmp_b4 = forceTerm_4_b_collide + xi_296 - xi_405 - xi_409 + xi_411 + xib_15_collide;
                  const float tmp_b5 = forceTerm_5_b_collide + xi_302 - xi_414 - xi_416 + xi_417 + xib_16_collide;
                  const float tmp_b6 = forceTerm_6_b_collide + xi_304 + xi_414 + xi_416 + xi_417 + xib_9_collide;
                  const float tmp_b7 = forceTerm_7_b_collide + xi_320 + xi_421 + xi_424 + xi_427 + xib_27_collide;
                  const float tmp_b8 = forceTerm_8_b_collide + xi_329 + xi_420 + xi_423 + xi_427 + xib_4_collide;
                  const float tmp_b9 = forceTerm_9_b_collide + xi_330 + xi_420 + xi_424 + xi_428 + xib_11_collide;
                  const float tmp_b10 = forceTerm_10_b_collide + xi_332 + xi_421 + xi_423 + xi_428 + xib_24_collide;
                  const float tmp_b11 = forceTerm_11_b_collide + xi_339 + xi_429 + xi_430 + xi_433 + xib_12_collide;
                  const float tmp_b12 = forceTerm_12_b_collide + xi_348 + xi_433 + xi_434 + xi_435 + xib_17_collide;
                  const float tmp_b13 = forceTerm_13_b_collide + xi_355 + xi_437 + xi_439 + xi_442 + xib_22_collide;
                  const float tmp_b14 = forceTerm_14_b_collide + xi_361 + xi_436 + xi_438 + xi_442 + xib_23_collide;
                  const float tmp_b15 = forceTerm_15_b_collide + xi_362 + xi_430 + xi_434 + xi_443 + xib_28_collide;
                  const float tmp_b16 = forceTerm_16_b_collide + xi_364 + xi_429 + xi_435 + xi_443 + xib_10_collide;
                  const float tmp_b17 = forceTerm_17_b_collide + xi_365 + xi_436 + xi_439 + xi_444 + xib_29_collide;
                  const float tmp_b18 = forceTerm_18_b_collide + xi_367 + xi_437 + xi_438 + xi_444 + xib_13_collide;
                  const float xirecolor_0 = tmp_a0 + tmp_b0;
                  const float xirecolor_1 = _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0] + _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0];
                  const float xirecolor_2 = ((1.0f) / (xirecolor_1));
                  const float xi_466 = xirecolor_2*_data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0];
                  const float xirecolor_3 = xirecolor_2*_data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0];
                  const float xirecolor_4 = tmp_a1 + tmp_b1;
                  const float xirecolor_5 = xi_222;
                  const float xirecolor_6 = xi_220;
                  const float xirecolor_7 = xi_240;
                  const float xi_445 = -xirecolor_7;
                  const float xirecolor_8 = xi_445;
                  const float xirecolor_9 = xi_227;
                  const float xirecolor_10 = -xirecolor_9;
                  const float xirecolor_11 = xi_229;
                  const float xi_446 = xirecolor_10 + xirecolor_11;
                  const float xirecolor_12 = xi_242;
                  const float xi_447 = xirecolor_12 + xirecolor_8;
                  const float xirecolor_13 = xi_446 + xi_447;
                  const float xirecolor_14 = xi_244;
                  const float xirecolor_15 = xi_245;
                  const float xi_448 = xirecolor_14 - xirecolor_15;
                  const float xirecolor_16 = -xi_448;
                  const float xirecolor_17 = xi_214;
                  const float xirecolor_18 = xi_215;
                  const float xirecolor_19 = xi_216;
                  const float xi_449 = xirecolor_18 + xirecolor_19;
                  const float xirecolor_20 = xi_217;
                  const float xirecolor_21 = xi_449 - xirecolor_17 + xirecolor_20;
                  const float xirecolor_22 = xi_213;
                  const float xirecolor_23 = xi_212;
                  const float xi_450 = xirecolor_22 - xirecolor_23;
                  const float xirecolor_24 = -xi_450;
                  const float xirecolor_25 = xi_248;
                  const float xirecolor_26 = -xirecolor_25;
                  const float xirecolor_27 = xi_234;
                  const float xirecolor_28 = xi_233;
                  const float xirecolor_29 = xi_250;
                  const float xi_451 = xirecolor_26 + xirecolor_29;
                  const float xirecolor_30 = xi_451 - xirecolor_27 + xirecolor_28;
                  const float xirecolor_31 = xirecolor_13 + xirecolor_16 + xirecolor_21 + xirecolor_24 + xirecolor_30 - xirecolor_5 + xirecolor_6;
                  const float xirecolor_32 = xi_219;
                  const float xirecolor_33 = xi_221;
                  const float xi_452 = xirecolor_32 + xirecolor_33;
                  const float xirecolor_34 = xi_452 + xirecolor_5 - xirecolor_6;
                  const float xirecolor_35 = xi_224;
                  const float xirecolor_36 = xi_225;
                  const float xi_453 = xirecolor_35 - xirecolor_36;
                  const float xirecolor_37 = -xi_453;
                  const float xirecolor_38 = xi_446;
                  const float xirecolor_39 = xi_231;
                  const float xirecolor_40 = xi_232;
                  const float xi_454 = xirecolor_39 + xirecolor_40;
                  const float xirecolor_41 = xirecolor_27 - xirecolor_28;
                  const float xirecolor_42 = xi_454 + xirecolor_41;
                  const float xirecolor_43 = xirecolor_38 + xirecolor_42;
                  const float xirecolor_44 = xi_450 + xirecolor_21 + xirecolor_34 + xirecolor_37 + xirecolor_43;
                  const float xirecolor_45 = xi_448;
                  const float xirecolor_46 = xi_447;
                  const float xirecolor_47 = xi_256;
                  const float xirecolor_48 = -xirecolor_47;
                  const float xirecolor_49 = xi_258;
                  const float xi_455 = xirecolor_48 + xirecolor_49;
                  const float xirecolor_50 = xi_453 + xi_455;
                  const float xirecolor_51 = xi_449 + xirecolor_17 - xirecolor_20 + xirecolor_24 + xirecolor_34 + xirecolor_45 + xirecolor_46 + xirecolor_50;
                  const float xirecolor_52 = powf((xirecolor_31*xirecolor_31) + (xirecolor_44*xirecolor_44) + (xirecolor_51*xirecolor_51), 0.5f);
                  const float xirecolor_53 = ((1.0f) / (xirecolor_52));
                  const float xi_456 = xirecolor_31*xirecolor_53;
                  const float xi_457 = xirecolor_51*xirecolor_53;
                  const float xi_458 = xirecolor_44*xirecolor_53;
                  const bool xirecolor_54 = xirecolor_52 > 0.0f;
                  const float xirecolor_55 = beta*((1.0f) / ((xirecolor_1*xirecolor_1)))*_data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0]*_data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0];
                  const float xirecolor_56 = xirecolor_55*(0.055555555555555552f*_data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0] + 0.055555555555555552f*_data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]);
                  const float xirecolor_57 = xirecolor_56*((xirecolor_54) ? (xi_456): (0.0f));
                  const float xirecolor_58 = tmp_a2 + tmp_b2;
                  const float xirecolor_59 = xirecolor_56*((xirecolor_54) ? (-xi_456): (0.0f));
                  const float xirecolor_60 = tmp_a3 + tmp_b3;
                  const float xirecolor_61 = xirecolor_56*((xirecolor_54) ? (-xi_457): (0.0f));
                  const float xirecolor_62 = tmp_a4 + tmp_b4;
                  const float xirecolor_63 = xirecolor_56*((xirecolor_54) ? (xi_457): (0.0f));
                  const float xirecolor_64 = tmp_a5 + tmp_b5;
                  const float xirecolor_65 = xirecolor_56*((xirecolor_54) ? (-xi_458): (0.0f));
                  const float xirecolor_66 = tmp_a6 + tmp_b6;
                  const float xirecolor_67 = xirecolor_56*((xirecolor_54) ? (xi_458): (0.0f));
                  const float xirecolor_68 = tmp_a7 + tmp_b7;
                  const float xirecolor_69 = xi_312;
                  const float xirecolor_70 = xi_313;
                  const float xirecolor_71 = xi_452;
                  const float xirecolor_72 = xirecolor_50 + xirecolor_71;
                  const float xirecolor_73 = xi_317;
                  const float xirecolor_74 = xirecolor_25 - xirecolor_29 + xirecolor_73;
                  const float xirecolor_75 = xi_314 - xirecolor_11 + xirecolor_41 - xirecolor_69 + xirecolor_70 + xirecolor_72 + xirecolor_74 + xirecolor_9;
                  const float xirecolor_76 = xirecolor_53*0.70710678118654757f;
                  const float xi_459 = xirecolor_75*xirecolor_76;
                  const float xirecolor_77 = xirecolor_55*(0.027777777777777776f*_data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0] + 0.027777777777777776f*_data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]);
                  const float xirecolor_78 = xirecolor_77*((xirecolor_54) ? (-xi_459): (0.0f));
                  const float xirecolor_79 = tmp_a8 + tmp_b8;
                  const float xirecolor_80 = xi_324;
                  const float xirecolor_81 = xi_325;
                  const float xi_460 = -xirecolor_80 + xirecolor_81;
                  const float xirecolor_82 = xi_326;
                  const float xirecolor_83 = xi_327 + xi_460 + xirecolor_30 + xirecolor_38 + xirecolor_72 + xirecolor_82;
                  const float xi_461 = xirecolor_76*xirecolor_83;
                  const float xirecolor_84 = xirecolor_77*((xirecolor_54) ? (xi_461): (0.0f));
                  const float xirecolor_85 = tmp_a9 + tmp_b9;
                  const float xirecolor_86 = xirecolor_77*((xirecolor_54) ? (-xi_461): (0.0f));
                  const float xirecolor_87 = tmp_a10 + tmp_b10;
                  const float xirecolor_88 = xirecolor_77*((xirecolor_54) ? (xi_459): (0.0f));
                  const float xirecolor_89 = tmp_a11 + tmp_b11;
                  const float xirecolor_90 = xi_454 + xirecolor_37 + xirecolor_71;
                  const float xirecolor_91 = -xi_445 - xi_460 - xirecolor_12;
                  const float xirecolor_92 = xi_335 + xirecolor_45 + xirecolor_74 + xirecolor_90 + xirecolor_91;
                  const float xi_462 = xirecolor_76*xirecolor_92;
                  const float xirecolor_93 = xirecolor_77*((xirecolor_54) ? (-xi_462): (0.0f));
                  const float xirecolor_94 = tmp_a12 + tmp_b12;
                  const float xirecolor_95 = xirecolor_16 + xirecolor_69 - xirecolor_70;
                  const float xirecolor_96 = xi_345 + xi_451 + xirecolor_46 + xirecolor_82 + xirecolor_90 + xirecolor_95;
                  const float xi_463 = xirecolor_76*xirecolor_96;
                  const float xirecolor_97 = xirecolor_77*((xirecolor_54) ? (-xi_463): (0.0f));
                  const float xirecolor_98 = tmp_a13 + tmp_b13;
                  const float xirecolor_99 = xi_353 + xi_455 + xirecolor_13 + xirecolor_42 + xirecolor_45 + xirecolor_73 + xirecolor_82;
                  const float xi_464 = xirecolor_76*xirecolor_99;
                  const float xirecolor_100 = xirecolor_77*((xirecolor_54) ? (-xi_464): (0.0f));
                  const float xirecolor_101 = tmp_a14 + tmp_b14;
                  const float xirecolor_102 = xi_359 + xirecolor_43 + xirecolor_47 - xirecolor_49 + xirecolor_91 + xirecolor_95;
                  const float xi_465 = xirecolor_102*xirecolor_76;
                  const float xirecolor_103 = xirecolor_77*((xirecolor_54) ? (-xi_465): (0.0f));
                  const float xirecolor_104 = tmp_a15 + tmp_b15;
                  const float xirecolor_105 = xirecolor_77*((xirecolor_54) ? (xi_463): (0.0f));
                  const float xirecolor_106 = tmp_a16 + tmp_b16;
                  const float xirecolor_107 = xirecolor_77*((xirecolor_54) ? (xi_462): (0.0f));
                  const float xirecolor_108 = tmp_a17 + tmp_b17;
                  const float xirecolor_109 = xirecolor_77*((xirecolor_54) ? (xi_465): (0.0f));
                  const float xirecolor_110 = tmp_a18 + tmp_b18;
                  const float xirecolor_111 = xirecolor_77*((xirecolor_54) ? (xi_464): (0.0f));
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + ctr_0] = xirecolor_0*xirecolor_3;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_4 + xirecolor_57;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_58 + xirecolor_59;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_60 + xirecolor_61;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_62 + xirecolor_63;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_64 + xirecolor_65;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_66 + xirecolor_67;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_68 + xirecolor_78;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_79 + xirecolor_84;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_85 + xirecolor_86;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_87 + xirecolor_88;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_89 + xirecolor_93;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_94 + xirecolor_97;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3 + ctr_0] = xirecolor_100 + xirecolor_3*xirecolor_98;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3 + ctr_0] = xirecolor_101*xirecolor_3 + xirecolor_103;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3 + ctr_0] = xirecolor_104*xirecolor_3 + xirecolor_105;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3 + ctr_0] = xirecolor_106*xirecolor_3 + xirecolor_107;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3 + ctr_0] = xirecolor_108*xirecolor_3 + xirecolor_109;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3 + ctr_0] = xirecolor_110*xirecolor_3 + xirecolor_111;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + ctr_0] = xi_466*xirecolor_0;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3 + ctr_0] = xi_466*xirecolor_4 - xirecolor_57;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3 + ctr_0] = xi_466*xirecolor_58 - xirecolor_59;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3 + ctr_0] = xi_466*xirecolor_60 - xirecolor_61;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3 + ctr_0] = xi_466*xirecolor_62 - xirecolor_63;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3 + ctr_0] = xi_466*xirecolor_64 - xirecolor_65;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3 + ctr_0] = xi_466*xirecolor_66 - xirecolor_67;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3 + ctr_0] = xi_466*xirecolor_68 - xirecolor_78;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3 + ctr_0] = xi_466*xirecolor_79 - xirecolor_84;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3 + ctr_0] = xi_466*xirecolor_85 - xirecolor_86;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3 + ctr_0] = xi_466*xirecolor_87 - xirecolor_88;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3 + ctr_0] = xi_466*xirecolor_89 - xirecolor_93;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3 + ctr_0] = xi_466*xirecolor_94 - xirecolor_97;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3 + ctr_0] = xi_466*xirecolor_98 - xirecolor_100;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3 + ctr_0] = xi_466*xirecolor_101 - xirecolor_103;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3 + ctr_0] = xi_466*xirecolor_104 - xirecolor_105;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3 + ctr_0] = xi_466*xirecolor_106 - xirecolor_107;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3 + ctr_0] = xi_466*xirecolor_108 - xirecolor_109;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3 + ctr_0] = xi_466*xirecolor_110 - xirecolor_111;
               }
            }
         }
      }
   }
}
}


void ColorGradientCollideSweepSinglePrecisionAVX::run(IBlock * block)
{
   
    auto phasefield = block->getData< field::GhostLayerField<float, 1> >(phasefieldID);
    auto force_b = block->getData< field::GhostLayerField<float, 3> >(force_bID);
    auto pdfs_a = block->getData< field::GhostLayerField<float, 19> >(pdfs_aID);
    auto force_a = block->getData< field::GhostLayerField<float, 3> >(force_aID);
    auto pdfs_b = block->getData< field::GhostLayerField<float, 19> >(pdfs_bID);
    auto rho_a = block->getData< field::GhostLayerField<float, 1> >(rho_aID);
    auto velocity = block->getData< field::GhostLayerField<float, 3> >(velocityID);
    auto rho_b = block->getData< field::GhostLayerField<float, 1> >(rho_bID);

    auto & beta = this->beta_;
    auto & sigma = this->sigma_;
    auto & omega_shear_b = this->omega_shear_b_;
    auto & omega_shear_a = this->omega_shear_a_;
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force_a->nrOfGhostLayers()))
    float * RESTRICT const _data_force_a = force_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force_b->nrOfGhostLayers()))
    float * RESTRICT const _data_force_b = force_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(force_b->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_a->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs_a = pdfs_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) pdfs_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_b->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs_b = pdfs_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) pdfs_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(phasefield->nrOfGhostLayers()))
    float * RESTRICT const _data_phasefield = phasefield->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) phasefield->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho_a->nrOfGhostLayers()))
    float * RESTRICT const _data_rho_a = rho_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) rho_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho_b->nrOfGhostLayers()))
    float * RESTRICT const _data_rho_b = rho_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) rho_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(velocity->nrOfGhostLayers()))
    float * RESTRICT const _data_velocity = velocity->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) velocity->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->xSizeWithGhostLayer(), int64_t(int64_c(force_a->xSize()) + 2))
    const int64_t _size_force_a_0 = int64_t(int64_c(force_a->xSize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->ySizeWithGhostLayer(), int64_t(int64_c(force_a->ySize()) + 2))
    const int64_t _size_force_a_1 = int64_t(int64_c(force_a->ySize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->zSizeWithGhostLayer(), int64_t(int64_c(force_a->zSize()) + 2))
    const int64_t _size_force_a_2 = int64_t(int64_c(force_a->zSize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    const int64_t _stride_force_a_1 = int64_t(force_a->yStride());
    const int64_t _stride_force_a_2 = int64_t(force_a->zStride());
    const int64_t _stride_force_a_3 = int64_t(1 * int64_t(force_a->fStride()));
    const int64_t _stride_force_b_1 = int64_t(force_b->yStride());
    const int64_t _stride_force_b_2 = int64_t(force_b->zStride());
    const int64_t _stride_force_b_3 = int64_t(1 * int64_t(force_b->fStride()));
    const int64_t _stride_pdfs_a_1 = int64_t(pdfs_a->yStride());
    const int64_t _stride_pdfs_a_2 = int64_t(pdfs_a->zStride());
    const int64_t _stride_pdfs_a_3 = int64_t(1 * int64_t(pdfs_a->fStride()));
    const int64_t _stride_pdfs_b_1 = int64_t(pdfs_b->yStride());
    const int64_t _stride_pdfs_b_2 = int64_t(pdfs_b->zStride());
    const int64_t _stride_pdfs_b_3 = int64_t(1 * int64_t(pdfs_b->fStride()));
    const int64_t _stride_phasefield_1 = int64_t(phasefield->yStride());
    const int64_t _stride_phasefield_2 = int64_t(phasefield->zStride());
    const int64_t _stride_rho_a_1 = int64_t(rho_a->yStride());
    const int64_t _stride_rho_a_2 = int64_t(rho_a->zStride());
    const int64_t _stride_rho_b_1 = int64_t(rho_b->yStride());
    const int64_t _stride_rho_b_2 = int64_t(rho_b->zStride());
    const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
    const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
    const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
    internal_196de4e60a60b4c0f3c29ed6234bf8bc::colorgradientcollidesweepsingleprecisionavx_colorgradientcollidesweepsingleprecisionavx(_data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_b, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_force_a_0, _size_force_a_1, _size_force_a_2, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3, beta, omega_shear_a, omega_shear_b, sigma);
    
}


void ColorGradientCollideSweepSinglePrecisionAVX::runOnCellInterval(const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers, IBlock * block)
{
   
    CellInterval ci = globalCellInterval;
    CellInterval blockBB = blocks->getBlockCellBB( *block);
    blockBB.expand( ghostLayers );
    ci.intersect( blockBB );
    blocks->transformGlobalToBlockLocalCellInterval( ci, *block );
    if( ci.empty() )
        return;

    auto phasefield = block->getData< field::GhostLayerField<float, 1> >(phasefieldID);
    auto force_b = block->getData< field::GhostLayerField<float, 3> >(force_bID);
    auto pdfs_a = block->getData< field::GhostLayerField<float, 19> >(pdfs_aID);
    auto force_a = block->getData< field::GhostLayerField<float, 3> >(force_aID);
    auto pdfs_b = block->getData< field::GhostLayerField<float, 19> >(pdfs_bID);
    auto rho_a = block->getData< field::GhostLayerField<float, 1> >(rho_aID);
    auto velocity = block->getData< field::GhostLayerField<float, 3> >(velocityID);
    auto rho_b = block->getData< field::GhostLayerField<float, 1> >(rho_bID);

    auto & beta = this->beta_;
    auto & sigma = this->sigma_;
    auto & omega_shear_b = this->omega_shear_b_;
    auto & omega_shear_a = this->omega_shear_a_;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(force_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(force_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(force_a->nrOfGhostLayers()))
    float * RESTRICT const _data_force_a = force_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    float * RESTRICT const _data_force_b = force_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(force_b->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs_a = pdfs_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) pdfs_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs_b = pdfs_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) pdfs_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    float * RESTRICT const _data_phasefield = phasefield->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) phasefield->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    float * RESTRICT const _data_rho_a = rho_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) rho_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    float * RESTRICT const _data_rho_b = rho_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) rho_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    float * RESTRICT const _data_velocity = velocity->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) velocity->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 2))
    const int64_t _size_force_a_0 = int64_t(int64_c(ci.xSize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 2))
    const int64_t _size_force_a_1 = int64_t(int64_c(ci.ySize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 2))
    const int64_t _size_force_a_2 = int64_t(int64_c(ci.zSize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    const int64_t _stride_force_a_1 = int64_t(force_a->yStride());
    const int64_t _stride_force_a_2 = int64_t(force_a->zStride());
    const int64_t _stride_force_a_3 = int64_t(1 * int64_t(force_a->fStride()));
    const int64_t _stride_force_b_1 = int64_t(force_b->yStride());
    const int64_t _stride_force_b_2 = int64_t(force_b->zStride());
    const int64_t _stride_force_b_3 = int64_t(1 * int64_t(force_b->fStride()));
    const int64_t _stride_pdfs_a_1 = int64_t(pdfs_a->yStride());
    const int64_t _stride_pdfs_a_2 = int64_t(pdfs_a->zStride());
    const int64_t _stride_pdfs_a_3 = int64_t(1 * int64_t(pdfs_a->fStride()));
    const int64_t _stride_pdfs_b_1 = int64_t(pdfs_b->yStride());
    const int64_t _stride_pdfs_b_2 = int64_t(pdfs_b->zStride());
    const int64_t _stride_pdfs_b_3 = int64_t(1 * int64_t(pdfs_b->fStride()));
    const int64_t _stride_phasefield_1 = int64_t(phasefield->yStride());
    const int64_t _stride_phasefield_2 = int64_t(phasefield->zStride());
    const int64_t _stride_rho_a_1 = int64_t(rho_a->yStride());
    const int64_t _stride_rho_a_2 = int64_t(rho_a->zStride());
    const int64_t _stride_rho_b_1 = int64_t(rho_b->yStride());
    const int64_t _stride_rho_b_2 = int64_t(rho_b->zStride());
    const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
    const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
    const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
    internal_196de4e60a60b4c0f3c29ed6234bf8bc::colorgradientcollidesweepsingleprecisionavx_colorgradientcollidesweepsingleprecisionavx(_data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_b, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_force_a_0, _size_force_a_1, _size_force_a_2, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3, beta, omega_shear_a, omega_shear_b, sigma);
    
}



} // namespace pystencils
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning pop
#endif
