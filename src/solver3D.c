/*
 * Copyright (c) 2018 Maxime Schmitt <max.schmitt@unistra.fr>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 * may be used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdlib.h>
#include <string.h>
#include <tgmath.h>

#include "solver3D.h"

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define max3(a, b, c) max(max(a, b), c)

enum border_condition {
  border_copy_value,
  border_negate_x,
  border_negate_y,
  border_negate_z,
};

static void add_source3D(size_t M, size_t N, size_t O,
                         float x[restrict M][N][O], float s[restrict M][N][O],
                         float dt) {
  for (size_t i = 0; i < M; ++i)
    for (size_t j = 0; j < N; ++j)
      for (size_t k = 0; k < O; ++k)
        x[i][j][k] += dt * s[i][j][k];
}

static void set_bnd3D(size_t M, size_t N, size_t O, enum border_condition b,
                      float x[M][N][O]) {

  if (b == border_negate_x) {
    for (size_t j = 1; j < N - 1; ++j) {
      for (size_t k = 1; k < O - 1; ++k) {
        x[0][j][k] = -x[1][j][k];
        x[M - 1][j][k] = -x[M - 2][j][k];
      }
    }
  } else {
    for (size_t j = 1; j < N - 1; ++j) {
      for (size_t k = 1; k < O - 1; ++k) {
        x[0][j][k] = x[1][j][k];
        x[M - 1][j][k] = x[M - 2][j][k];
      }
    }
  }
  if (b == border_negate_y) {
    for (size_t i = 1; i < M - 1; ++i) {
      for (size_t k = 1; k < O - 1; ++k) {
        x[i][0][k] = -x[i][1][k];
        x[i][N - 1][k] = -x[i][N - 2][k];
      }
    }
  } else {
    for (size_t i = 1; i < M - 1; ++i) {
      for (size_t k = 1; k < O - 1; ++k) {
        x[i][0][k] = x[i][1][k];
        x[i][N - 1][k] = x[i][N - 2][k];
      }
    }
  }
  if (b == border_negate_z) {
    for (size_t i = 1; i < M - 1; ++i) {
      for (size_t j = 1; j < N - 1; ++j) {
        x[i][j][0] = -x[i][j][1];
        x[i][j][O - 1] = -x[i][j][O - 2];
      }
    }
  } else {
    for (size_t i = 1; i < M - 1; ++i) {
      for (size_t j = 1; j < N - 1; ++j) {
        x[i][j][0] = x[i][j][1];
        x[i][j][O - 1] = x[i][j][O - 2];
      }
    }
  }

  const float AThird = 1.f / 3.f;

  x[0][0][0] = AThird * (x[0][0][1] + x[0][1][0] + x[1][0][0]);
  x[0][0][O - 1] = AThird * (x[0][0][O - 2] + x[0][1][O - 1] + x[1][0][O - 1]);
  x[0][N - 1][0] = AThird * (x[0][N - 1][1] + x[0][N - 2][0] + x[1][N - 1][0]);
  x[0][N - 1][O - 1] =
      AThird * (x[0][N - 1][O - 2] + x[0][N - 2][O - 1] + x[1][N - 1][O - 1]);
  x[M - 1][0][0] = AThird * (x[M - 1][0][1] + x[M - 1][1][0] + x[M - 2][0][0]);
  x[M - 1][0][O - 1] =
      AThird * (x[M - 1][0][O - 2] + x[M - 1][1][O - 1] + x[M - 2][0][O - 1]);
  x[M - 1][N - 1][0] =
      AThird * (x[M - 1][N - 1][1] + x[M - 1][N - 2][0] + x[M - 2][N - 1][0]);
  x[M - 1][N - 1][O - 1] =
      AThird * (x[M - 1][N - 1][O - 2] + x[M - 1][N - 2][O - 1] +
                x[M - 2][N - 1][O - 1]);
}

extern double rand_skip_percent;

static void lin_solve3D(size_t M, size_t N, size_t O, enum border_condition b,
                        float x[restrict M][N][O], float x0[restrict M][N][O],
                        float a, float c) {

  for (size_t P = 0; P < 10; P++) {
    for (size_t i = 1; i < M - 1; ++i) {
      for (size_t j = 1; j < N - 1; ++j) {
        for (size_t k = 1; k < O - 1; ++k) {
          if (drand48() < rand_skip_percent)
            continue;
          x[i][j][k] = (x0[i][j][k] + a * (x[i - 1][j][k] + x[i + 1][j][k] +
                                           x[i][j - 1][k] + x[i][j + 1][k] +
                                           x[i][j][k - 1] + x[i][j][k + 1])) /
                       c;
#ifndef NDEBUG
          if (isnan(x[i][j][k]))
            exit(1);
#endif
        }
      }
    }
    set_bnd3D(M, N, O, b, x);
  }
}

static void diffuse3D(size_t M, size_t N, size_t O, enum border_condition b,
                      float x[restrict M][N][O], float x0[restrict M][N][O],
                      float diff, float dt) {
  size_t max_bnd = max3(M - 2, N - 2, O - 2);
  float a = dt * diff * max_bnd * max_bnd * max_bnd;
  lin_solve3D(M, N, O, b, x, x0, a, 1 + 6 * a);
}

static void advect3D(size_t M, size_t N, size_t O, enum border_condition b,
                     float d[restrict M][N][O], float d0[M][N][O],
                     float u[M][N][O], float v[M][N][O], float w[M][N][O],
                     float dt) {

  const float dtx = dt * (M - 2);
  const float dty = dt * (N - 2);
  const float dtz = dt * (O - 2);
  const float AHalf = 1.f / 2.f;

  for (size_t i = 1; i < M - 1; ++i) {
    for (size_t j = 1; j < N - 1; ++j) {
      for (size_t k = 1; k < O - 1; ++k) {

        float x = i - dtx * u[i][j][k];
        if (x < AHalf)
          x = AHalf;
        if (x > (M - 2) + AHalf)
          x = (M - 2) + AHalf;
        size_t i0 = (size_t)x;
        size_t i1 = i0 + 1;

        float y = j - dty * v[i][j][k];
        if (y < AHalf)
          y = AHalf;
        if (y > (N - 2) + AHalf)
          y = (N - 2) + AHalf;
        size_t j0 = (size_t)y;
        size_t j1 = j0 + 1;

        float z = k - dtz * w[i][j][k];
        if (z < AHalf)
          z = AHalf;
        if (z > (O - 2) + AHalf)
          z = (O - 2) + AHalf;
        size_t z0 = (size_t)z;
        size_t z1 = z0 + 1;

        float s1 = x - i0;
        float s0 = 1.f - s1;
        float t1 = y - j0;
        float t0 = 1.f - t1;
        float u1 = z - z0;
        float u0 = 1.f - u1;
        d[i][j][k] = s0 * t0 * (u0 * d0[i0][j0][z0] + u1 * d0[i0][j0][z1]) +
                     s0 * t1 * (u0 * d0[i0][j1][z0] + u1 * d0[i0][j1][z1]) +
                     s1 * t0 * (u0 * d0[i1][j0][z0] + u1 * d0[i1][j0][z1]) +
                     s1 * t1 * (u0 * d0[i1][j1][z0] + u1 * d0[i1][j1][z1]);
      }
    }
  }
  set_bnd3D(M, N, O, b, d);
}

static void project3D(size_t M, size_t N, size_t O, float u[restrict M][N][O],
                      float v[restrict M][N][O], float w[restrict M][N][O],
                      float p[restrict M][N][O], float div[restrict M][N][O]) {

  const float MinusAHalfDividedByM = -1.f / (2.f * (M - 2));
  const float MinusAHalfDividedByN = -1.f / (2.f * (N - 2));
  const float MinusAHalfDividedByO = -1.f / (2.f * (O - 2));
  for (size_t i = 1; i < M - 1; ++i) {
    for (size_t j = 1; j < N - 1; ++j) {
      for (size_t k = 1; k < O - 1; ++k) {
        div[i][j][k] =
            MinusAHalfDividedByM * (u[i + 1][j][k] - u[i - 1][j][k]) +
            MinusAHalfDividedByN * (v[i][j + 1][k] - v[i][j - 1][k]) +
            MinusAHalfDividedByO * (w[i][j][k + 1] - w[i][j][k - 1]);
      }
    }
  }
  set_bnd3D(M, N, O, border_copy_value, div);
  memset(p, 0, sizeof(float[M][N][O]));

  lin_solve3D(M, N, O, border_copy_value, p, div, 1, 6);

  const float MinusAHalfTimesM = (M - 2) / -2.;
  const float MinusAHalfTimesN = (N - 2) / -2.;
  const float MinusAHalfTimesO = (O - 2) / -2.;
  for (size_t i = 1; i < M - 1; ++i) {
    for (size_t j = 1; j < N - 1; ++j) {
      for (size_t k = 1; k < O - 1; ++k) {
        u[i][j][k] += MinusAHalfTimesM * (p[i + 1][j][k] - p[i - 1][j][k]);
        v[i][j][k] += MinusAHalfTimesN * (p[i][j + 1][k] - p[i][j - 1][k]);
        w[i][j][k] += MinusAHalfTimesO * (p[i][j][k + 1] - p[i][j][k - 1]);
      }
    }
  }
  set_bnd3D(M, N, O, border_negate_x, u);
  set_bnd3D(M, N, O, border_negate_y, v);
  set_bnd3D(M, N, O, border_negate_z, w);
}

void dens_step3D(size_t M, size_t N, size_t O, float x[M][N][O],
                 float x0[M][N][O], float u[M][N][O], float v[M][N][O],
                 float w[M][N][O], float diff, float dt) {

  add_source3D(M, N, O, x, x0, dt);

  diffuse3D(M, N, O, border_copy_value, x0, x, diff, dt);

  advect3D(M, N, O, border_copy_value, x, x0, u, v, w, dt);
}

void vel_step3D(size_t M, size_t N, size_t O, float u[M][N][O],
                float v[M][N][O], float w[M][N][O], float u0[M][N][O],
                float v0[M][N][O], float w0[M][N][O], float visc, float dt) {

  add_source3D(M, N, O, u, u0, dt);
  add_source3D(M, N, O, v, v0, dt);
  add_source3D(M, N, O, w, w0, dt);

  diffuse3D(M, N, O, border_negate_x, u0, u, visc, dt);
  diffuse3D(M, N, O, border_negate_y, v0, v, visc, dt);
  diffuse3D(M, N, O, border_negate_z, w0, w, visc, dt);

  project3D(M, N, O, u0, v0, w0, u, v);

  advect3D(M, N, O, border_negate_x, u, u0, u0, v0, w0, dt);
  advect3D(M, N, O, border_negate_y, v, v0, u0, v0, w0, dt);
  advect3D(M, N, O, border_negate_z, w, w0, u0, v0, w0, dt);

  project3D(M, N, O, u, v, w, u0, v0);
}

void add_force_source3D(
    size_t sim_step, size_t M, size_t N, size_t O, float d[restrict M][N][O],
    float u[restrict M][N][O], float v[restrict M][N][O],
    float w[restrict M][N][O], size_t num_dens,
    float dens_add[restrict num_dens], size_t dens_sourceX[restrict num_dens],
    size_t dens_sourceY[restrict num_dens],
    size_t dens_sourceZ[restrict num_dens],
    size_t dens_frequency[restrict num_dens], size_t num_vel,
    float vel_force[restrict num_vel], size_t vel_sourceX[restrict num_vel],
    size_t vel_sourceY[restrict num_vel], size_t vel_sourceZ[restrict num_vel],
    float vel_angleXY[restrict num_vel], float vel_angleXZ[restrict num_vel],
    size_t vel_frequency[restrict num_vel]) {

  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      for (size_t k = 0; k < O; ++k) {
        u[i][j][k] = v[i][j][k] = w[i][j][k] = d[i][j][k] = 0.0f;
      }
    }
  }

  for (size_t i = 0; i < num_dens; ++i) {
    if (sim_step % dens_frequency[i] == 0) {
      d[dens_sourceX[i]][dens_sourceY[i]][dens_sourceZ[i]] = dens_add[i];
    }
  }
  for (size_t i = 0; i < num_vel; ++i) {
    if (sim_step % vel_frequency[i] == 0) {
      float cosValXY = cos(vel_angleXY[i]);
      float sinValXY = sin(vel_angleXY[i]);
      float sinValXZ = sin(vel_angleXZ[i]);

      u[vel_sourceX[i]][vel_sourceY[i]][vel_sourceZ[i]] =
          cosValXY * vel_force[i];
      v[vel_sourceX[i]][vel_sourceY[i]][vel_sourceZ[i]] =
          sinValXY * vel_force[i];
      w[vel_sourceX[i]][vel_sourceY[i]][vel_sourceZ[i]] =
          sinValXZ * vel_force[i];
    }
  }
}
