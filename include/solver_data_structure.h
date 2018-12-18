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

#ifndef FLUID_DATA_STRUCTURE_H
#define FLUID_DATA_STRUCTURE_H

#include <stddef.h>

struct twoDCons {
  size_t M;
  size_t N;
  void *velX;
  void *velY;
  void *dens;
  void *velX_prev;
  void *velY_prev;
  void *dens_prev;
};

struct threeDCons {
  size_t M;
  size_t N;
  size_t O;
  void *velX;
  void *velY;
  void *velZ;
  void *dens;
  void *velX_prev;
  void *velY_prev;
  void *velZ_prev;
  void *dens_prev;
};

struct twoDSim {
  size_t num_dens, num_vel;
  FLOAT_TYPE *dens_add;
  size_t *dens_sourceX;
  size_t *dens_sourceY;
  size_t *dens_frequency;
  size_t *vel_sourceX;
  size_t *vel_sourceY;
  FLOAT_TYPE *vel_angleXY;
  FLOAT_TYPE *vel_force;
  size_t *vel_frequency;
};

struct threeDSim {
  size_t num_dens, num_vel;
  FLOAT_TYPE *dens_add;
  size_t *dens_sourceX;
  size_t *dens_sourceY;
  size_t *dens_sourceZ;
  size_t *dens_frequency;
  size_t *vel_sourceX;
  size_t *vel_sourceY;
  size_t *vel_sourceZ;
  FLOAT_TYPE *vel_angleXY;
  FLOAT_TYPE *vel_angleXZ;
  FLOAT_TYPE *vel_force;
  size_t *vel_frequency;
};

enum solverDimensions {
  solver_dimension_twoD,
  solver_dimension_threeD,
};

struct solverData {
  enum solverDimensions dim;
  FLOAT_TYPE dt, diff, visc;
  union {
    struct twoDCons twoD;
    struct threeDCons threeD;
  } cons_var;
  union {
    struct twoDSim twoD;
    struct threeDSim threeD;
  } sim_decl;
};

#define arrayCast2DFloat(arr, name, X, Y)                                      \
  arrayCast2D(arr, name, X, Y, FLOAT_TYPE)

#define arrayCast2D(arr, name, X, Y, type) type(*name)[Y] = ((type(*)[Y])arr);

#define arrayCast3DFloat(arr, name, X, Y, Z)                                   \
  arrayCast3D(arr, name, X, Y, Z, FLOAT_TYPE)

#define arrayCast3D(arr, name, X, Y, Z, type)                                  \
  type(*name)[Y][Z] = ((type(*)[Y][Z])arr);

#define consVarLocals2D(solverData)                                            \
  size_t M = solverData->cons_var.twoD.M;                                      \
  size_t N = solverData->cons_var.twoD.N;                                      \
  arrayCast2DFloat(solverData->cons_var.twoD.velX, velX,                       \
                   solverData->cons_var.twoD.M, solverData->cons_var.twoD.N);  \
  arrayCast2DFloat(solverData->cons_var.twoD.velY, velY,                       \
                   solverData->cons_var.twoD.M, solverData->cons_var.twoD.N);  \
  arrayCast2DFloat(solverData->cons_var.twoD.dens, dens,                       \
                   solverData->cons_var.twoD.M, solverData->cons_var.twoD.N);  \
  arrayCast2DFloat(solverData->cons_var.twoD.velX_prev, velX_prev,             \
                   solverData->cons_var.twoD.M, solverData->cons_var.twoD.N);  \
  arrayCast2DFloat(solverData->cons_var.twoD.velY_prev, velY_prev,             \
                   solverData->cons_var.twoD.M, solverData->cons_var.twoD.N);  \
  arrayCast2DFloat(solverData->cons_var.twoD.dens_prev, dens_prev,             \
                   solverData->cons_var.twoD.M, solverData->cons_var.twoD.N);  \
  hideUnusedConsValLocals2D()

#define hideUnusedConsValLocals2D()                                            \
  (void)M;                                                                     \
  (void)N;                                                                     \
  (void)velX;                                                                  \
  (void)velY;                                                                  \
  (void)dens;                                                                  \
  (void)velX_prev;                                                             \
  (void)velY_prev;                                                             \
  (void)dens_prev;

#define consVarLocals3D(solverData)                                            \
  size_t M = solverData->cons_var.threeD.M;                                    \
  size_t N = solverData->cons_var.threeD.N;                                    \
  size_t O = solverData->cons_var.threeD.O;                                    \
  arrayCast3DFloat(                                                            \
      solverData->cons_var.threeD.velX, velX, solverData->cons_var.threeD.M,   \
      solverData->cons_var.threeD.N, solverData->cons_var.threeD.O);           \
  arrayCast3DFloat(                                                            \
      solverData->cons_var.threeD.velY, velY, solverData->cons_var.threeD.M,   \
      solverData->cons_var.threeD.N, solverData->cons_var.threeD.O);           \
  arrayCast3DFloat(                                                            \
      solverData->cons_var.threeD.velZ, velZ, solverData->cons_var.threeD.M,   \
      solverData->cons_var.threeD.N, solverData->cons_var.threeD.O);           \
  arrayCast3DFloat(                                                            \
      solverData->cons_var.threeD.dens, dens, solverData->cons_var.threeD.M,   \
      solverData->cons_var.threeD.N, solverData->cons_var.threeD.O);           \
  arrayCast3DFloat(solverData->cons_var.threeD.velX_prev, velX_prev,           \
                   solverData->cons_var.threeD.M,                              \
                   solverData->cons_var.threeD.N,                              \
                   solverData->cons_var.threeD.O);                             \
  arrayCast3DFloat(solverData->cons_var.threeD.velY_prev, velY_prev,           \
                   solverData->cons_var.threeD.M,                              \
                   solverData->cons_var.threeD.N,                              \
                   solverData->cons_var.threeD.O);                             \
  arrayCast3DFloat(solverData->cons_var.threeD.velZ_prev, velZ_prev,           \
                   solverData->cons_var.threeD.M,                              \
                   solverData->cons_var.threeD.N,                              \
                   solverData->cons_var.threeD.O);                             \
  arrayCast3DFloat(solverData->cons_var.threeD.dens_prev, dens_prev,           \
                   solverData->cons_var.threeD.M,                              \
                   solverData->cons_var.threeD.N,                              \
                   solverData->cons_var.threeD.O);                             \
  hideUnusedConsValLocals3D()

#define hideUnusedConsValLocals3D()                                            \
  (void)M;                                                                     \
  (void)N;                                                                     \
  (void)O;                                                                     \
  (void)velX;                                                                  \
  (void)velY;                                                                  \
  (void)velZ;                                                                  \
  (void)dens;                                                                  \
  (void)velX_prev;                                                             \
  (void)velY_prev;                                                             \
  (void)velZ_prev;                                                             \
  (void)dens_prev;

#define simDeclLocals2D(solverData)                                            \
  size_t num_dens = solverData->sim_decl.twoD.num_dens;                        \
  size_t num_vel = solverData->sim_decl.twoD.num_vel;                          \
  FLOAT_TYPE *dens_add = solverData->sim_decl.twoD.dens_add;                   \
  size_t *dens_sourceX = solverData->sim_decl.twoD.dens_sourceX;               \
  size_t *dens_sourceY = solverData->sim_decl.twoD.dens_sourceY;               \
  size_t *dens_frequency = solverData->sim_decl.twoD.dens_frequency;           \
  size_t *vel_sourceX = solverData->sim_decl.twoD.vel_sourceX;                 \
  size_t *vel_sourceY = solverData->sim_decl.twoD.vel_sourceY;                 \
  FLOAT_TYPE *vel_angleXY = solverData->sim_decl.twoD.vel_angleXY;             \
  FLOAT_TYPE *vel_force = solverData->sim_decl.twoD.vel_force;                 \
  size_t *vel_frequency = solverData->sim_decl.twoD.vel_frequency;             \
  hideUnusedSimDecl2D()

#define hideUnusedSimDecl2D()                                                  \
  (void)num_dens;                                                              \
  (void)num_vel;                                                               \
  (void)dens_sourceX;                                                          \
  (void)dens_sourceY;                                                          \
  (void)dens_frequency;                                                        \
  (void)vel_sourceX;                                                           \
  (void)vel_sourceY;                                                           \
  (void)vel_angleXY;                                                           \
  (void)vel_force;                                                             \
  (void)vel_frequency;                                                         \
  (void)dens_add;

#define simDeclLocals3D(solverData)                                            \
  size_t num_dens = solverData->sim_decl.threeD.num_dens;                      \
  size_t num_vel = solverData->sim_decl.threeD.num_vel;                        \
  FLOAT_TYPE *dens_add = solverData->sim_decl.threeD.dens_add;                 \
  size_t *dens_sourceX = solverData->sim_decl.threeD.dens_sourceX;             \
  size_t *dens_sourceY = solverData->sim_decl.threeD.dens_sourceY;             \
  size_t *dens_sourceZ = solverData->sim_decl.threeD.dens_sourceZ;             \
  size_t *dens_frequency = solverData->sim_decl.threeD.dens_frequency;         \
  size_t *vel_sourceX = solverData->sim_decl.threeD.vel_sourceX;               \
  size_t *vel_sourceY = solverData->sim_decl.threeD.vel_sourceY;               \
  size_t *vel_sourceZ = solverData->sim_decl.threeD.vel_sourceZ;               \
  FLOAT_TYPE *vel_angleXY = solverData->sim_decl.threeD.vel_angleXY;           \
  FLOAT_TYPE *vel_angleXZ = solverData->sim_decl.threeD.vel_angleXZ;           \
  FLOAT_TYPE *vel_force = solverData->sim_decl.threeD.vel_force;               \
  size_t *vel_frequency = solverData->sim_decl.threeD.vel_frequency;           \
  hideUnusedSimDecl3D()

#define hideUnusedSimDecl3D()                                                  \
  (void)num_dens;                                                              \
  (void)num_vel;                                                               \
  (void)dens_sourceX;                                                          \
  (void)dens_sourceY;                                                          \
  (void)dens_sourceZ;                                                          \
  (void)dens_frequency;                                                        \
  (void)vel_sourceX;                                                           \
  (void)vel_sourceY;                                                           \
  (void)vel_sourceZ;                                                           \
  (void)vel_angleXY;                                                           \
  (void)vel_angleXZ;                                                           \
  (void)vel_force;                                                             \
  (void)vel_frequency;                                                         \
  (void)dens_add;

#endif // FLUID_DATA_STRUCTURE_H
