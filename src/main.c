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

#include <getopt.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "simulation_setup.h"
#include "solver.h"
#include "solver_data_structure.h"
#include "time_measurement.h"

static const struct option opt_options[] = {
    {"verbose", no_argument, 0, 'v'},
    {"help", no_argument, 0, 'h'},
    {"3d", no_argument, 0, 'p'},
    {"setup-file", required_argument, 0, 's'},
    {"output-file", required_argument, 0, 'o'},
    {"rand-skip", required_argument, 0, 'R'},
    {"sort-skip", no_argument, 0, 'S'},
    {0, 0, 0, 0}};

static const char help_string[] =
    "Options:"
    "\n  -h --help        : Print this help"
    "\n  -p --3d          : Use the 3D solver"
    "\n  -v --verbose     : Print simulation avancement percentage"
    "\n  -s --setup-file  : Load a simulation setup file"
    "\n  -o --output-file : Simulation output file";

static const char options[] = ":ps:o:vhR:S";

double rand_skip_percent = 0.;
bool sort_skip = false;

int main(int argc, char **argv) {
  struct solverData solverData = {.dim = solver_dimension_twoD};
  char *config_path = NULL;
  char *output_filename = NULL;
  bool verbose = false;

  while (true) {
    int optchar = getopt_long(argc, argv, options, opt_options, NULL);
    if (optchar == -1)
      break;
    switch (optchar) {
    case 'p':
      solverData.dim = solver_dimension_threeD;
      break;
    case 's':
      config_path = optarg;
      break;
    case 'o':
      output_filename = optarg;
      break;
    case 'v':
      verbose = true;
      break;
    case 'h':
      printf("Usage: %s <options>\n%s\n", argv[0], help_string);
      exit(EXIT_SUCCESS);
    case ':':
      switch (optopt) {
      case 'o':
        output_filename = "fluidData.dat";
        break;
      default:
        fprintf(stderr, "Option %c requires an argument\n", optopt);
        exit(EXIT_FAILURE);
      }
      break;
    case 'R': {
      int sscanf_return = sscanf(optarg, "%lf", &rand_skip_percent);
      if (sscanf_return == EOF || sscanf_return == 0 ||
          rand_skip_percent < 0. || rand_skip_percent > 1.) {
        fprintf(stderr,
                "Please enter a floating point number between 0 and 1 for the "
                "random skip value "
                "instead of \"-%c %s\"\n",
                optchar, optarg);
      }
    } break;
    case 'S':
      sort_skip = true;
      break;
    default:
      fprintf(stderr, "Unrecognized option %c\n", optopt);
      exit(EXIT_FAILURE);
      break;
    }
  }

  if (config_path == NULL) {
    fprintf(stderr, "Error: No configuration file selected\n");
    exit(EXIT_FAILURE);
  }

  srand48(time(NULL));

  size_t num_steps;
  initialize_solver(config_path, &num_steps, &solverData);
  float percentage;
  size_t verbose_step;
  if (num_steps >= 20) {
    verbose_step = num_steps / 20;
    percentage = 100. * verbose_step / (float)num_steps;
  } else {
    percentage = 100. / (float)num_steps;
    verbose_step = 1;
  }
  time_measure startTime, endTime;
  get_current_time(&startTime);
  for (size_t i = 0; i < num_steps; ++i) {
    if (verbose && i % verbose_step == 0) {
      printf("\rSolver avancement %.0f%%", i / verbose_step * percentage);
      fflush(stdout);
    }
    solver_step(i, &solverData);
  }
  if (verbose)
    printf("\rSolver avancement 100%%\n");
  get_current_time(&endTime);
  double elapsed_time = measuring_difftime(startTime, endTime);

  printf("Kernel time %.4fs\n", elapsed_time);

  if (output_filename) {
    print_data_to_file(output_filename, &solverData);
  }
  free_solver_data(&solverData);
}
