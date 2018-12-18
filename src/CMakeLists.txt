add_executable(fluidsolver main.c solver.c solver2D.c solver3D.c simulation_setup.c)
target_include_directories(fluidsolver PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(fluidsolver PRIVATE m)
set_property(TARGET fluidsolver
  PROPERTY C_STANDARD 11)

target_compile_definitions(fluidsolver
  PRIVATE "-DFLOAT_TYPE=${FLUIDSOLVER_FLOAT_TYPE}"
  PRIVATE "-DFLOAT_PRINT_PRECISION=\"${FLUIDSOLVER_PRINT_PRECISION}\"")

# Compile Options
include(compile-flags-helpers)
include(${PROJECT_SOURCE_DIR}/optimization_flags.cmake)

if (DEFINED ADDITIONAL_BENCHMARK_COMPILE_OPTIONS)
  add_compiler_option_to_target_type(fluidsolver Benchmark PRIVATE ${ADDITIONAL_BENCHMARK_COMPILE_OPTIONS})
endif()

foreach(compile_type IN ITEMS Release RelWithDebInfo)
  add_compiler_option_to_target_type(fluidsolver ${compile_type} PRIVATE ${ADDITIONAL_RELEASE_COMPILE_OPTIONS})
  add_linker_option_to_target_type(fluidsolver ${compile_type} PRIVATE ${ADDITIONAL_RELEASE_LINK_OPTIONS})
endforeach()

add_compiler_option_to_target_type(fluidsolver Debug PRIVATE ${ADDITIONAL_DEBUG_COMPILE_OPTIONS})

# Linker Options

if (DEFINED ADDITIONAL_BENCHMARK_LINK_OPTIONS)
  add_linker_option_to_target_type(fluidsolver Benchmark PRIVATE ${ADDITIONAL_BENCHMARK_LINK_OPTIONS})
endif()

add_sanitizers_to_target(fluidsolver Debug PRIVATE address undefined)

include(CheckIPOSupported)
check_ipo_supported(RESULT result)
if((result) AND USE_IPO)
  set_property(TARGET fluidsolver PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

install(TARGETS fluidsolver RUNTIME DESTINATION bin)
