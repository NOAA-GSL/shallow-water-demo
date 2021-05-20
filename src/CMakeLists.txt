
list( APPEND shallow_water_model_src_files
module_netcdf_utils.f90
module_shallow_water_kind.f90
module_shallow_water_geometry_config.f90
module_shallow_water_model_config.f90
module_shallow_water_geometry.f90
module_shallow_water_state.f90
module_shallow_water_model.f90
module_shallow_water_tl.f90
module_shallow_water_adj.f90
)

include(GNUInstallDirs)

#ecbuild_add_library( TARGET     shallow_water
#                     SOURCES    ${shallow_water_model_src_files}
#                     INSTALL_HEADERS LISTED
#                     HEADER_DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME} )

add_library( shallow_water
             ${shallow_water_model_src_files} )

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
    target_compile_options(shallow_water PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:-fallow-argument-mismatch>)
endif()

target_link_libraries(${PROJECT_NAME} PUBLIC NetCDF::NetCDF_C)
target_link_libraries(${PROJECT_NAME} PUBLIC NetCDF::NetCDF_Fortran)
if(OpenMP_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_Fortran)
endif()
target_link_libraries(${PROJECT_NAME} PUBLIC MPI::MPI_Fortran)

#Fortran module output directory for build interface
set(MODULE_DIR module/${PROJECT_NAME}/${CMAKE_Fortran_COMPILER_ID}/${CMAKE_Fortran_COMPILER_VERSION})
set_target_properties(${PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR})
install(DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR}/ DESTINATION ${MODULE_DIR})
target_include_directories(${PROJECT_NAME} INTERFACE
                                $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${MODULE_DIR}>
                                $<INSTALL_INTERFACE:${MODULE_DIR}>)

#ecbuild_add_executable( TARGET  shallow_water.x
#                        SOURCES shallow_water.f90
#                        LIBS    shallow_water )

add_executable(shallow_water.x shallow_water.f90)
target_link_libraries(shallow_water.x ${PROJECT_NAME})
