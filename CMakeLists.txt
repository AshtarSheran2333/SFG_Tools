#made with help of https://www.atomwitch.net/2022/04/02/cmake-fortran.html

cmake_minimum_required(VERSION 3.1)

project(SFG_TOOLS VERSION 2.0.0
		DESCRIPTION "A program suite to calculate SFG spectra from MD trajectories"
		LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_Fortran_EXTENSIONS OFF)

#TODO it would be nice to have a standard code for F18 or something...
#compiler setup
if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
# O2 optimalization (should be default anyway)
# -qopenmp (to enable openmp)
# -fpp (to enable preprocessor)
set(FVERSION "-O2 -qopenmp -fpp")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
# O2 optimalization (should be default anyway)
# -fopenmp (to enable openmp)
# -cpp (to enable preprocessor)
# -ffree-line-length-none (unlimited line lengths)
set(FVERSION "-O2 -fopenmp -cpp -std=gnu -ffree-line-length-none")
endif()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FVERSION}")

#pack all modules into a library
add_library(sfg_lib STATIC
	lib/BINDER_FILE.f90
	lib/BOXDATA.f90
	lib/DENSITY.f90
	lib/DXDRZ_DB.f90
	lib/FRAME_READERS.f90
	lib/GRO_READER.f90
	lib/INSTANTANEOUS_SURFACE.f90
	lib/PRETTY_PRINT.f90
	lib/SFG_CORRELATION.f90
	lib/SFG_STRUCTURE.f90
	lib/SFG_UTILS.f90
	lib/switches.f90
	lib/TRR_READER.f90
	lib/UTILS_ERROR.f90
	lib/XYZ_READER.f90
)

target_include_directories(sfg_lib PUBLIC ${PROJECT_SOURCE_DIR}/lib)

#the executables to be compiled
add_executable( Interface.x
	Interface/INTERFACE_UI.f90
	Interface/Interface.f90
)

add_executable( Binder.x
	Binder/BINDER_UI.f90
	Binder/Binder.f90
)

add_executable( SFG_CORR.x
	SFG_CORR/CORR_UI.f90
	SFG_CORR/SFG_CORR.f90
)

#link the library
target_link_libraries(Interface.x PRIVATE sfg_lib)
target_link_libraries(Binder.x    PRIVATE sfg_lib)
target_link_libraries(SFG_CORR.x  PRIVATE sfg_lib)

#this is a directive for install option
install (TARGETS 
			Interface.x
			Binder.x
			SFG_CORR.x
		 DESTINATION
			"bin"
)
