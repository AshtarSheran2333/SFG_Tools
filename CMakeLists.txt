#made with help of https://www.atomwitch.net/2022/04/02/cmake-fortran.html

cmake_minimum_required(VERSION 3.1)

project(MD_SFG VERSION 1.0.0
		DESCRIPTION "A program suite to calculate SFG spectra from MD trajectories"
		LANGUAGES Fortran)

enable_language(Fortran)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f90)

#compiler setup
#this does not work for some reason, and the compiler needs to be specifiedn in the compile script
set(CMAKE_Fortran_COMPILER ifx)
# O2 optimalization (should be default anyway)
# -qopenmp (to enable openmp)
# -fpp (to enable preprocessor)
set(FVERSION "-O2 -qopenmp -fpp")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FVERSION}")

#pack all includes as a library
add_library(sfg_lib STATIC
	lib/bdata.f90
	lib/SFG_TYPES.f90
	lib/SFG_FRAMES.f90
	lib/switches.f90
)

target_include_directories(sfg_lib
	PUBLIC
		${PROJECT_SOURCE_DIR}/lib
)

#the executables to be compiled
add_executable(Interface.x
	Interface_P/Interface_P.f90
)

add_executable(Binder.x
	Binder_P/Binder_P.f90
)

add_executable(SSP_CORR.x
	SSP_CORR_P/SSP_CORR_P.f90
)

add_executable(Corr_to_spectrum.x
	Corr_to_spectrum/Corr_to_spectrum.f90
)

#link the library
target_link_libraries(Interface.x
	PRIVATE
		sfg_lib
)

#link the library
target_link_libraries(Binder.x
	PRIVATE
		sfg_lib
)

#link the library
target_link_libraries(SSP_CORR.x
	PRIVATE
		sfg_lib
)

#link the library
target_link_libraries(Corr_to_spectrum.x
	PRIVATE
		sfg_lib
)

#this is a directive for install option
install (TARGETS 
			Interface.x
			Binder.x
			SSP_CORR.x
			Corr_to_spectrum.x			
		 DESTINATION
			"bin"
)


