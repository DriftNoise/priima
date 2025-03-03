# Copyright 2021, Kai Wah Chan, Mitja Echim, Andreas Folkers

# This file is part of PRIIMA.
# PRIIMA is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
# PRIIMA is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with
# PRIIMA. If not, see https://www.gnu.org/licenses/gpl-3.0.html.

set(LIB_NAME ${PROJECT_NAME})

find_package(OpenCV REQUIRED COMPONENTS core imgproc calib3d)
find_package(OpenMP REQUIRED)

if (${OpenMP_FOUND})
	message(STATUS "Found OpenMP")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(
	${OpenCV_INCLUDE_DIRS}
	${CMAKE_CURRENT_SOURCE_DIR} # local opencv2 adaptations
	${PYTHON3_INCLUDE_DIR}
	${PROJECT_SOURCE_DIR}/pyboostcvconverter/include
)
message(${PROJECT_SOURCE_DIR})

add_library(${LIB_NAME} SHARED
	tps_trans.cpp
)

#target_include_directories(${LIB_NAME} PRIVATE )

target_link_libraries(${LIB_NAME} LINK_PUBLIC
	pbcvt
	${OpenCV_LIBS}
)

get_target_property(TEST_OUT_PATH ${LIB_NAME} LIBRARY_OUTPUT_NAME)

set_target_properties(${LIB_NAME} 
	PROPERTIES
		ARCHIVE_OUTPUT_NAME ${LIB_NAME}${ARCHIVE_OUTPUT_SUFFIX}  # prevent name conflict for python2/3 outputs
		PREFIX ""
		OUTPUT_NAME ${LIB_NAME}
)