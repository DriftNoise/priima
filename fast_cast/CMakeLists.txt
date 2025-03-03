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

cmake_minimum_required(VERSION 3.5)

set(PROJECT_NAME FastCast)

project(${PROJECT_NAME})

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS FALSE) # we want to use -std=c++** (FALSE) instead of -std=gnu++** (TRUE)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Update if necessary
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-long-long -pedantic")
endif()

set(-DPYTHON_DESIRED_VERSION=3.X)
add_compile_definitions(NUMPY_IMPORT_ARRAY_RETVAL=NULL)

set(UNSUPPORTED_OS "Your OS does not seem to be supported - You may face problems compiling Fast Cast.")

if (UNIX AND NOT APPLE)
    file(READ "/etc/issue" OS_Name)
    string(REGEX MATCH "Ubuntu *([0-9]+[.][0-9]+)[.]" Version "${OS_Name}")

    if (${CMAKE_MATCH_COUNT} EQUAL 1)

	if (${CMAKE_MATCH_1} VERSION_EQUAL "18.04")

	    find_package(PythonInterp 3 REQUIRED)
	    find_package(PythonLibs 3 REQUIRED)

	    set(PYTHON3_INCLUDE_DIR ${PYTHON_INCLUDE_DIRS})
	    set(PYTHON3_LIBRARY ${PYTHON_LIBRARIES})

	    set(PYTHON3_VERSION_MINOR 6)
	else()

	    MESSAGE(WARNING ${UNSUPPORTED_OS})
	endif()

    else()

	MESSAGE(WARNING ${UNSUPPORTED_OS})
    endif()

else()

    MESSAGE(WARNING ${UNSUPPORTED_OS})
endif()

add_subdirectory(pyboostcvconverter)
add_subdirectory(src)
