###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s):<christophe.prudhomme@feelpp.org>
#       Date: 06 Feb 2015
#
#  Copyright (C) 2015-2017 Feel++ Consortium
#
# Distributed under the GPL(GNU Public License):
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#
cmake_minimum_required(VERSION 3.3)

project(projet_eads)
set(PROJECT_SHORTNAME "pw")

if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
  #find_package(Feel++ PATHS $ENV{FEELPP_DIR}/share/feelpp/feel/cmake/modules /usr/share/feelpp/feel/cmake/modules /usr/local/share/feelpp/feel/cmake/modules )
  find_package(Feel++ PATHS $ENV{FEELPP_DIR}/share/feel/cmake/modules /usr/share/feel/cmake/modules /usr/local/share/feel/cmake/modules )
  if(NOT FEELPP_FOUND)
      message(FATAL_ERROR "Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  endif()
endif()

enable_testing()

if ( NOT TARGET check )
  add_custom_target(check)
endif()
 
if ( NOT TARGET projet )
  add_custom_target(projet)
endif()


#
# projet M1 CSMI 
#
feelpp_add_application( eads_2d SRCS projet_eads.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 EXEC PROJET_EADS_2D INSTALL )
#add_dependencies(projet ${PW_EADS_2D}  ${PW_EADS_2D} )

add_custom_target(install-projet
  DEPENDS
  projet
  COMMAND
  "${CMAKE_COMMAND}"
  -DCMAKE_INSTALL_COMPONENT=testcases  -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Quickstart -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Bin -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_install.cmake"
  )

if ( NOT TARGET install-feelpp-base )
  add_custom_target( install-feelpp-base DEPENDS install-projet)
  else()
  add_dependencies(install-feelpp-base install-projet)
endif()

