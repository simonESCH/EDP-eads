###  MakeLists.txt; coding: utf-8 --- 

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

project(projet_CSMI)
set(PROJECT_SHORTNAME "csmi")

message("\n")
message(FEELPP_DIR=$ENV{FEELPP_DIR})
message("\n")

if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
    find_package(Feel++ PATHS $ENV{FEELPP_DIR}/share/feelpp/feel/cmake/modules /usr/share/feelpp/feel/cmake/modules /usr/local/share/feelpp/feel/cmake/modules )
    if(NOT FEELPP_FOUND)
        message(FATAL_ERROR "Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
    endif()
endif()

enable_testing()

if ( NOT TARGET check )
    add_custom_target(check)
endif()

if ( NOT TARGET projet_CSMI )
    add_custom_target(projet_CSMI)
endif()

#
# Laplacian
#
#feelpp_add_application( laplacian_2d SRCS qs_laplacian.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 EXEC QS_LAPLACIAN_2D INSTALL TESTS)
#feelpp_add_application( laplacian_3d SRCS qs_laplacian.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 EXEC QS_LAPLACIAN_3D INSTALL TESTS)
#feelpp_add_testcase( laplacian PREFIX qs CATEGORY projet_CSMI )
#add_dependencies(projet_CSMI ${QS_LAPLACIAN_2D} ${QS_LAPLACIAN_3D}  )
#add_dependencies(check ${QS_LAPLACIAN_2D} ${QS_LAPLACIAN_3D} )


#
# EADS
#
feelpp_add_application( eads_2d SRCS eads/projet_eads.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 EXEC CSMI_EADS_2D INSTALL TESTS)
#feelpp_add_application( eads_3d SRCS eads/projet_eads.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 EXEC CSMI_EADS_3D INSTALL TESTS)
feelpp_add_testcase( eads PREFIX csmi CATEGORY projet_CSMI )
add_dependencies(projet_CSMI ${CSMI_EADS_2D}  )
add_dependencies(check ${CSMI_EADS_2D} )

add_custom_target(install-projet_CSMI
DEPENDS
projet_CSMI
COMMAND
"${CMAKE_COMMAND}"
-DCMAKE_INSTALL_COMPONENT=testcases  -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_install.cmake"
-DCMAKE_INSTALL_COMPONENT=Projet_CSMI -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_install.cmake"
-DCMAKE_INSTALL_COMPONENT=Bin -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_install.cmake"
)

if ( NOT TARGET install-feelpp-base )
    add_custom_target( install-feelpp-base DEPENDS install-projet_CSMI)
else()
    add_dependencies(install-feelpp-base install-projet_CSMI)
endif()

