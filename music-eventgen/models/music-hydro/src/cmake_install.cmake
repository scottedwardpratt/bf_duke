# Install script for directory: /Users/scottpratt/git/bf_duke/music-eventgen/models/music-hydro/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/Users/scottpratt/git/bf_duke/music-eventgen/models/music-hydro")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/scottpratt/git/bf_duke/music-eventgen/models/music-hydro/libmusic.dylib")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/Users/scottpratt/git/bf_duke/music-eventgen/models/music-hydro" TYPE SHARED_LIBRARY FILES "/Users/scottpratt/git/bf_duke/music-eventgen/models/music-hydro/src/libmusic.dylib")
  if(EXISTS "$ENV{DESTDIR}/Users/scottpratt/git/bf_duke/music-eventgen/models/music-hydro/libmusic.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/Users/scottpratt/git/bf_duke/music-eventgen/models/music-hydro/libmusic.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/scottpratt/git/bf_duke/music-eventgen/models/music-hydro"
      "$ENV{DESTDIR}/Users/scottpratt/git/bf_duke/music-eventgen/models/music-hydro/libmusic.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -x "$ENV{DESTDIR}/Users/scottpratt/git/bf_duke/music-eventgen/models/music-hydro/libmusic.dylib")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/scottpratt/git/bf_duke/music-eventgen/models/music-hydro/src/CMakeFiles/music.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/scottpratt/git/bf_duke/music-eventgen/models/music-hydro/MUSIChydro")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/Users/scottpratt/git/bf_duke/music-eventgen/models/music-hydro" TYPE EXECUTABLE FILES "/Users/scottpratt/git/bf_duke/music-eventgen/models/music-hydro/src/MUSIChydro")
  if(EXISTS "$ENV{DESTDIR}/Users/scottpratt/git/bf_duke/music-eventgen/models/music-hydro/MUSIChydro" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/Users/scottpratt/git/bf_duke/music-eventgen/models/music-hydro/MUSIChydro")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/scottpratt/git/bf_duke/music-eventgen/models/music-hydro/src"
      -add_rpath "/Users/scottpratt/git/bf_duke/music-eventgen/models/music-hydro"
      "$ENV{DESTDIR}/Users/scottpratt/git/bf_duke/music-eventgen/models/music-hydro/MUSIChydro")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}/Users/scottpratt/git/bf_duke/music-eventgen/models/music-hydro/MUSIChydro")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/scottpratt/git/bf_duke/music-eventgen/models/music-hydro/src/CMakeFiles/MUSIChydro.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

