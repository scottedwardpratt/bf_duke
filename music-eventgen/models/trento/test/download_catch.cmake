message(STATUS "Downloading Catch test header")

file(DOWNLOAD
  "https://raw.github.com/catchorg/Catch2/master/single_include/catch.hpp"
  "/Users/scottpratt/git/bf_duke/music-eventgen/models/trento/test/catch.hpp"
  TIMEOUT 20
  STATUS status
  TLS_VERIFY ON)

list(GET status 0 status_code)
list(GET status 1 status_string)

if(NOT status_code EQUAL 0)
  message(FATAL_ERROR
    "download failed: code ${status_code} ${status_string}")
endif()