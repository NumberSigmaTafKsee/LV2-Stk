# Install script for directory: /home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
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

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/quake/Projects/Current/SoundWave/c++/Aquila/build/lib/unittestpp/libUnitTest++.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/UnitTest++" TYPE FILE FILES
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/AssertException.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/CheckMacros.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/Checks.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/CompositeTestReporter.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/Config.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/CurrentTest.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/DeferredTestReporter.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/DeferredTestResult.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/ExceptionMacros.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/ExecuteTest.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/HelperMacros.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/MemoryOutStream.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/ReportAssert.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/ReportAssertImpl.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/Test.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/TestDetails.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/TestList.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/TestMacros.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/TestReporter.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/TestReporterStdout.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/TestResults.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/TestRunner.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/TestSuite.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/TimeConstraint.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/TimeHelpers.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/UnitTest++.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/UnitTestPP.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/XmlTestReporter.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/UnitTest++/Posix" TYPE FILE FILES
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/Posix/SignalTranslator.h"
    "/home/quake/Projects/Current/SoundWave/c++/Aquila/lib/unittestpp/UnitTest++/Posix/TimeHelpers.h"
    )
endif()

