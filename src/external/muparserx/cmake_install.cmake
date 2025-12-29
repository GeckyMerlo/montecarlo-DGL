# Install script for directory: /home/domenico/amsc/montecarlo_1/src/external/muparserx

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

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/domenico/amsc/montecarlo_1/src/external/muparserx/libmuparserx.a")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "/home/domenico/amsc/montecarlo_1/src/external/muparserx/muparserx.pc")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/cmake/muparserx" TYPE FILE FILES
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/cmake/muparserxConfig.cmake"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/muparserxConfigVersion.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/muparserx" TYPE FILE FILES
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpDefines.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpError.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpFuncCmplx.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpFuncCommon.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpFuncMatrix.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpFuncNonCmplx.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpFuncStr.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpFwdDecl.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpICallback.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpIOprt.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpIOprtBinShortcut.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpIPackage.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpIPrecedence.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpIToken.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpIValReader.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpIValue.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpIfThenElse.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpMatrix.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpMatrixError.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpOprtBinAssign.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpOprtBinCommon.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpOprtBinShortcut.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpOprtCmplx.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpOprtIndex.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpOprtMatrix.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpOprtNonCmplx.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpOprtPostfixCommon.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpPackageCmplx.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpPackageCommon.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpPackageMatrix.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpPackageNonCmplx.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpPackageStr.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpPackageUnit.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpParser.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpParserBase.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpParserMessageProvider.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpRPN.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpScriptTokens.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpStack.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpStringConversionHelper.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpTest.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpTokenReader.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpTypes.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpValReader.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpValue.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpValueCache.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/mpVariable.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/suSortPred.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/suStringTokens.h"
    "/home/domenico/amsc/montecarlo_1/src/external/muparserx/parser/utGeneric.h"
    )
endif()

