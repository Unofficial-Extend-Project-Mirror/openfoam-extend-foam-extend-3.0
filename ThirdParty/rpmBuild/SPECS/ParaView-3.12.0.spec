#------------------------------------------------------------------------------
# =========                 |
# \\      /  F ield         | foam-extend: Open Source CFD
#  \\    /   O peration     |
#   \\  /    A nd           | For copyright notice see file Copyright
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
# License
#     This file is part of foam-extend.
#
#     foam-extend is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     foam-extend is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.
#
# Script
#     RPM spec file for ParaView-3.12.0
#
# Description
#     RPM spec file for creating a relocatable RPM
#
# Author:
#     Martin Beaudoin, Hydro-Quebec, (2012)
#
#------------------------------------------------------------------------------

# We grab the value of WM_THIRD_PARTY and WM_OPTIONS from the environment variable
%{expand:%%define _WM_THIRD_PARTY_DIR %(echo $WM_THIRD_PARTY_DIR)}
%{expand:%%define _WM_OPTIONS         %(echo $WM_OPTIONS)}

# Disable the generation of debuginfo packages
%define debug_package %{nil}

# The topdir needs to point to the $WM_THIRD_PARTY/rpmbuild directory
%define _topdir	 	%{_WM_THIRD_PARTY_DIR}/rpmBuild
%define _tmppath	%{_topdir}/tmp

# Will install the package directly $WM_THIRD_PARTY_DIR
#   Some comments about package relocation:
#   By using this prefix for the Prefix:  parameter in this file, you will make this 
#   package relocatable. 
#
#   This is fine, as long as your software is itself relocatable.
#
#   Simply take note that libraries built with libtool are not relocatable because the
#   prefix we specify will be hard-coded in the library .la files.
#   Ref: http://sourceware.org/autobook/autobook/autobook_80.html
#
#   In that case, if you ever change the value of the $WM_THIRD_PARTY_DIR, you will
#   not be able to reutilize this RPM, even though it is relocatable. You will need to 
#   regenerate the RPM.
#
%define _prefix         %{_WM_THIRD_PARTY_DIR}

%define name		ParaView
%define release		%{_WM_OPTIONS}
%define version 	3.12.0

%define buildroot       %{_topdir}/BUILD/%{name}-%{version}-root

BuildRoot:	        %{buildroot}
Summary: 		ParaView
License: 		Unkown
Name: 			%{name}
Version: 		%{version}
Release: 		%{release}
URL:            http://www.paraview.org/files/v3.12/
Source: 		%url/%{name}-%{version}.tar.gz
Prefix: 		%{_prefix}
Group: 			Development/Tools
Patch0:         ParaView-3.12.0.patch_darwin
Patch1:         paraview-gcc47.patch

%define _installPrefix  %{_prefix}/packages/%{name}-%{version}/platforms/%{_WM_OPTIONS}

#--------------------------------------------------------------------------
#
# Here, we define default compiling options for cmake
#
# One can override the option on the commande line : --define='MACRO EXPR'
#
%{!?_withVerbose:     %define _withVerbose     false}
%{!?_withMesa:        %define _withMesa        false}
%{!?_withMPI:         %define _withMPI         true}
%{!?_withPython:      %define _withPython      true}
%{!?_withQt:          %define _withQt          true}
%{!?_qmakePath:       %define _qmakePath       Undefined}
%{!?_mesaIncludePath: %define _mesaIncludePath Undefined}
%{!?_mesaLibPath:     %define _mesaLibPath     Undefined}
%{!?_pythonLibPath:   %define _pythonLibPath   Undefined}

#--------------------------------------------------------------------------

%description
%{summary}

%prep
%setup -q 

%ifos darwin
%patch0 -p1
%endif

%patch1 -p1

%build
#
# set CMake cache variables
#
    addCMakeVariable()
    {
        while [ -n "$1" ]
        do
            CMAKE_VARIABLES="$CMAKE_VARIABLES -D$1"
            shift
        done
   }

    # export WM settings in a form that GNU configure recognizes
    [ -n "$WM_CC" ]         &&  export CC="$WM_CC"
    [ -n "$WM_CXX" ]        &&  export CXX="$WM_CXX"
    [ -n "$WM_CFLAGS" ]     &&  export CFLAGS="$WM_CFLAGS"
    [ -n "$WM_CXXFLAGS" ]   &&  export CXXFLAGS="$WM_CXXFLAGS"
    [ -n "$WM_LDFLAGS" ]    &&  export LDFLAGS="$WM_LDFLAGS"

    set +x
    echo ""
    echo "Compilation options:"
    echo "     _withVerbose     : %{_withVerbose}"
    echo "     _withMesa        : %{_withMesa}"
    echo "     _withMPI         : %{_withMPI}"
    echo "     _withPython      : %{_withPython}"
    echo "     _withQt          : %{_withQt}"
    echo "     _qmakePath       : %{_qmakePath}"
    echo "     _mesaIncludePath : %{_mesaIncludePath}"
    echo "     _mesaLibPath     : %{_mesaLibPath}"
    echo "     _pythonLibPath   : %{_pythonLibPath}"
    echo ""
    set -x

    # start with these general settings
    addCMakeVariable  VTK_USE_TK:BOOL=OFF
    addCMakeVariable  BUILD_SHARED_LIBS:BOOL=ON  VTK_USE_RPATH:BOOL=OFF
    addCMakeVariable  CMAKE_BUILD_TYPE:STRING=Release

    # include development files in "make install"
    addCMakeVariable  PARAVIEW_INSTALL_DEVELOPMENT:BOOL=ON

    # new alternative to "make HTMLDocumentation"
    addCMakeVariable  PARAVIEW_GENERATE_PROXY_DOCUMENTATION:BOOL=ON

 %ifos darwin
    # Additional installation rules for Mac OS X 
    addCMakeVariable  PARAVIEW_EXTRA_INSTALL_RULES_FILE:FILEPATH=%{_topdir}/BUILD/%{name}-%{version}/Applications/ParaView-3.12.0_extra_install_Darwin.cmake
%endif

    # Add the value of _qmakePath for QT_QMAKE_EXECUTABLE
    addCMakeVariable  QT_QMAKE_EXECUTABLE:FILEPATH=%{_qmakePath}

    echo "CMAKE_VARIABLES: $CMAKE_VARIABLES"

%ifos darwin
    # For Mac OSX:
    # The configuration of Paraview will be using the environment variable MACOSX_DEPLOYMENT_TARGET.
    # This variable was initialized using 'sw_vers -productVersion' in etc/bashrc.
    # We need to get rid of the revision number from this string. eg turn "10.7.5" into "10.7"
    v=( ${MACOSX_DEPLOYMENT_TARGET//./ } )
    export MACOSX_DEPLOYMENT_TARGET="${v[0]}.${v[1]}"
    echo "Resetting MACOSX_DEPLOYMENT_TARGET to ${MACOSX_DEPLOYMENT_TARGET}"
%endif

    mkdir -p ./buildObj
    cd ./buildObj

    cmake \
        -DCMAKE_INSTALL_PREFIX:PATH=%{_installPrefix} \
        $CMAKE_VARIABLES \
	..

    [ -z "$WM_NCOMPPROCS" ] && WM_NCOMPPROCS=1
    make -j $WM_NCOMPPROCS

%install
    # On OpenSUSE, rpmbuild, will choke when detecting unreferenced symlinks
    # created during installation.
    # Qt version 4.6.3 will generate some unreferenced symlinks when
    # ParaView is compiled and installed. By enabling the following
    # environment variable, the command brp-symlink will still complain
    # about missing link targets, but it won't stop rpmbuild from generating
    # the final rpm.
    # For all other Unix distros, this is a no-op.
    export NO_BRP_STALE_LINK_ERROR=yes

    cd buildObj
    make install DESTDIR=$RPM_BUILD_ROOT

    # Creation of foam-extend specific .csh and .sh files"

    echo ""
    echo "Generating foam-extend specific .csh and .sh files for the package %{name}-%{version}"
    echo ""
    #
    # Generate package specific .sh file for foam-extend
    #
mkdir -p $RPM_BUILD_ROOT/%{_installPrefix}/etc
cat << DOT_SH_EOF > $RPM_BUILD_ROOT/%{_installPrefix}/etc/%{name}-%{version}.sh
# Load %{name}-%{version} libraries and binaries if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

export PARAVIEW_DIR=\$WM_THIRD_PARTY_DIR/packages/%{name}-%{version}/platforms/\$WM_OPTIONS
export PARAVIEW_BIN_DIR=\$PARAVIEW_DIR/bin
export PARAVIEW_LIB_DIR=\$PARAVIEW_DIR/lib
export PARAVIEW_INCLUDE_DIR=\$PARAVIEW_DIR/include

export PARAVIEW_VERSION=%{version}

# NB: It is important to set the PV_PLUGIN_PATH location to a directory containing only the ParaView plugins.
#     Otherwise, paraview will try to automatically autoload each and every dynamic library it can find in the
#     specified directory to see if a given library is a paraview plugin. 
#     In the case of \$FOAM_LIBBIN, with over 80 libraries, this is a total waste of time that will slow down the
#     startup of paraview or even make paraview crash on startup.
export PV_PLUGIN_PATH=\$FOAM_LIBBIN/paraview_plugins

[ -d \$PARAVIEW_LIB_DIR/paraview-3.12 ] && _foamAddLib \$PARAVIEW_LIB_DIR/paraview-3.12

# Enable access to the package applications if present
[ -d \$PARAVIEW_BIN_DIR ] && _foamAddPath \$PARAVIEW_BIN_DIR

# Additional binary path if running on Mac OS X
[ -d \$PARAVIEW_BIN_DIR/paraview.app/Contents/MacOS ] && _foamAddPath \$PARAVIEW_BIN_DIR/paraview.app/Contents/MacOS

DOT_SH_EOF

    #
    # Generate package specific .csh file for foam-extend
    #
cat << DOT_CSH_EOF > $RPM_BUILD_ROOT/%{_installPrefix}/etc/%{name}-%{version}.csh
# Load %{name}-%{version} libraries and binaries if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setenv PARAVIEW_DIR \$WM_THIRD_PARTY_DIR/packages/%{name}-%{version}/platforms/\$WM_OPTIONS
setenv PARAVIEW_BIN_DIR \$PARAVIEW_DIR/bin
setenv PARAVIEW_LIB_DIR \$PARAVIEW_DIR/lib
setenv PARAVIEW_INCLUDE_DIR \$PARAVIEW_DIR/include

setenv PARAVIEW_VERSION %{version}

# NB: It is important to set the PV_PLUGIN_PATH location to a directory containing only the ParaView plugins.
#     Otherwise, paraview will try to automatically autoload each and every dynamic library it can find in the
#     specified directory to see if a given library is a paraview plugin. 
#     In the case of \$FOAM_LIBBIN, with over 80 libraries, this is a total waste of time that will slow down the
#     startup of paraview or even make paraview crash on startup.
setenv PV_PLUGIN_PATH \$FOAM_LIBBIN/paraview_plugins

if ( -e \$PARAVIEW_BIN_DIR ) then
    _foamAddPath \$PARAVIEW_BIN_DIR
endif

if ( -e \$PARAVIEW_LIB_DIR/paraview-3.12 ) then
    _foamAddLib \$PARAVIEW_LIB_DIR/paraview-3.12
endif


# Additional binary path if running on Mac OS X
if ( -e \$PARAVIEW_BIN_DIR/paraview.app/Contents/MacOS ) then
    _foamAddPath \$PARAVIEW_BIN_DIR/paraview.app/Contents/MacOS
endif
DOT_CSH_EOF

    #finally, generate a .tgz file for systems where using rpm for installing packages
    # as a non-root user might be a problem.
    (mkdir -p  %{_topdir}/TGZS/%{_target_cpu}; cd $RPM_BUILD_ROOT/%{_prefix}; tar -zcvf %{_topdir}/TGZS/%{_target_cpu}/%{name}-%{version}.tgz  packages/%{name}-%{version})

    # Make sure we clean the installation directory. Somehow, make install generates files there. Need to revisit this eventually.
    (rm -rf %{_prefix}/packages/%{name}-%{version}    )

%clean
rm -rf $RPM_BUILD_ROOT    
rm -rf %{_prefix}/packages/%{name}-%{version}   

%files
%defattr(-,root,root)
%{_installPrefix}

