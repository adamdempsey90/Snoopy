#	This file is part of the Snoopy code.
#
#    Snoopy code is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    Snoopy code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with the Snoopy code.  If not, see <http://www.gnu.org/licenses/>.



AC_INIT([snoopy], [2.0])
AC_PROG_CC

# Save debug flags
AC_SUBST(CFLAGS_DEBUG)
AC_SUBST(LDFLAGS_DEBUG)
CFLAGS_DEBUG=$CFLAGS
LDFLAGS_DEBUG=$LDFLAGS

#---------------------------------------------------
# Check standard libraries

AC_CHECK_HEADERS([stdlib.h])
AC_CHECK_HEADERS([math.h])
AC_CHECK_HEADERS([complex.h],[ok=yes],[
			echo "******************************************"
	              	echo "Error! Complex type support not found."
			echo "******************************************"
                  	exit -1
                  ])
#----------------------------------------------------
#Check best optimisation
AX_CC_MAXOPT 

#----------------------------------------------------
#if FFTW is not in the standard PATH...

# FFT Path: Path where to find the fft libraries (e.g. ~/usr)

if test "$FFTPATH" != ""; then
	CPPFLAGS="$CPPFLAGS -I$FFTPATH/include"
	CFLAGS="$CFLAGS -I$FFTPATH/include"
	CFLAGS_DEBUG="$CFLAGS_DEBUG -I$FFTPATH/include"
	LDFLAGS="$LDFLAGS -L$FFTPATH/lib"
	LDFLAGS_DEBUG="$LDFLAGS_DEBUG -L$FFTPATH/lib"
	
fi

#----------------------------------------------------
# Find Math library
AC_CHECK_LIB(m,cexp,[],[
			echo "******************************************"
	              	echo "Error! Math library not found."
			echo "******************************************"
					exit -1
				  ])

#----------------------------------------------------
#Check we have FFTW

AC_CHECK_HEADERS([fftw3.h],[ok=yes],[
				  echo "******************************************"
	              echo "Error! FFTW library not found."
				  echo "******************************************"
                  exit -1
                  ])


AC_CHECK_LIB(fftw3,fftw_execute, [], )


#----------------------------------------------------
# MPI Feature (--enable-mpi)

AC_SUBST(MPI_ACTIVE)

MPI_ACTIVE=notconfigured

AC_ARG_ENABLE(mpi, [  --enable-mpi     enable mpi parallelisation], ok=yes, ok=no)
if test "$ok" = "yes"; then
	ACX_MPI(MPI_ACTIVE=yes,[
				  echo "******************************************"
	              echo "Error! MPI libraries not found."
				  echo "******************************************"
                  exit -1
                  ])
fi

#----------------------------------------------------
#FFTW_MPI Feature (--enable-fftw-mpi)

AC_SUBST(FFTW_MPI_ACTIVE)

FFTW_MPI_ACTIVE=notconfigured

AC_ARG_ENABLE(fftw-mpi, [  --enable-fftw-mpi     enable fftw3 experimental mpi support], ok=yes, ok=no)
if test "$ok" = "yes"; then
	CCSAVE=$CC	## Save the current compiler, since we need mpicc for the following test
	CC=$MPICC
	AC_SEARCH_LIBS(fftw_mpi_init, fftw3_mpi fftw3 fftw3_threads,[ ],[
				  echo "******************************************"
	              echo "Error! FFTW MPI libraries not found."
				  echo "******************************************"
                  exit -1
                  ])
	FFTW_MPI_ACTIVE=yes
	CC=$CCSAVE
fi

#----------------------------------------------------
# OpenMP feature (--enable-openmp)

AC_SUBST(OPENMP_CFLAGS)
AC_SUBST(OPENMP_ACTIVE)

OPENMP_ACTIVE=notconfigured

AC_ARG_ENABLE(openmp, [  --enable-openmp     enable OpenMP support], ok=yes, ok=no)
if test "$ok" = "yes"; then
	#check if openMP is supported by CC
	AX_OPENMP(ok=yes,[
				  echo "******************************************"
	              echo "Error! OpenMP not supported by compiler."
				  echo "******************************************"
                  exit -1
                  ])

	#check if Threads or OpenMP is supported by fftw
	#save CFLAGS and add openmp flag for next test

	CFLAGS_SAVE=$CFLAGS
	CFLAGS="$CFLAGS $OPENMP_CFLAGS"
	
	AC_SEARCH_LIBS(fftw_init_threads, fftw3_threads, OPENMP_ACTIVE=yes,[
				  echo "************************************************"
	              echo "Error! Threads or OpenMP not supported by FFTW. "
				  echo "************************************************"
                  exit -1
                  ])

	#Restore CFLAGS
	CFLAGS=$CFLAGS_SAVE
				  
fi

#----------------------------------------------------
# Debug options (--enable-debug)

AC_SUBST(DEBUG_ACTIVE)

DEBUG_ACTIVE=no

AC_ARG_ENABLE(debug, [  --enable-debug     enable Debug options], ok=yes, ok=no)
if test "$ok" = "yes"; then
	DEBUG_ACTIVE=yes
fi

#----------------------------------------------------
# if a problem is specified, create a gvars.h and snoopy.cfg
# corresponding to this problem

AC_ARG_WITH(problem, [  --with-problem=PBM	Problem generator from cfg/*], problem=$withval, problem=notset)
if test ! "$problem" = "notset"; then
	# initialize the problem $problem
	if test -e src/problem/$problem/gvars.h; then
		cp src/problem/$problem/gvars.h src/
		cp src/problem/$problem/snoopy.cfg .
	else
		echo "********************************************************"
		echo "Invalid problem directory, valid names listed below:"
		ls src/problem/
		echo "********************************************************"
		exit -1
	fi
fi

#----------------------------------------------------
#Create a code configuration if none exist
if test ! -e src/gvars.h; then
	echo "Using default problem"
	problem=Default
	cp -f src/problem/default/gvars.h src/gvars.h
fi
if test ! -e snoopy.cfg; then
	problem=Default
	cp -f src/problem/default/snoopy.cfg .
fi

#----------------------------------------------------
# Warnings in the relevant files

AC_SUBST(WARNING1)
AC_SUBST(WARNING2)
WARNING1="WARNING! This file has been automatically generated by configure."
WARNING2="Any changes to it will be overwritten the next time configure is run."



AC_CONFIG_FILES([Makefile])
AC_OUTPUT
echo "*******************************************************"
echo "I will use the following configuration:"
echo " Problem:     $problem"
echo " Compiler:    $CC"
echo " MPI:         $MPI_ACTIVE"
echo " FFTW MPI:    $FFTW_MPI_ACTIVE"
echo " OpenMP:      $OPENMP_ACTIVE"
echo " Debug:       $DEBUG_ACTIVE"
echo " CFlags:      $CFLAGS"
echo " LDFlags:     $LDFLAGS"
echo ""
echo "Type make or make check to compile"
echo "*******************************************************"
