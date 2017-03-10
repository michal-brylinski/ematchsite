#===============================================================================
#          ______                  _        _    _           
#         |  ___ \       _        | |      | |  (_)_         
#     ____| | _ | | ____| |_  ____| | _     \ \  _| |_  ____ 
#    / _  ) || || |/ _  |  _)/ ___) || \     \ \| |  _)/ _  )
#   ( (/ /| || || ( ( | | |_( (___| | | |_____) ) | |_( (/ / 
#    \____)_||_||_|\_||_|\___)____)_| |_(______/|_|\___)____)
#
#
#   eMatchSite - sequence order independent binding site alignment
#
#   Computational Systems Biology Group
#   Department of Biological Sciences
#   Center for Computation & Technology
#   Louisiana State University
#   407 Choppin Hall, Baton Rouge, LA 70803, USA
#
#   http://www.brylinski.org
#
#   Report bugs to michal@brylinski.org
#
#   Copyright 2013 Michal Brylinski
#
#   This file is part of eMatchSite.
#
#   eMatchSite is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   eMatchSite is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with eMatchSite. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

CXX = g++

CC = gcc

FC = gfortran

SH = sh

EXE = ematchsite ematchsite_kcombu

SRC_EMATCHSITE = src-ematchsite
SRC_GZSTREAM   = src-gzstream
SRC_LIBSVM     = src-libsvm
SRC_MUNKRES    = src-munkres

CPPFLAGS = -O2 -Wall -Wno-write-strings -std=c++11 -fPIC -I. -I$(SRC_EMATCHSITE) -I$(SRC_GZSTREAM) -I$(SRC_LIBSVM) -I$(SRC_MUNKRES)

FFLAGS = -O2

LDFLAGS = -lz -lm -L.

OBJ_EMATCHSITE = $(SRC_EMATCHSITE)/coords.o \
                 $(SRC_EMATCHSITE)/data.o \
                 $(SRC_EMATCHSITE)/ematchsite.o \
                 $(SRC_EMATCHSITE)/fisherpitman.o \
                 $(SRC_EMATCHSITE)/protein.o \
                 $(SRC_EMATCHSITE)/rmsd.o \
                 $(SRC_EMATCHSITE)/runsvm.o \
                 $(SRC_EMATCHSITE)/screen.o \
                 $(SRC_EMATCHSITE)/sitepair.o \
                 $(SRC_EMATCHSITE)/tanimoto.o \
                 $(SRC_EMATCHSITE)/walltime.o \
                 $(SRC_GZSTREAM)/gzstream.o \
                 $(SRC_LIBSVM)/svm.o \
                 $(SRC_MUNKRES)/munkres.o

default: $(EXE)

all: $(EXE)

ematchsite: $(OBJ_EMATCHSITE)
	$(CXX) -o $@ $(OBJ_EMATCHSITE) $(LDFLAGS)
	@mkdir -p bin/
	@mv ematchsite bin/

ematchsite_kcombu:
	$(SH) $(SRC_EMATCHSITE)/ematchsite_kcombu.shar
	@chmod +x ematchsite_kcombu
	@mkdir -p bin/
	@mv ematchsite_kcombu bin/

#=== eMatchSite ================================================================

ematchsite.o: ematchsite.C
	$(CXX) $(CPPFLAGS) -c -o ematchsite.o ematchsite.C

coords.o: coords.C
	$(CXX) $(CPPFLAGS) -c -o coords.o coords.C

data.o: data.C
	$(CXX) $(CPPFLAGS) -c -o data.o data.C

fisherpitman.o: fisherpitman.C
	$(CXX) $(CPPFLAGS) -c -o fisherpitman.o fisherpitman.C

protein.o: protein.C
	$(CXX) $(CPPFLAGS) -c -o protein.o protein.C

rmsd.o: rmsd.f
	$(FC) $(FFLAGS) -c -o rmsd.o rmsd.f

runsvm.o: runsvm.C
	$(CXX) $(CPPFLAGS) -c -o runsvm.o runsvm.C

screen.o: screen.C
	$(CXX) $(CPPFLAGS) -c -o screen.o screen.C

sitepair.o: sitepair.C
	$(CXX) $(CPPFLAGS) -c -o sitepair.o sitepair.C

tanimoto.o: tanimoto.C
	$(CXX) $(CPPFLAGS) -c -o tanimoto.o tanimoto.C

walltime.o: walltime.C
	$(CXX) $(CPPFLAGS) -c -o walltime.o walltime.C

#=== gzstream ==================================================================

gzstream.o: gzstream.C
	$(CXX) $(CPPFLAGS) -c -o gzstream.o gzstream.C

#=== libsvm ====================================================================

svm.o: svm.C
	$(CXX) $(CPPFLAGS) -c -o svm.o svm.C

#=== munkres ===================================================================

munkres.o: munkres.cpp
	$(CXX) $(CPPFLAGS) -c -o munkres.o munkres.cpp

#=== clean =====================================================================

clean:
	@(rm -f ${EXE} bin/ematchsite bin/ematchsite_kcombu ${OBJ_EMATCHSITE})
