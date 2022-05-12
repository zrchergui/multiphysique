#!/bin/sh

### Create a folder for the dune and dumux modules
### Go into the folder and execute this script

#if [ -d dune-common ]; then
#  echo "error: A directory named dune-common already exists."
#  echo "Aborting."
#  exit 1
#fi

### Clone the necessary modules
git clone https://gitlab.dune-project.org/core/dune-common.git
git clone https://gitlab.dune-project.org/core/dune-geometry.git
git clone https://gitlab.dune-project.org/core/dune-grid.git
git clone https://gitlab.dune-project.org/core/dune-istl.git
git clone https://gitlab.dune-project.org/core/dune-localfunctions.git
git clone https://gitlab.dune-project.org/staging/dune-uggrid.git
git clone https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git dumux
git clone https://github.com/francesco-patacchini/vertical-equilibrium.git Becker2018a
#git clone https://git.iws.uni-stuttgart.de/dumux-pub/Becker2018a.git

### Go to specific branches
cd dune-common && git checkout releases/2.5 && cd ..
cd dune-geometry && git checkout releases/2.5 && cd ..
cd dune-grid && git checkout releases/2.5 && cd ..
cd dune-istl && git checkout releases/2.5 && cd ..
cd dune-localfunctions && git checkout releases/2.5 && cd ..
cd dune-uggrid && git checkout releases/2.5 && cd ..
cd dumux && git checkout releases/2.11 && cd ..

### Run dunecontrol
#./dune-common/bin/dunecontrol --opts=./dumux/optim.opts all
