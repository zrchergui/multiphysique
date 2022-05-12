Summary
=======
This is the DuMuX module containing the code for producing the results
published in:

B. Becker, B. Guo, K. Bandilla, M. A. Celia, B. Flemisch, R. Helmig (2018)<br>
An adaptive multiphysics model coupling vertical equilibrium and full multidimensions for multiphase flow in porous media.<br>
Water Resources Research, 54, 4347–4360. https://doi.org/10.1029/2017WR022303

Installation
============

You can build the module just like any other DUNE
module.
The easiest way to install this module is to create a new folder and to execute
the file [installBecker2018a.sh](https://git.iws.uni-stuttgart.de/dumux-pub/Becker2018a/raw/master/installBecker2018a.sh)
in this folder.
You can copy the following to a terminal:
```bash
mkdir -p Becker2018a && cd Becker2018a
wget -q https://git.iws.uni-stuttgart.de/dumux-pub/Becker2018a/raw/master/installBecker2018a.sh
chmod +x installBecker2018a.sh && ./installBecker2018a.sh
```
For more detailed informations on installation, have a look at the
[DuMuX installation guide](https://dumux.org/installation/)
or use the [DuMuX handbook](https://dumux.org/docs/).

Applications
============

The applications can be found in the folder `appl` and `test`, particularly in the
following subfolders. Each program is related to a corresponding section in the
publication.

* `energyStorage`: program using a full multidimensional model for calculating reference solutions.
* `decoupled/2pve`: program using a single vertical equilibrium model.
* `modelcoupling/2pve2pfullmono`: program using the multiphysics model coupling vertical equilibrium and full multidimensions.


Used Versions and Software
==========================

For an overview on the used versions of the DUNE and DuMuX modules as well as
grid managers, have a look at [installBecker2018a.sh](https://git.iws.uni-stuttgart.de/dumux-pub/Becker2018a/raw/master/installBecker2018a.sh).

In addition, the linear solver SuperLU has to be installed.

The module has been tested successfully with GCC 5.0.
The autotools-based DUNE buildsystem has been employed.