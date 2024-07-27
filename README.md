# four-top-with-w-boson-production
This repository contains all the analysis codes of the four top in association with w boson.
First, the data is analyzed at the MadGraph level.
The detector effects are added using the Delphes modular framework. 
Multivariate analysis is performed using BDT.

## Installing MadGraph5
```
wget http://madgraph.physics.illinois.edu/Downloads/MG5_aMC_v2.x.y.tar.gz
tar xf MG5_aMC_v2.x.y.tar.gz
rm MG5_aMC_v2.x.y.tar.gz
cd MG5_aMC
./bin/mg5_aMC
install lhapdf5
install pythia8
install mg5amc_py8_interface
install Delphes #it is recommended to install Delphes externally
```
Installing Delphes:
```
wget http://cp3.irmp.ucl.ac.be/downloads/Delphes-3.2.0.tar.gz
tar -xf Delphes-3.2.0.tar.gz
cd Delphes-3.2.0/
```
for lxplus users, run the following commands before make
```
source /afs/cern.ch/sw/lcg/external/gcc/4.9.3/x86_64-slc6/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.00/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh
````
```
make
```
## Signal Generation at MadGraph Level
```
MG5_aMC> generate p p > t t~ t t~ w- [NLO]
MG5_aMC> output fourtop_wminus
MG5_aMC> launch
shower = Pythia8
madspin = ON
```

1) run_card.dat changes:
```
   lhapdf = pdlabel
   260000 = lhaid
   True = fixed_ren_scale
   True = fixed_fac_scale
   344.00  = scale! fixed ren scale
   344.00  = dsqrt_q2fact1  ! fixed fact scale for pdf1
   344.00  = dsqrt_q2fact2  ! fixed fact scale for pdf2
   False  = use_syst
```

2) madspin_card.dat changes:
```
   decay t > w+ b, w+ > all all
   decay t~ > w- b~, w- > all all
   decay w+ > all all
   decay w- > all all
   decay z > all all
   #select your final state instead of "all".
```
## Detector Simulation
Detector simulation is performed using the Delphes Modular framework.
Go to the Delphes directory
To reconstruct different high-level objects (electrons, muons, jets, b-tag jets, missing energy) in a given experiment of the collider use:

```
./DelphesHepMC cards/delphes_card_CMS.tcl path/to/directory/output_file.root path/to/directory/input_file.hepmc
```
###### Write a code (c++/python) with the help of examples in the example directory to analyze the output of the Delphes.

