# four-top-with-w-boson-production
# This repository contains all the analysis codes of the four top in association with w boson.
# First, the data is analyzed at the MadGraph level.
# The detector effects are added using the Delphes modular framework. 
# Multivariate analysis is performed using BDT.

**Signal Generation at MadGraph Level:**

MG5_aMC> generate p p > t t~ t t~ w- [NLO]
MG5_aMC> output fourtop_wminus
MG5_aMC> launch
shower = Pythia8
madspin = ON
1) run_card.dat changes:
   lhapdf = pdlabel
   260000 = lhaid
   True = fixed_ren_scale
   True = fixed_fac_scale
   344.00  = scale! fixed ren scale
   344.00  = dsqrt_q2fact1  ! fixed fact scale for pdf1
   344.00  = dsqrt_q2fact2  ! fixed fact scale for pdf2
   False  = use_syst

3) madspin_card.dat changes:
   decay t > w+ b, w+ > all all
   decay t~ > w- b~, w- > all all
   decay w+ > all all
   decay w- > all all
   decay z > all all
   #select your final state instead of all.

**Detector Simulation**
# Detector simulation is performed using the Delphes Modular framework.
# Go to the Delphes directory
# To reconstruct different high-level objects (electrons, muons, jets, b-tag jets, missing energy) in a given experiment of the collider use:

./DelphesHepMC cards/delphes_card_CMS.tcl path/to/directory/output_file.root path/to/directory/pythia8_events.hepmc
# Write a code (c++/python) with the help of examples in the example directory to analyze the output of the Delphes.

