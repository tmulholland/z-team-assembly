# z-team-assembly
Text files and common code for RA2b Zinv analysers

# Instructions for Bill's integration script

To run the integration code:

cd src/

root -l RA2bin_driver.C

# Instructions for generating new dat files (this takes time)

cd python/

./getExtrapDatFiles.py sig >&! ../datFiles/DY_signal.dat

./getExtrapDatFiles.py hdp >&! ../datFiles/DY_hdp.dat

./getExtrapDatFiles.py ldp >&! ../datFiles/DY_ldp.dat

./getDoubleRatioDatFiles.py sig >&! ../datFiles/DR_signal.dat

./getDoubleRatioDatFiles.py hdp >&! ../datFiles/DR_hdp.dat

./getDoubleRatioDatFiles.py ldp >&! ../datFiles/DR_ldp.dat

# Instructions for writing the efficiencies to the effHists.root file

cd python/

./setPhotonEffs.py

./setZllEffs.py


# code to make histograms

To go here

# syncronization stuff

To go here
