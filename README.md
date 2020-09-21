# nanoAOD producer customized for R(J/\psi) analysis (forked from RK analisis)
See at the end for UL

## Getting started

```shell
cmsrel CMSSW_10_2_15
cd CMSSW_10_2_15/src
cmsenv
git cms-init
```

## Add the latest code and model (2019Aug07) for the electron ID 
( skip if you don't use electrons)

```shell
git cms-addpkg RecoEgamma/EgammaElectronProducers
git cms-merge-topic CMSBParking:from-CMSSW_10_2_15_2019Aug07
git cms-addpkg RecoEgamma/ElectronIdentification
scram b

# Check $CMSSW_BASE/external exists before this step (e.g. run 'scram b' to create it)
git clone --single-branch --branch 102X_LowPtElectrons_2019Aug07 git@github.com:CMSBParking/RecoEgamma-ElectronIdentification.git $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data

# The following step is required if running on CRAB
mv $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data/LowPtElectrons $CMSSW_BASE/src/RecoEgamma/ElectronIdentification/data 
```

## Add the modification needed to use post-fit quantities for electrons  

```shell
git cms-addpkg TrackingTools/TransientTrack
git cms-merge-topic -u CMSBParking:GsfTransientTracks
```

## Add the modification needed to use the KinematicParticleVertexFitter  

```shell
git cms-merge-topic -u CMSBParking:fixKinParticleVtxFitter
```

## Add the BParkingNano package and build everything

```shell
git clone git@github.com:friti/BParkingNANO.git  ./PhysicsTools
git cms-addpkg PhysicsTools/NanoAOD
scram b
```

## To run on a test file

```shell
cd PhysicsTools/BParkingNano/test/
cmsenv 
cmsRun run_nano_cfg.py
```

***

# For UL
```shell
cmsrel CMSSW_10_6_14
cd CMSSW_10_6_14/src
cmsenv
git cms-init
git cms-merge-topic -u friti:TransientTracks
git cms-merge-topic -u friti:KinParticleVtxFitter
git clone git@github.com:friti/BParkingNANO.git  ./PhysicsTools
git cms-addpkg PhysicsTools/NanoAOD
scram b
```
