# nanoAOD producer customized for RJPSI analysis (from RK/K*/phi analysis)

# For UL
```shell
cmsrel CMSSW_10_6_14
cd CMSSW_10_6_14/src
cmsenv
git cms-init
git cms-merge-topic -u friti:TransientTracks
git cms-merge-topic -u friti:KinParticleVtxFitter
git clone --single-branch --branch ul git@github.com:friti/BParkingNANO.git  ./PhysicsTools
git cms-addpkg PhysicsTools/NanoAOD
scram b
```
