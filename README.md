# TauRecoCMSPOS

```
# setup the CMSSW release
cmsrel CMSSW_9_2_10
cd CMSSW_9_2_10/src
cmsenv

# checkout the necessary packages
git-cms-addpkg /DataFormats/PatCandidates/
git-cms-addpkg /DataFormats/TauReco/
git-cms-addpkg /PhysicsTools/PatAlgos/
git-cms-addpkg /RecoTauTag/RecoTau/

# add Jan’s repository
git remote add jan https://github.com/steggema/cmssw.git
git fetch jan

# create a new branch and make it point to Jan’s
git checkout CMSSW_9_2_X_TauRecoMiniAOD

# now compile
scram b -j16

# move to the RecoTau package
cd RecoTauTag/RecoTau/test/

# get the cfg file and run it
wget https://gist.githubusercontent.com/steggema/f5b65e09ebee723f5f27b1bef53dfa03/raw/05e064b65d0fe4c3aeeb65b028ccbda857b27814/tau_miniaod.py
voms-proxy-init --voms cms
cmsRun tau_miniaod.py
```


