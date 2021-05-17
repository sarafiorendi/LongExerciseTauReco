### 1. Tau Reco

**Set up code:**
```bash
source /cvmfs/grid.cern.ch/emi3ui-latest/etc/profile.d/setup-ui-example.sh
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

export SCRAM_ARCH=slc7_amd64_gcc820
cmsrel CMSSW_10_2_26
cd CMSSW_10_2_26/src
cmsenv

# set up cmssw git repo and check out tau reco package
git cms-init
git cms-addpkg RecoTauTag/RecoTau/

# compile
scram b -j 8

# move to testing area
cd RecoTauTag/RecoTau/test/

# get code examples
git clone git@github.com:sarafiorendi/LongExerciseTauReco.git

scp -r lxplus.cern.ch:/afs/cern.ch/work/j/jbechtel/public/CMSPOS_2019/TauRECO .

cp -rs TauRECO/*root .

```
**Available Scripts**
> We advice you to always check the scripts before running them to get a feeling on what to script will actually do
>  
***Re-run the tau reconstruction on miniAOD:***

```bash
cmsRun LongExerciseTauReco/customise_tau_reco_miniaod_cfg.py
```
First we rerun the Tau Reconstruction on the provided Samples ( a ZTT Sample and a QCD Sample). During this step, different Parameters of the Tau Reconstruction can be modified, and their impact on e.g. the efficiency of the tau reconstruction can be checked. 

***Ntupelize the output file:***

```bash
python LongExerciseTauReco/read_taus.py --file (ZTT|QCD)
```
After the Tau Reconstruction is rerun, you can convert the MiniAOD file to a flat ntuple. These flat ntuples can then be used to plot and calculate various quantities. 

The script will result in two nuples:

1. containing one entry for each generated Tau
2. containing one entry for each reconstructed Jet

***Plot the efficiency curve and smearing of tau pt:***
```bash
python LongExerciseTauReco/plotting.py  --file (ZTT|QCD)
```
The plotting script is used to generate plots from the flat ntuples.

**Tasks:**

Perform the following tasks on the miniAOD Z->tau tau sample and the QCD sample provided:

1. Re-run the tau reconstruction on the both samples

2. Ntupelize the two output miniAOD files. Make sure to set the correct input file!

3. Calculate the overall tau reconstuction efficiency using the ZTT sample by applying the following cuts. You can add the calculation in the script `plotting.py`.

    * gen. tau $p_{T} > 18$ GeV
    *  gen. tau $|\eta| < 2.4$ 

    and calculating $\frac{\text{Reconstructed taus from genuine tau}}{\text{Number genuine taus}}$ (Efficiency)

4. Caclulate the overall misidentification propability using the QCD sample by applying the following cuts. You can add the calculation in the script `plotting.py`. 

    * jet $p_{T} > 18$ GeV
    * jet $|\eta| < 2.4$

    and calculating $\frac{\text{Reconstructed taus from jet}}{\text{Total jets}}$ (Misidentification probability)

5. Plot the reconstruction efficiency for for genuine taus as a function of the tau $p_{T}$, as well as the smearing of the reconstructed tau $p_{T}$ using `plotting.py`. This is already implemented.

6. **Adapt `plotting.py`** to additionally plot the reconstruction efficiency of jets misidentified as tau leptons for the QCD sample.

10. **Adapt `plotting.py`** to plot the distribution of the reconstructed and true decay modes for tau leptons passing the $p_{T}$  and $\eta$ cuts of part 3 in the same plot. Use different colors and a ROOT TLegend to label the two histograms. The variable decay_mode is an integer with `decay mode = 5 * (n_charged - 1) + n_pi0`. The most common decay modes are 0 (1 charged hadron, no pi0), 1 (1 charged hadron, 1 pi0), and 10 (three charged hadrons, no pi0)

12. Read out the tau-isolation discriminator `byLooseIsolationMVArun2v1DBoldDMwLT`. It can be accessed for reconstructed taus by requesting `tau.tauID("name of ID")` when ntupelizing the miniAOD file. 

9. Repeat steps 3.-7., however now only use reconstructed tau leptons passing the isolation discriminator.

10. How did the reconstruction efficiencies and (mainly) the misidentification probability change?

11. *For early finishers only*:
Use ```LongExerciseTauReco/customise_tau_reco_miniaod_cfg.py``` to rerun the tau reconstruction using custom settings, and see how it changes the quality of your tau reconstruction. You could change e.g.:
    * Decay modes to be considered for valid tau leptons
    * dR isolation cone size
    * Minimum quality requirements on tracks to be considered 



