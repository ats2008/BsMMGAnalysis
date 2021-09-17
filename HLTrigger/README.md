# HLT Path Development for RUN3

## Setting up the environment
```
$ cmsrel CMSSW_X_Y_Z
$ cd CMSSW_X_Y_Z/src/
$ git clone git@github.com:ats2008/BsMMGAnalysis.git
$ git checkout hltDev
$ cd BsMMGAnalysis
$ scram b -j 2
$ cd HLTrigger/test

```
## HLT Menu Development
### Helper Scripts
bs2mmgParts.py  bsmmg_getHltMenu1130Data.sh  bsmmg_getHltMenu1130.sh  customizeForBs2MMG.py  customizeForPath.py  GetMenus  menudev  putCustIntoFiles.py
  * bsmmg_getHltMenu1130Data.sh , bsmmg_getHltMenu1130.sh
    * use these to obtain configs from confDB
      ```
      $ ./bsmmg_getHltMenu1130XX.sh
      ```
  * putCustIntoFiles.py
    * use it to customize the config file with POG suggested correction as well as BMMG customizations. The original file will be backed up with .bak extension and file will be replaced with the customized version
      ```
      $ ./putCustIntoFiles.py <pathToConfigFile>
      ```
  * customizeForBs2MMG.py
    * contains various edits required for the running of the path
    * has the customizer function for production of reduced edmTree for hlt Ntuplizing step

  * customizeForPath.py
    * produses a customizer file that has prescribed parameter set

## Rate Measurement

```
$ cd test/RateAndEfficiency
```
* Make the header and shared objects corresponding to the HLTNTuple tree, and edit `rate_production.cc` and `efficiency_production.cc` appropriately
* Use rate_production.cc for making the rate/eficiency measurements  .. configuration file required (rateConfig.cfg)
* Use efficiency_production.cc for producing root file with efficiency/distribution histograms  [ for variables  of interest ] 
* Use helper functions in parameterSetup.py to setup condor jobs using HLTNTuples
    * see jobFileMaker_efficMeasureFullWindow.py

