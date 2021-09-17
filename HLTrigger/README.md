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
## Helper Scripts
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
    * contains various eduts required for the running of the path

  * customizeForPath.py
    * produses a customizer file that has prescribed parameter set
