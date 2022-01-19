## PhotonID Dev Workflow

Photon ID Trainer can run on any trees provided that the tree contains all the variables required by the MVA defenition.

#### Step 1
  Prepare the trianing/testing trees from the BMMG ntuples. One can craft new variables from the existing data and add to the tree.
  Configure the cfg file to piont to proper files and setup the vaules properly
  ```
  cmsenv
  make mvaDataMaker
  ./mvaDataMaker configs/abc.config
  ```
#### Step 2
  Edit the Trainer configs in the `config/` folder. Add the methods, files and variables of your interest.
  ```
  cmsenv
  make mvaDataTrainer
  ./mvaDataTrainer configs/pqr.config
  ```
  ___
  
  NOTE :
  Funtions like `AddSCTree("subLeadPi0Gamma_SCTree")` and `AddSCHistos("mergedPi0_")` provides more extendability of the codebase, so please add similar function or extend the same funtions while adding new vars.
