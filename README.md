# EMShower Studies at FASER

Hola!! This repository contains scripts for EM Shower studies with testbeam 2023 at the FASER experiment. Below is an overview of the script structure and their functionalities.

## Reports and Talks
* [Feb 19, 2025](https://indico.cern.ch/event/1516827/contributions/6383034/attachments/3017347/5322117/testbeamFeb19.pdf): Statistics checks and higher statistics
* [Jan 29, 2025](https://indico.cern.ch/event/1507491/contributions/6344011/attachments/3004641/5295953/testbeamJan29.pdf): MC checks and a few comparisons between MC & data; updated & original methods
* [Dec 10, 2024](https://indico.cern.ch/event/1487433/contributions/6270288/attachments/2983941/5254816/testbeamDec10.pdf): DBSCAN clustering applied
* [Oct 18, 2024](https://indico.cern.ch/event/1467312/contributions/6177681/attachments/2950416/5186242/testbeamOct18.pdf): MC analysis, and background analysis
* [Oct 8, 2024](https://indico.cern.ch/event/1462870/contributions/6158954/attachments/2943129/5171474/testbeamOct7.pdf): Analysis methodology introduced in detail at _FASER Physics Meeting_
* [Aug 27, 2024](https://indico.cern.ch/event/1449445/contributions/6102055/attachments/2916517/5118392/testbeamAug27.pdf): Paper draft, timeline and analysis sketch
* [June 26, 2024](https://indico.cern.ch/event/1365995/contributions/5973369/attachments/2884695/5055518/testbeamFASERColMeet.pdf): Analysis methodology introduced in detail at _FASER Collaboration Meeting_
* [June 11, 2024](https://indico.cern.ch/event/1425413/contributions/5995955/attachments/2875222/5035023/testbeam24June11.pdf): Updated methodology and variable distributions
* [May 28, 2024](https://indico.cern.ch/event/1420209/contributions/5971632/attachments/2865639/5015734/testbeam24May28.pdf): Efficiency, Density, and Angular Distributions
* [Apr 16, 2024](https://indico.cern.ch/event/1406120/): Variable distributions check
* [Mar 19, 2024](https://indico.cern.ch/event/1396062/): Original methodology and event display
* [Feb 6, 2024](https://indico.cern.ch/event/1379177/contributions/5798290/attachments/2794640/4874672/testbeamFeb6.pdf): Track density checks

## Coding

### Utils

- **`emshower.py`** - using DBSCAN + clustering methods to reconstruct shower directions and EM showers.
- **`showerrec.py`** - define segment class; count cylinder
- **`robustfit.py`** - robust fitting module

### Scripts

#### MC Studies

- **`truth_MC_study.ipynb`** – Plot the variable distributions using the truth MC information
- **`dbscan_MC.ipynb`** - Study the DBSCAN Clustering Method on electron MC, and compare its performance with our old method.

#### Testbeam Data Studies

- **`data_preprocess.ipynb`** - Prepare the csv data file we need.
- **`testbeam_study.ipynb`** - Apply the cylinder & kinematic selections to testbeam data, and get the variable distributions.

#### Technique Demonstrations

- **`dbscan.ipynb`** - DBSCAN demonstration.
- **`hdbscan.ipynb`** - HDBSCAN demonstration.

### Data (*.csv)

#### **`data/testbeam`** - Testbeam data
- **`data/testbeam/select`** - Original scanned data with a basic preliminary cut `nseg>=3`
- **`data/testbeam/cylinder`** - testbeam data after cylinder selections. The segments and tracks in the same cylinder are regarded as the same event, and are assigned the same event ID.
- **`data/testbeam/shower`** - testbeam data shower information. Including important shower variables.

#### **`data/MC`** - MC data
- **`data/MC/select.csv`**: MC sample with a basic preliminary cut `nseg>=3``
- **`data/MC/axis_old.csv`**: axis information obtained by our "old method"
- **`data/MC/nue_truth.csv`**: truth information in nueCC 
- **`data/MC/cylinder_truth.csv`**: MC sample after the **truth** cylinder selection. Which means the cylinder axis is the truth primary electron direction recorded in nueCC samples.
- **`data/MC/shower_truth.csv`**: MC sample truth shower information after the **truth** cylinder selection. Including important shower variables.
