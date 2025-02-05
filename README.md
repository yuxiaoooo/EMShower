# EMShower Studies at FASER

Hola!! This repository contains scripts for EM Shower studies with testbeam 2023 at the FASER experiment. Below is an overview of the script structure and their functionalities.

## Utils

- **`emshower.py`** - using DBSCAN + clustering methods to reconstruct shower directions and EM showers.
- **`showerrec.py`** - define segment class; count cylinder
- **`robustfit.py`** - robust fitting module

## Scripts

- **`truth_MC_study.ipynb`** – Plot the variable distributions using the truth MC information
- **`dbscan_MC.ipynb`** - Study the DBSCAN Clustering Method on electron MC, and compare its performance with our old method.

## Data (*.csv)

- **`data/testbeam`** - Testbeam data
- `data/testbeam/select` - Original scanned data with a basic preliminary cut `nseg>=3`
- `data/testbeam/cylinder` - testbeam data after cylinder selections. The segments and tracks in the same cylinder are regarded as the same event, and are assigned the same event ID.
- `data/testbeam/shower` - testbeam data shower information. Including important shower variables.

- **`data/MC`** - MC data
- `data/MC/select.csv`: MC sample with a basic preliminary cut `nseg>=3``
- `data/MC/axis_old.csv`: axis information obtained by our "old method"
- `data/MC/nue_truth.csv`: truth information in nueCC 
- `data/MC/cylinder_truth.csv`: MC sample after the **truth** cylinder selection. Which means the cylinder axis is the truth primary electron direction recorded in nueCC samples.
- `data/MC/shower_truth.csv`: MC sample truth shower information after the **truth** cylinder selection. Including important shower variables.
