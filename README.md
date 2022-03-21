[![DOI](https://zenodo.org/badge/304295763.svg)](https://zenodo.org/badge/latestdoi/304295763)

# BirdGAN

This repository includes data and code for the paper "Using Generative Adversarial Networks to Simulate Central-Place Foraging Trajectories"

# Data

Trajectory datasets are contained in the **data** folder. These trajectories consist in foraging trips of distinct breeding seabirds species. During breeding period, seabirds need to feed and protect their chick. Therefore, they perform relatively short central-place foraging trips, where they leave their nest to get some food from the ocean and then they return to their nest to feed and protect their offspring. 

| Species           | Country   | Nb of trips |
|-------------------|-----------|-------------|
| *Sula variegata*  | Peru      | 78          |
| *Sula dactylatra* | Brazil    | 50          | 
| *Sula sula*       | Brazil    | 30          |


# Code

Code for the paper are available in the **code** folder : 

- **1_gan_selection_20_steps_SS.ipynb** : consist in the comparison of different GAN architecture (CNN or LSTM) on a simplified dataset
- **2_gan_vs_hmm_200_steps_SD.ipynb** : consist in the training of DCGAN with spectral regularization on the **Sula dactylatra** dataset 
- **2_gan_vs_hmm_200_steps_SV.ipynb** : consist in the training of DCGAN with spectral regularization on the **Sula variegata** dataset 

# Tutorial 

A tutorial notebook is also proposed **tutorial.ipynb**. This notebook can be directly run on Google Colab by clicking [here](https://colab.research.google.com/github/AmedeeRoy/BirdGAN/blob/master/tutorial.ipynb)

