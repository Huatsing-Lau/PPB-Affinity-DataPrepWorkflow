# PPB-Affinity-DataPrepWorkflow

**This is the code for data preprocessing of PPB-Affinity**

PPB-Affinity: Protein-Protein Binding Affinity dataset for AI-based protein drug discovery

Prediction of protein-protein binding (PPB) affinity plays an important role in large-molecular drug discovery. Deep learning (DL) has been adopted to predict the change of PPB binding affinity upon mutation, but there was a scarcity of studies predicting the PPB affinity itself. The major reason is the paucity of open-source dataset concerning PPB affinity. Therefore, the current study aimed to introduce and disclose a PPB affinity dataset (PPB-Affinity), which will definitely benefit the development of applicable DL to predict the PPB affinity. The PPB-Affinity dataset contains key information such as crystal structures of protein-protein complexes (with or without protein mutation patterns), PPB affinity, receptor protein chain, ligand protein chain, etc. To the best of our knowledge, this is the largest and publicly available PPB-Affinity dataset, which may finally help the industry in improving the screening efficiency of discovering new large-molecular drugs. We also developed a deep-learning benchmark model with this dataset to predict the PPB affinity, providing a foundational comparison for the research community.

## How do we preprocess the source data-set to get PPB-Affinity？

### 1.Preprocess each source data-set separately
We preprocess the source dataset by running the notebook files separately, there are:

> process_Affinity Benchmark v5.5
> process_ATLAS.ipynb
> process_PDBbind v2020.ipynb
> process_SAbDab.ipynb
> process_SKEMPI v2.0.ipynb

### 2. Merge the preprocessed source data-sets
Run **merge.ipynb** to merge the preprocessed source data-sets.


## How to develop an affinity prediction AI utilizing PPB-Affinity？

We have provided a demonstration case to illustrate how to develop an affinity prediction AI using this dataset. For more detailed information and access to the code, please visit **https://github.com/ChenPy00/PPB-Affinity**.
