# brca-NOVUS
Protocol for classifying variants linked with breast cancer in BRCA1 and BRCA2 genes. 


## About
The brca-NOVUS ML tool has been created to facilitate the classification of variants of uncertain significance. We have trained two models on a one-of-its kind database comprising of gold-standard variants classified by ACMG/AMP Guidelines.  The algorithm used is XGBoost, since traditionally tree ensemble models are relied on for tabular data classification. This README document illustrates the installation, preprocessing and procedure for running both the models.


## Installation
Download and install the latest version of [Anaconda](https://docs.anaconda.com/anaconda/install/linux/). Next, download and install [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/user-guide/download/) and [Loftee](https://github.com/konradjk/loftee) for hg38 through this [link](https://clingen.igib.res.in/ML_install.tar.gz). The unzipped folder will install both the tools through the instructions given below.
**Around 100GB free space would be required to download the file, and around 400GB would be needed for the unzipped folder.**
Note: Our pipeline has been standardized on Ubuntu 18.04.6 LTS. 

Unzip the folder and install the tools using the following commands:
```
tar -xvzf ML_install.tar.gz
cd ML_install
bash INSTALL.sh
```

Next, place this repository inside the unzipped folder:
```
git clone https://github.com/aastha-v/brca-NOVUS.git
```

Now create and activate the brcaNOVUS conda environment:
```
conda env create -f brca-NOVUS/brca-NOVUS.yml
conda activate brcaNOVUS
```

## Using the brca-NOVUS tool
Prepare your VCF file as shown in the brca-NOVUS/input_folder/brca1.vcf file. Please do not Add "chr" in front of the chromosome number. To run your own file, place your vcf in brca-NOVUS/input_folder.
To process the VCF into the appropriate input format and generating predictions for your BRCA1 variants using the brca-NOVUS tool, run the first command below. This will generate a file called 'predictions_brca1.csv', and place it in the brca-NOVUS/BRCA1 folder. Similarly the second command will save 'predictions_brca2.csv' to the brca-NOVUS/BRCA2 folder. The predictions files contain pathogenicity predictions where '1' indicates that the variant is pathogenic and '0' indicates it is benign.
```
bash brca-NOVUS/scripts/predictions_brca1.sh
bash brca-NOVUS/scripts/predictions_brca2.sh
```
In case the command does not work, please try using absolute paths in the vep command in preprocessing_brca1/2.sh. If that does not work, or if the VEP plugin does not work, please run the following:
```
bash brca-NOVUS/troubeshooting.sh
```


To process a VCF with a different name than "sample.vcf", modify the "input_filename" in the predictions_brca1.sh 1/2.sh script and then run it:
```
input_filename='myvcf.vcf'
```

## Generating a New Model
In order to train a new model based on your own data, the following commands can be run. Place your VCF in the brca-NOVUS/input_folder (example: train_brca1.vcf) along with a file in which the first column represents your variant, and the second its pathogenicity (example: op_outcome_train_brca1). The first command will create and save your model as brca-NOVUS/model_creation/selfmodel_brca1.txt, and generate and save the model's performance metrics in the brca-NOVUS/model_creation/brca1 folder.
```
brca-NOVUS/scripts/model_creation_brca1.sh
brca-NOVUS/scripts/model_creation_brca2.sh
```

In order to generate predictions using your models, run the following commands. The predictions will be saved to the brca-NOVUS/model_creation/brca1 folder as selfpredictions_brca1.csv
```
python3 brca-NOVUS/scripts/self_predictions_brca1.py 
python3 brca-NOVUS/scripts/self_predictions_brca2.py 
```
