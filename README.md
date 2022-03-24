## OpenSNP to PLINK

This script is my attempt at converting raw OpenSNP data to PLINK (.bed, .bim, .fam) format for research purposes.

### Usage
1. Download OpenSNP data from [https://opensnp.org](http://opensnp.org/)
2. Clone the git
3. Extract the OpenSNP raw data to the RAW folder
4. Run convert.py
5. I ran into a plink error that the bim file still had a split chromosome? Don't know what that is, but it was easily fixed with `plink.exe --bfile opensnp --make-bed --out opensnp2`

Note from older tutorials I've seen it seems like the file-naming convention of opensnp has changed?
This script assumes they are in the following format: `user1_file9_yearofbirth_1985_sex_XY.23andme.txt` should that change this code will likely break.

### How it works
The code first verifies which files it can/can't parse and creates a fam file. This is also required to assess the total number of participants in the eventual Plink file.
Then it goes through the files it knows it can parse and converts the genotypes to PLINK format.

The current parser is far from perfect, skips quite a few files, removes triallelic SNPs, cannot handle >1 file per participant, etc. So use this at your own risk.

### Acknowledgements
Many thanks to [superbobry's snpy](https://github.com/superbobry/snpy) which is used here. 