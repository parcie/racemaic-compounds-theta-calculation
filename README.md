# racemaic-compounds-theta-calculation
this project is written for calculating the theta angle that defined by our article in racemic compound crystals. it is mostly aimed at crystals whose Z prime value is 1.

This project uses:
	•	Python 3.9.0
	•	numpy 1.26.4
	•	pandas 2.2.3
	•	gemmi 0.6.7
	•	ccdc 3.3.1
	•	openpyxl 3.1.5

The CSD Python API can be installed on all supported platforms via conda from the CCDC conda channel. The  detailed instructions can be found in the CCDC official tutorials at 
https://downloads.ccdc.ac.uk/documentation/API/ installation_notes
Or, refer to code in GitHub:
https://github.com/ccdc-opensource/csd-python-api-scripts.
The ccdc module need official license that given by CCDC

To recreate the environment:
conda env create -f environment.yml
conda activate csd
