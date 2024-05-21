# TMEscore

[![PyPI version](https://badge.fury.io/py/TMEscore.svg)](https://badge.fury.io/py/TMEscore)
[![License: GPL-2]

`TMEscore` is a Python package that implements the ESTIMATE algorithm, ISTMEscore method, and SIA method. The ESTIMATE and ISTMEscore methods were originally available only in R, and we've ported them to Python for broader accessibility. Additionally, the SIA method, which did not have any existing packages, has been manually implemented in Python.

## Features

- Implementation of the ESTIMATE algorithm for estimating stromal, immune and estimate scores in tumor samples. Estimate tumor purity for Affymetrix platform data. 
- Implementation of the ISTMEscore method for improved tumor microenvironment (TME) immune and stromal scoring. The ISTME TME subtypes are also provided.
- Novel implementation of the SIA method for comprehensive TME analysis.

## Installation

You can install the package via pip:

```bash
pip install TMEscore
```

Or install from the source code:
```
git clone https://github.com/yourusername/TMEscore.git
cd TMEscore
pip install .
```

## Usage

Here are some basic usage examples:

```
# Example 1: Using the ESTIMATE algorithm
from TMEscore import estimateScore
estimate_scores = estimateScore.ESTIMATEscore(data)
print(estimate_scores)

# Example 2: Using the ISTMEscore method
from TMEscore import ISTMEscore
istmescore = ISTMEscore.ISTMEscore(data)
print(istmescore)

# Example 3: Using the SIA method
from TMEscore import SIAscore
sia_score = SIAscore.sia_score(data)
print(sia_score)
```

## Documentation
Detailed documentation is available at [your documentation link].

## Contributing
Contributions are welcome! 

## License
This project is licensed under the MIT License. See the LICENSE file for more details.

## Contact
If you have any questions or feedback, feel free to open an issue on GitHub Issues.

## Acknowledgements
The ESTIMATE algorithm from Yoshihara et al.
The ISTMEscore method from Zeng et al.
The SIA method from Mezheyeuski et al.

## Citations

If you use TMEscore in your research, please cite the following papers:

Yoshihara, K., Shahmoradgoli, M., Martínez, E. et al. Inferring tumour purity and stromal and immune cell admixture from expression data. Nat Commun 4, 2612 (2013). https://doi.org/10.1038/ncomms3612

Zeng, Z., Li, J., Zhang, J. et al. Immune and stromal scoring system associated with tumor microenvironment and prognosis: a gene-based multi-cancer analysis. J Transl Med 19, 330 (2021). https://doi.org/10.1186/s12967-021-03002-1

Mezheyeuski, A., Backman, M., Mattsson, J., Martín-Bernabé, A., Larsson, C., Hrynchyk, I., Hammarström, K., Ström, S., Ekström, J., Mauchanski, S., Khelashvili, S., Lindberg, A., Agnarsdóttir, M., Edqvist, P. H., Huvila, J., Segersten, U., Malmström, P. U., Botling, J., Nodin, B., Hedner, C., … Sjöblom, T. (2023). An immune score reflecting pro- and anti-tumoural balance of tumour microenvironment has major prognostic impact and predicts immunotherapy response in solid cancers. EBioMedicine, 88, 104452. https://doi.org/10.1016/j.ebiom.2023.104452





