# SpacePhenotyper
SpacePhenotyper: A computational method introduces a framework for translating information between different (but related) types of data

## Installing required Python packages

#### OS console/terminal

```
python -m pip install numpy pandas seaborn matplotlib scikit-survival numba
```

#### Jupyter

```
%pip install numpy pandas seaborn matplotlib scikit-survival numba
```

#### Conda

```
conda install -c conda-forge numpy pandas seaborn matplotlib scikit-survival numba
```

#### From inside the Python console itself

```python
import subprocess
import sys

subprocess.check_call([
    sys.executable, "-m", "pip", "install",
    "numpy", "pandas", "seaborn", "matplotlib",
    "scikit-survival", "numba"
])
```


## Usage in Python

#### Example

```python
import SpacePhenotyper

dfPredictions = SpacePhenotyper.SpacePhenotyper_TrainAndPredict(dfBulk, dfPhi, dfSRT, dfSRTCoords, phenotype="continuous", CV=5, CVNumShuffles=100, percentile=0.10, numPCs=10, maxNumPCs=30)[0]
```

#### Parameters

- `dfBulk` (`pandas.DataFrame`): DataFrame of bulk gene expression. Rows are samples, and columns are genes.
- `dfPhi` (`pandas.DataFrame`): DataFrame of phenotype quantities for bulk samples. Rows are samples. In the case of continuous phenotypes, the single column is the phenotype quantity for each sample. In the case of hazard, there need to be two columns named 'vital_status' and 'survival_time', containing patient status (integer 1 if dead, 0 if alive) and survival/censoring time.
- `dfSRT` (`pandas.DataFrame`): DataFrame of SRT gene expression. Rows are spatial spots, and columns are genes.
- `dfSRTCoords` (`pandas.DataFrame`): DataFrame of coordinates for SRT spots. Rows are spatial spots, and columns 'array_row' and 'array_col' are 2-dimensional y and x coordinates.
- `CV` (`int`): The number of folds for cross-validation. If the number is not greater than 1, no cross-validation will be done, and the provided parameter `numPCs` will be used. Otherwise, the value of `numPCs` will be automatically determined by cross-validation.
- `CVNumShuffles` (`int`): The number of random shuffles of samples before splitting them into folds for cross-validation. The total number of pairs of cross-validation training and test datasets will be `CV` * `CVNumShuffles`.
- `percentile` (`float`): The top percentage of genes that will be retained after sorting them by the absolute value of their association with the phenotype.
- `numPCs` (`int`): The number of top principal directions per gene set that will be used in the model (see Methods of the SpacePhenotyper paper). This parameter will be decided automatically by cross-validation if `CV`>1, and need not be supplied in that case.
- `maxNumPCs` (`int`): The maximum number of top principal directions per gene set that can be considered by cross-validation. The actual maximum will be the minimum of this parameter and the number of top principal directions that suffice to reach 99% of variance in the data (see Methods of the SpacePhenotyper paper).

#### Output

- `dfPredictions` (`pandas.DataFrame`): DataFrame of predictions of phenotype quantities in SRT spots. Columns `array_row` and `array_col` are the y and x coordinates, and the column `SpacePhenotyper` contains the predicted phenotype quantity.


## How to run on the data used in the paper

The `./SpacePhenotyper-code/main.py` file contains a short number of lines needed to reproduce the results of the SpacePhenotyper paper, and it can be run directly from console via `python main.py` or its contents can be imported into a Jupyter notebook.

All the used real and synthetic data is found in the `Data` folder.

It makes use of the `SpacePhenotyper.SpacePhenotyper_IOWrapper()` function that handles some of the IO details for the data used in the paper, and all the output data is stored in the folder represented by the provided `outputRoot` parameter.

## Reproducing the synthetic data generation

See `README.txt` in the Simulator folder.
