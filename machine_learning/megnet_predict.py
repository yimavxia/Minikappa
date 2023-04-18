from pathlib import Path

import numpy as np
import pandas as pd
import pymatgen as pmg
import tensorflow as tf
from megnet.config import set_global_dtypes
from megnet.models import MEGNetModel
from pymatgen.io.vasp.inputs import Poscar

"""
author: @dgaines2
An example script to perform predictions of kL_min for structures at 600 K
"""

set_global_dtypes("32")
tf.compat.v1.disable_eager_execution()

TEMPERATURE = 600.0
NUMBER_OF_MEGNET_PREDICTIONS = 10

# Folder for containing folder with {unitcell_filename} for each compound to predict
compound_directories = list(sorted(Path("structures").glob("*")))
print(f"{len(compound_directories)} structures found")

# Set name of {poscar} that can be found in each folder
unitcell_filename = "POSCAR-unitcell"

# Set the temperature, state vecotr is list of lists for MEGNet
state = [[TEMPERATURE]]
print(f"Temperature = {TEMPERATURE}")
print()

# Create list of all structures for prediction
names = []
structures = []
for compound_directory in compound_directories:
    dir_name = compound_directory.name
    names.append(dir_name)

    poscar_path = compound_directory / unitcell_filename
    poscar = Poscar.from_file(
        poscar_path, check_for_POTCAR=False, read_velocities=False
    )
    structure = poscar.structure.copy()
    structure.state = state
    structures.append(structure)

# Load MEGNet Model
model = MEGNetModel.from_file("models/model.hdf5")

# Run predictions and average
y_pred = []
for i in range(NUMBER_OF_MEGNET_PREDICTIONS):
    y_pred.append([model.predict_structure(x)[0] for x in structures])
y_pred = np.mean(y_pred, axis=0)

# Save predictions to csv
megnet_pred_df = pd.DataFrame.from_dict(
    {
        "name": names,
        "kL_min_pred": y_pred,
    }
)
print(megnet_pred_df)

filename = "megnet_predictions.csv"
megnet_pred_df.to_csv(filename, index=False)
print(f"Saved predictions to {filename}")
