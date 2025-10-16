# Chlamydomonas photosynthesis model
This model captures the basic photosynthetic electron transport chain in chlamydomonas rheinhaardtii with state transitions.
For a more detailed description, check out the model gitlab at: https://git.rwth-aachen.de/computational-life-science/chlamydomonas-photosynthesis-model 

To run a simulation simulating low light conditions, run the runmodel.ipynb file. It gets the model from the params_and_vars.py file where the basic parameters and variables, rates and derived values are added.
The description of the rate functions is in the rate.py file and the description of the derived values in the derived.py file. cfunctions.py contains helper functions. The get_model() function defined the params_and_variables.py file has an argument to use a static or dynamic photosystem II description, which is why the description of PSII is seperate. By default, a dynamic description is used.

In future versions pixi will be used. For now, you have to install dependecies manually, but there are currently only three packages necessary.

requires-python = ">=3.12"
dependencies = [
    "jupyter>=1.1.1",
    "mxlpy>=0.22.0",
    "ruff>=0.12.1",
]

For review, please give feedback on wether the general file structure is understandable. Files to review in more detail could be rates.py and derived_quantities.py and the other files could be reviewed in the next session.

Happy roasting

Cheers, Matthias