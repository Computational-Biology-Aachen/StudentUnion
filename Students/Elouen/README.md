# Elouen's Code

## GreenSloth

In this directory you will find code that is part of Elouen's [GreenSloth](https://github.com/ElouenCorvest/GreenSloth) project. To run any code found in this directory, please use the `pixi` environment provided. To do so, install [`pixi`](http://pixi.sh/latest/) and run the following inside of this directory:

```bash
pixi install
```

If you wish to use another way to install python and other packages, you can find the requirements inside the [`pixi.toml`](GreenSloth/pixi.toml).

### Docs

The files containing actual project code are: [`model_validation.py`](GreenSloth/model_validation.py)

Each function **should** have an intiutive docstring that explains what each function does. However, an easier way to see the docs is by opening the generated [HTML](GreenSloth/docs/index.html) through `pdoc`.

## Runners

In the [`runner.ipynb`](GreenSloth/runner.ipynb) you can run two different models and see their validation process. The models used for testing were taken from [`mxlbricks`](https://pypi.org/project/mxlbricks/). The tests are not succesful to 100%, however appropriate `TODO`s have written for future reference. This runner will create a markdown file of the saadat2021 model.