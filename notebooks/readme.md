# Notebooks

## Run the notebooks locally

In order to run the notebooks locally, you will first need a local copy of the repository. You can clone the repository in a local folder:

```bash
git clone https://github.com/jmiramont/benchmark-test.git
```

After this, you can use  [```poetry```](https://python-poetry.org/docs/) to easily install the dependencies of the project in a virtual environment.

```bash
poetry install -E notebooks
```

The ```-E notebooks``` option will install the Jupyter dependencies needed for running the notebooks.
