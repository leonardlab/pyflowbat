# PyFlowBAT

<img src="./resources/pyflowbat-logo.png" width=40% height=40%>

PyFlowBAT is an open-source Python package for performing high-throughput batch analysis of flow cytometry data.

PyFlowBAT can be used to

---

## Requirements

PyFlowBAT needs Python version `>=3.10` and poetry installed to work.
After both have been installed, running the following will install the necessary Python packages.

```text
poetry lock
poetry install
```

### Linux installation note

On some versions of Linux, the Python tkinter package may not be installed.
Since this is a standard library package, it cannot be installed using pip.
Instead, to install run the following command, or the equivalent on your Linux distro:

```text
sudo apt get install python3-tkinter
```


### Poetry installation

#### Linux

On Ubuntu, run:

```text
curl -sSL https://install.python-poetry.org | python3 -
```

## (add documentation here)

<!-- [ ] MAKE THIS AN ACTUAL README -->

## Advanced usage

### Custom R functions

<!-- [ ] make sure below instructions are actually CORRECT -->

Custom R functions (i.e., R functions not included in the `pyflowbat.r_gating` module) can be called from a PyFlowBAT workspace after `pyflobat.Workspace.init_r()` has been called if the proper dependencies are installed.

To call custom R functions, two things are needed:

1. a Python function that uses the `rpy2` library to call the R function and
2. an R function that is called by the Python function.

The R function can itself call other R functions; it just needs to be set up in such a way that it can accept the arguments supplied by the Python call and returns data that can be used by Python.

Generally, this means that the R function will take the following form:

```R
library(flowCore)

FUNCTION_NAME <- function(data, nr, nc,
    channel_names, gating_channels, ARGUMENTS) {
    data <- matrix(unlist(data), ncol = nc, nrow = nr)
    channel_names <- unlist(channel_names)
    colnames(data) <- channel_names
    fr <- flowFrame(exprs = data)
    gating_channels <- unlist(gating_channels)
    g <- FUNCTION_CALL(ARGUMENTS)
    res <- Subset(fr, g)
    return(exprs(res))
}
```

Note that `ARGUMENTS` can be either a single object (as must be the case if the built-in `pyflowbat.r_gating_general.r_gate()` function is used to call it) or multiple parameters.

To call the custom R function from Python, there are two methods:

1. use the `pyflowbat.r_gating_general.r_gate()` function or
2. write your own Python function for calling a R functions using the `rpy2` library.

If you choose to write a new Python function, it should take the form of the following with `KEY_WORD_ARGUMENTS` being the necessary key word arguments to call your function.
Recall that like all gating functions called by PyFlowBAT, including `**kwargs` is required to ensure that PyFlowBAT can pass potentially unused keyword arguments to the function without error.

```python
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
import numpy as np
rpy2.robjects.numpy2ri.activate()

r = ro.r

def FUNCTION_NAME_HERE(data, gating_channels, KEY_WORD_ARGUMENTS, **kwargs):
    r['source'](R_FUNCTION_FILE_NAME)

    fcs_data = data.copy()
    fcs_data_nparr = np.asarray(data)

    r_func = ro.globalenv[R_FUNCTION_NAME]

    nr,nc = fcs_data_nparr.shape

    r_result = r_func(data = fcs_data_nparr.T.tolist(), nr = nr, nc = nc, channel_names = list(fcs_data.channels), gating_channels = gating_channels, KEY_WORD_ARGUMENTS)

    nr, nc = r_result.shape
    fcs_data = fcs_data[0:nr, :]
    fcs_data[:, :] = r_result
    return fcs_data
```

For detailed examples of how to write the Python and R functions, look at the functions defined in `pyflowbat.gating.r` and `pyflowbat.r_gating.py`.
