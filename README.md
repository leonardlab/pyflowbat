# PyFlowBAT: An Open-Source Software Package for Performing High-Throughput Batch Analysis of Flow Cytometry Data

<img src="./pyflowbat/resources/pyflowbat-logo.png" width=40% height=40%>

## Requirements

PyFlowBAT requires Python version `>=3.10` and the following packages:

```text
numpy >= 1.24.2
matplotlib >= 3.7.1
FlowCal >= 1.3.0
statsmodels >= 0.13.5
scipy >= 1.10.1
pandas >= 2.0.0
```

To install PyFlowBAT and all required dependencies, run

```text
pip install -e pyflowbat
```

This will install `pyflowbat` to your Python `site-packages` like a regular `pip` package.

### Linux installation note

On some versions of Linux, the Python tkinter package may not be installed.
Since this is a standard library package, it cannot be installed using pip.
To install tkinter, run the following&mdash;or equivalent for your Linux distro&mdash;command:

```text
sudo apt get install python3-tkinter
```
