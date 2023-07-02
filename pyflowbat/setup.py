from setuptools import setup, find_packages

setup(
    author='Eleftheria Beres <elli.beres@u.northwestern.edu>',
    description='PyFlowBAT: An Open-Source Software Package for Performing High-Throughput Batch Analysis of Flow Cytometry Data ',
    name='pyflowbat',
    version='0.0.1',
    packages=find_packages(include=['pyflowbat','pyflowbat.*']),
    install_requires=[
         "numpy == 1.24.2",
         "matplotlib == 3.7.1",
         "FlowCal == 1.3.0",
         "statsmodels == 0.13.5",
         "scipy == 1.10.1",
         "pandas == 2.0.0"
    ],
    python_requires='>=3.10'
)
