# LLWP - Luis' Loomis-Wood Program

LLWP allows you to efficiently and confidently assign (typically rotational or rovibrational) spectra by relying on Loomis-Wood plots.

A quickstart guide is given down below. For more information see its [website](https://llwp.astro.uni-koeln.de).

If you want to acknowledge LLWP, please cite the paper [LLWP - A new Loomis-Wood software at the example of Acetone-13C1](https://doi.org/10.1016/j.jms.2022.111674).


## Quickstart Guide

The preferred way to install LLWP is via Python's package manager pip.
Run the following command in a terminal to install LLWP:

```bash
pip install llwp
```

After installing LLWP via pip you can run it from any terminal by simply running

```bash
llwp
```

To see and assign your first series

1. open your spectrum and prediction files via drag and drop or *Files > Add Files*
2. specify the correct reference series in the *Reference Series* window
3. choose the fitfunction under *Fit > Choose Fit Function*
4. select the area around the experimental peak with the mouse to fit the data

### ASAP Mode

To start the [ASAP](https://doi.org/10.1016/j.jms.2015.02.014) mode of LLWP run

```bash
lasap
```

To see and assign your first cross-correlation peaks

1. open your spectrum, \*.egy, and \*.cat file via drag and drop or *Files > Add Files*
2. specify the correct energy levels in the *ASAP Settings* window
3. specify the correct unit conversion factor for the \*.cat file in the *Units Cat File* field (e.g. 3.335641e-05 for \*.cat file in MHz and \*.egy file in wavenumbers)
3. press *Calculate Cross Correlation*
3. choose the fitfunction under *Fit > Choose Fit Function*
4. select the area around the experimental peak with the mouse to fit the data