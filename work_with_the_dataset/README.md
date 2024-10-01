# SKAMPI DATA

Here s an example on how to deal and work with the data.

```
python EXAMPLE_TO_ACCESS_AND_PLOT.py -h
```

== Example how to access and plot SKAMPI Data == 

Usage: EXAMPLE_TO_ACCESS_AND_PLOT.py [options]

Options:
  -h, --help            show this help message and exit
  --DATA_FILE=DATAFILE  DATA - HDF5 file of the Prototyp
  --DOPLOT_FINAL_SPEC   Plot the final spectrum after Flagging
  --DOPLOT_FINAL_WATERFALL
                        Plot the final waterfall after Flagging
  --DOPLOT_ON_SKY       Plot the final flags in ra and dec
  --DOSAVEPLOT          Save the plots as figures
  --DOUSEMASK           Use the default mask to plot the data
  --HELP                Show info on input



```
python EXAMPLE_TO_ACCESS_AND_PLOT.py --DATA_FILE=EDD_2022-10-18T11:07:55.459719UTC_jvpgs.hdf5 --DOPLOT_FINAL_SPEC --DOPLOT_FINAL_WATERFALL --DOPLOT_ON_SKY --DOSAVEPLOT
```

Generates for each polarisation and noise diode combination a file.

As examples:

- 1d spectrum

- waterfall spectrum

- sky coverage