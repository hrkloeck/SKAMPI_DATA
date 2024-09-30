# Obtain SKAMPI DATA 

Example datasets are stored in the [MPG keeper
service](https://keeper.mpdl.mpg.de/d/26112717a86f48b8a30d/).

Datasets are stored into subfolders OBS_<OBSERVATION_ID>, containing a
preprocessed hdf5 file and the automatically generated plots at the
time of data processing. 

The hdf5 file may have been split into several parts of 1 GiB size as e.g. 

foo.hdf5_part000
foo.hdf5_part001
foo.hdf5_part002
...

these can be joined on the commandline via $cat *_part* > foo.hdf5. Note that the single parts can not be read individually and have to be merged into a single file after download. They can be deleted after merging.

We downloaded the
[OBS_18373](https://keeper.mpdl.mpg.de/d/26112717a86f48b8a30d/?p=%2FOBS_18373&mode=list)

and the directory contains the following files:

```
EDD_2022-10-18T11:07:55.459719UTC_jvpgs.hdf5_part000
EDD_2022-10-18T11:07:55.459719UTC_jvpgs.hdf5_part001
EDD_2022-10-18T11:07:55.459719UTC_jvpgs.hdf5_part002
EDD_2022-10-18T11:07:55.459719UTC_jvpgs.hdf5_part003
EDD_2022-10-18T11:07:55.459719UTC_jvpgs.hdf5_part004
EDD_2022-10-18T11:07:55.459719UTC_jvpgs.hdf5_part005
EDD_2022-10-18T11:07:55.459719UTC_jvpgs.hdf5_part006
EDD_2022-10-18T11:07:55_120_scan-equatorial_000.pdf
EDD_2022-10-18T11:07:55_120_scan_local.pdf
EDD_2022-10-18T11:07:55_200_spectrum.pdf
EDD_2022-10-18T11:07:55_350_dynamic-spectrum.pdf
```

- the next step is to merge all the parts (EDD_2022-10-18T11:07:55.459719UTC_jvpgs.hdf5_part000 ... part006)
  into a single file.
  ```
  cat *_part* > EDD_2022-10-18T11:07:55.459719UTC_jvpgs.hdf5
  ```

- for completness the pdf files are stored here.

  [EDD_2022-10-18T11:07:55_120_scan-equatorial_000.pdf](https://github.com/hrkloeck/SKAMPI_DATA/blob/main/obtain_example_dataset/OBS_18373/EDD_2022-10-18T11%3A07%3A55_120_scan-equatorial_000.pdf)
  
  [EDD_2022-10-18T11:07:55_120_scan_local.pdf](https://github.com/hrkloeck/SKAMPI_DATA/blob/main/obtain_example_dataset/OBS_18373/EDD_2022-10-18T11%3A07%3A55_120_scan_local.pdf)

Please note that SKAMPI is measureinig 
  [EDD_2022-10-18T11:07:55_120_spectrum.pdf](https://github.com/hrkloeck/SKAMPI_DATA/blob/main/obtain_example_dataset/OBS_18373/EDD_2022-10-18T11%3A07%3A55_200_spectrum.pdf)

  [EDD_2022-10-18T11:07:55_120_dynamic-spectrum.pdf](https://github.com/hrkloeck/SKAMPI_DATA/blob/main/obtain_example_dataset/OBS_18373/EDD_2022-10-18T11%3A07%3A55_350_dynamic-spectrum.pdf)


