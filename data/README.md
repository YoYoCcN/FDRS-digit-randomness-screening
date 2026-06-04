# Data directory

Place input data files here or in the working directory specified in the R scripts.

## Input format

Each input file should be a plain-text file with one numerical value per line and no header.

Example:

```text
0.3167
0.4281
0.5029
0.3915
```

## Default filenames used by the scripts

- `RawData.txt`
- `ErrData.txt`
- `RealRawData2.txt`
- `RealRawData3.txt`

## Data-sharing note

Raw numerical research data may contain unpublished or sensitive experimental information. Before uploading any data to GitHub, confirm that the data are approved for public release.

The `.gitignore` file is configured to avoid accidentally committing common data files in this directory.
