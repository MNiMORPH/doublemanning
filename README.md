[![DOI](https://zenodo.org/badge/388610692.svg)](https://zenodo.org/badge/latestdoi/388610692)

# doublemanning

## Package contents

This package contains two command-line utilities
* `doublemanning-fit` inverts stage&mdash;discharge data to generate a rating curve and the associated double-Manning-equation parameters:
  * A Manning's-equation relationships for an approximately rectangular channel
  * A generalized Manning's equation (power-law) relationship for flows across the floodplain
* `doublemanning-calc` uses this fit to perform forward computations:
  * Stage &rarr; discharge
  * Flow depth &rarr; discharge
  * Discharge &rarr; stage
  * Discharge &rarr; flow depth

## Installation

### From PyPI using Pip

This command will install the most recent stable release of `doublemanning`.

```bash
pip install doublemanning
```

### Editable, from a local directory

These instructions allow you to use the most recent version of `doublemanning` and to make your own edits.

These instructions assume that you have the GitHub CLI installed. If you do not, just change the repository-cloning line to a standard `git` command.

```bash
gh repo clone MNiMORPH/doublemanning
cd doublemanning
pip install -e .
```

## Running the double-Manning software

You should be able to run both commands by simply typing their names on the command line.
Here we provide the outputs form the "`-h`" help flag.
Such outputs are also provided if you enter the commands with no arguments.

### `doublemanning-fit`

**Note:** Although many command-line options exist to fit the data, I recommend using the `YAML` configuration-file option alone. This will allow you to access the full functionality of `doublemanning-fit` (the command-line options include only a smaller subset) and will self-document your work.

```bash
>> doublemanning-fit -h
usage: doublemanning-fit [-h] [-y CONFIGFILE] [-f DATAFILE]
                         [--delimiter DELIMITER] [-b CHANNEL_WIDTH]
                         [-H CHANNEL_DEPTH] [-s SLOPE] [-o OUTFILE]
                         [--use_depth] [--us_units] [--plot] [-v]

Pass channel and flow characteristics to obtain a "Double Manning" --
Manning\'s Equation (channel) + generic power-law (floodplain) stage--discharge
-- relationship.

options:
  -h, --help            show this help message and exit
  -y CONFIGFILE, --configfile CONFIGFILE
                        YAML file from which all inputs are read.
  -f DATAFILE, --datafile DATAFILE
                        file with two columns: Discharge, Stage
  --delimiter DELIMITER
                        "tab", "comma", or "semicolon"
  -b CHANNEL_WIDTH, --channel_width CHANNEL_WIDTH
                        river-channel width
  -H CHANNEL_DEPTH, --channel_depth CHANNEL_DEPTH
                        river-channel depth (not flow depth)
  -s SLOPE, --slope SLOPE
                        channel slope
  -o OUTFILE, --outfile OUTFILE
                        Stores fit parameters.
  --use_depth           Use flow depth instead of hydraulic radius.
  --us_units            Convert imported data from cfs and feet
  --plot                Plot stage-discharge relationship
  -v, --verbose         Plot stage-discharge relationship
```

### `doublemanning-calc`

The `doublemanning-calc` program returns a scalar value to `stdout`.

```bash
>> doublemanning-calc -h
usage: doublemanning-calc [-h] [-p PARAMFILE] [-zQ STAGE_DISCHARGE]
                          [-hQ DEPTH_DISCHARGE] [-Qz DISCHARGE_STAGE]
                          [-Qh DISCHARGE_DEPTH]

Return stage or discharge based on a double-Manning fit. All values are SI
(mks).

options:
  -h, --help            show this help message and exit
  -p PARAMFILE, --paramfile PARAMFILE
                        CSV file for double-Manning parameters.
  -zQ STAGE_DISCHARGE, --stage_discharge STAGE_DISCHARGE
                        Calculate discharge from this stage.
  -hQ DEPTH_DISCHARGE, --depth_discharge DEPTH_DISCHARGE
                        Calculate discharge from this flow depth.
  -Qz DISCHARGE_STAGE, --discharge_stage DISCHARGE_STAGE
                        Calculate stage from this discharge.
  -Qh DISCHARGE_DEPTH, --discharge_depth DISCHARGE_DEPTH
                        Calculate flow depth from this discharge.
```

## Physical and mathematical basis

### Core equation

The double-Manning approach applies the following combination of Manning's equation for in-channel flows (left of the $+$ sign) and a power-law equation for overbank flow (right of the $+$ sign):

$$Q = \frac{b}{n_\mathrm{ch}} h R_h^{2/3} S^{1/2} + k_\mathrm{fp} \left(h - h_\beta \right)^{P_\mathrm{fp}}$$

### Variables

| **Variable**    | **Description**                                                                                        | **Units [SI]**                                  |
|-----------------|--------------------------------------------------------------------------------------------------------|-------------------------------------------------|
| $Q$             | Discharge                                                                                              | m$\mathrm{m}^3 \text{ s}^{-1}$                  |
| $b$             | Channel width                                                                                          | m                                               |
| $B$             | Valley-bottom width                                                                                    | m                                               |
| $B-b$           | Floodplain width                                                                                       | m                                               |
| $z_b$           | River-bed elevation (compared to an arbitrary datum)                                                   | m                                               |
| $z_s$           | River stage: water-surface elevation (compared to the same arbitrary datum)                            | m                                               |
| $h$             | Flow depth: $h = z_s - z_b$                                                                            | m                                               |
| $h_b$           | Channel-bank height                                                                                    | m                                               |
| $R_h$           | Hydraulic radius; for the assumed rectangular channel, $R_h = b \cdot h / (b + 2 (h \wedge h_\beta) )$ | m                                               |
| $n_\mathrm{ch}$ | Manning's roughness coefficient within the channel                                                     | m                                               |
| $S$             | River-channel slope                                                                                    | &mdash;                                         |
| $k_\mathrm{fp}$ | Floodplain-flow coefficient                                                                            | $\mathrm{m}^{3 - P_\mathrm{fp}} \text{ s}^{-1}$ |
| $P_\mathrm{fp}$ | Floodplain-flow exponent                                                                               | &mdash;                                         |
