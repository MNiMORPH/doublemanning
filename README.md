[![DOI](https://zenodo.org/badge/388610692.svg)](https://zenodo.org/badge/latestdoi/388610692)

# rating-curve-2x-manning

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

# Installation

## From PyPI using Pip

This command will install the most recent stable release of `doublemanning`.

```bash
pip install doublemanning
```

## Editable, from a local directory

These instructions allow you to use the most recent version of `doublemanning` and to make your own edits.

These instructions assume that you have the GitHub CLI installed. If you do not, just change the repository-cloning line to a standard `git` command.

```bash
gh repo clone MNiMORPH/doublemanning
cd doublemanning
pip install -e .
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
