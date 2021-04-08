# PowerModelsDistribution.jl

<img src="https://lanl-ansi.github.io/PowerModelsDistribution.jl/dev/assets/logo.svg" align="left" width="200" alt="PowerModelsDistribution logo">

[![CI](https://github.com/lanl-ansi/PowerModelsDistribution.jl/workflows/CI/badge.svg)](https://github.com/lanl-ansi/PowerModelsDistribution.jl/actions?query=workflow%3ACI) [![Documentation](https://github.com/lanl-ansi/PowerModelsDistribution.jl/workflows/Documentation/badge.svg)](https://lanl-ansi.github.io/PowerModelsDistribution.jl/stable/)

PowerModelsDistribution.jl is an extention package of PowerModels.jl for Steady-State Power Distribution Network Optimization.  It is designed to enable computational evaluation of emerging power network formulations and algorithms in a common platform.  The code is engineered to decouple problem specifications (e.g. Power Flow, Optimal Power Flow, ...) from the power network formulations (e.g. AC, linear-approximation, SOC-relaxation, ...).
This enables the definition of a wide variety of power network formulations and their comparison on common problem specifications.

## Core Problem Specifications

- Power Flow (pf)
  - ACP, ACR, IVR, LinDist3Flow, NFA, DCP
- Optimal Power Flow (opf)
  - ACP, ACR, IVR, LinDist3Flow, NFA, DCP
- Continuous load shed, minimum load delta (mld)
  - ACP, LinDist3Flow, NFA
- Optimal Power Flow with on-load tap-changer (opf_oltc)
  - ACP

**Note: The documentation is somewhat lagging behind development and the parings of network features with problem specifications with formulations has not been enumerated. We are working to correct this for you.**

## Core Network Formulations

- Nonlinear
  - ACP
  - ACR
  - IVR
- Relaxations
  - SDP BFM
  - SOC BFM and BIM relaxation (W-space)
- Linear Approximations
  - LinDist3Flow
  - NFA
  - DCP

## Network Data Formats

- OpenDSS ".dss" files

## Documentation

Please see our [online documentation](https://lanl-ansi.github.io/PowerModels.jl/stable/) for information about how to install and use PowerModelsDistribution.

## Examples

Examples of how to use PowerModelsDistribution can be found in the main documentation and in Jupyter Notebooks inside the `/examples` directory

## Development

Community-driven development and enhancement of PowerModelsDistribution are welcome and encouraged. Please fork this repository and share your contributions to the master with pull requests.

## Acknowledgments

This code has been developed as part of the Advanced Network Science Initiative at Los Alamos National Laboratory.  The primary developers are David Fobes(@pseudocubic) and Carleton Coffrin(@ccoffrin) with support from the following contributors:

- Sander Claeys (@sanderclaeys) KU Leuven, transformer models and ACR formulation
- Frederik Geth (@frederikgeth) CSIRO, Distribution modeling advise

## Citing PowerModelsDistribution

If you find PowerModelsDistribution useful for your work, we kindly request that you cite the following [publication](https://doi.org/10.1016/j.epsr.2020.106664):

```bibtex
@article{FOBES2020106664,
title = "PowerModelsDistribution.jl: An open-source framework for exploring distribution power flow formulations",
journal = "Electric Power Systems Research",
volume = "189",
pages = "106664",
year = "2020",
issn = "0378-7796",
doi = "https://doi.org/10.1016/j.epsr.2020.106664",
url = "http://www.sciencedirect.com/science/article/pii/S0378779620304673",
author = "David M. Fobes and Sander Claeys and Frederik Geth and Carleton Coffrin",
keywords = "Nonlinear optimization, Convex optimization, AC optimal power flow, Julia language, Open-source",
abstract = "In this work we introduce PowerModelsDistribution, a free, open-source toolkit for distribution power network optimization, whose primary focus is establishing a baseline implementation of steady-state multi-conductor unbalanced distribution network optimization problems, which includes implementations of Power Flow and Optimal Power Flow problem types. Currently implemented power flow formulations for these problem types include AC (polar and rectangular), a second-order conic relaxation of the Branch Flow Model (BFM) and Bus Injection Model (BIM), a semi-definite relaxation of BFM, and several linear approximations, such as the simplified unbalanced BFM. The results of AC power flow have been validated against OpenDSS, an open-source “electric power distribution system simulator”, using IEEE distribution test feeders (13, 34, 123 bus and LVTestCase), all parsed using a built-in OpenDSS parser. This includes support for standard distribution system components as well as novel resource models such as generic energy storage (multi-period) and photovoltaic systems, with the intention to add support for additional components in the future."
}
```

The associated Power Systems Computation Conference talk can be found on [YouTube](https://youtu.be/S7ouz2OP0tE)

## License

This code is provided under a BSD license as part of the Multi-Infrastructure Control and Optimization Toolkit (MICOT) project, LA-CC-13-108.
