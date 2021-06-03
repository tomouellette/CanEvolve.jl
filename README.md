# CanEvolve <img align="right" width="208" height="129" src="icon.svg">

[![Build Status](https://travis-ci.com/tomouellette/CanEvolve.jl.svg?branch=master)](https://travis-ci.com/tomouellette/CanEvolve.jl)

#### Overview

This package is used to generate realistic variant allele frequency (VAF) distributions that are observed in bulk-sequenced single tumour biopsies. We use a stochastic branching process, under a multiplicative fitness framework, to grow tumour populations under positive selection. To achieve speed and scale, we utilize a small N approximation for synthetic tumour populations that still fully recapitulates the observed VAFs in patient tumours. We note that because we use small N in our stochastic simulations, neutrally evolving synthetic tumours have a higher probability of observing spurious subclones due to stronger genetic drift. To avoid the generation of spurious subclones, we implement a generative synthetic data step for neutral tumours by sampling VAF distributions from Pareto distributions parametrized by empirically relevant shape and scape values. 

We also integrate feature engineering functions to convert the VAF distribution into deep learning friendly data structures for evolutionary inference. For deep learning applications using this synthetic data, please see *ADD GITHUB REPO HERE*


#### Acknowledgements

This package is an extension of the [CancerSeqSim.jl](https://github.com/marcjwilliams1/CancerSeqSim.jl) framework developed by Marc Williams. In addition, we generate neutral synthetic data using a Pareto distribution with shape and scale parameters retrieved from PCAWG fits by [Caravagna et al., 2020](https://www.nature.com/articles/s41588-020-0675-5).
