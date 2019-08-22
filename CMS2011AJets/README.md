# CMS2011AJets

The **CMS2011AJets** collection uses the Jet Primary Dataset from Run 2011 A of the CMS Open Data as well as associated Monte Carlo datasets.

## MODProducer

Contained in the `producer` folder, this component contains code used to extract events from the CMS-provided AOD files into portable, text-based MOD files.

## MODAnalyzer

Contained in the `analyzer` folder, this component contains C++ code used to extract jets from the MOD files as well as python code and jupyter notebooks used for the rest of the analysis.

## JetPrimaryDataset

Contains auxilliary files relevant for the Jet Primary Dataset including Run 2011 luminosity information, jet energy corrections, and file manifests.

## QCDSimDatasets

Contains auxilliary files relevant for the Monte Carlo datasets including jet energy corrections and file manifests.