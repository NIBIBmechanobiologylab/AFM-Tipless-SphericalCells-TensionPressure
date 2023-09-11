# "Actomyosin Cortical Mechanical Properties in Nonadherent Cells Determined by Atomic Force Microscopy" 

## **What is this analysis?**

This repository is based on the methodology developed by Cartagena-Rivera et al.: “Actomyosin Cortical Mechanical Properties in Nonadherent Cells Determined by Atomic Force Microscopy” which measures mechanics such as surface tension and hydrostatic pressure of soft spherical specimens.

![image](https://github.com/mechanobiologylab/test/assets/104796244/54e2ce79-6a35-4e41-b911-cec7e9b08b1f)

Specifically, this is a Matlab-based analysis of Bruker AFM data output files (Bruker Catalyst) with accompanying GUI that enables easier visualization and user-defined fitting of linear regions of force curves. This analysis is meant to analyze force curves taken using soft tipless AFM cantilevers, as presented in the [main document](https://www.sciencedirect.com/science/article/pii/S0006349516302375?via%3Dihub#sec2).

The software package was created by Dr. Cameron Parvini (2021).

### Matlab requirements:
Matlab (supports R2020b)
  - Curve Fitting Toolbox

### Installation

### Getting started

Run the GUI "TiplessGUI.mlapp" through Matlab and specify a candidate directory containing your Bruker files using the "Select Path" option. Furthermore, you can specify optional arguments for stiffness (nN/nm) (already stored in individual Bruker raw output file) and cell radius. Fitting regions of a curve can be adjusted by changing the start and end index values, which can be visualized using "Re-Plot". Once a region has been fitted, hit "Accept" to proceed to the next curve. After the batch of curves are completed, .xlsx and .txt summary files are generated.

### Citations
[[1]](https://www.sciencedirect.com/science/article/pii/S0006349516302375?via%3Dihub#sec2)	Cartagena-Rivera, A. X., et al. (2016). "Actomyosin Cortical Mechanical Properties in Nonadherent Cells Determined by Atomic Force Microscopy." Biophys J 110(11): 2528-2539.
