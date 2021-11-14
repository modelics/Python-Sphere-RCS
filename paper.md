---
title: 'Python Sphere RCS: a Library to Compute Mie Series Scattering of Lossy Dielectric Spheres'
tags:
  - Python
  - computational electromagnetics
  - mie series
  - scattering
  - dielectric spheres
  - Radar Cross Section
authors:
  - name: Iliya Shofman # note this makes a footnote saying 'co-first author'
    orcid: 0000-0003-1798-8337
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Damian Marek # note this makes a footnote saying 'co-first author'
    affiliation: 1
  - name: Shashwat Sharma # note this makes a footnote saying 'co-first author'
    affiliation: 1
  - name: Piero Triverio # note this makes a footnote saying 'co-first author'
    affiliation: 1
affiliations:
 - name: Dept. of Electrical & Computer Engineering, University of Toronto
   index: 1

date: 11 November 2021
bibliography: paper.bib


---

# Summary

Mie series provide an analytic solution to the problem of scattering of a plane electromagnetic wave by a homogenous sphere. The aggregate result of this calculation is known as the Radar Cross Section (RCS), which is a pattern of the relative intensity of the scattered electromagnetic wave. Knowledge of the Radar Cross Section, and its dependence on frequency, angle, sphere material, etc., are useful in many engineering and scientific applications ranging from Radar Imaging to atmospheric science. Furthermore, RCS values obtained from the analytic results from the Mie series are used to verify the correctness of numerical results from computational electromagnetic solvers. 

![Pattern of scattered electric field intensity for different angles of observation.](compare_bistatic_materials.png)

# Statement of need

While several open-source libraries for Mie series RCS calculations exist, those codes suffer from a lack of robustness for highly conductive objects. A literature review revealed that several existing libraries were unable to yield results for lossy dielectrics and at high frequencies – a significant problem considering the importance of these regimes for a wide variety of electromagnetic engineering problems.

`Python Sphere RCS` can perform calculations in the high-frequency high-conductivity regime by incorporating arbitrary precision mathematics in the critical portions of the calculation, while solving several computational issues that beset the existing publicly available solvers. Results from this solver were cross validated with numerical solvers to good visual agreement (see plot below).

![Validation Plot against Results obtained with Method of Moments Solver.](bistatic_validation.png)

The focus of this code is producing and plotting RCS through easy-to-use interface functions. Examples of various plots that can be made are provided along with the code used to generate them. The `Python Sphere RCS` library can benefit the broader scientific community to enable research and educational endeavors. 


# Acknowledgements
We acknowledge the support of the Natural Sciences and Engineering Research Council of Canada (NSERC).
