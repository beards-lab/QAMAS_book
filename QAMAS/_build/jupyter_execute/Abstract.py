#!/usr/bin/env python
# coding: utf-8

# # Title Page
# 
# ## Authors
# 
# E. Benjamin Randall$^{1}$, Marcus Hock$^{2}$, Rachel Lopez$^{1}$, Bahador Marzban$^{1}$, Collin Marshall$^{1}$, Daniel A. Beard$^{1*}$
# 
# $^{1}$ *Department of Molecular and Integrative Physiology, University of Michigan, Ann Arbor, MI*
# 
# $^{2}$ *Department of Bioengineering, University of California at San Diego, San Diego, CA*
# 
# *Corresponding author
# 
# *Email addresses*: ebrandal@umich.edu (E.B. Randall), m1hock@eng.ucsd.edu (M. Hock), ralopez@umich.edu (R. Lopez), bmarzban@umich.edu (B. Marzban), colmar@umich.edu (C. Marshall), beardda@umich.edu (D.A. Beard). 
# 
# 
# ## Abstract 
# 
# We present a computational framework for analyzing and simulating mitochondrial ATP synthesis using basic thermodynamic and kinetic principles. The framework invokes detailed descriptions of the thermodynamic driving forces associated with the processes of the electron transport chain, mitochondrial ATP synthetase, and phosphate and adenine nucleotide transporters. Assembling models of these discrete processes into an integrated model of mitochondrial ATP synthesis, we illustrate how to analyze and simulate in vitro respirometry experiments and how models identified from in vitro experimental data effectively explain cardiac respiratory control in vivo. Computer codes for these analyses are embedded as Python scripts in a Jupyter Book to facilitate easy adoption and modification of the concepts developed here. This accessible framework may also prove useful in supporting educational applications. All source codes are available on at <a href="https://beards-lab.github.io/QAMAS_book/">https://beards-lab.github.io/QAMAS_book/</a>. 
# 
# ## Highlights 
# 
# -   A kinetic and thermodynamic framework for mitochondrial energetics is developed.
# -   The framework is applied to simulate ATP synthesis and respiratory control.
# -   We illustrate how respiratory control in vitro translates to energetics in vivo.
# -   Computer codes are available at DOI: 10.5281/zenodo.4915567. 
# 
# 
# ## Funding 
# 
# This work supported by NIH grant HL144657.
# 

# In[ ]:





# 
# ```{toctree}
# :hidden:
# :titlesonly:
# 
# 
# Abbreviations
# Introduction
# Principles
# BuildingModel
# InVitroModel
# InVivoModel
# Summary
# References
# ```
# 
