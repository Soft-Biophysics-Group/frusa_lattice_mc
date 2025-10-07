---
title: Testing frusa_mc
author: Vincent Ouazan-Reboul
date: 2025/04/07
---

## Abstract

`frusa_mc` is giving odd results for the cubic self-limiting aggregates ("camemberts") I was
trying to make.

I am now not 100% trustful in the code I have written: something could easily be wrong on the
Python side of things. To make sure all is well, I am writing some simple tests that should
diagnose all the issues.

## Contents

- `00_plot_all_cubes`: Make sure the orientations in the python code and the Blender plotting
  match.
