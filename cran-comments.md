# This is an update of the NPflow package  

## Test environments  
 * local R installation, R 4.4.3 on macOS 15.5
 * Linux (Ubuntu 24.04), macOS (14.7) and Windows (Server 2022 10.0), R devel and release (through GitHub Actions)
 * Rhub
 * Win-builder

## R CMD check results  
0 error | 0 warning | 0 note

On some architectures, the large volume of compiled file may trigger a note.
These compiled functions are necessary to ensure reasonable computation time
for NPflow.
I updated the Suggested packages as MASS will not be imported by ggplot2 anymore

## Reverse dependencies 
This update does not introduce new notes, warnings or errors in the reverse 
dependent package bqror


Thanks, Boris Hejblum

---
