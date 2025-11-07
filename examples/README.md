# Run Tlab

Workflow summary. Details of each tool and corresponding input data to be read in the corresponding source file.

## Preprocessing Tools

| binary    | inputs                           | outputs |
| --------- | -------------------------------- | ------- |
|inigrid.x  | tlab.ini                         | grid |
|inirand.x  | tlab.ini, grid                   | [flow,scal].rand.? |
|iniflow.x  | tlab.ini, grid [,flow.rand.?]    | flow.ics.?
|iniscal.x  | tlab.ini, grid [,scal.rand.?]    | scal.ics.?

## Simulation Tools

| binary    | inputs                            | outputs |
| --------- | --------------------------------- | ------- |
|dns.x      |tlab.ini, grid, flow.*.?, scal.*.? |   flow.*.?, scal.*.? |

## Postprocessing Tools

| binary    | inputs                            | outputs |
| --------- | --------------------------------- | ------- |
|visuals.x  | tlab.ini, grid, flow.*.?, scal.*.?| *variable files* |
|averages.x | tlab.ini, grid, flow.*.?, scal.*.?| avg*
|pdfs.x     | tlab.ini, grid, flow.*.?, scal.*.?| pdf*
|spectra.x  | tlab.ini, grid, flow.*.?, scal.*.?| xsp*, zsp*

## List of examples

### Boundary-free flows: mixing layers and gravity waves

* Case01. Shear layer with broadband ICs. Uniform grid. Kelvin-Helmholtz.  
* Case02. Same as Case01, but with stretched grid.  
* Case03. Same as Case02, but 2 scalars with different Schmidt numbers instead of one.

* Case06. Stably stratified density interface with discrete ICs: oscillating inversion.  
* Case07. Unstable stratified density interface: Rayleigh-Taylor.  
* Case08. Stably stratified shear layer.  

* Case11. Wave maker in incompressible case.

### Wall-bounded flows: channel flows and Rayleigh-Benard convection

* Case21. Rayleigh-Benard convection with Dirichlet boundary conditions.

* Case25. Channel flow.  
* Case26. Stably stratified channel flow.
* Case27. Half Channel flow.  
* Case28. Same as 27, implicit solver. TBD.
* Case29. Rotating channel flow (rotation term for turbulent transition).  

### Boundary layers

* Case31. Heated plate.  
* Case32. Convective boundary layer.

* Case41. 1D perturbed laminar Ekman layer.
* Case42. 1D perturbed laminar Ekman layer, implicit solver. TBD.
* Case43. 3D Neutral Ekman layer with pasive scalar.
* Case44. 3D Stable Ekman layer.
* Case45. 3D Neutral Ekman layer with interactive BC at the bottom.  

### Anelastic flows

* Case51. Air anelastic formulation of CBL.  
* Case52. Airvapor anelastic formulation of CBL.  

* Case61. Airwater equilibrium anelastic formulation of stratocumulus-topped boundary layer.
* Case62. Same as Case61, but adding subsidence and sedimentation.  
* Case63. Same as Case62, but with dimensions.  
* Case64. Same as Case63, but with gray radiation model.  
* Case65. Same as Case63, but with band radiation model.  

* Case71. Airwater non-equilibrium (1-moment bulk scheme) anelastic formulation of stratocumulus-topped boundary layer.
          Same as Case61, but with explicit non-equilibrium formulation for mono-disperse droplets.
* Case72. Same as Case62, but with semi-implicit non-equilibrium formulation for log-normal DSD droplets.
* Case73. Same as Case62, but with fully implicit non-equilibrium formulation for log-normal DSD droplets.

<!-- make checkrl/checkdb runs the check.sh bash-script inside each directory -->
