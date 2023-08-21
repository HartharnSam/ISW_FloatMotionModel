# ISW Float Motion Model
Code to implement Float Motion Model described in _"Interactions between Internal Solitary Waves and Sea Ice"_
Takes an input of horizontal (1-D) velocity through time (x, t) and particle properties and outputs the path of the particle. 


## How to use:

### Float Motion Model
- ```FloatMotionModel.m```  : Function code computing particle tracks, using either the basic model (as described in the paper), or a drag-based "advanced" model.
- ```run_lab_FMM```: Parsing script for Lab Data into the FMM
- ```run_DJL_FMM```: Parsing script for DJL data into the FMM
- ```run_model_FMM``` : Parsing script for SPINS model data into the FMM

### Plotting and explaining scripts
- ```DJLPlotsAnimation.m``` : Produces an animation of the DJL Speed Plots function
- ```DJL_HistAnimation``` :   Produces an animated Histogram representing the model for DJL
- ```DJL_SpeedvsFloatLength``` : Plots calculated float motion according to DJL waves and the FMM model
- ```DJL_speeds_plot``` : Plots the fluid speeds at point A and point B for a float moving with the fluid according to the FMM
- ```LabFMM_speeds_plot``` : Plots the fluid speeds at point A and point B for lab measured PIV data
- ```plot_FMM``` :  [Early version] of a general plotting script, left in for potential utility 




## Other:

Some of the dependencies in this toolbox require tools contained within the ISWLabToolkit repository, at [GitHub](https://github.com/HartharnSam/ISWLabToolkit)

## Development:
The FMM and basic parsing scripts work, and have been used in work due to be published. However the various animation and plotting scripts, provided in the hope they will prove useful, but with no guarantee to this utility, their documentation or functionality. There are no current plans to document or improve those elements further, unless I have time

## Installation:
This can be installed from github at : https://github.com/HartharnSam/ISW_FloatMotionModel

## See Also:
- ISWLabToolkit : In particular ISWLabToolkit/DigiFlow_read/particletracking
- [DJLES](https://github.com/mdunphy/DJLES/)
