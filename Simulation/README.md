This repository will be for the full beamline simulation. The process is as follows:
1. Opal-T for gun to Yag4, because nothing in that section changes; optimized to minimize $\epsilon_x$, allowing the beam to pass through the small structure, while trying to reduce $s_{rms}$, creating a higher peak current. Increasing peak current (desired for creating wakefields) has the potential to blow up emittance, hence the inclusion of both. 
2. Elegant simulation from Yag4 to the structure to optimize quadrupoles to create beam waist in center of structure
3. Take optimized quadrupole information and put into Opal-T for Yag4 to start of structure
4. WarpX to propagate beam (uploaded from Opal-T) through the structure
5. Elegant to optimize quadrupoles from end of structure to create beam waist at final Yag screen 
6. Opal-T from end of structure to final Yag screen using optimized quadrupole settings

The final simulation used in my dissertation was completed on the NIU cluster Metis. Inside each directory is the list of files and analysis for the single-beam case with another directory of the two-beam case. Please see directory for more information. 

## Opal-T
Currently, I do not convert an external beam description to an Opal-T input file. Instead, I set the beamline conditions and run an optimizer for Opal-T to achieve specific results at the end of the simulation. The optimization was created in Python my Dr. Philippe Piot. It allows the user to specify certain ranges for the gun injection amplitude and phase, the focusing, main, and bucking solenoid currents, the radius of the laser, and the amplitude and phase of specific linacs. The kinetic energy, number of particles, emittance, rms size, and rms_ps are recorded and the specific fitness parameters desired are taken into account. 

The Output units are:
[$\beta_{x,y,z} \gamma$] 

## Elegant
Input is SDDS file, but you also need the numerical value of the Central momentum of the beamline [MeV/c]. For individual particles' information, x, y, z are in [m], compared to where the reference particle are and the momenta are in $\beta \gamma$. These are the same units as used in Opal-T. 
Input SDDS for "elegant" requires (x [m], xp [dimensionless], y [m], yp [dimensionless], t [s], p [$\beta \gamma$])




## WarpX
Another transfer script was created to convert the Elegant output to a WarpX input file. Output units are: [$\gamma v$] for momenta. This file uses a specific function to create the various glass tubing structures, as seen in its epsilon and sigma functions. We start the initialization plane 10 mm behind where the structure starts (**Note: this is also where we take the particles from in the Elegant simulation**). This is so that the particles are able to create their fields before entering the structure. The structure lengths were all set to 90 mm with the wall thickness set to the glass dimensions and the conductor thickness set to be at least a couple mesh points. The glass relative dielectric constant was 4.26 and the conductance of the coating is 5.8e7. More information regarding the pixel size, boundaries, etc is in the ReadMe of WarpX directory.









