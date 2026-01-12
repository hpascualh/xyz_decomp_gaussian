XYZ decomposition of IR or Raman in Gaussian
- we need to display the polarizabilities in the freq job
"#P PW91PW91/genecp nosymm freq=(raman,savenormalmodes) iop(7/33=1)"
- important to have the #P - it displays more info; and the iop keyword. 
- important to not have any atom fixed (-1) in the system. Otherwise, the python code will crash

After running the freq job, copy the gauss_vib_xyz_raman_hacky_v2.py file in the same folder as the out.log file. 
Open a python terminal in that directory and run: 
> python gauss_vib_xyz_raman_hacky_v2.py


If it runs successfully, you will find the following new files in your folder: 
- gaus_ir, gaus_irx, gaus_iry, gaus_irz, gaus_raman
These will contain the IR activities of each component. To generate the spectra-looking data, open each gaus_ir_i in GabeEdit (read from an ASCII XY file) and save the data from there. 


