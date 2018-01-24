Tutorial
=========

The steps required to quickly and efficiently generate source files that are used to simulate a GRB on a given potential BurstCube arrangement with a given strength and sky position are as follows:

1.	Creating source files

config.yaml allows you to select the desired geometry, desired energy range, and desired source position respectively. 

Using python, this yaml file is read in using the configurator, which is contained in the  simGenerator module. 

The simGenerator module, along with the rest of the modules required for this process are located within the bin. If it is the case where these functions need to be done in a new or different directory, make sure to export them.
 
Next, create an object with the configurator by pointing to the yaml file like follows:

conf = configurator(‘~/config.yaml’)

and now with this object, run the function conf.createSourceFiles(), which outputs all of the source files into the desired directory. 

2.	Running simulations 

runSims.py is a function which runs the source files and outputs the .sim files to be further analyzed. 

In order to run runSim.py, simply run 

python <location>/runsims.py –runCosima TRUE –sourcepath ./ x 
where x is the desired amount of source files to run at once. 


After runSims has completed its operation, the sim files will be outputted in the designated directory for further analysis. 

3.	Analyzing results

With the sim files generated, these can be plotted using python. 

Once uploading the sim files to python, they can be plotted using the module plotSim, which provides many options as to how the various simulated detector responses can be analyzed. 

