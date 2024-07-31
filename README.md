# RamanMatch
Application for Raman Spectroscopy Analysis

A manuscript detailing the project is submitted for publication and may be cited as below: 

Berrada, M., McFall, A., and Chen, B. (Submitted in October 2023) Raman Match: Application for Raman Spectroscopy Analysis. 
Journal of Open Research Software.

Abstract

Raman spectroscopy is a powerful analytical technique used in various scientific disciplines, including chemistry, materials science, and biology. Analyzing Raman spectra and identifying the presence of specific substances within a sample can be a complex task due to the large available database. The Raman Match application was developed to simplify the sample identification process. The application integrates the well-established RRUFF Raman database with the python programming language and provides a user-friendly graphical interface for loading Raman spectra, identifying peaks, matching peaks to reference libraries, visualizing the results, and produce publication ready figures. 





# Description of Repository content: 

Window EXE: Contains a link to the executable .exe file for the application to run on Windows. 

Databases: For users who wish to work on the code, the Jupyter Notebook for the creation of the two databases is available in the repository. The files may be easily converted to .py.
The databases are also available here: https://drive.google.com/file/d/1ZoIzt2MXQP5-1vqUuQQsuP0d3J7I2G1e/view?usp=sharing and https://drive.google.com/file/d/1A6PjJOPmuCa-wqkDWG6O5zD4smgfvkV8/view?usp=sharing 

Main.py: The main python code of the Raman Match application. 

Forsterite.txt: This is an example file you may use of the Raman spectrum of a non-oriented Forsterite crystal. 
Basalt.txt: This is an example file you may use of the Raman spectrum of a rock sample, presumaly basalt, collected in Oahu.

icon.ico: The icon for the application is provided in case the user wishes to compile a full version on their machine. 

For MacOS/Linux system users, the database codes should be ran once in your machine so the .db files are created (a link to the .db files is also available above). Then you may simply run the python code in your terminal to use the software, without having to compile it. You may also use it through Jupyer Notebook. 

# INSTRUCTIONS
Windows: Download .exe file or run the python code (GUI) on your preferred platform.

MacOS/Linus: Download the python code (GUI) on your preferred platform. Run the database codes first. 

For any questions of inquiries, please contact Meryem Berrada at berrada@hawaii.edu

# UPDATES
July 31, 2024: Databases updated to include all of the RRUFF files. Previously, the code only gathered the "excellent" and "fair" files. 
