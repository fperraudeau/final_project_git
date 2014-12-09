stats215a
=========

Shared repository for the Berkeley's STAT215A Final Project

#####Author : Fanny Perraudeau

Please follow the instructions below to reproduce my project. 
I assume the user is comfortable using the command line and github. 
In case of doubts, please email me fperraudeau@berkeley.edu

#####BASIC REQUIREMENTS

    R, git, latex

#####INSTRUCTIONS FOR UNIX SYSTEMS:

Make sure that all the following files are in the same folder. You should at least have:
1. 
2

You will also have tables in .csv format that I have generated on the SCF cluster. 
You may want to backup these tables in a different folder so they won't be overwritten.

Open initialize.R and change "work.dir" to your working directory where all the files are located.
Also, make sure that all the libraries listed in initialize.R are installed.

Don't forget to put fMRIdata.RData in the same folder as well.


Follow the steps to reproduce the lab
1. Open the command prompt on a SCF computer
2. Create folder to store this repository (mkdir fannyFinalProject)
3. Cd into the folder, cd ~/fannyFinalProject
4. Clone this repository ( git clone https://github.com/fannyperraudeau/blablabla )
5. Cd into the repository folder ( cd stats215a )
6. Change permissions for the script that compiles the pdf
7. Run with the command line qsub -pe smp 6 shell_example.sh
7. Run the script to compile the pdf
9. Open the pdf ( .pdf )