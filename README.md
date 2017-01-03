#DataProcSandbox

***
##Summary
This package contains a collection C++ classes and scripts designed for processing and cleaning datasets. It is designed to read input ROOT files. There are tools for producing ROOT plots and threading jobs to batch clusters. Furthermore, there are analysis wrappers for processing Delphes output files.

####Structure
The source executables are in the "src" directory. Headers are in the DataProcSandbox directory. Output plots and tables are saved to the "output" directory. Running make will compile executables in the "run" directory. Objects are compiled in the "bin" directory. The Makefile should be called along with the name of the executable that you want to build. Here is an example:
     # Compile the source code in src/DelphesTtbarReconstruction ...
     make DelphesTtbarReconstruction
     ./run/DelphesTtbarReconstruction

Source the "setup.sh" script each time you log in before running in order to ensure that your environment is set up for using ROOT and Delphes libraries.

***
##Git Refresher (for Terminal Work)

####Clone the repository on your local machine:
    git clone https://github.com/jwebste2/DataProcSandbox.git DataProcSandbox

####Basic Pulling:
    git pull

####Basic Pushing:
    # To include changes to files that are already tracked
    git add -u

    # To add new files that are not yet tracked
    git add filenameA filenameB ...

    # Then commit the changes    
    git commit -m "Add a log message"

    # Then push the commit
    git push

####Merging your branch into master
    # Pull updates from master to your branch
    git pull

    # Commit and push your branch updates
    git commit -m "log message"
    git push

    # Switch to the master branch and then merge in your changes from your_branch
    git checkout master
    git pull
    git merge your_branch

    # Go back to editing your branch if you want
    git checkout your_branch
    git rebase master


####Checking Logs:
    # Format ==> git log -NEntriesToPrint, e.g.
    git log -5
