### feMorph Application

### Installation

First, clone the git repository in a local folder:

```
git clone git@marsdenfileserver1.stanford.edu:schiopaz/femorph-application.git
```

Create a directory for out-of-source build (e.g., feBin):

```
mkdir feBin
```

Enter this folder 

```
cd feBin
```

run cmake with the following command:

```
cmake ../femorph-application/
```

Run make:

```
make (or make -jn for parallel compilation with n processes)
```

#### Installing the documentation

The smpFilter documentation is written using Sphinx. To compile the documentation go to the docs folder:

```
cd docs
```

and run the make utility and specify the html target

```
make html
```

An **index.html** document will be created in the **docs/build/html** folder.

#### Running the code

To run the Poisson Pressure Equation solver and export the results in VTK format, use the following command:

```
feMorph -z -f input.dat -o output.vtk
```

