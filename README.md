### feMorph Application

### How to Install

#### Prerequisites

To successfully compile the feMoprh application, the following software is required:

- Boost Libraries (release 1.40 or later).
- Sphinx (this documentation).
- One implementation of the MPI libraries (optional).

One of the following linear algebra packages:

- CSparse.
- Trilinos (the MPI library needs to be also installed).
- Armadillo.

First, clone the git repository in a local folder. Then create a directory for out-of-source build (e.g., bin)

```
mkdir bin
```

Enter this folder

```
cd bin
```

run ccmake with the following command (you can alternative use the cmake GUI): 

```
ccmake ../feMorph/
```

and set the following options:

- **useMPI**. Use the MPI libraries.
- **useCSparse**. Use the CSparse library.
- **useTrilinos**. Use the Trilinos library.
- **useArmadillo**. Use the Armadillo library.

Run make 

```
make (or make -jn for parallel compilation with n processes)
```

#### Installing the documentation

The smpFilter documentation is written using Sphinx. To compile the documentation go to the docs folder

```
cd docs
```

and run the make utility and specify the html target

```
make html
```

An **index.html** document will be created in the **docs/build/html** folder.

#### Running the code

To run the Poisson Pressure Equation solver and export the results in VTK format, use the following command

```
feMorph -z -f input.dat -o output.vtk
```
