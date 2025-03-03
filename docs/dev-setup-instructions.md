<!-- 
Copyright 2025, Drift+Noise GmbH

This file is part of PRIIMA.
PRIIMA is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.
PRIIMA is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with
PRIIMA. If not, see https://www.gnu.org/licenses/gpl-3.0.html.
-->

# Local development environment setup instructions

## Basic setup

Clone the repository to your local computer and enter the project directory.

```
$ git clone <repository>
$ cd priima
```

Install the basic Linux packages:

```
$ sudo apt install gdal-bin libgdal-dev python3-dev build-essential virtualenv \
    nco cdo cmake libopencv-dev libboost-all-dev libfreetype-dev docker-compose
```

Copy landmask shapefiles for the Arctic and Antarctic into the project to ensure
that PRIIMA has correct landmask information to avoid that land parts of Sentinel-1
images are transformed but only sea-ice.

Then run the `install-deps` script to set up the project and install the
required Python dependencies:

```
$ ./devops/bin/install-deps
```

Activate virtual environment to use locally installed Python libraries.

```
$ source venv/bin/activate
````

## Fast Cast

You can compile the Fast Cast Library for improved performance by the following steps.

Using the Fast Cast module requires a submodule called `pyboostcvconverter` which communicates between opencv (C++) and Python.
This communication is necessary because the FastCast code is written in C++. To get this library run the following command:

```
$ git submodule update --init --recursive
```

And finally build the Fast Cast Library:

```
$ make fast_cast
```

Check the Fast Cast installation by running the tests. The `fast_cast`
components should not be skipped - this would be indicated by an orange "s"!

```
$ make test
```

## Troubleshooting

### OpenCV library not available as OS package

In some cases it might not be possible to install the OpenCV library via the
operating system packages.  In this case, one can install it by following
the [instructions to install it from
source](https://docs.opencv.org/3.4/d2/de6/tutorial_py_setup_in_ubuntu.html).

### OpenCV Python wrapper not working

In case that the OpenCV Python wrapper is not working, try to install the
`opencv-python-headless` package into the virtual environment:

```
$ pip install opencv-python-headless
```
