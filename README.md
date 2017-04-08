# PDB2CDextrapolator
A (relatively) simple extrapolation script for PDB2CD web server data, to allow it to be compatible with Dichroweb and other Circular Dichroism tools.

**PDB2CD** is a webserver (http://pdb2cd.cryst.bbk.ac.uk/) which predicts/simulates Circular Dichroism spectral data, however it only does this up to 240 nm. This makes it incompatible with Dichroweb (http://dichroweb.cryst.bbk.ac.uk/html/home.shtml), the webserver which analyses the same kind of data, from experiments. Many of the reference sets require 260 nm worth of data though, even if they dont use it within the analysis.

# PDB2CDextrapolate.py
This python script takes in the plain text output of PDB2CD spectra, handles all the file manipulation and returns some extended spectral data from simple curve fitting via a Finite Difference Taylor Series. PDB2CD data looks like this:


    #################################
    # Created By PDB2CD             #
    # 1/01/70 00:00:00              #
    # Version 1.0                   #
    # Author: Dr. Lazaros Mavridis  #
    # E-mail: l.mavridis@qmul.ac.uk #
    #################################

    Data for the run: 1abc
    ----------

    Wavelength	CD
    175	-2.065
    176	-1.974
    177	-1.929
    ...
    238	-0.203
    239	-0.155
    240	-0.123
    
# Installation
None required, just clone this repository and run the script. The exception is if you wish to use the simple `matplotlib` functionality embedded in one of the methods to visualise the extrapolation (`-p|--plot`) - if so you'll need to make sure you've installed the module before.

    $ python -m pip install matplotlib
    
is usually sufficient.

# Execution
There are a few different funtionalities build in to the script. The script takes the following options:

`-s | --spectrum` - The input file from PDB2CD (plain text format).
`-b | --bibliography` - Display some salient references and exit (you cannot use this option with other options).
`-f | --forecast` - The number of periods (x axis increments) to extend the data by. Default of 20 nm in 1 nm increments.
`-o | --outfile` - Have the script write a separate new file for the extended spectrum. Prints to screen by default.
`-v | --verbose` - Output some useful info about variable values and progress messages. Supports 2 levels of verbosity at present (-v,-vv).
`-p | --plot` - Invoke `matplotlib` at the end of the script to visualise the extended spectrum. Requires `matplotlib` to be preinstalled - see above.
`-d | --delta` - The span of data from which gradients are calculated, this is a default of 3 x-increments. **IMPORTANT - see the note below about usage**.

A basic script invocation would be:

    $ python PDB2CDextrapolate.py -s 1abc.txt
    
To view the plot and output a new file:

    $ python PDB2CDextrapolate.py -s 1abc.txt -o ~/path/to/outputfile.txt -p


# The Delta parameter and the Taylor expansion:
Delta, the interval along X to perform the finite difference calculations on has to have a minimum value of 3 so that the average gradient and difference calculations can work. The script by default will use 3 data points, and calculate the first 2 derivative terms of the Taylor series (the first term is just itself):

http://mathworld.wolfram.com/images/equations/TaylorSeries/NumberedEquation1.gif

Therefore order of the series used is always *delta* - 1 .

This is because each successive derivative calculation, relies on the differences of the previous, which reduces the dimensions of the array you calculate from by 1 each time (eventually you'd have an empty array, and the average calculations which divide by the `len(array)`, would give a divide-by-zero error.

The most the script will use is the 3rd derivative, which requires at least 5 data points. If the default is used, the third derivative is simply set at 0.

# Disclaimers:
The script is not fully tested, and I'm a moderately proficient user of python, but a pretty useless user of calculus. Some stuff may be broken/untested and the calculations may not be perfect (and may straight up not work for some weird CD spectra - I can't test against all scenarios). That said, I'd love to hear suggestions for improvements etc.
