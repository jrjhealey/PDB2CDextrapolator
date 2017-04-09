#!/usr/bin/env python

"""
This script takes the simulated output spectrum from PDB2CD.
It extrapolates further wavelengths beyond the returned 240nm
to allow the spectra to be compatible with Dichroweb.
"""


# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.


import traceback
import sys

__author__ = "Joe R. J. Healey"
# Contributions to Taylor series calculation from Lewis Baker

__version__ = "1.0.0"
__title__ = "PDB2CDextraploate"
__license__ = "GPLv3"
__author_email__ = "J.R.J.Healey@warwick.ac.uk"


def parseArgs():
    """Parse command line arguments"""

    import argparse

    try:
        parser = argparse.ArgumentParser(
            description='This script takes in the spectral output of PDB2CD (http://pdb2cd.cryst.bbk.ac.uk/) and extrapolates the spectra to allow compatibility with Dichroweb')

        parser.add_argument('-s',
                            '--spectrum',
                            action='store',
                            default='None',
                            help='The spectra output by PDB2CD (plain text, 2 column tab separated after header (Wavelength | CD).')
        parser.add_argument('-b',
                            '--bibliography',
                            action='store_true',
                            help='Print out relevant references.')
        parser.add_argument('-f',
                            '--forecast',
                            action='store',
                            type=int,
                            default=20,
                            help='The length of the data (x axis) to extrapolate the gradient until [Default = 20 nm].')
        parser.add_argument('-o',
                            '--outfile',
                            action='store',
                            default='None',
                            help='File to store the new spectrum in. Prints to screen by default.')
        parser.add_argument('-v',
                            '--verbose',
                            action='count',
                            default=0,
                            help='Verbose behaviour, printing parameters of the script.')
        parser.add_argument('-p',
                            '--plot',
                            action='store_true',
                            help='Output a plot to examine the new spectrum (requires Matplotlib installed).')
        parser.add_argument('-d',
                            '--delta',
                            action='store',
                            type=int,
                            default=3,
                            help='The step distance to calculate the gradient over. [Default & minimum = 3]')

    except:
        print "An exception occurred with argument parsing. Check your provided options."
        traceback.print_exc()

    return parser.parse_args()


def displayRefs():
    """Display relevant references"""
    print('''
    - PDB2CD has been previously published here:
        Mavridis and Janes, Bioinformatics, (2017), 33(1): 56-63

    - Dichroweb has been previously published here:
        Whitmore, L. and Wallace, B.A. (2008) Biopolymers 89: 392-400.
        Whitmore, L. and Wallace, B.A. (2004) Nucleic Acids Research 32: W668-673.
    ''')
    sys.exit(1)


def getSpectraValues(spectrum):
    """Extract arrays of values from the input data file"""
    import re

    with open(spectrum, 'r') as ifh:
        Xvals = []
        Yvals = []
        # Skip extra info at top of file
        header = []
        noheader = []

        for linenumber, line in enumerate(ifh):
            if linenumber < 12:
                header.append(line)
            elif linenumber > 12:
                noheader.append(line)

        for line in noheader:
            # Ensure we're dealing with spectral data with wavelength in column 1
            if re.match(r"^\d\d\d", line):
                x, y = line.split('\t', 1)
                Xvals.append(float(x))
                Yvals.append(float(y))
            else:
                print("Could not find spectral values. Ensure the input spectra is formed correctly.")
                sys.exit(1)

    return Xvals, Yvals, header

def arrayAvg(array):
    return sum(array)/float(len(array))

def arrayDiff(array):
    return [b - a for a, b in zip(array[:-1], array[1:])]


def forecastSpectra(Xvals,Yvals, forecast, delta, verbose):
    """Compute the gradient (average difference) of the last delta values of each X and Y array"""
    import math

    for a in range(forecast):
        Xsub = Xvals[-delta:]
        Ysub = Yvals[-delta:]
        if verbose > 1:
            print("Extrapolating from Xvalues:")
            print(Xsub)
            print("Extrapoliting from Yvalues:")
            print(Ysub)

        # Compute first x and y differences as arrays
        x1diffs = arrayDiff(Xsub)
        y1diffs = arrayDiff(Ysub)
        x1avg = arrayAvg(arrayDiff(Xsub))

        # Average first order difference
        d1y = arrayAvg(arrayDiff(Ysub))/arrayAvg(arrayDiff(Xsub))

        if verbose > 1: print("First order Taylor term = " + str(d1y))

        # Make array of adjacent gradients
        diffarray1 = []
        for x, y in zip(x1diffs, y1diffs):
            diffarray1.append(y/x)

        # Average second order difference
        diffarray2 = arrayDiff(diffarray1)
        print(diffarray2)
        if verbose > 1: print("Second order Taylor term = " + str(d2y))

        d2y = arrayAvg(diffarray2)/x1avg

        if delta >= 5:
        # Third order terms if delta is at least 5 else not enough data
        # and zero divison error will occur.
            diffarray3 = arrayDiff(diffarray2)
            print(diffarray3)
            d3y = arrayAvg(diffarray3)/x1avg
            if verbose > 1: print("Third order Taylor term = " + str(d3y))
        else:
            d3y = 0

        # Simple to extend in the case of x as this is just linear
        Xvals.append(Xvals[-1] + x1avg)

        # For y, a finite difference partial Taylor expansion approximates the curve
        # Maximum number of terms = number of points fitted..

        Yvals.append(Yvals[-1] +
                     Yvals[-1] * d1y +
                     (((Yvals[-1]**2) * d2y) * (1/math.factorial(2))) +
                     (((Yvals[-1]**3) * d3y) * (1/math.factorial(3)))
                     )

    # y (first term
    # y x First derivative (second term)
    # 0.5*y^2 x Second derivative (third term)
    return Xvals, Yvals


def plot(Xvals, Yvals, verbose):
    """"Create a quick plot via matplotlib to visualise the extended spectrum"""

    if verbose > 0: print("Plotting data")

    import matplotlib.pyplot as plt
    plt.plot(Xvals, Yvals)
    plt.show()


def main():
    """Extrapolate from PDB2CD spectra"""

    # Parse arguments
    args = parseArgs()

    # Convert object elements to standard variables for functions
    spectrum = args.spectrum
    forecast = args.forecast
    outfile = args.outfile
    verbose = args.verbose
    makeplot = args.plot
    delta = args.delta

    if args.bibliography is True:
        displayRefs()
        sys.exit(0)

    if spectrum == 'None':
        print('You have not specified a spectrum file (use the -s or --spectrum flags)')
        sys.exit(1)


    # Recover list of each variable:
    Xvals, Yvals, header = getSpectraValues(spectrum)

    # Project new data to output
    Xvals, Yvals = forecastSpectra(Xvals, Yvals, forecast, delta, verbose)

    if outfile is not 'None':
        with open(outfile, 'w') as ofh:
            ofh.write(header)
            for c1, c2 in zip(Xvals, Yvals):
                ofh.write(str(c1) + '\t' + str(round(c2,3)))
    else:
        for line in header:
            print(line.rstrip("\n"))
        for c1, c2 in zip(Xvals, Yvals):
            print(str(c1) + '\t' + str(round(c2,3)))


    if makeplot is True:
        plot(Xvals,Yvals, verbose)

    if verbose > 0: print("All done.")


if __name__ == '__main__':
    main()
