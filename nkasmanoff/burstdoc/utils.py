#!/usr/bin/env python

import os


def setPath(desiredPath = "", desiredPathName = 'BURSTCUBE'):

    '''Checks for desiredPathName path.  Returns 0 if ok, 1 if bad.'''

    if desiredPath:
        if '~' in desiredPath:
            os.environ['desiredPathName'] = os.path.expanduser(desiredPath)
            return 0
        else:
            os.environ['desiredPathName'] = desiredPath
            return 0
    elif not (desiredPathName in os.environ):
        print('Set or provide ' + desiredPathName)
        return 1
    else:
        return 0


def getFilenameFromDetails(details):

    '''Takes a dictionary of details and makes a machine readible filename
    out of it.'''

    filename = "{}_{:.3f}keV_Cos{:.3f}".format(details['base'],
                                               details['keV'],
                                               details['Cos'])

    return filename


def getDetailsFromFilename(filename):

    '''Function to get the energy and angle from a filename.
    Really should be meta data.'''

    details = {}
    info = filename.split('_')

    details['base'] = info[0]

    details['keV'] = float(info[1][:-3])

    angle = info[2].split('.')
    details['Cos'] = float("{}.{}".format(angle[0][3:], angle[1]))

    return details

