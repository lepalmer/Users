#!/usr/bin/env python

from astropy.tests.helper import pytest

try:
    from simGenerator import configurator
except ImportError:
    pass

try:
    from simGenerator import createSourceString
except ImportError:
    pass


@pytest.fixture(scope='module')
def create_configurator(request):

    from os import path

    testdir = path.expandvars('$BURSTCUBE/Simulation/MEGAlib/test/')
    conf = configurator(testdir+'config.yaml')

    return conf


def test_configurator_setup(create_configurator):

    conf = create_configurator
    print(conf.config['run']['basename'])


def test_createSourceString(create_configurator):

    refstr = 'Version 1\nGeometry '
    refstr += '$BURSTCUBE/Simulation/MEGAlib/test/'
    refstr += 'BurstCube_1Cylinder.geo.setup\n'
    refstr += 'CheckForOverlaps 1000 0.01\nPhysicsListEM Livermore\n'
    refstr += 'StoreCalibrate True\nStoreSimulationInfo True\n'
    refstr += 'StoreOnlyEventsWithEnergyLoss True\nDiscretizeHits True\n\n'
    refstr += 'Run FFPS\nFFPS.Filename $BURSTCUBE/Simulation/MEGAlib/test/'
    refstr += 'test_100.000keV_Cos0.100\n'
    refstr += 'FFPS.NTriggers 1000\nFFPS.Source One\nOne.ParticleType 1\n'
    refstr += 'One.Beam FarfieldPointSource 5.73 0\n'
    refstr += 'One.Spectrum Mono 100.0\nOne.Flux 1000.0\n'

    conf = create_configurator
    sourcestr = createSourceString(conf.config, 100., 0.1)

    assert(sourcestr == refstr)


def test_createSourceFiles(create_configurator, tmpdir_factory):

    from os.path import isfile

    srcdir = tmpdir_factory.mktemp('source')
    files = ('test_100.000keV_Cos0.500.source',
             'test_100.000keV_Cos0.750.source',
             'test_100.000keV_Cos1.000.source',
             'test_173.205keV_Cos0.500.source',
             'test_173.205keV_Cos0.750.source',
             'test_173.205keV_Cos1.000.source',
             'test_300.000keV_Cos0.500.source',
             'test_300.000keV_Cos0.750.source',
             'test_300.000keV_Cos1.000.source')

    conf = create_configurator
    conf.createSourceFiles(dir=srcdir.__str__())

    for f in files:
        srcfile = srcdir.__str__() + '/' + f
        assert(isfile(srcfile))
