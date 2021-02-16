#!/usr/bin/env python
## a class to process single pulse ephemeris data and single pulse file path.
## -- Time-stamp: <2020-05-12T10:23:27.577136-04:00 hedfp> --
##
## based on getarrtime from Jodrell Bank.
## 2018-11 /CS/ Carl Schmiedekamp

from __future__ import absolute_import, division
from __future__ import print_function, unicode_literals

import sys

## redefine input() for Python 2
if sys.version_info[0] < 3:
    def input( str):
        return raw_input( str)

## Q: should this class read the file or
##    should that be done in separate function??
## A: yes. If all ephemeris data is stored in class in later version,
##    Then *.spe file only needs the epoch of the observation and
##    The single pulse text file data.
##
##  *SPE.txt file format version 1:
##    <file format version number>  <Epoch of observation>
##    <line from JB ephemeris closest to epoch of observation>
##    Following lines are 'single pulse text file' lines.
##    That is, two lines are added at the front of the usual single
##      pulse text file.

    
class spephemeris:

    ## construct instance based on SPE file.
    @staticmethod
    def readspe( spe_file):
        with open(spe_file) as f:
            line1 = f.readline().split()
            if( len( line1) < 2):  ## after item 2 is comment
                print( 'SPE file format error.')
            else:
                fileversion = int( line1[0])
                obsepoch = float( line1[1])
#                print( 'DBug: fileversion, obsepoch:', fileversion, obsepoch)
                assert fileversion == 1

            SPFile = spe_file
            
            line3 = f.readline().split()
#            print( 'DBug: line3:', line3)
            
            MJDd = float( line3[3])
            ##MJDs = float( line3[4])  ### Use MIT value
            MJDs = float( line3[5])  ### Use JPL value
#            MJD = MJDd + MJDs/86400
#            print( 'DBug: MJDd, MJDs, MJD:', MJDd, MJDs, MJD)
    
            nu = float( line3[6])
            nudot = float( line3[8])
#            print( 'DBug: nu, nudot:', nu, nudot)
            
            f.close()

            return spephemeris( obsepoch, MJDd, MJDs, nu, nudot,
                                SPFile, fileversion)
    

    def __init__( self, obsepoch, mjdd, mjds, nu, nudot,
                  SPFile=None, fileversion=1):
        """ Set up power series coefficients using Jodrell Bank ephemeris
        epoch: MJD (with fractional days) of observation
        mjdd: modified julian date (integer part), from ephemeris (days)
        mjds: seconds after mjd, from ephemeris (s)
        nu: pulsar rotation frequency (Hz)
        nudot: rate of change of nu times (in 1e-15 Hz/s)
        """


        self.obsepoch = obsepoch
        self.mjdd = mjdd
        self.mjds = mjds
        self.epoch = mjdd + mjds/86400
        self.nu = nu
        self.nudot = nudot = nudot*1e-15
        self.p = p = 1.0/nu                  ## period in seconds
        self.pdot = pdot = -nudot/(nu*nu)    ## period derivative in s/s (dimensionless)
        self.nuddot = 2.0*pdot*pdot/(p*p*p)  ## second derivative of frequency
        self.SPFile = SPFile                 ## single pulse text file
        self.fileversion = fileversion       ## format version for SP file


        
    def period( self, treq):
        """ Return period at specified MJD with fractional days.
          time: MJD (with fractional days) of time when period is wanted.
        """
        nu = self.nu
        epoch = self.epoch
        nudot = self.nudot
#        p = self.p            ## period in seconds
#        pdot = self.pdot      ## period derivative in s/s (dimensionless)
        nuddot = self.nuddot  ## second derivative of frequency in Hz/s
        
        ti = (treq-epoch)*86400    ## difference from ephemeris epoch in seconds

        ##print('DBug: ti, p:', ti, p)
        ##print('DBug: pdot, nuddot:', pdot, nuddot)
        freq = nu + ti*nudot + (ti*ti/2)*nuddot  ## freqency at required time
        preq = 1.0/freq                          ## period at required time
        ##print('DBug: freq, preq:', freq, preq)
        return preq

    def arrtime( self, treq):
        """ Return arrival time of first pulse after time requested.
        """
        nu = self.nu
        nudot = self.nudot
        nuddot = self.nuddot
        epoch = self.epoch

        ti = (treq-epoch)*86400    ## difference from ephemeris epoch in seconds
        freq = nu + ti*nudot + (ti*ti/2)*nuddot  ## freqency at required time
        preq = 1.0/freq                          ## period at required time
        ti = (treq-epoch)*86400    ## difference from ephemeris epoch in seconds
#        p = 1.0/nu                 ## period in seconds
#        pdot = -nudot/(nu*nu)      ## period derivative in s/s (dimensionless)
        
        turns = nu*ti + nudot*ti*ti/2 +nuddot*ti*ti*ti/6.0  ## Num. turns
        ##print( "\nti = {:18.15f}, p = {:20.15f}, pdot = {:18.15e},
            ## nuddot = {:15.9e}, freq = {:15.9e}".format(
        ##    ti, p, pdot, nuddot, freq))
        ##print('DBug: turns:', turns)

        turns = turns - int( turns)   ## fractional number of turns.
        tarr = treq + (1.0 - turns)*preq/86400.0

        ##print('DBug: tarr:', tarr)
        ##print( "\nFirst pulse arrival time after required epoch (mjd & secs):{:d} {:.8f}\n".
        ##   format( int( tarr), (tarr - int( tarr))*86400.0))
        ##print( 'DBug: arrival MJD: {:.11f}'.format( tarr))
        ##print( "Preq(secs) and nu(Hz): {:.11f} {:.11f}\n\n".
        ##   format( preq, freq))

        return tarr

    def phase( self, treq):
        """ Return fraction of turn (the phase 0 .. 1) of the time requested.
        """
        nu = self.nu
        nudot = self.nudot
        nuddot = self.nuddot
        epoch = self.epoch

        ti = (treq-epoch)*86400    ## difference from ephemeris epoch in seconds
        freq = nu + ti*nudot + (ti*ti/2)*nuddot  ## freqency at required time
        preq = 1.0/freq                          ## period at required time
        ti = (treq-epoch)*86400    ## difference from ephemeris epoch in seconds
        p = 1.0/nu                 ## period in seconds
        pdot = -nudot/(nu*nu)      ## period derivative in s/s (dimensionless)
        
        turns = nu*ti + nudot*ti*ti/2 +nuddot*ti*ti*ti/6.0  ## Num. turns
        ##print( "\nti = {:18.15f}, p = {:20.15f}, pdot = {:18.15e},
           ## nuddot = {:15.9e}, freq = {:15.9e}".format(
        ##    ti, p, pdot, nuddot, freq))
        ##print('DBug: turns:', turns)

        turns = turns - int( turns)   ## fractional number of turns.
        if( turns < 0):
            turns = 1 + turns
        ##print('DBug: turns:', turns)
        return turns
        
##################################################

if( __name__=="__main__"):
    print( 'Testing spephemeris class.')

    ## define comparison function with tollerance
    def close( value1, value2, tollerance):
        return abs( value1 - value2) < tollerance

##Test case: GBO 20 m observation: 26186_26260
## observation epoch: 57794.019103417333099    (from inf file)
## nearest ephemeris line is:
##    15 FEB 17  57799  0.012400 100 29.6451772650  3 -368879.38 0.39   56.7983  -0.0228   100
## first giant pulse is 12.202279 s into observation
 #   (MJD = 57794.01924464741750853).
##
    ck1 = spephemeris( obsepoch=57794.019103417333099,
                       mjdd=57799,
                       mjds=0.012400,
                       nu=29.6451772650,
                       nudot=-368879.38)
    print( '\n--------')
    print( 'DBug: ck1 type info:', type( ck1))
    print( 'DBug: ck1.obsepoch:', ck1.obsepoch)
    print( 'DBug: ck1.SPFile, ck1.epoch, ck1.nu:',
           ck1.SPFile, ck1.epoch, ck1.nu)
    pulsetime = 12.202279/86400 + 57794.019103417333099
    testph = ck1.phase( pulsetime)
    print( 'DBug: ck1 phase at first pulse:', testph)
    assert close( testph, 0.046407, 0.0001)
    testat = ck1.arrtime( pulsetime)
    print( 'DBug: ck1 arrival time "after" first pulse:', testat)
    assert close( testat, 57794.0192454101, 1e-7)
    testp = ck1.period( pulsetime)
    print( 'DBug: ck1 period at first pulse:', testp)
    assert close( testp, 0.033732119, 0.000000001)
    print( '---------------------\n')

    print( '\n---------------------')
    testfile = '26186_26260.SPE.txt'
    print( '    Requires file {} in same directory.'.format( testfile))
    ck2 = spephemeris.readspe( testfile)
    print( '\nck2 type info:', type( ck2))
    print( 'DBug: ck2.obsepoch:', ck2.obsepoch)
    print( 'DBug: ck2.SPFile, ck2.epoch, ck2.nu, p:',
           ck2.SPFile, ck2.epoch, ck2.nu, 1/ck2.nu) 
    print( '----------\n')

    newp = ck2.period( ck2.epoch+1)
    print( 'DBug: ck2.period():', newp)
    
    print( 'DBug: ck2.arrtime():', '{:.17f}'.format(ck2.arrtime( ck2.epoch+1)))

    print( 'DBug: ck2.phase():', ck2.phase( ck2.epoch+1))

    for i in range( 1,6):
        f = i/5
        treq = newp*f/86400 + ck2.epoch+1
        treqrev = -newp*f/86400 + ck2.epoch+1
        print( 'DBug: +{:4.0f}% of period'.format( f*100), ck2.phase( treq))
        print( 'DBug: -{:4.0f}% of period'.format( f*100), ck2.phase( treqrev))

    ## should return zero phase:
    testphz = ck2.phase( ck2.arrtime( treq))
    print( 'DBug: zero phase??:', testphz)
    assert close( testphz, 0.0, 0.00001)

    print( 'DBug: ck2.obsepoch:', ck2.obsepoch)
    newp = ck2.period( ck2.obsepoch)
    print( 'DBug: ck2.period():', newp)
    
    print( 'DBug: ck2.arrtime():', '{:.17f}'.format(ck2.arrtime( ck2.obsepoch)))

    print( 'DBug: ck2.phase():', ck2.phase( ck2.obsepoch))

