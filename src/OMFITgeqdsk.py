############################
# G-FILE CLASS OMFITgeqdsk #
############################
from __future__ import print_function
import numpy as np
import re

class OMFITgeqdsk(dict):#(SortedDict, OMFITascii):
    """
    class used to interface G files generated by EFIT
    :param filename: filename passed to OMFITascii class
    :param \**kw: keyword dictionary passed to OMFITascii class
    """
    def __init__(self, filename, **kw):
        """
        Method used to read g-files
        :param raw: load gEQDSK exactly as it's on file, regardless of COCOS
        """
        def splitter(inv, step=16):
            value = []
            for k in range(int(len(inv) / step)):
                value.append(inv[step * k:step * (k + 1)])
            return value

        def merge(inv):
            return ''.join(inv)
        
        self.clear()

        # clean lines from the carriage returns
        with open(filename, 'r') as f:
            EQDSK = f.read().splitlines()

        # first line is description and sizes
        self['CASE'] = np.array(splitter(EQDSK[0][0:48], 8))
        try:
            tmp = list([_f for _f in EQDSK[0][48:].split(' ') if _f])
            [IDUM, self['NW'], self['NH']] = list(map(int, tmp[:3]))
        except ValueError: # Can happen if no space between numbers, such as 10231023
            IDUM = int(EQDSK[0][48:52])
            self['NW'] = int(EQDSK[0][52:56])
            self['NH'] = int(EQDSK[0][56:60])
            tmp = []
#            printd('IDUM, NW, NH',IDUM,self['NW'],self['NH'],topic='OMFITgeqdsk.load')
        if len(tmp) > 3:
            self['EXTRA_HEADER'] = EQDSK[0][49 + len(re.findall('%d +%d +%d ' % (IDUM, self['NW'], self['NH']), EQDSK[0][49:])[0]) + 2:]
        offset = 1

        # now, the next 20 numbers (5 per row)
        [self['RDIM'], self['ZDIM'], self['RCENTR'], self['RLEFT'], self['ZMID'], \
         self['RMAXIS'], self['ZMAXIS'], self['SIMAG'], self['SIBRY'], self['BCENTR'], \
         self['CURRENT'], self['SIMAG'], XDUM, self['RMAXIS'], XDUM, \
         self['ZMAXIS'], XDUM, self['SIBRY'], XDUM, XDUM] = list(map(eval, splitter(merge(EQDSK[offset:offset + 4]))))
        offset = offset + 4

        # now I have to read NW elements
        nlNW = int(np.ceil(self['NW'] / 5.))
        self['FPOL'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlNW])))))
        offset = offset + nlNW
        self['PRES'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlNW])))))
        offset = offset + nlNW
        self['FFPRIM'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlNW])))))
        offset = offset + nlNW
        self['PPRIME'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlNW])))))
        offset = offset + nlNW
        try:
            # official gEQDSK file format saves PSIRZ as a single flat array of size rowsXcols
            nlNWNH = int(np.ceil(self['NW'] * self['NH'] / 5.))
            self['PSIRZ'] = np.reshape(np.fromiter(splitter(''.join(EQDSK[offset:offset + nlNWNH])), dtype=np.float),(self['NH'], self['NW']))
            offset = offset + nlNWNH
        except ValueError:
            # sometimes gEQDSK files save row by row of the PSIRZ grid (eg. FIESTA code)
            nlNWNH = self['NH'] * nlNW
            self['PSIRZ'] = np.reshape(np.fromiter(splitter(''.join(EQDSK[offset:offset + nlNWNH])), dtype=np.float),(self['NH'], self['NW']))
            offset = offset + nlNWNH
        self['QPSI'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlNW])))))
        offset = offset + nlNW

        # now vacuum vessel and limiters
        [self['NBBBS'], self['LIMITR']] = list(map(int, [_f for _f in EQDSK[offset:offset + 1][0].split(' ') if _f]))
        offset = offset + 1

        nlNBBBS = int(np.ceil(self['NBBBS'] * 2 / 5.))
        self['RBBBS'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlNBBBS]))))[0::2])
        self['ZBBBS'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlNBBBS]))))[1::2])
        offset = offset + max(nlNBBBS, 1)

        try:
            # this try/except is to handle some gEQDSK files written by older versions of ONETWO
            nlLIMITR = int(np.ceil(self['LIMITR'] * 2 / 5.))
            self['RLIM'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlLIMITR]))))[0::2])
            self['ZLIM'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlLIMITR]))))[1::2])
            offset = offset + nlLIMITR
        except ValueError:
            # if it fails make the limiter as a rectangle around the plasma boundary that does not exceed the computational domain
            self['LIMITR'] = 5
            dd = self['RDIM'] / 10.
            R = np.linspace(0, self['RDIM'], 2) + self['RLEFT']
            Z = np.linspace(0, self['ZDIM'], 2) - self['ZDIM'] / 2. + self['ZMID']
            self['RLIM'] = np.array([max([R[0], np.min(self['RBBBS']) - dd]),
                                        min([R[1], np.max(self['RBBBS']) + dd]),
                                        min([R[1], np.max(self['RBBBS']) + dd]),
                                        max([R[0], np.min(self['RBBBS']) - dd]),
                                        max([R[0], np.min(self['RBBBS']) - dd])])
            self['ZLIM'] = np.array([max([Z[0], np.min(self['ZBBBS']) - dd]),
                                        max([Z[0], np.min(self['ZBBBS']) - dd]),
                                        min([Z[1], np.max(self['ZBBBS']) + dd]),
                                        min([Z[1], np.max(self['ZBBBS']) + dd]),
                                        max([Z[0], np.min(self['ZBBBS']) - dd])])