#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 14:42:09 2020

@author: jguterl
"""
import os
try:
    import omfit
    from omfit.classes.omfit_eqdsk import OMFITgeqdsk
    from omfit.classes.omfit_eqdsk import OMFITaeqdsk
except:
    pass


def ImportEqFiles(Shot: int, Time: float, Path: str = os.getcwd(), Verbose: bool = True, Device: str = 'D3D', SNAPFile: str = 'EFIT01', ForceReload: bool = False)->tuple:

    if not os.path.isdir(Path):
        print('Cannot access {}'.format(Path))
        return None

    if type(Time) is not list: Time = [Time]
    Output = []
    for T in Time:
        if Verbose: print("Looking for a-file and g-file for the shot #{} at t={} at {} with SNAPFile:{}".format(Shot, T, Device, SNAPFile))
        gfilename = 'gfile.' + '{}.{}'.format(Shot, T)
        afilename = 'afile.' + '{}.{}'.format(Shot, T)
        afile = os.path.join(Path, afilename)
        gfile = os.path.join(Path, gfilename)
        print('gfile: {}, afile: {}'.format(gfile, afile))

        if not os.path.isfile(gfile) or ForceReload or not os.path.getsize(gfile) > 0:

            if Verbose: print(' > Loading gfile from Mdsplus ')

            g = OMFITgeqdsk(gfilename)
            g.from_mdsplus(device=Device, shot=Shot, time=T, exact=False, SNAPfile=SNAPFile, time_diff_warning_threshold=10, fail_if_out_of_range=True, show_missing_data_warnings=True, quiet=not Verbose)
            g.saveas(gfile)
            g.close()
            if os.path.getsize(gfile) == 0:
                print(" >  gfile {}  is empty .... Deleting it  ".format(gfile))
                os.remove(gfile)
                gfile = None
            else:
                print(" > gfile {} saved in {} ".format(gfile, Path))
        else:

            if Verbose: print(' > Found gfile:{} '.format(gfile))

        if not os.path.isfile(afile) or ForceReload or not os.path.getsize(afile) > 0:

            if Verbose: print(' > Loading afile from Mdsplus ')

            a = OMFITaeqdsk(afilename)
            a.from_mdsplus(device=Device, shot=Shot, time=T, exact=False, SNAPfile=SNAPFile, time_diff_warning_threshold=10, fail_if_out_of_range=True, show_missing_data_warnings=True, quiet=not Verbose)
            a.saveas(afile)
            a.close()

            if os.path.getsize(afile) == 0:
                    print(" >  afile {}  is empty .... Deleting it  ".format(afile))
                    os.remove(afile)
                    afile = None
            else:
                    print(" > afile {} saved in {} ".format(afile, Path))
        else:
            if Verbose: print(' > Found afile:{} '.format(afile))

        Output.append((afile, gfile))
    if len(Output) == 1:
                Output = Output[0]

    return Output
