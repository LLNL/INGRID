from Ingrid import Ingrid,LSN
import e2dgrid as e


i=Ingrid()
i.yaml['eqdsk']='/Users/holma2/Dropbox (Aalto)/UEDGE/INGRID/INGRID/seq#1/eqdsk'
i.OMFIT_read_psi()
i.psi_norm=0
i.eq=0
i.plate_data=0
l=LSN(i)
l.set_gridue_manual(*e.grid('/Users/holma2/Dropbox (Aalto)/UEDGE/E2D-EIR_integration/exporte2dgrid/aholm_mar1719_11'))
i.write_gridue(l.gridue_params,fname = 'gridue_e2d')

