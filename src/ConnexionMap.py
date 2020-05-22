#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 15:06:10 2020

@author: jguterl
"""
List=['IDL', 'IPF', 'ISB', 'ICB', 'IST', 'ICT', 'OST', 'OCT', 'OSB', 'OCB', 'ODL', 'OPF']
ConnexionMap={}
for p in List:
    if 'C' in p:
        ConnexionMap[p]={'N':p[0]+'S'+p[2]+'_S'}

  ConnexionMap['IDL']={'N':'IDL_S'}
  ConnexionMap['IPF']={'N':'IDL_S'}
  ConnexionMap['ISB']={'N':'IDL_S'}
  ConnexionMap['ICB']={'N':'IDL_S'}
  ConnexionMap['IST']={'N':'IDL_S'}
  ConnexionMap['ICT']={'N':'IDL_S'}
  ConnexionMap['OST']={'N':'IDL_S'}
  ConnexionMap['OCT']={'N':'IDL_S'}
  ConnexionMap['OSB']={'N':'IDL_S'}
  ConnexionMap['OCB']={'N':'IDL_S'}
  ConnexionMap['ODL']={'N':'IDL_S'}
  ConnexionMap['OPF']={'N':'IDL_S'}
  
  