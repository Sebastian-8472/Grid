#-------------------------------------------------------------------------------
# Name:         Programma
# Purpose:      Questo script serve per poter eseguire tutti i moduli senza doverli avviare singolarmente.
#
# Author:      Sebastiano
#
# Created:     23/02/2016
# Copyright:   (c) Sebastiano 2016
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import os
Cartella_corrente = os.getcwd()
import sys
sys.path.append(Cartella_corrente)
if __name__ == '__main__':
    from M0_Creazione_Griglia_di_punti_sparsi import*
    from M1_Geometria_base import*
    from M2_Figure_piane import*
else:
    print('Questo script non e\' un modulo')


