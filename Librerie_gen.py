# -------------------------------------------------------------------------------
# Name:        Moduli comuni
# Purpose:     Serve per limitare il nome di importazioni di moduli generici e
#              anche per omogneizzare i nomi dei moduli importati nei vari file
#
# Author:      Sebastiano
#
# Created:     15/02/2016
# Copyright:   (c) Sebastiano 2016
# Licence:     <your licence>
# -------------------------------------------------------------------------------
import os   # Gestione Cartelle e sistema operativo.
import sys  # Gestione Cartelle e sistema operativo.
import sympy  # Gestione delle espressioni simboliche per aumentare la precisione
# from numba import jit, vectorize
from PIL import Image  # Gestione immagine grafici.
import numpy as np  # Gestione Vettori e altre funzioni matematiche.
import random as Rnd  # CCreazione di numeri casuali.
from Sebastiano import*  # Funzioni Base che uso in tutti i moduli
#  Numero dei decimali per i vari tipi di dati.
precisione = int(3)  # Numero dei decimali per valori generici.
Cord_dec = precisione  # Numero dei decimali per le coordinate.
Angl_dec = precisione + 1  # Numero dei decimali per gli angoli.
Vers_dec = int((2 * precisione) + 1)  # Numero dei decimali per i versori e vettori.

Tipi_di_numeri = [type(int()), type(float()), type(np.float16()),
                  type(np.float32()), type(np.float64())]


def main():
    print(os)
    print(sys)
    print(sympy)
    print(Image)
    print(np)
    print(Rnd)
    print("precisione Generale {} decimali \n precisione Coordinate {} "
          "decimali \n precisione angolare  {} decimali \n precisione Vettori "
          "{} decimali".format(precisione, Cord_dec, Angl_dec, Vers_dec))


if __name__ == '__main__':
    main()
