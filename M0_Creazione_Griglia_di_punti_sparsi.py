# -------------------------------------------------------------------------------
# Name:        CReazione serie di punti in posizione sparsa nello spazio
# Purpose:     Crea una serie di punti Sparsi allo scopo di testare il programma
# di creazione modello terreno da punti generici.
#
# Author:      SEBASTIANO
#
# Created:     27 Luglio 2015
# Copyright:   (c) SEBASTIANO 2015
# Licence:     <your licence>
# -------------------------------------------------------------------------------
from Librerie_gen import *


def Creatore_Punto(max_x, max_y, max_z):
    """
Questa funzione genera un array con le coordinate di un singolo punto.
Le coordinate sono comprese tra zero e il valore massimo assegnato
=================================================================================
ARGOMENTI:
---------------------------------------------------------------------------------
    :param max_x:       ==>     TIPO: numero.
                                Valore massimo della x

    :param max_y:       ==>     TIPO: numero
                                Valore massimo della y

    :param max_z:       ==>     TIPO:  Numero
                                Valore massimo della z
=================================================================================
RETURN:
    Punto       ==>     Tipo: numpy.array
                        Array

    """
    Punto = np.round(np.random.random((1, 3)), precisione) * np.array([max_x, max_y, max_z])
    # print(Punto)  stampa delle coordinate per il controllo della funzione
    return Punto


def Lista_punti(Num_punti, max_x=None, max_y=None, max_z=None):
    """
    Questa funzione crea una lista di punti dalle coordinate casuali.
    Num_punti = numero dei punti nella lista
    max_x = Valore massimo della x
    max_y = Valore massimo della y
    max_z = Valore massimo della z
    """
    if max_x is None:
        max_x = Rnd.randrange(Num_punti)
    if max_y is None:
        max_y = Rnd.randrange(Num_punti)
    if max_z is None:
        max_z = Rnd.randrange(Num_punti)
    Elenco_punti = []  # Creazione lista vuota dove salvare i punti
    for n in range(Num_punti):
        P = Creatore_Punto(max_x, max_y, max_z)  # Creazione dei punti
        Elenco_punti.append(P)  # aggiunta punto alla lista dei punti
    return Elenco_punti


def Save_lista_to_File(Lista, Cartella=None, Nome_file=None):
    """
    Questa funzione esporta la lista di punti in un file di testo.
    VARIABILI
    Lista      ==>  lista di punti da esportare
    Cartella   ==>  Stringa contenente il percorso per la cartella
                    di salvataggio del file
    """
    # Selezione cartella nel quale salvare il file di testo
    if Cartella is None:
        Percorso_file = Choose_folder()
    else:
        Percorso_file = Cartella
    # Scelta Nome del File.
    if Nome_file is None:
        file_name = input('Scegliere nome file')
    else:
        file_name = Nome_file
    # Creazione File-
    Posizione_file = Percorso_file + '\\ '+ file_name +'.txt'
    File_X_lista = open(Posizione_file, 'w')
    # Scrittura dei dati.
    for e in range(len(Lista)):
        File_X_lista.writelines(str(Lista[e][0])+'\t'+str(Lista[e][1])+'\t'+str(Lista[e][2])+'\n')
    File_X_lista.close()
    return Posizione_file


### TEst ###
def main():
    Punti = Lista_punti(20, 2, 2000, 2000)
    for ele in Punti:
        print(ele)
    # Save_lista_to_File(Punti, Cartella='D:\\Documenti\\Phyton Modules')


if __name__ == '__main__':
    main()
