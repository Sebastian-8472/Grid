# -------------------------------------------------------------------------------
# Name:         T.I.M.
#
# Purpose:      Create a grid-like list of points coordinate for a terrain model starting from zero or a text file contaning the coordinates of each point.
#               This model could be used to could be used to simulate a computer model of an existing site, or to generate a computer model of an imaaginary site to use in architectural design.
#
# Author:      SEBASTIANO
#
# Created:     29/10/2014
# -------------------------------------------------------------------------------

"""
The module allows to create a computer model of the terrain in a site as long
as the grid as the same numbers of rows and colons. This number has to be even .
this module also allows to manipulate the Terrain computer module and also
trasform the list of points to an image in grey scale or color.
It also allows to increase the number of points or change from a list a
points describing the vertex of the grid to a description based on a point and
the normal of the face in that point.
it is been design to help the construction of terrain models in programs like
Autocad and Lumion; so most of is out puts are design to work for them.
"""

### Importazione dei moduli necessari.

##   module to manage folders and files.
from Librerie_gen import *  # Libreria generale
import Sebastiano as Seba  # Modulo con le funzioni piu' comuni
from PIL import Image  # Modulo per la gestione e creazione di immagini
##   module for vectorized calculation of mathematical functions.
import numpy as NP
import random as rnd


## Importing a List from a File
def Import_List(File_name=None, Dir=str(os.getcwd())):
    """
    This function allows to import text file list frome the files created by the
    matcad files.

    Import_List(File_name,Dir = str(os.getcwd()))  =>   array[]

    The resulting Array is in the form of a list of this kind:

        [[ X , Y , Z ]
         [ X , Y , Z ]
         [ X , Y , Z ]
             .....   ]]

    This function requires the name's file if it is in the same directory  of the
     Python file you are running; you can
    also set a diferent directory writing the path in a string in the parameter Dir.
    """
    if File_name is None:  # Setting the file name and opening it.
        File_raw = Seba.Carica()
        inizio = File_raw.find('name') + 6
        fine = File_raw.find('mode') - 2
        File = File_raw[inizio:fine]
    else:
        File = Dir + "\\" + File_name
    input_file = open(File, "r")
    C = []  # Transfering the data in the file to a provisory list.
    for line in input_file:
        a = line.strip()
        b = a.split()
        for i in range(len(b)):
            C.append(float(b[i]))
    input_file.close()
    A = []  # Creating a list for the Coordinates of the points.
    B = None
    for i in range(0, len(C), 3):
        X_punto = round(C[i], Cord_dec)
        Y_punto = round(C[i + 1], Cord_dec)
        Z_punto = round(C[i + 2], Cord_dec)
        c = [X_punto, Y_punto, Z_punto]
        print((len(C) - i))
        A.append(c)
        B = NP.array(A)
    return B


def Media_su_diag_principale(matrice):
    """
    questa funzione permette di calcolare la media tra i valori lungo il
    verso della diagonale principale.
    =============================================================================
    ARGOMENTI:
    -----------------------------------------------------------------------------
        :param matrice:         ==>         Tipo: Numpy.array
                                            Matrice su cui si vuole eseguire
                                            il calcolo.
    =============================================================================
    :return:        Tipo: numpy.Array
                    Matrice co le medie dei verice eseguita nella direzione
                    della matrice principale.
    """
    vertici_1 = matrice[:-1]
    vertici_1 = vertici_1.transpose()[:-1].transpose()
    vertici_2 = matrice[1:]
    vertici_2 = vertici_2.transpose()[1:].transpose()
    media_vertici = (vertici_1 + vertici_2)/2
    return NP.round(media_vertici, Cord_dec)


def Media_su_diag_secondaria(matrice):
    """
    questa funzione permette di calcolare la media tra i valori lungo il
    verso della diagonale secondaria.
    =============================================================================
    ARGOMENTI:
    -----------------------------------------------------------------------------
        :param matrice:         ==>         Tipo: Numpy.array
                                            Matrice su cui si vuole eseguire
                                            il calcolo.
    =============================================================================
    :return:        Tipo: numpy.Array
                    Matrice co le medie dei verice eseguita nella direzione
                    della matrice secondaria.
    """
    vertici_1 = matrice[:-1]
    vertici_1 = vertici_1.transpose()[1:].transpose()
    vertici_2 = matrice[1:]
    vertici_2 = vertici_2.transpose()[:-1].transpose()
    media_vertici = (vertici_1 + vertici_2)/2
    return NP.round(media_vertici, Cord_dec)


def Media_su_riga(matrice):
    """
    questa funzione permette di calcolare la media tra i valori lungo il
    verso della riga.
    =============================================================================
    ARGOMENTI:
    -----------------------------------------------------------------------------
        :param matrice:         ==>         Tipo: Numpy.array
                                            Matrice su cui si vuole eseguire
                                            il calcolo.
    =============================================================================
    :return:        Tipo: numpy.Array
                    Matrice co le medie dei verice eseguita nella direzione
                    della riga.
    """
    vertici_1 = matrice[:]
    vertici_1 = vertici_1.transpose()[:-1].transpose()
    vertici_2 = matrice[:]
    vertici_2 = vertici_2.transpose()[1:].transpose()
    media_vertici = (vertici_1 + vertici_2)/2
    return NP.round(media_vertici, Cord_dec)


def Media_su_colonna(matrice):
    """
    questa funzione permette di calcolare la media tra i valori lungo il
    verso della riga.
    =============================================================================
    ARGOMENTI:
    -----------------------------------------------------------------------------
        :param matrice:         ==>         Tipo: Numpy.array
                                            Matrice su cui si vuole eseguire
                                            il calcolo.
    =============================================================================
    :return:        Tipo: numpy.Array
                    Matrice co le medie dei verice eseguita nella direzione
                    della riga.
    """
    vertici_1 = matrice[:-1]
    vertici_2 = matrice[1:]
    media_vertici = (vertici_1 + vertici_2)/2
    return NP.round(media_vertici, Cord_dec)


class Griglia:
    """
    Questa classe permette di creare griglie con strutture regolari per la
    rappresentazione attraverso punti del terreno
    """
    def __init__(self, col, row, Lato, d_z, Q=0):
        """
        Questa funzione permette di creare un oggetto base in cui sono contenute
        le matrice che descrivono la griglia
        dei punti nel terreno.
        =========================================================================
        PARAMETRI:
        -------------------------------------------------------------------------
            :parameter col:     ==>     tipo: numero intero
                                        Numero delle colonne della griglia.

            :parameter row:     ==>     Tipo: numero intero
                                        Numero delle righe della griglia

            :parameter Lato:    ==>     Tipo: Numero
                                        Lunghezza del lato della cella
                                        della griglia

            :parameter d_z:     ==>     Tipo: numero float
                                        Variazione massima della quota da un
                                        punto all'altro

            :parameter Q:       ==>     Tipo: Numero float
                                        Quota di partenza del vertice dal quale
                                        si sviluppa la griglia

        """
        self.lato_cella = Lato
        self.righe = row
        self.colonne = col
        self.X = NP.zeros((row, col))  # Setting the empty matrix for the X
        self.Y = NP.zeros((row, col))  # Setting the empty matrix for the Y
        self.Z = NP.zeros((row, col))  # Setting the empty matrix for the Z
        self.Z[:] = Q
        for r in range(row):
            for c in range(col):
                self.X[r][c] = np.float64(c * Lato)
                self.Y[r][c] = np.float64(r * Lato)
                bump = d_z * float((rnd.random() - rnd.random()))
                Q_r = self.Z[r-1][c]  # riferimento quota riga precedente
                Q_c = self.Z[r][c-1]  # riferimento quota colonna precedente
                Q_z = round(bump + NP.mean(np.array([Q_r, Q_c])), Cord_dec)
                self.Z[r][c] = Q_z

    @property
    def Area_Piana(self):
        """
        Questa funzione calcola l'area di estensione sul piano xy della griglia
        del terreno.
        =========================================================================
        Argomenti
        -------------------------------------------------------------------------
        self: oggetto stesso
        :return:        Tipo: Numero
                        Valore dell'area piana
        """
        return pow(self.lato_cella, 2) * self.righe * self.colonne

    @property
    def Lunghezza_tot(self):
        """
        Questa funzione calcola la lunghezza totale dei due lati della griglia.
        =========================================================================
        Argomenti
        -------------------------------------------------------------------------
        self: oggetto stesso
        =========================================================================
        :return:        Tipo: numpy.array
                        Array contenente i valori della lunghezza totale per
                        le righe e per le colonne.
        """
        return NP.array([self.righe, self.colonne]) * self.lato_cella

    @property
    def Points(self):
        """
        Questa funzione restituisce la matrice delle coordinate dei punti.
        in questa matrice i singoli elementi sono i vettori posizionale di
        ciascun punto
        :return:    numpy.array
        """
        return np.dstack((self.X, self.Y, self.Z))

    def Grid_to_List(self):
        """
        Questa funzione permette di trasformare l'ogetto griglia in una lista
        delle coordinate dei punti dell modello TIM.
        =========================================================================
        ARGOMENTI:
        -------------------------------------------------------------------------
        self:
        =========================================================================
        :return:        Tipo: Numpy.array
                        Array della lista delle coordinate dei punti

        """

        return np.dstack((self.X.flatten(), self.Y.flatten(), self.Z.flatten()))

    @property
    def Highest_Point(self):
        """
        Questa funzione permette di ottenere le coordinate del vertice piu'
        alto della Griglia.
        =========================================================================
        ARGOMENTI:
        -------------------------------------------------------------------------
            :parameter self:        ==>         Tipo: object
                                            Oggetto Griglia nella quale si
                                            vuole effettuare la ricerca.
        =========================================================================
            :return:        Tipo: Numpy.array
                            Array contenente le coordinate del vertice piu'
                            alto nella griglia.
        """
        z_max = NP.max(self.Z)  # Valore massimo della zeta
        b = self.Z == z_max  # matrice boleana per la posizione nelle matrici
        x = self.X[b][0]  # Coordinata x del vertice
        y = self.Y[b][0]  # Coordinata y del vertice
        return NP.array([x, y, z_max])

    def Find_highest_position(self):
        """
        Questa funzione permette di trovare la posizione nella griglia del
        punto piu' alto.
        =========================================================================
        ARGOMENTI:
        -------------------------------------------------------------------------
            :parameter self:        ==>         Tipo: object
                                            Oggetto Griglia nella quale si
                                            vuole effettuare la ricerca.
        =========================================================================
            :return:        Tipo: tuple
                            Coppia dei valori di riga e colonna della griglia
                            dove si trova il vertice piu' alto.
        """
        z_max = NP.max(self.Z)  # Valore minimo della zeta
        b = self.Z == z_max  # matrice boleana per la posizione nelle matrici
        riga = NP.nonzero(b)[0][0]
        colonna = NP.nonzero(b)[1][0]
        return riga, colonna

    @property
    def Lowest_Point(self):
        """
        Questa funzione permette di ottenere le coordinate del vertice piu'
        basso della Griglia.
        =========================================================================
        ARGOMENTI:
        -------------------------------------------------------------------------
            :parameter self:        ==>         Tipo: object
                                            Oggetto Griglia nella quale si
                                            vuole effettuare la ricerca.
        =========================================================================
            :return:        Tipo: Numpy.array
                            Array contenente le coordinate del vertice piu'
                            basso nella griglia.
        """
        z_min = NP.min(self.Z)  # Valore minimo della zeta
        b = self.Z == z_min  # matrice booleana per la posizione nelle matrici
        x = self.X[b][0]  # Coordinata x del vertice
        y = self.Y[b][0]  # Coordinata y del vertice
        return NP.array([x, y, z_min])

    def Find_lowest_position(self):
        """
        Questa funzione permette di trovare la posizione nella griglia del
        punto piu' basso.
        =========================================================================
        ARGOMENTI:
        -------------------------------------------------------------------------
            :parameter self:        ==>         Tipo: object
                                            Oggetto Griglia nella quale si
                                            vuole effettuare la ricerca.
        =========================================================================
            :return:        Tipo: tuple
                            Coppia dei valori di riga e colonna della griglia
                            dove si trova il vertice piu' basso
        """
        z_min = NP.min(self.Z)  # Valore minimo della zeta
        b = self.Z == z_min  # matrice booleana per la posizione nelle matrici
        riga = NP.nonzero(b)[0][0]
        colonna = NP.nonzero(b)[1][0]
        return riga, colonna

    def Estrai_colonne(self, colonna):
        """
        Questa funzione serve per estrarre una determinata colonna o un gruppo
        di colonne.
        =========================================================================
        ARGOMENTI:
        -------------------------------------------------------------------------
            :param colonna:     ==>     Tipo:Lista
                                        Colonna che si vuole estrarre dalla
                                        griglia, l'indice deve essere inserito
                                        come lista di numeri anche se si vuole
                                        avere una colonna singola.
        =========================================================================
        :return:        Tipo: Griglia.object
                            Griglia con le nuove colonne
        """
        X, Y, Z = [], [], []
        if not isinstance(colonna, list):
            colonna = [colonna]
        for c in colonna:
            X.append(self.X.transpose()[c])
            Y.append(self.Y.transpose()[c])
            Z.append(self.Z.transpose()[c])
        return NP.dstack((X, Y, Z))[0]

    def Estrai_righe(self, righe):
        """
        Questa funzione serve per estrarre una determinata colonna o un gruppo
        di colonne.
        =========================================================================
        ARGOMENTI:
        -------------------------------------------------------------------------
            :param righe:     ==>     Tipo:Lista
                                        riga che si vuole estrarre dalla
                                        griglia, l'indice deve essere inserito
                                        come lista di numeri anche se si vuole
                                        avere una colonna singola.
        =========================================================================
        :return:        Tipo: Griglia.object
                            Griglia con le nuove righe.
        """
        X, Y, Z = [], [], []
        if not isinstance(righe, list):
            righe = [righe]
        for c in righe:
            X.append(self.X[c])
            Y.append(self.Y[c])
            Z.append(self.Z[c])
        return NP.dstack((X, Y, Z))[0]

    def Punti_centrali(self, diagonale='P'):
        """
        Questa funzione permette di calcolare le coordinate dei dei punti
        centrali per ogni singola cella della griglia. La funzione restituisce
        una griglia dei punti centrali.
        =========================================================================
        ARGOMENTI:
        -------------------------------------------------------------------------
            :param diagonale:       ==>     Tipo: Stringa.
                                        Identifica la diagonale da usare come
                                        riferimento.
                                        I possibili valori sono:
        "P" => La funzione segue la diagonale principale e divide la griglia
        in triangoli cosi:
                            +---+---+---+
                            |\  |\  |\  |
                            | + | + | + |
                            |  \|  \|  \|
                            +---+---+---+
        "S" => La funzione segue la diagonale secondaria e divide la griglia
        in triangoli cosi:
                            +---+---+---+
                            |  /|  /|  /|
                            | + | + | + |
                            |/  |/  |/  |
                            +---+---+---+

        "E" => a funzione segue entrambe le e divide la griglia in triangoli
        cosi:
                            +-----+-----+-----+
                            |\   /|\   /|\   /|
                            | \ / | \ / | \ / |
                            |  X  |  X  |  X  |
                            | / \ | / \ | / \ |
                            |/   \|/   \|/   \|
                            +-----+-----+-----+
        =========================================================================
            :return:        Tipo: Griglia
                            Griglia dei punti centrali della griglia originale
        """
        n_r = self.righe - 1
        n_c = self.colonne -1
        g = Griglia(n_r, n_c, self.lato_cella, 0, 0)
        if diagonale == 'P':
            g.X = Media_su_diag_principale(self.X)
            g.Y = Media_su_diag_principale(self.Y)
            g.Z = Media_su_diag_principale(self.Z)
        elif diagonale == 'S':
            g.X = Media_su_diag_secondaria(self.X)
            g.Y = Media_su_diag_secondaria(self.Y)
            g.Z = Media_su_diag_secondaria(self.Z)
        elif diagonale == 'E':
            x_p = Media_su_diag_principale(self.X)
            x_s = Media_su_diag_secondaria(self.X)
            g.X = NP.round((x_p+x_s)/2, Cord_dec)
            y_p = Media_su_diag_principale(self.Y)
            y_s = Media_su_diag_secondaria(self.Y)
            g.Y = NP.round((y_p + y_s) / 2, Cord_dec)
            z_p = Media_su_diag_principale(self.Z)
            z_s = Media_su_diag_secondaria(self.Z)
            g.Z = NP.round((z_p + z_s) / 2, Cord_dec)
        else:
            raise ValueError('Parametro non valido per diagonale.\nI parametri validi sono P, S, E')
        return g

    def Punti_in_Anello(self, indice):
        """
        Questa funzione permette di calcolare il numero dei punti presente in
        un particolare Anello.
        =========================================================================
        ARGOMENTI:
        -------------------------------------------------------------------------
            :param indice:      ==>     Tipo: Numero
                                        Indie dell'anello di cui si vuole
                                        sapere il numero di punti
        =========================================================================
        :return: Numero dei punti presenti in un particolare anello
        """
        num_in_righe = self.colonne - 2*indice
        num_in_colon = self.righe - 2*(1+indice)
        num_punti = 2*num_in_colon + 2*num_in_righe
        if num_punti < 1:
            return 1
        else:
            return num_punti

    def Numero_di_anelli(self):
        """
        Questa funzione permette di determinare il numero di anelli che e
        possibile creare per la griglia.
        :return: Numero
        """
        base = min(self.righe, self.colonne)
        if NP.mod(base, 2) == 0:
            num_= int(base/2)
        else:
            num_ = int(base/2)+1
        return num_

    def Anello(self, indice):
        """
        Questa funzione permette di costruire uno dei possibili anelli dalla
        matrice della griglia.
        La funzione restituisce l'array delle coordinate dei singoli punti
        che fanno parte dell'anello.
        =========================================================================
        ARGOMENTI:
        -------------------------------------------------------------------------
            :param indice:       ==>     Tipo: Numero intero
                                    Indice dell'anello che si vuole ottenere
                                    dalla griglia.
        =========================================================================
        :return:        Tipo: array
                        Array degli anelli
        """
        L = []
        num_punti_riga = (self.colonne-(2*indice))
        num_punti_colo = (self.righe-(2*indice))-2
        num_punti = 2*num_punti_colo + 2*num_punti_riga
        if indice > self.Numero_di_anelli():
            raise IndexError(indice)
        else:
            r_sup = self.Estrai_righe(indice)[indice: self.colonne-indice]
            r_inf = self.Estrai_righe(self.righe-indice-1)[indice: self.colonne-indice]
            c_ini = self.Estrai_colonne(indice)[indice+1:-1-indice]
            c_fin = self.Estrai_colonne(self.colonne-indice-1)[
                    indice+1:-1-indice]
            if num_punti < 2:
                return r_sup[0]
            else:
                anello = [r_sup, c_fin, NP.flip(r_inf, axis=0), NP.flip(c_ini, axis=0)]
                for ele in anello:
                    for p in ele:
                        L.append(p)
                return L

    def Anelli(self, indici=None):
        """
        Questa funzione permette di trasformare la griglia di punti in un
        array degli anelli dei punti. gli anelli vengono contanti partendo
        dall'anello esetrno; il primo anello ha inidce zero.
        =========================================================================
        ARGOMENTI:
        -------------------------------------------------------------------------
            :param indici:       ==>     Tipo: Numero intero/ lista
                                    Indice dell'anello che si vuole ottenere
                                    dalla griglia.
                                    Il valore predefinito e' None, e permette
                                    di estrarre tutti gli anelli
        =========================================================================
        :return:        Tipo: array
                        Array degli anelli
        """
        L = []
        # Determinazione del numero di anelli possibili
        if indici is None:
            for a in range(self.Numero_di_anelli()+1):
              L.append(self.Anello(a))
        elif isinstance(indici, tuple) or isinstance(indici, list):
            for a in indici:
                L.append(self.Anello(a))
        else:
            index = [indici]
            self.Anelli(index)
        return L

    def Turn_into_pict(self, Mod="B|N"):
        """
        This function transforms the site terrain grid model into a picture.
        The parameters for this function are:

            Lista   =>      Lista or grid you want to trasform into a picture.

            Mod     =>      This serves to decide if you want the picure in grey's scale or in color.
                            The default value is set "B|N"
                            The possible values for  this parameter are:

                                "B|N" => set the picture in grey scale.

                                "Col" => set the picture in colors:

                                    blue for the water.

                                    green and marron for the graund

        """
        q = self.Z.shape
        r = q[0]
        c = q[1]
        #   Checks the initial input data and if transform it in a grid like format if its needed.
        Z_min = NP.min(self.Z)
        Z_max = NP.max(self.Z)
        I_array = np.array(self.Z.shape)
        if Mod == "B|N":  # GRey scale picture.
            D_Q = Z_max - Z_min
            coef = 255 / D_Q
            for x in range(c):
                for y in range(r):
                    Z = (self.Z[y][r]) * coef
                    I_array[r][c] = Z
        else:  # Color picture
            coef_blu = 255 / abs(Z_min)
            coef_red = 255 / (Z_max / 2)
            coef_gre = coef_red
            for x in range(c):
                for y in range(r):
                    q = self.Z[r][c]
                    if q < 0:
                        Z = abs(q) * coef_blu
                        pixel = np.array([Z, 0, 0])
                    elif q > Z_max / 2:
                        Z = (q - (Z_max / 2)) * coef_red
                        pixel = np.array([0, Z, 0])
                    else:
                        Z = q * coef_gre
                        pixel = np.array([0, 0, Z])
                    I_array[r][c] = pixel
        picture_ = Image.fromarray(I_array)
        picture_.show()  # showing the picture.

    def Aumenta_risoluzione(self, n):
        """
        Thi function increase the resolution of the grid.
        Add_points(Lista,n)     =>      List[((r-1)*n)+1][((c-1)*n)+1]
        The parameters are:
            Lista => the lista you want to increase the resolution.
            n     => the factor of resolution increase.
        """
        q = self.Z.shape
        r = q[0]
        c = q[1]
        # Setting the new grid parameters.
        New_r = ((r - 1) * n) + 1
        New_c = ((c - 1) * n) + 1
        New_L = self.lato_cella /n
        # Creating the new matrix.
        New_Grid = Griglia(New_c, New_r, New_L, 0)
        for r in New_Grid.righe:
            for c in New_Grid.colonne:
                if NP.mod(r, n) == 0:  # indice riga divisibile per n
                    if NP.mod(c, n) == 0:  # indice colonna divisibile per n
                        New_Grid.Z[r][c] = self.Z[r/n][c/n]
                    else:
                        r_ = int(r/n)  # indice riga riferimento in originale
                        c_p = int(c/n)  # colonna precedente griglia originale
                        c_s = c_p +1  # colonna successiva griglia originale
                        d_z = self.Z[r_][c_p] - self.Z[r_][c_s]
                        z_ =self.Z[r_][c_p] + ((d_z/n) * NP.mod(c, n))
                        New_Grid.Z[r][c] = z_ + (NP.random.rand()/(n*10))
                else:  # indice non divisibile per n
                    if NP.mod(c, n) == 0:  # indice colonna divisibile per n
                        r_p = int(r / n)  # riga precedente in originale
                        r_s = r_p + 1  # riga successiva griglia originale
                        c_ = int(c / n)  # colonna griglia originale
                        d_z = self.Z[r_p][c_] - self.Z[r_s][c_]
                        z_ = self.Z[r_p][c_] + ((d_z / n) * NP.mod(r, n))
                        New_Grid.Z[r][c] = z_ + (NP.random.rand() / (n * 10))
                    else:
                        # indici dei punti matrici originali
                        r_p = int(r/n)
                        r_s = r_p +1
                        c_p = int(c/n)
                        c_s = c_p +1
                        # resto degli indici
                        resto_r = NP.mod(r, n)
                        resto_c = NP.mod(c, n)
                        z_0 = self.Z[r_p][c_p]
                        if resto_r == resto_c:  # delta principale
                            delta_dp = self.Z[r_s][c_s] - self.Z[r_p][c_p]
                            z = z_0 + delta_dp * (resto_r/n)
                            New_Grid.Z[r][c] = z
                        elif resto_r > resto_c:  # triangolare inferiore
                            delta_c1 = self.Z[r_s][c_p] - self.Z[r_p][c_p]
                            delta_r2 = self.Z[r_s][c_s] - self.Z[r_s][c_p]
                            delta = NP.mean([(delta_c1*(resto_c/n)),
                                             (delta_r2*(resto_c/n))])
                            New_Grid.Z[r][c] = z_0 + delta
                        else:  # triangolare superiore
                            delta_r1 = self.Z[r_p][c_s] - self.Z[r_p][c_p]
                            delta_c2 = self.Z[r_s][c_s] - self.Z[r_p][c_s]
                            delta = NP.mean([(delta_r1*(resto_c/n)),
                                             (delta_c2*(resto_c/n))])
                            New_Grid.Z[r][c] = z_0 + delta
        self.Z = New_Grid.Z
        self.X = New_Grid.Z
        self.Y = New_Grid.Y
        self.colonne = New_Grid.colonne
        self.righe = New_Grid.righe
        self.lato_cella = New_Grid.lato_cella

    def Move_Grid(self, V):
        """
        This function allow to move the all Grid in all 3 direction.
        Lista => grid of points.
        V => movement you want to impose to the grid.
        """
        self.X += V[0]
        self.Y += V[1]
        self.Z += V[2]

    def Scale_Grid(self, S, Node=None, Pos=None):
        """
        This function allow to scalare all Grid in all 3 direction.

        Lista => grid of points.

        Scala => Factor for scaling the grid.

                [x,y,z]

                x = Scale the x values
                y = Scale the y values
                z = Scale the z values

        Node => node around wich you want to scale the grid.

        Pos  => coordinate of the center.
        """
        if isinstance(S, list) or isinstance(S, tuple):
            S = S
        elif type(S) in Tipi_di_numeri:
            S = [S, S, S]
        else:
            S = None
        if all((Node is None, Pos is None)):
            self.X *= S[0]
            self.Y *= S[1]
            self.Z *= S[2]
        else:
            if Node is not None:
                r, c = Node[0], Node[1]
                V = NP.array([self.X[r][c], self.Y[r][c], self.Z[r][c]])
                self.Move_Grid(V*(0-1))
                self.X *= S[0]
                self.Y *= S[1]
                self.Z *= S[2]
                self.Move_Grid(V)
            else:
                self.Move_Grid(Pos*(0-1))
                self.X *= S[0]
                self.Y *= S[1]
                self.Z *= S[2]
                self.Move_Grid(Pos)

    def Invert_Quote(self):
        """
        Inverts the quotes of the grid.
                negatives quotes to positive quotes
                positive quotes to negative quotes.
        """
        self.Z *= -1

    def Set_to_zero(self):
        """
        This function sets the negative quotes to zero.
        """
        b = self.Z < 0
        self.Z[b] = 0




# TODO Griglia - Creare una funzione per ruotare la griglia - ne eseste una
#  tra gli appunti ma non se se abbia senso fa ruotare la griglia.

# TODO Griglia - creare funzione per trovari punti con z minore/ maggiore/ di
#  una certo valore  (per istruzioni su come articolare la funzione 6/36)

# TODO Griglia - Modificare la funzione per aumentare la risoluzione in modo
#  do poter scegliere una delle tre varianti per i punti lungo la diagonale

'''
# Rotate the  Grid
def Rot_Grid(Lista, Ang, Node=None, Pos=None, Ang_unit="Dec"):
    """
    This function allow to rotate all Grid in all 3 direction.

    Lista => grid of points.

    Ang => angle you want to rotate the grid.

            [a,b,g]

            a = rotazione intorno all'asse z
            b = rotazione intorno all'asse y
            g = rotazione intorno all'asse x

    Node => node around wich you want to spin the grid.

    Pos  => coordinate of the center of rotation.

    Ang_unit => Specify the angle unit.
                The default value is set to Decimals
                The possible values are:

                    "Dec"   =>  Decimals

                    "Rad"   =>  Radiants

    """
    q = Shape_Grid(Lista)
    r = q[1]
    c = q[0]
    #   Checks the initial input data and if transform it in a grid like format if its needed.
    if q[3] == "List":
        Griglia = List_to_Grid(Lista)
    elif q[3] == "Grid":
        Griglia = Lista
    else:
        Griglia = None
    # Setting the indentity matrix of rank 3 as a base for the rotation matrix
    Rot_mat = NP.identity(3)
    # Setting the angles to Radinants unit
    if Ang_unit == "Dec":
        for i in range(3):
            Ang[i] = NP.deg2rad(Ang[i])
    #   Assigning the angle to the variable.
    a = Ang[0]
    b = Ang[1]
    g = Ang[2]
    #   Creating the Rotation Matrix.
    Rot_mat[0][0] = cos(a) * cos(b)
    Rot_mat[0][1] = (cos(g) * sin(a)) - (sin(g) * sin(b) * cos(a))
    Rot_mat[0][2] = (sin(g) * sin(a)) + (cos(g) * sin(b) * cos(a))
    Rot_mat[1][0] = sin(a) * cos(b) * (0 - 1)
    Rot_mat[1][1] = (cos(g) * cos(a)) - (sin(g) * sin(b) * sin(a))
    Rot_mat[1][2] = (sin(g) * sin(a)) + (cos(g) * sin(b) * sin(a))
    Rot_mat[2][0] = (0 - 1) * sin(b)
    Rot_mat[2][1] = (0 - 1) * sin(g) * cos(b)
    Rot_mat[2][2] = cos(g) * cos(b)
    #   Setting the movement in case of the rotation around a specified node or a specified position.
    Vector_M = NP.array([0, 0, 0])
    if Pos is not None:
        for i in range(len(Pos)):
            Vector_M[i] = Pos[i]
    if Node is not None:
        x = Node[0] * q[2]
        y = Node[1] * q[2]
        Vector_M = NP.array([x, y, 0])
    #   Changing the origin of the oxix to the specified point.
    mov_1 = Vector_M * (0 - 1)
    GRID_1 = Move_Grid(Griglia, mov_1)
    #   Rotating the Grid
    Rot = Griglia
    for col in range(c):
        for row in range(r):
            P = NP.array(GRID_1[col][row])
            N_P = NP.dot(Rot_mat, P)
            for i in range(3):
                Rot[col][row][i] = round(N_P[i], Dec_prec)
    #   Moving back the the original axes system.
    Grid = Move_Grid(Rot, Vector_M)
    #   Returning the rotate grid.
    return Grid
    

# Esports the rings to a files.
def Esport_Rings_to_files(Lista, n=None, f_type=".txt", separator=","):
    """
    This function allows to extract the rows or the columns to a file.
        Extract_line_to_file(Lista,Line = "Rows",n = [],form = ".dat",Dir = os.getcwd(),sep = ",") => File/nested list.

    The input for this function are:

        Lista   =>      List or grid of the site terrain computer model.

        n       =>      Indicates the lines to extract.
                        The default value is set to [] meaning it extracts all linees in the grid.
                        This value has to have always the form of a list

        f_type    => indicates the format of text file you want to write.
                    The default value is set to ".txt"


        separator     =>  Set the separator for the values  in the files.
                    The default value is set to ","

    """
    q = Shape_Grid(Lista)
    c = q[0]
    #   Checks the initial input data and if transform it in a grid like format if its needed.
    if q[3] == "List":
        M_Cols = List_to_Grid(Lista)
    elif q[3] == "Grid":
        M_Cols = Lista
    else:
        M_Cols = None
    M_Rows = Ext_Rows(Lista)
    Num_rings = int(c / 2)
    #   Setting the rings to extract.
    if n is None:
        n = list(range(Num_rings))
    else:
        n = n
    Folder = Seba.Choose_folder()
    for i in n:
        fin = abs(i - c)
        b = list(range(4))
        b[0] = list(M_Cols[i][i:fin])  # C_0
        b[1] = list(M_Rows[fin - 1][i + 1:fin])  # R_0
        b[2] = list(M_Cols[fin - 1][i:fin - 1])  # C_1
        b[2].reverse()
        b[3] = list(M_Rows[i][i:fin - 1])  # R_1
        b[3].reverse()
        Ring = []
        for e in range(4):
            for s in range(len(b[e])):
                Ring.append(b[e][s])
        F_N = "Ring" + str(n[i])
        Esport_List_to_Files(Ring, File_name=F_N, form=f_type, Dir=Folder, sep=separator)
    print(Folder)


def Rnd_Flat_cells(Lista, N_F):
    """
    This funcion creates a List of N_F random faces for the function Add_Flat planes;

    """
    q = Shape_Grid(Lista)
    #   Setting the Data imput in a grid like format.
    Griglia = List_to_Grid(Lista)
    #   Creation of the face to be changed selection list.
    P = []
    for i in range(N_F):
        #   it is subctract 2 to the max value to easy the handling of the cpy of face tha could accour usinf arandom function.
        x = rnd.randint(0, q[0] - 2)
        y = rnd.randint(0, q[1] - 2)
        Z = (Griglia[x][y][2] + Griglia[x][y + 1][2] + Griglia[x + 1][y + 1][2] + Griglia[x + 1][y][2] +
             Griglia[x][y + 1][2] + Griglia[x + 1][y + 1][2]) * (1.0 / 4.0)
        Q = round(Z, 3)
        P.append([x, y, Q])
    #   Sortinf the List of faces to change ( sortin by columns).
    P = sorted(P)
    #   Subviding the lisst of faces to change according to the columns in a grid like list.
    I = list(range(q[0]))
    for o in range(len(I)):
        I[o] = []
    c = 0
    for e in range(1, N_F):
        if e != 0:
            if P[e][0] != P[e - 1][0]:
                c = c + 1
                I[c].append(P[e])
            else:
                I[c].append(P[e])
        else:
            I[0] = P[0]

    #   Checking the near faces to make sure they have the same elevation.
    for k in range(4000):
        print(k)
        #   First in the direction from the beggining zero to the end
        for col in range(len(I)):
            for row in range(len(I[col])):
                Q = I[col][row][2]
                Row = I[col][row][1]
                #   Check on same column.
                for r_i in range(row + 1, len(I[col])):  #
                    r_0 = I[col][r_i][1]
                    z_0 = I[col][r_i][2]
                    if Row == r_0:
                        I[col][r_i][1] = I[col][row][1] + 1
                        I[col][r_i][2] = 0.5 * (Q + z_0)
                    if Row == r_0 - 1:
                        I[col][r_i][2] = 0.5 * (Q + z_0)
                    if Row == r_0 + 1:
                        I[col][r_i][2] = 0.5 * (Q + z_0)
                #   Check on next column.
                for r_n in range(len(I[col + 1])):
                    Z_1 = I[col + 1][r_n][1]
                    Q_R = I[col + 1][r_n][2]
                    if Row == Z_1:
                        I[col + 1][r_n][2] = 0.5 * (Q + Q_R)
                    if Row == 1 + Z_1:
                        I[col + 1][r_n][2] = 0.5 * (Q + Q_R)
                    if Row == Z_1 - 1:
                        I[col + 1][r_n][2] = 0.5 * (Q + Q_R)
        for cl in range(len(I) - 1, 0, -1):
            for rw in range(len(I[cl]) - 1, 0, -1):
                Q = I[cl][rw][2]
                Row = I[cl][rw][1]
                #   Check on same column.
                for r_i in range(len(I[cl]) - 1, rw + 1, -1):  #
                    r_0 = I[cl][r_i][1]
                    z_0 = I[cl][r_i][2]
                    if Row == r_0:
                        I[cl][r_i][1] = I[cl][rw][1] + 1
                        I[cl][r_i][2] = 0.5 * (Q + z_0)
                    if Row == r_0 - 1:
                        I[cl][r_i][2] = 0.5 * (Q + z_0)
                    if Row == r_0 + 1:
                        I[cl][r_i][2] = 0.5 * (Q + z_0)
                #   Check on next column.
                for r_n in range(len(I[cl - 1]) - 1, 0, -1):
                    Z_1 = I[cl - 1][r_n][1]
                    Q_R = I[cl - 1][r_n][2]
                    if Row == Z_1:
                        I[cl - 1][r_n][2] = 0.5 * (Q + Q_R)
                    if Row == 1 + Z_1:
                        I[cl - 1][r_n][2] = 0.5 * (Q + Q_R)
                    if Row == Z_1 - 1:
                        I[cl - 1][r_n][2] = 0.5 * (Q + Q_R)
    # changin the Index from agrid like format to a list-like format
    Index = []
    for n in range(len(I)):
        for m in range(len(I[n])):
            Index.append(I[n][m])
    return Index


##   Adding flat planes.
def Add_Flat_plane(Lista, N_flat=1, mod=None, P=None):
    """
    This function adds flat planes. It is possible to choose a point in the surface, o nod or to put them random in the grid.

    Lista   =>  The list/grid of the points of the terrain.

    N_Flat  =>  numbers of plans to add.

    Q       =>  Set the elevation of the new flat plane.

    Mod     =>  Allows to substitue a determinated plane.

                The possible values are:

                    "x,y"  =>  Changes the plane containing the point with the coordinates indicated.

                    "Face"      =>  Changes the plane that are in the cross  of r and c.

    P       =>  Indicate the list of position/coordinates.
    """
    q = Shape_Grid(Lista)
    #   Checks the initial input data and if transform it in a grid like format if its needed.
    if q[3] == "List":
        Griglia = List_to_Grid(Lista)
    elif q[3] == "Grid":
        Griglia = Lista
    else:
        Griglia = None
    if mod == "Face":  # setting the parametres to add the flat plane knowing the face
        N_flat = len(P)
        Index = []
        for i in range(N_flat):
            x = int(P[i][0])
            y = int(P[i][1])
            if len(P) == 3:
                z = P[1][2]
            else:
                z = (Griglia[x][y][2] + Griglia[x][y + 1][2] + Griglia[x + 1][y + 1][2] + Griglia[x + 1][y][2] +
                     Griglia[x][y + 1][2] + Griglia[x + 1][y + 1][2]) * (1.0 / 4.0)
            Index.append([x, y, z])
    elif mod == "x,y":  # Setting the parameters to sobstitute a face with the point in it.
        N_flat = len(P)
        Index = []
        for i in range(N_flat):
            x = int(P[i][0] / q[2])
            y = int(P[i][1] / q[2])
            if len(P) == 3:
                z = P[1][2]
            else:
                z = (Griglia[x][y][2] + Griglia[x][y + 1][2] + Griglia[x + 1][y + 1][2] + Griglia[x + 1][y][2] +
                     Griglia[x][y + 1][2] + Griglia[x + 1][y + 1][2]) * (1.0 / 4.0)
            Index.append([x, y, z])
    else:  # Setting the parameters to substitute a flat plane at random.
        Index = Rnd_Flat_cells(Lista, N_flat)
    for i in range(N_flat):  # Adding the faces
        if i == 0:
            c_P = Index[i][0]
            r_P = Index[i][1]
            z_P = Index[i][2]
            Griglia[c_P][r_P][2] = z_P
            Griglia[c_P + 1][r_P][2] = z_P
            Griglia[c_P + 1][r_P + 1][2] = z_P
            Griglia[c_P][r_P + 1][2] = z_P
        else:
            c_P = Index[i][0]
            r_P = Index[i][1]
            z_P = Index[i][2]
            Griglia[c_P][r_P][2] = z_P
            Griglia[c_P + 1][r_P][2] = z_P
            Griglia[c_P + 1][r_P + 1][2] = z_P
            Griglia[c_P][r_P + 1][2] = z_P
    return Griglia



## Create a list to create a snake like polyline.
def Snake_Poly(Lista):
    q = Shape_Grid(Lista)
    r = q[1]
    c = q[0]
    #   Checks the initial input data and if transform it in a grid like format if its needed.
    if q[3] == "List":
        Griglia = List_to_Grid(Lista)
    elif q[3] == "Grid":
        Griglia = Lista
    else:
        Griglia = None
    Grid = list(range(c))
    for col in range(c):
        b = list(range(r))
        for row in range(r):
            b[row] = [0, 0, 0]
            for i in range(3):
                b[row][i] = Griglia[col][row][i]
        Grid[col] = b
    for col in range(1, c, 2):
        Grid[col].reverse()
    return Grid


# Create an flat external border.
def New_Border(Lista, z=0):
    """
    Thi funcion create a flat external border by adding two rinf of points the the original Grid-like Array of terrain points.

    """
    q = Shape_Grid(Lista)
    #   making sure the imput data is in grid like format.
    if q[3] == "List":
        Griglia = List_to_Grid(Lista)
    elif q[3] == "Grid":
        Griglia = Lista
    else:
        Griglia = None
    Ori_col = q[0]
    Ori_row = q[1]
    #   Setting up columns and rows size of the new matrix.
    New_col = Ori_col + 4
    New_row = Ori_row + 4
    #   Creating the new matrix.
    Border_Grid = Create_Grid(New_col, New_row, q[2], 0, Q=z)
    #   Movement calculation.
    Move = [(2 - 4) * q[2], (2 - 4) * q[2], 0]
    #   Moving the border matrix
    B_Grid = Move_Grid(Border_Grid, Move)
    #   Creating the border.
    for c in range(Ori_col):
        for r in range(Ori_col):
            B_Grid[c + 2][r + 2] = Griglia[c][r]
    return B_Grid


# Sovrapone a flat plane.
def Sovrapone_plane(Lista, z=0):
    """
    Thi funcion create a flat external border by adding two rinf of points the the original Grid-like Array of terrain points.

    """
    q = Shape_Grid(Lista)
    #   making sure the imput data is in grid like format.
    if q[3] == "List":
        Griglia = List_to_Grid(Lista)
    elif q[3] == "Grid":
        Griglia = Lista
    else:
        Griglia = None
    #   Moving down the grid.
    B_Grid = Move_Grid(Griglia, [0, 0, 0 - z])
    #   flattening the negatives vlues.
    C_Grid = Set_to_zero(B_Grid)
    #   Moving up the Grid.
    D_Grid = Move_Grid(C_Grid, [0, 0, z])
    return D_Grid


def Emersione(Lista, Ele=0):
    """
    This funcion counts the number of points andure a certain elevation and those above .
    """
    S = Shape_Grid(Lista)
    if S[3] == "List":
        L = Lista
    elif S[3] == "Grid":
        L = Grid_to_List(Lista)
    else:
        L = None
    sup = 0.00
    inf = 0.00
    ugu = 0.00
    for e in range(len(L)):
        if L[e][2] < Ele:
            inf += 1
        elif L[e][2] > Ele:
            sup += 1
        else:
            ugu += 1
    p_sup = round((sup / (sup + inf + ugu)), 3) * 100
    p_inf = round((inf / (sup + inf + ugu)), 3) * 100
    p_ugu = round((ugu / (sup + inf + ugu)), 3) * 100
    Data = [[sup, p_sup, "above "], [ugu, p_ugu, "equal to "], [inf, p_inf, "below "]]
    for i in range(3):
        a = str(Data[i][1]) + " % of points are " + str(Data[i][2]) + str(Ele)
        print(a)
    return Data


## Extract the columns to file to import in CAD.
def Esp_Cols_for_CAD(Lista, n=None, f_type=".txt", separator=","):
    """
    This function allows to extract the rows or the columns to a file to use in AutoCAd.
        Extract_line_to_file(Lista,Line = "Rows",n = [],form = ".dat",Dir = os.getcwd(),sep = " ") => File

    The input for this function are:

        Lista   =>      List or grid of the site terrain computer model.

        n       =>      Indicates the lines to extract.
                        The default value is set to [] meaning it extracts all linees in the grid.

        f_type    => indicates the format of text file you want to write.
                    The default value is set to ".txt"

        Separator     =>  Set the separator for the values  in the files.
                    The default value is set to ","

    """
    q = Shape_Grid(Lista)
    #   Checks the initial input data and if transform it in a grid like format if its needed.
    if q[3] == "List":
        Griglia = List_to_Grid(Lista)
    elif q[3] == "Grid":
        Griglia = Lista
    else:
        Griglia = None
    #   Setting the cols to extract.
    D = Ext_Cols(Griglia, n)
    Data = Create_Grid(q[0], q[1] + 1, 0, 0)
    for c in range(q[0]):
        for r in range(len(Data[c])):
            if r == 0:
                Data[c][r] = D[c][r]
            else:
                Data[c][r] = D[c][r - 1]
    #   Creating the file.
    if n is None:
        n = list(range(len(Data)))
    Folder = media.choose_folder()
    for i in range(len(Data)):
        F_N = "Colonna " + str(n[i])
        Esport_to_File(Data[i], File_name=F_N, form=f_type, Dir=Folder, sep=separator)
    print(Folder)


## exctracting the rows to a file to import in CAD.
def Esp_Rows_for_CAD(Lista, n=None, f_type=".txt", separator=","):
    """
    This function allows to extract the rows or the columns to a file  to Import in AutoCad.
        Extract_line_to_file(Lista,Line = "Rows",n = [],form = ".dat",Dir = os.getcwd(),sep = " ") => File/nested list.

    The input for this function are:

        Lista   =>      List or grid of the site terrain computer model.

        n       =>      Indicates the lines to extract.
                        The default value is set to [] meaning it extracts all linees in the grid.
                        This value has to have always the form of a list

        f_type    => indicates the format of text file you want to write.
                    The default value is set to ".txt"

        separator     =>  Set the separator for the values  in the files.
                    The default value is set to ","

    """
    q = Shape_Grid(Lista)
    #   Checks the initial input data and if transform it in a grid like format if its needed.
    if q[3] == "List":
        Griglia = List_to_Grid(Lista)
    elif q[3] == "Grid":
        Griglia = Lista
    else:
        Griglia = None
    #   Setting the cols to extract.
    #   Setting the cols to extract.
    D = Ext_Rows(Griglia, n)
    Data = Create_Grid(q[0], q[1] + 1, 0, 0)
    for c in range(q[0]):
        for r in range(len(Data[c])):
            if r == 0:
                Data[c][r] = D[c][r]
            else:
                Data[c][r] = D[c][r - 1]
    #   Creating the file.
    if n is None:
        n = list(range(len(Data)))
    Folder = media.choose_folder()
    for i in range(len(Data)):
        F_N = "Riga " + str(n[i])
        Esport_to_File(Data[i], File_name=F_N, form=f_type, Dir=Folder, sep=separator)
    print(Folder)


#### Test
##
## Lista = Import_List("noti 1 lati.dat" ,"D:\\Documenti\\Modello Terreno" )
##
## F = Rnd_Flat_cells(Lista,100)




## Esport the data in to a text file.
def Esport_to_File(Lista, File_name=None, form=".txt", Dir=os.getcwd(), sep=","):
    """
    This function allows you to esport your lista or grid to a file.
    The imput are:
        Lista = > The Data matrix you want to esport to file in the form of a list.
                    This data can be a simple list or an Grid like list.
        File_name => The name of the file in the form of a string.
                     The deffault value is set to None, so it prompts  a window wher you can choose folder file name and format.

        form    =>  Indicates the format of the file you want to write.
                    The default value is set to ".tx.t"

        Dir     =>  Set the directory where the file is saved.
                    The default value is set to the direcory where the python file is saved.

        sep     =>  Set the separator for the values  in the files.
                    The default value is set to ","


    """
    #   Set parameters for calculations.
    Data = None
    q = Shape_Grid(Lista)
    #   Checks the initial input data and transforms it in a list format if its needed.
    if q[3] == "List":
        Data = Lista
    elif q[3] == "Grid":
        Data = Grid_to_List(Lista)
    #   Set the path to save the file.
    if File_name is None:  # if you d'ont have specified a file name and a directory.
        File_path = Seba.Choose_file_save_path()
    else:  # In case you have specified a file name and a directory.
        File_path = Dir + "\\" + File_name + form
    #   Writes the files.
    File = open(File_path, "w")  # create the file.
    for n in range(len(Data)):
        a = str(Data[n][0])
        for i in range(1, len(Data[n])):
            a = a + sep + str(Data[n][i])
        a = a + "\n"
        File.write(a)  # Writes the data in the file
    print((str(File_path)))  # Prints the file path to find the file :)


## Extract the columns to file.
def Esport_Cols_to_file(Lista, n=None, f_type=".txt", separator=","):
    """
    This function allows to extract the rows or the columns to a file.
        Extract_line_to_file(Lista,Line = "Rows",n = [],form = ".dat",Dir = os.getcwd(),sep = ",") => File

    The input for this function are:

        Lista   =>      List or grid of the site terrain computer model.

        n       =>      Indicates the lines to extract.
                        The default value is set to [] meaning it extracts all linees in the grid.

        f_type    => indicates the format of text file you want to write.
                    The default value is set to ".txt"

        Separator     =>  Set the separator for the values  in the files.
                    The default value is set to ","

    """
    q = Shape_Grid(Lista)
    #   Checks the initial input data and if transform it in a grid like format if its needed.
    if q[3] == "List":
        Griglia = List_to_Grid(Lista)
    elif q[3] == "Grid":
        Griglia = Lista
    else:
        Griglia = None
    #   Setting the cols to extract.
    Data = Ext_Cols(Griglia, n)
    #   Creating the file.
    if n is None:
        n = list(range(len(Data)))
    Folder = Seba.Choose_folder()
    for i in range(len(Data)):
        F_N = "Colonna " + str(n[i])
        Esport_to_File(Data[i], File_name=F_N, form=f_type, Dir=Folder, sep=separator)
    print(Folder)


## Esport the List-like data in to a text file.
def Esport_List_to_Files(Lista, File_name=None, form=".txt", Dir=os.getcwd(), sep=","):
    """
    This functio allows you to esport your lista or grid to a file.
    The imput are:
        Lista = > The Data matrix you want to esport to file in the form of a list.
                    This data can be a simple list or an Grid like list.
        File_name => The name of the file in the form of a string.
                     The deffault value is set to None, so it prompts  a window wher you can choose folder file name and format.

        form    =>  Indicates the format of the file you want to write.
                    The default value is set to ".tx.t"

        Dir     =>  Set the directory where the file is saved.
                    The default value is set to the direcory where the python file is saved.

        sep     =>  Set the separator for the values  in the files.
                    The default value is set to ","


    """
    #   Set the path to save the file.
    if File_name is None:  # if you d'ont have specified a file name and a directory.
        File_path = Seba.Choose_file_save_path()
    else:  # In case you have specified a file name and a directory.
        File_path = Dir + "\\" + File_name + form
    #   Writes the files.
    File = open(File_path, "w")  # create the file.
    for n in range(len(Lista)):
        a = str(Lista[n][0])
        for i in range(1, len(Lista[n])):
            a = a + sep + str(Lista[n][i])
        a = a + "\n"
        File.write(a)  # Writes the data in the file
    print((str(File_path)))  # Prints the file path to find the file :)



## Turning list-like setup into grid-like setup.
def List_to_Grid(Lista):
    """
    Thi function allows to trasform an array of point's coordinate from a list-like format into a grid-like format.

    Lista   =>      List-like array data points.

    """
    Cols = Shape_Grid(Lista)[0]  # Setting the number of columns.
    rows = Shape_Grid(Lista)[1]  # Setting the number of rows.
    if type(Lista[0][0]) != NP.ndarray:  # Case 1: the data input is a list-like array
        Grid = Create_Grid(Cols, rows, Shape_Grid(Lista)[2], 0)
        i = 0
        for c in range(Cols):
            for r in range(rows):
                Grid[c][r] = Lista[i]
                i += 1
    else:  # Case 2: the data input is a grid like case.
        Grid = Lista
    return Grid


## extracting the rows to a file.
def Esport_Rows_to_files(Lista, n=None, f_type=".txt", separator=","):
    """
    This function allows to extract the rows or the columns to a file.
        Extract_line_to_file(Lista,Line = "Rows",n = [],form = ".dat",Dir = os.getcwd(),sep = ",") => File/nested list.

    The input for this function are:

        Lista   =>      List or grid of the site terrain computer model.

        n       =>      Indicates the lines to extract.
                        The default value is set to [] meaning it extracts all linees in the grid.
                        This value has to have always the form of a list

        f_type    => indicates the format of text file you want to write.
                    The default value is set to ".txt"

        separator     =>  Set the separator for the values  in the files.
                    The default value is set to ","

    """
    q = Shape_Grid(Lista)
    #   Checks the initial input data and if transform it in a grid like format if its needed.
    if q[3] == "List":
        Griglia = List_to_Grid(Lista)
    elif q[3] == "Grid":
        Griglia = Lista
    else:
        Griglia = None
    #   Setting the cols to extract.
    Data = Ext_Rows(Griglia, n)
    #   Creating the file.
    if n is None:
        n = list(range(len(Data)))
    Folder = Seba.Choose_folder()
    for i in range(len(Data)):
        F_N = "Riga " + str(n[i])
        Esport_to_File(Data[i], File_name=F_N, form=f_type, Dir=Folder, sep=separator)
    print(Folder)


    def Center_Normal(self, Diag="P"):
        """
        This function returns a list whit the coordinate of a center point and
        the normal veros of the terrain in that point.
        Center_Normal(Lista, Diag = "P") => Array[[x,y,z,Nx,Ny,Nz],[....],[...]]
        The input data are:
            Lista => Lista of points or a grid of points that defines a the site terrain model.
            Diag => Idifies the diagonal used to divide the cells in triangules.
                    The defaoult value is set to "P"
                    The possible values are:

                        "P" => Thi allows to follow the principal diagonal, and divides the grid in triagules like this:

                        -------------
                        |+ /|+ /|+ /|
                        | / | / | / |
                        |/ +|/ +|/ +|
                        -------------

                        "S" => Thi allows to follow the second diagonal, and divides the grid in triagules like this:

                        -------------
                        |\ +|\ +|\ +|
                        | \ | \ | \ |
                        |+ \|+ \|+ \|
                        -------------

                        "E" => Thi allows to follow both diagonals, and divides the grid in triagules like this:

                        +-----+-----+-----+
                        |\ + /|\ + /|\ + /|
                        | \ / | \ / | \ / |
                        |+ X +|+ X +|+ X +|
                        | / \ | / \ | / \ |
                        |/ + \|/ + \|/ + \|
                        +-----+-----+-----+
        The Out pun of this function is a list.
        """
        #   Set parametres for calculations.
        q = self.Z.shape
        r = q[1]
        c = q[0]

        if Diag == "P":  # Set the parameters for following the principal diagonal.
            NUM = (r - 1) * (c - 1) * 2
            st = -1
            end = 3
            pas = 2
        elif Diag == "S":  # Set the parameters for following the secondary diagonal.
            NUM = (r - 1) * (c - 1) * 2
            st = 0
            end = 4
            pas = 2
        elif Diag == "E":  # Set the parameters for following both diagonals.
            NUM = (r - 1) * (c - 1) * 4
            st = 0
            end = 4
            pas = 1
        else:
            NUM = None
            st = None
            end = None
            pas = None
        Cen_Sup = NP.zeros((NUM, 3), float)  # Provisory list to store the center point's coordinate
        Ve_norm = NP.zeros((NUM, 3), float)  # Provisory list to store the normal versor
        e = 0
        A = []
        for N_c in range(c - 1):
            for N_r in range(r - 1):
                # Separating the vertex of the cells
                a = NP.array(Griglia[N_c][N_r], float)
                b = NP.array(Griglia[N_c + 1][N_r], float)
                c = NP.array(Griglia[N_c + 1][N_r + 1], float)
                d = NP.array(Griglia[N_c][N_r + 1], float)
                Cella = NP.array([a, b, c, d], float)
                B = NP.vstack((Cella, Cella))
                for num in range(st, end, pas):
                    # Central point of the triangolar surface
                    Cen_Sup[e] = (1.0 / 3.0) * (NP.array((B[num] + B[num + 1] + B[num + 2])))
                    v = num + 1
                    # Normal versor calculation
                    V_1 = (2 - 3) * (NP.array((B[v] - B[v + 1]), float))
                    V_2 = (2 - 3) * (NP.array((B[v] - B[v - 1]), float))
                    V_n = NP.cross(V_1, V_2)
                    E = NP.sqrt(NP.sum(NP.square(V_n)))
                    Ve_norm[e] = (NP.true_divide(V_n, E))  # Normal versor
                    A.append([Cen_Sup[e][0], Cen_Sup[e][1], Cen_Sup[e][2], Ve_norm[e][0], Ve_norm[e][1], Ve_norm[e][2]])
        B = NP.array(A)
        return B






     SAND BOX
 #   Checking the near faces to make sure they have the same elevation.
    #   Second in the direction from end to the beggining.
        for col in range(len(I),0):
            for row in range(len(I[col],0)):
                Q = I[col][row][2]
                Row = I[col][row][1]
    #   Check on same column.
                for r_i in range(len(I[col]),row+1): #
                    Z_0 = I[col][r_i][1]
                    if Row == Z_0:
                        I[col][r_i][1] = I[col][row][1]+1
                        I[col][r_i][2] = 0.5 * (Q + I[col][r_i][2])
                    if Row == Z_0 +1:
                        I[col][r_i][2] = 0.5 * (Q + I[col][r_i][2])
    #   Check on next column.
                for r_n in range(len(I[col-1]),0):
                    Z_1 = I[col-1][r_n][1]
                    Q_R = I[col-1][r_n][2]
                    if Row == Z_1:
                         I[col-1][r_n][2] = 0.5 * (Q + Q_R)
                    if Row == 1 + Z_1:
                         I[col-1][r_n][2] = 0.5 * (Q + Q_R)
                    if Row == Z_1 -1 :
                         I[col-1][r_n][2] = 0.5 * (Q + Q_R)


    return(Index)
'''
