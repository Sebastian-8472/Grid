# -------------------------------------------------------------------------------
# Name:        Geometria Base
# Purpose:      Avere un modulo che permette di gestire piani e rette usando i
#               metodi e le propriet? dell'algebra lineare
#
# Author:      Sebastiano Pc fisso
#
# Created:     10/08/2015
# Copyright:   (c) Sebastiano Pc fisso 2015
# Licence:     <your licence>
# -------------------------------------------------------------------------------
### Importazione dei moduli
import M0_Creazione_Griglia_di_punti_sparsi as M0
# Modulo per creare lista di punti per eseguire i test
from Librerie_gen import *
from numba import jit

### Origine del sistema di riferimento
O = np.zeros(3)
# list che raccoglie tutti i tipi di numeri


def Conversione_in_stringhe(punto):
    """
    Questa funzione permette di convertire ciascuna delle coordinate del
    punto da numero in stringa.
=================================================================================
    ARGOMENTI
    -----------------------------------------------------------------------------
        :param punto:       ==>     TIPO: numpy.array
                                    punto di cui si vuole trasformare le
                                    coordinate danumero in stringa
=================================================================================
    :return:
        :return P:      ==>     Tipo: numpy.array
                                array del punto dove le coordinate sono delle
                                stringhe invece che numeri
    """
    P = []
    for c in punto:
        stringa = str(c)
        s = stringa[:stringa.find('.') + Cord_dec]
        P.append(s)
    return P


def Importa_lista_punti(percorso_file=''):
    """
Funzione che permette di importare una lista di punti da un file di testo.
=================================================================================
ARGOMENTI:
---------------------------------------------------------------------------------
    :param percorso_file:       ==>     Tipo:   Stringa
                                        Permette l'inserimento del percorso
                                        del file di coordinate da importare.
                                        Il percorso va inserito in forma di stringa.
                                        Se non specificato si apre un finestra
                                        per la scelta del file da importare.
=================================================================================
RETURN:
    Lista_Coordinate_Punti      ==>     Tipo: List
                                        Lista dei punti che si trovano nel file.
    """
    if percorso_file == '':
        Lista_punti_file_percorso = Choose_file_path()  # Permette la scelta del file
    else:
        # Permette l'inserimento del percorso del file
        Lista_punti_file_percorso = percorso_file
    punti = open(Lista_punti_file_percorso)  # Aperture del file
    # Creazione Lista vuota in cui inserire le coordinate
    Lista_Coordinate_Punti = []
    for line in punti:
        # Legge e raccoglie i dati in maniera grezza, con parentesi ecc. ...
        Raw_data_cordinate = line.split()
        # Trasformazione in array
        Data = []  # Variabile dove raccogliere i valori delle coordinate
        for element in Raw_data_cordinate:
            num = np.trunc(np.round(np.float(element) * (10 ** Cord_dec),
                                    0)) / 10 ** Cord_dec
            # Trasformazione in numeri tipo float e creazione punto
            Data.append(round(num, Cord_dec))
            # Aggiunta punto alla lista dei punti come array
        Lista_Coordinate_Punti.append(np.array(Data))
    return Lista_Coordinate_Punti


def Grandezza_vettore(Vettore):
    """
Calcola la grandezza di un vettore.
=================================================================================
ARGOMENTI
---------------------------------------------------------------------------------
    :param Vettore:     ==>     Tipo: numpy.array
                                Vettore, deve essere in formato array del
                                modulo Numpy.
=================================================================================
RETURN
    :return Grand:      ==>     Tipo: Numero
                                Grandezza del vettore.
    """
    Grand = round(np.sqrt(np.dot(Vettore, Vettore)), Vers_dec)
    return Grand


def Ordina_punti(Lista_punti):
    """
Questa funzione ordina i punti all'interno della lista. I punti vengono
ordinati in ordine crescente a meno che non si indichi diversamente.
=================================================================================
ARGOMENTI:
---------------------------------------------------------------------------------
    :param Lista_punti:     ==>     Tipo: List
                                    Lista dei punti da ordinare.
=================================================================================
RETURN
    l_       ==>         Tipo: List
                        Lista dei punti nell'ordine indicato.
    """
    # Creazione del database
    diz = {}
    for v in Lista_punti:
        if len(v) == 1:
            vet = Conversione_in_stringhe(v[0])
        else:
            vet = Conversione_in_stringhe(v)
        if vet[0] in diz.keys():
            if vet[1] in diz[vet[0]].keys():
                diz[vet[0]][vet[1]].append(vet[2])
            else:
                diz[vet[0]][vet[1]] = [vet[2]]
        else:
            print(v, vet, len(vet))
            diz.__setitem__(vet[0], {vet[1]: [vet[2]]})
    L_ = []
    x = sorted(diz.keys())
    for x_ in x:
        y = sorted(diz[x_])
        for y_ in y:
            z = sorted(diz[x_][y_])
            for z_ in z:
                p = np.array([np.float64(x_), np.float64(y_), np.float64(z_)])
                L_.append(p)
    return L_


def Dist(punto_1, punto_2, Plane=''):
    """
Questa funzione Restituisce la distanza tra i punti P_1 e P_2.
=================================================================================
VARIABILI
---------------------------------------------------------------------------------
    :param punto_1:     ==>     Tipo: numpy.array
                                Array delle coordinate del punto 1.

    :param punto_2:     ==>     Tipo: numpy.array.
                                Array delle coordinate del punto 2.

    :param Plane:       ==>     Tipo: stringa:
                                Indica un piano di preferenza sul quale
                                calcolare la distanza.
                                'xy' ==> Calcola la distanza sul piano xy.
                                'zy' ==> Calcola la distanza sul piano zy.
                                'xz' ==> Calcola la distanza sul piano xz.
=================================================================================
RETURN
    :return dist:       ==>     Tipo: Numero
                                Distanza tra i due punti indicati.
                                Se si e' indicato un piano e' la distanza
                                misurata su quel piano.
    """
    piano = Plane.lower()  # Trasforma in minuscolo le lettere inserite per
    # semplificare il confronto.
    if type(punto_1) == type(list()):
        punto_1 = np.array(punto_1)
    elif type(punto_1) != type(np.zeros((3, 1))):
        raise TypeError('non e\' possibile calcolare la distanza per %s' % (type(punto_1)))
    if type(punto_2) == type(list()):
        punto_2 = np.array(punto_2)
    elif type(punto_2) != type(np.zeros((3, 1))):
        raise TypeError('non e\' possibile calcolare la distanza per %s' % (type(punto_2)))
    v = (punto_1 - punto_2)
    Dist_ = 0.00
    if piano == '':
        Dist_ = Grandezza_vettore(v)
    elif piano == 'xy' or piano == 'yx':
        Dist_ = Grandezza_vettore(v[0:2])
    elif piano == 'zy' or piano == 'yz':
        Dist_ = Grandezza_vettore(v[1:3])
    elif piano == 'xz' or piano == 'zx':
        Dist_ = Grandezza_vettore(v[1:3:2])
    dist = round(Dist_, precisione)
    return dist


def Versori(Vettore):
    """
Questa funzione restituisce il versore del vettore che viene passato come
argomento.
=================================================================================
VARIABILI
---------------------------------------------------------------------------------
    :param Vettore:     ==>     Tipo: numpy array
                                Vettore di cui si vuole ottenere i versori.
=================================================================================
RETURN
    :return vers:       ==>     Tipo: numpy.array
                                Versore del vettore indicato
    """
    Dim_vet = Grandezza_vettore(Vettore)
    if Dim_vet == 0.0:
        vers = Vettore
    else:
        vers = np.array(Vettore / Dim_vet)
    vers.round(Vers_dec)
    return vers


def Azimut_oriz(Vettore):
    """
Questa funzione l'azimut tra due punti sul piano xy, partendo dalle
coordinate dei punti.
_L_'angolo restituito ? quello mostrato in figura.
       y^
        | ang  /
        |    /
        |  /
        |/
        -------------> x
=================================================================================
VARIABILI
---------------------------------------------------------------------------------
    :param Vettore:      ==>     Tipo: numpy array.
                                Vettore di cui si vuole calcolare l'angolo di
                                inclinazione rispetto all'asse y nel piano xy.
=================================================================================
RETURNS
    :return Alpha:      ==>     Tipo: numero.
                                Angolo cercato.
    """
    Delta_X = round(Vettore[0], Vers_dec)
    Delta_Y = round(Vettore[1], Vers_dec)
    Alpha = 0.00
    if Delta_Y == 0.0:
        if Delta_X > 0.0:
            Alpha = np.pi * 0.5
        else:
            Alpha = np.pi * (3.0 / 2.0)
    elif Delta_X == 0.0:
        if Delta_Y > 0.0:
            Alpha = 0.0
        else:
            Alpha = np.pi
    elif Delta_X > 0.0:
        if Delta_Y > 0.0:
            Alpha = np.arctan((Delta_X / Delta_Y))
        else:
            Alpha = np.pi + np.arctan((Delta_X / Delta_Y))
    elif Delta_X < 0.0:
        if Delta_Y > 0.0:
            Alpha = (2 * np.pi) + np.arctan((Delta_X / Delta_Y))
        else:
            Alpha = np.pi + np.arctan((Delta_X / Delta_Y))
    Alpha = round(Alpha, Angl_dec)
    return Alpha


def Azimut_vert(Vettore):
    """
Questa funzione l'inclinazione verticale tra due punti partendo dalle coordinate dei
punti.
    _L_'angolo restituito ? quello mostrato in figura.
       z^
        |  /
        | /
        |/ angolo
        -------------> Plane xy
=================================================================================
VARIABILI
---------------------------------------------------------------------------------
    :param Vettore:      ==>     Tipo: numpy array.
                                Vettore di cui si vuole calcolare l'angolo di
                                inclinazione rispetto all'asse y nel piano xy.
=================================================================================
RETURNS
    :return Alpha:      ==>     Tipo: numero.
                                Angolo cercato.
    """
    origine = np.zeros((1, 3))
    radius = Dist(origine, Vettore)
    if radius == 0.0:
        Ang = 0.0
    else:
        cos_ang = Vettore[2] / radius
        Ang = ((np.pi / 2) - (np.arccos(cos_ang)))
    return round(Ang, Angl_dec)


def Polar_comp(punto):
    """
    Trasforma le coordinate cartesiane in coordinate polari.
    =================================================================================
    ARGOMENTI
    ---------------------------------------------------------------------------------
        :param punto:       ==>     Tipo: numpy.array
                                    Vettore in coordinate cartesiane che si vuole
                                    trasformare in vettore di coordinate polari.
    =================================================================================
    RETURNING
    TIPO: numpy.array
    Restituisce le coordinate polari partendo dalle coordinate cartesiane
    nella forma: numpy.array([ radius , Ang_oriz , Ang_vert ])
    """
    origine = np.zeros((1, 3))
    # Calcola il raggio distanza dall'origine.
    radius = Dist(origine, punto)
    # Calcola l'inclinazione rispetto all'asse delle "Y"
    ang_oriz = Azimut_oriz(punto)
    # Calcola l'inclinazione rispetto all'asse delle "Z"
    ang_vert = Azimut_vert(punto)
    return np.array([radius, ang_oriz, ang_vert])


def Matrice_rotazione_asse_x(angolo):
    """
Crea una matrice di rotazione rispetto all'asse x.
=================================================================================
ARGOMENTI:
---------------------------------------------------------------------------------
    :param angolo:      ==>     Tipo: Numero
                                Angolo di cui si vuole eseguire la rotazione
                                rispetto all'asse x.
                                _L_'angolo deve essere espresso in Radianti.
                                I valori positivi dell'angolo esprimono rotazioni
                                antiorarie.
=================================================================================
RETURNING
    type numpy.array   Matrice di rotazione rispetto all,asse x
    """
    R_1 = [1.0, 0.0, 0.0]  # Valori della prima riga
    R_2 = [0.0, np.cos(angolo), -np.sin(angolo)]  # Valori della seconda riga
    R_3 = [0.0, np.sin(angolo), np.cos(angolo)]  # Valori della terza riga.
    Matrix = np.array([R_1, R_2, R_3])
    return Matrix


def Matrice_rotazione_asse_y(angolo):
    """
Crea una matrice di rotazione rispetto all'asse y.
=================================================================================
ARGOMENTI:
---------------------------------------------------------------------------------
    :param angolo:      ==>     Tipo: Numero
                                Angolo di cui si vuole eseguire la rotazione
                                rispetto all'asse y.
                                _L_'angolo deve essere espresso in Radianti.
                                I valori positivi dell'angolo esprimono rotazioni
                                antiorarie.
=================================================================================
RETURNING
    :type numpy.array:   Matrice di rotazione rispetto all'asse y
    """
    R_1 = [np.cos(angolo), 0.0, np.sin(angolo)]  # Valori della prima riga
    R_2 = [0.0, 1.0, 0.0]  # Valori della seconda riga
    R_3 = [-np.sin(angolo), 0.0, np.cos(angolo)]  # Valori della terza riga.
    Matrix = np.array([R_1, R_2, R_3])
    return Matrix


def Matrice_rotazione_asse_z(angolo):
    """
Crea una matrice di rotazione rispetto all'asse z.
=================================================================================
ARGOMENTI:
---------------------------------------------------------------------------------
    :param angolo:      ==>     Tipo: Numero
                                Angolo di cui si vuole eseguire la rotazione
                                rispetto all'asse z.
                                _L_'angolo deve essere espresso in Radianti.
                                I valori positivi dell'angolo esprimono rotazioni
                                antiorarie.
=================================================================================
RETURNING
    :type numpy.array:   Matrice di rotazione rispetto all'asse z
    """
    R_1 = [np.cos(angolo), np.sin(angolo), 0.0]  # Valori della prima riga
    R_2 = [-np.sin(angolo), np.cos(angolo), 0.0]  # Valori della seconda riga
    R_3 = [0.0, 0.0, 1.0]  # Valori della terza riga.
    Matrix = np.array([R_1, R_2, R_3])
    return Matrix


@jit
def Rotazione(vettore, Angoli, Centro_Rotazione=None, asse_rotazione=None,
              unit='rad'):
    """
Questa funzione permette di eseguire rotazioni nello spazio. La rotazione puo'
essere eseguita attorni agli assi di una base canonica, rispetto a un punto e
rispetto a un asse particolare.
=================================================================================
ARGOMENTI:
---------------------------------------------------------------------------------
    :param vettore:     ==>     Tipo: numpy.array
                                punto o vettore da ruotare
    :param Angoli:      ==>     Tipo: numpy.array
                                Angoli di rotazione rispetto agli assi canonici
                                a = angolo di rotazione intorno all'asse x
                                b = angolo di rotazione intorno all'asse y
                                c = angolo di rotazione intorno all'asse z
                                Se viene assegnato un solo valore insieme all'asse
                                di rotazione si intende come angolo di rotazione
                                intorno all'asse indicato.

    :param Centro_Rotazione:    ==>     Tipo: numpy.array
                                Punto intorno a cui avviene la rotazione

    :param asse_rotazione:      ==>      Tipo:  numpy.array
                                Asse intorno a cui avviene la rotazione

    :param unit:        ==>     Tipo: stringa
                                unita di misura degli angoli il valore e'
                                impostato su radianti
                                'dec' per decimali
                                'cen' per centesimali
                                'rad' per radianti
=================================================================================
RISULTATI:
    :type numpy.array:      vettore ruotato in forma di numpy array.
    """
    if type(Angoli) == type(np.zeros((1, 3))[0]):
        ang_r = Angoli
    elif type(Angoli) == type([]):
        ang_r = np.array(Angoli)
    else:  # valore numerico singolo
        ang_r = np.array([Angoli, 0, 0])  # Rotazione intorno all'asse x
    # Cambio unita di misura
    if unit == 'dec':
        ang_r = (ang_r * np.pi) / 180.00
    elif unit == 'cen':
        ang_r = (ang_r * np.pi) / 200.00
    if type(asse_rotazione) == type(None):  # vettore degli angoli di rotazione intorno agli assi canonici
        R_x = Matrice_rotazione_asse_x(ang_r[0])
        R_y = Matrice_rotazione_asse_y(ang_r[1])
        R_z = Matrice_rotazione_asse_z(ang_r[2])
        Rot_mat = R_x * R_y * R_z  # * equals multiplication row * column
        if type(Centro_Rotazione) == type(None):
            V = np.dot(vettore, Rot_mat)  # Rotazione Semplice
            V = np.array(V[0])[0]
        else:
            V = np.dot((vettore - Centro_Rotazione), Rot_mat) + \
                Centro_Rotazione  # Rotazione rispetto al punto
            V = np.array(V[0])[0]
    else:
        asse = Versori(asse_rotazione)  # asse attorno a cui eseguire la rotazione.
        # Calcoli degli angoli direzionali.
        pro_asse_su_xy = Proiezione_Punto_su_Piano(asse, piano_xy)  # Riduzione
        # del vettore da tre a due dimensioni
        l_asse_xY = Grandezza_vettore(pro_asse_su_xy)  # Grandezza vettore a due dimensioni.
        alpha = np.arccos(asse[0] / l_asse_xY)  # Angolo dell'asse rispetto all'asse x
        if np.isnan(alpha):
            alpha = np.nan_to_num(alpha)
        beta = np.arccos(asse[2])  # Angolo dell'asse rispetto all'asse z
        if np.isnan(beta):
            beta = np.nan_to_num(beta)
        # Calcolo della matrice cambiamento coordinate in modo da avere asse
        # z = asse vettore
        M_1 = Matrice_rotazione_asse_z(alpha) * Matrice_rotazione_asse_y(beta)
        # Matrice di rotazione
        M_R = Matrice_rotazione_asse_z(Angoli)
        # Ritorno al sistema di riferimento originario
        M_1_inv = Matrice_rotazione_asse_y(-beta) * Matrice_rotazione_asse_z(-alpha)
        # Matrice di rotazione def
        Rot_mat = M_1 * M_R * M_1_inv
        if Centro_Rotazione is None:
            V = np.dot(vettore, Rot_mat)  # Rotazione Semplice
            V = np.array(V)[0]
        else:
            # Rotazione rispetto al punto
            V = np.dot((vettore - Centro_Rotazione), Rot_mat) + Centro_Rotazione
            V = np.array(V)[0]
    return V


@jit(parallel=True)
def Ordinazione_polare(Lista):
    """
Questa funzione ordina i punti secondo le loro cordinate polari, partendo dal
raggio, poi l'angolo orizzontale e l'ultimo parametro preso in considerazione
e' l'angolo verticale.
=================================================================================
ARGOMENTI:
---------------------------------------------------------------------------------
    :param Lista:       ==>         :type list:
                                    Lista dei punti che si vuole ordinare.
=================================================================================
RETURNING
    :type list:     Lista ordinate secondo i criteri indicati.
    """
    # Creazione della lista con le co ordinate cartesiane e polari del punto.
    M_polare = []
    for punto in Lista:
        P_polar = Polar_comp(punto)
        v = punto
        for value in P_polar:
            v = np.append(v, value)
        M_polare.append(v)
    # Ordinazione rispetto al Raggio
    Lista_1 = [M_polare[0]]  # Lista dei punti ordinati secondo il raggio
    for punto in M_polare[1:]:
        for n in range(len(Lista_1)):
            if punto[3] < Lista_1[n][3]:
                Lista_1.insert(n, punto)
                break
            elif punto[3] == Lista_1[n][3]:  # Ordinamento secondo angolo
                # orizzontale
                if punto[4] < Lista_1[n][4]:
                    Lista_1.insert(n, punto)
                    break
                elif punto[4] == Lista_1[n][4]:  # Ordinamento secondo angolo verticale
                    if punto[5] <= Lista_1[n][5]:
                        Lista_1.insert(n, punto)
                        break
                    else:
                        if n == len(Lista_1) - 1:
                            Lista_1.append(punto)
                            break
                        elif n < len(Lista_1) - 1:
                            if punto[4] == Lista_1[n + 1][4]:
                                pass
                            else:
                                Lista_1.append(punto)
                                break

                elif punto[4] > Lista_1[n][4]:
                    if n == len(Lista_1) - 1:
                        Lista_1.append(punto)
                        break
                    elif n < len(Lista_1) - 1:
                        if punto[3] == Lista_1[n + 1][3]:
                            pass
                        else:
                            Lista_1.append(punto)
                            break
        if punto[3] > Lista_1[-1][3]:
            Lista_1.append(punto)  # Creazione della lista con solo le coordinate
    Lista = []  # Lista dove verrano salvate solo le coordinate cartesiane del punto.
    for punto in Lista_1:
        Lista.append(punto[0:3])
    return Lista


def Angoli_tra_vettori(Vettore_1, Vettore_2, metodo=None):
    """
Questa funzione permette di calcolare l'angolo tra i due vettori.
Questa funzione usa sia il teorema del coseno che la definizione di prodotto
scalare. La funzione sotto intende che i vettori assegnati siano incidenti,
cioe' le rette su cui giacciono abbiano un punto in comune, o i piani di cui sono
i versori siano incidenti.
=================================================================================
ARGOMENTI
---------------------------------------------------------------------------------
    :PARAM Vettore_1:       ==>     Tipo: numpy.array
                                    Primo vettore.

    :param Vettore_2:       ==>     Tipo: numpy.array
                                    Secondo vettore

    :param metodo:          ==>     Tipo: numero intero
                                    Metodo con cui viene calcolato l'angolo:
                                    None    => Vengono usati tutti i metodi e valutati i risultati.
                                    1       => Viene usato solo il metodo 1 col prodotto scalare
                                    2       => Viene usato solo il metodo 2 col teorema di carnot
=================================================================================
RISULTATI
    :type numpy.array/Numero:         ==>     Angolo trai due vettori
    """
    l_1 = Grandezza_vettore(Vettore_1)  # Grandezza/lunghezza del primo vettore
    l_2 = Grandezza_vettore(Vettore_2)  # Grandezza/lunghezza del secondo vettore
    # Vettori opposti:
    if Grandezza_vettore(Vettore_1 + Vettore_2) <= pow(10, -6):
        return round(np.pi, Angl_dec)
    elif Grandezza_vettore(Vettore_1 - Vettore_2) <= pow(10, -6):
        return round(0.0, Angl_dec)
    else:
        # Metodo 1: prodotto scalare
        Prod_scalare = np.dot(Vettore_1, Vettore_2)
        cos_1 = Prod_scalare / (l_1 * l_2)
        ang_1 = round(np.arccos(cos_1), Angl_dec)  # Angolo cercato
        # Correzione in caso di valore nan
        if np.isnan(ang_1):
            ang_1 = np.nan_to_num(ang_1)
        # Metodo 2: teorema del coseno a^2 = b^2 + c^2 + 2*b*c*cos(alpha)
        v_3 = Vettore_2 - Vettore_1  # Vettore che chiude il triangolo
        l_3 = Grandezza_vettore(v_3)  # Grandezza/lunghezza del terzo vettore
        Numeratore = pow(l_3, 2) - pow(l_2, 2) - pow(l_1, 2)
        Denominatore = -2 * l_1 * l_2
        if Denominatore == 0.0:
            ang_2 = round(np.pi, Angl_dec)
        else:
            cos_2 = Numeratore / Denominatore
            ang_2 = round(np.arccos(cos_2), Angl_dec)  # Angolo cercato
            # Correzione in caso di valore nan
        if np.isnan(ang_2):
            ang_2 = np.nan_to_num(ang_2)
        if metodo is None:
            if ang_1 == ang_2:  # se i valori dei due metodi coincidono restituisci la media dei valori
                ang = np.mean((ang_1, ang_2))
                return ang
            else:  # Se i valori non sono uguali
                somma = round(ang_1 + ang_2, Angl_dec)
                # print(somma,round(np.pi,3),somma-round(np.pi,3))
                # se la loro somma e' un angolo giro restituisce il valore calcolato con il secondo metodo
                if abs(round(somma, 3) - round(np.pi, 3)) < pow(10,
                                                                0 - (Angl_dec + 1)):

                    return ang_2
                else:
                    raise Exception("Errore nel calcolo dell'angolo; i valori "
                                    "calcolati con i due metodi sono diversi e la loro somma non e' un angolo piatto.\n metodo 1 : %s\nmetodo 2 %s" % (
                                    ang_1, ang_2))
        elif metodo == 1:
            return ang_1
        elif metodo == 2:
            return ang_2
        else:
            raise TypeError


class Piano(object):
    """
    Questa classe crea un oggetto "piano"
    """

    def __init__(self, P_1, P_2, P_3):
        """
_L_'iniziazione di questa classe parte dal equazione parametrica del piano,
restituendo i punti e il versore direzionale e poi determina il vettore normale.
=================================================================================
VARIABILI
---------------------------------------------------------------------------------
    :param P_1:         ==>         Tipo: numpy.array
                                    primo punto del piano.

    :param P_2:         ==>         Tipo: numpy.array
                                    secondo punto del piano.

    :param P_3:         ==>         Tipo: numpy.array
                                    terzo punto del piano
=================================================================================
RISULTATO
    :parameter P_1:     ==>     Tipo: numpy.array
                                primo punto del piano.

    :parameter P_2:     ==>     Tipo: numpy.array
                                secondo punto del piano.

    :parameter P_3:     ==>     Tipo: numpy.array
                                terzo punto del piano

    :parameter Vers_1:      ==>         Tipo: numpy.array
                                Primo versore parametrico del piano.

    :parameter Vers_2:      ==>         Tipo: numpy.array
                                Secondo versore parametrico del piano.

    :parameter Vers_N:      ==>         Tipo: numpy.array
                                Versore direzionale del piano.

    :parameter d:           ==>         Tipo: Numero
                                Discriminante del piano
        """
        self.P1 = P_1  # Primo punto del piano. ( usato anche come origine locale del piano)
        self.P2 = P_2  # Secondo punto del piano.
        self.P3 = P_3  # Terzo punto del piano.
        self.Vers_1 = Versori(self.P2 - self.P1)  # Primo versore parametrico del piano.
        self.Vers_2 = Versori(self.P3 - self.P1)  # Primo versore parametrico del piano.
        # Versore direzionale del piano.
        _N = np.cross(self.Vers_1, self.Vers_2)
        if _N[2] < 0:
            self.Vers_N = (0 - 1) * Versori(_N)
        else:
            self.Vers_N = Versori(_N)
        # Discriminante del piano
        d = round(((self.Vers_N[0] * self.P1[0]) + (self.Vers_N[1] * self.P1[1]) + (
                self.Vers_N[2] * self.P1[2])), Vers_dec)
        if abs(d) <= (pow(10, (0 - 1) * 10)):
            d = 0.0
        self.d = d

    def Basi_ortonormali(self, Base=1):
        """
    Funzione per avere le basi orto normali di un piano.
    =================================================================================
    ARGOMENTI
    ---------------------------------------------------------------------------------
        :parameter Base:    ==>     Tipo: Numero intero
                                    Indica il versore parametrico che si vuole usare
                                    come partenza.
                                    I valori 1 e 0 indicano il primo versore.
                                    I valori 2 e superiori indicano il secondo versore.
    =================================================================================
    RISULTATI
    :return:        ==>     Tipo: numpy.array
                            La funzione restituisce un dizionario nel quale
                            sono salvate le basi nel modo seguente:
                            Basi = np.matrix() = matrix([[riga1],[riga 2], ... ])
        """
        if Base > 1:
            B_1 = self.Vers_2
            if np.dot(B_1, self.Vers_1) == 0:
                B_2 = self.Vers_1
            else:
                # Assicura che la coppia di vettori sia destrogira
                B_2 = np.cross(B_1, self.Vers_N)
                if Angoli_tra_vettori(B_1, B_2) > np.round(np.pi / 2.0,
                                                           Angl_dec):
                    B_2 = np.cross(self.Vers_N, B_1)
            B_2 = np.round(B_2, Vers_dec)
        else:
            B_1 = self.Vers_1
            if np.dot(B_1, self.Vers_2) == 0:
                B_2 = self.Vers_2
            else:
                # Assicura che la coppia di vettori sia destrogira.
                B_2 = np.cross(B_1, self.Vers_N)
                if Angoli_tra_vettori(B_1, B_2) > np.round(np.pi / 2.0,
                                                           Angl_dec):
                    B_2 = np.cross(self.Vers_N, B_1)
            B_2 = np.round(B_2, Vers_dec)
        B = []
        for ele in range(len(B_2)):
            B.append([B_1[ele], B_2[ele]])
        return np.array(B)

    def Punto_by_Param(self, t=None, s=None, base=None):
        """
    Funzione per calcolare le coordinate di un punto generico appartenente alla retta;
    usando l'equazione parametrica della retta.
    =================================================================================
    ARGOMENTI
    ---------------------------------------------------------------------------------
        :param t:       ==>     Tipo: Numero
                                Valore del parametro.

        :param s:       ==>     Tipo: Numero
                                Valore del parametro.

        :param base:      ==>   Tipo: numpy.array
                                Indica la base che si vuole usare per calcolare i
                                parametri.
                                La base puo' essere come indicata come una matrice
                                le cui colonne sono i vettori oppure come lista dei
                                vettori che costituiscono la base.
    =================================================================================
    RETURNING
        :type numpy.array:      Vettore che descrive il punto.
        """
        # determinazione Parametro.
        if t is None:
            par_t = Rnd.randrange(-100, 100)
        else:
            par_t = t
        if s is None:
            par_s = Rnd.randrange(-100, 100)
        else:
            par_s = s
        # Scelta della Base
        if type(base) == type([]):
            v_1 = base[0]
            v_2 = base[1]
        elif type(base) == type(np.array([0])):
            # il doppio indice serve ad ottenere un array semplice.
            v_1 = np.array(base.T[0])[0]
            v_2 = np.array(base.T[1])[0]
        else:
            v_1 = self.Vers_1
            v_2 = self.Vers_2
        # Calcolo del punto
        P = np.round(self.P1 + (par_t * v_1) + (par_s * v_2), Cord_dec)
        return P

    def Punto_by_Cord(self, Cord_val, Cord=''):
        """
Restituisce un punto partendo da una determinata coppia di coordinate.
Se non si sceglie un determinato piano, il punto viene calcolato usando la
valutazione parametrica.
=================================================================================
ARGOMENTI
---------------------------------------------------------------------------------
    :param Cord:        ==>     Tipo: Stringa
                                Stringa indicante la coordinata scelta; se non
                                assegnata il valore predefinito sara' XY
                                valori accettati: XY, YZ ,ZX.
                                le altre lettere verrano trattatee come
                                parametro generico

    :param Cord_val:   ==>      Tipo: numpy.array
                                lista dei valori delle coordinate scelte.
=================================================================================
RETURNING
    :returns numpy.array:       Vettore che descrive la posizione del punto
        """
        Cord = str(Cord).lower()
        if Cord == 'yx' or Cord == 'xy':
            # Calcolo Coordinata Z
            z = ((0 - 1) * (self.d + (self.Vers_N[0] * Cord_val[0]) + (self.Vers_N[1] * Cord_val[1]))) / self.Vers_N[2]
            z = round(z, Cord_dec)
            return np.array([Cord_val[0], Cord_val[1], z])
        elif Cord == 'yz' or Cord == 'zy':
            # Calcolo Coordinata x
            x = ((0 - 1) * (self.d + (self.Vers_N[1] * Cord_val[0]) + (self.Vers_N[2] * Cord_val[1]))) / self.Vers_N[0]
            x = round(x, Cord_dec)
            return np.array([x, Cord_val[0], Cord_val[1]])
        elif Cord == 'zx' or Cord == 'xz':
            # Calcolo Coordinata Y
            y = ((0 - 1) * (self.d + (self.Vers_N[0] * Cord_val[0]) + (self.Vers_N[2] * Cord_val[1]))) / self.Vers_N[1]
            y = round(y, Cord_dec)
            return np.array([Cord_val[0], y, Cord_val[1]])
        else:  # valutazione parametrica pura
            t = Cord_val[0]
            s = Cord_val[1]
            return self.Punto_by_Param(t, s, self.Basi_ortonormali())

    def Parametri_punto(self, punto, base=None):
        """
        Questa funzione permette di calcolare i parametri che corrispondono al punto secondo una deteminata
         base.
        ARGOMENTI:
            punto:              punto di cui si vogliono calcolare i parametri rispetto la piano.
                                Il tipo deve essere un array o una lista colonna.

            base:               Lista dei vettori che costituiscono la base su cui calcolare i parametri.
                                I vettori sono array; se non viene indicata nessuna base viene usata una generica base ortonormale.
        RISULTATI:
            Lista dei paramtri.
        """
        # Decisione su quale base usare.
        if type(base) == type(None):
            Base = self.Basi_ortonormali()
        else:
            B = np.array(base)
            Base = B.T
        # Calcolare i parametri.
        Base_xy = np.array(Base.A[:2])  # Matrice che considera la prima e la seconda equazione.
        v_noto = punto - self.P1  # Vettore dei termini noti.
        if np.linalg.det(Base_xy) != 0:
            parametri = np.linalg.solve(Base_xy, v_noto[:2])
        else:
            Base_yz = np.array(Base.A[1:3])
            if np.linalg.det(Base_yz) != 0:
                parametri = np.linalg.solve(Base_yz, v_noto[1:3])
            else:
                Base_zx = np.array(Base.A[:3:2])
                if np.linalg.det(Base_zx) != 0:
                    parametri = np.linalg.solve(Base_zx, v_noto[:3:2])
                else:
                    raise Exception
        return parametri

    def Inclinazione(self, V_ref=None):
        """
        Questa funzione calcola l'inclinazione del piano rispetto a un piano
        orizontale se non viene specificato il vettore di un piano particolare.
        Variabili:
            self        ==>     piano [ oggetto ]
            V_ref       ==>     Vettore del piano di riferimento [ np.array() ]
        Risultati:
            Inc         ==>     Angolo [ numero ]
        """
        if V_ref is None:
            V_z = [0.0, 0.0, 1.0]  # Versore asse z/ Versore direzionale asse orizontale
        else:
            V_z = V_ref
        V_p = self.Vers_N  # Versore faccia
        cos = np.dot(V_z, V_p) / (Grandezza_vettore(V_z) * Grandezza_vettore(V_p))  # coseno dell'angolo di inclinazione
        Ang_c = round(np.arccos(cos), Angl_dec)
        return Ang_c

    def Sposta(self, Vettore_spostamento):
        """
        Questa funzione sposta il piano secondo il vettore direzione assegnato.
        Argomenti:
            Vettore_spostamento =>      Il vettore che definisce lo spostamento lungo gli assi x,y e z.
                                        Puo' essere sia una array del modulo numpy sia una semplice lista.
        Risultato:
            Plane
        """
        P1 = self.P1 + Vettore_spostamento
        P2 = self.P2 + Vettore_spostamento
        P3 = self.P3 + Vettore_spostamento
        piano_spostato = Piano(P1, P2, P3)
        return piano_spostato

    def Ruota(self, Angl, Cen_Rot=None, asse_=None, un='rad'):
        """
        Questa funzione permette di eseguire rotazioni nello spazio. La rotazione
        puo' essere eseguita attorni agli assi di una base canonica, rispetto a un
        punto e rispetto a un asse particolare.
        ARGOMENTI:
            Angl      =>      Angoli di rotazione rispetto agli assi canonici
                                                np.array() = [a,b,c]
                                    a = angolo di rotazione intorno all'asse x
                                    b = angolo di rotazione intorno all'asse y
                                    c = angolo di rotazione intorno all'asse z
                                se viene assegnato un solo valore insieme all'asse
                                di rotazione si intende come angolo di rotazione
                                intorno all'asse indicato.
        Cen_Rot    => Punto intorno a cui avviene la rotazione [ np.arra() ]
        asse_       =>      asse intorno a cui avviene la rotazione
        unit        =>      unita di misura degli angoli
                                    il valore e' impostato su radianti
                                    'dec' per decimali
                                    'cen' per centesimali
                                    'rad' per radianti
        """
        P_1 = Rotazione(self.P1, Angl, Centro_Rotazione=Cen_Rot, asse_rotazione=asse_, unit=un)
        P_2 = Rotazione(self.P2, Angl, Centro_Rotazione=Cen_Rot, asse_rotazione=asse_, unit=un)
        P_3 = Rotazione(self.P3, Angl, Centro_Rotazione=Cen_Rot, asse_rotazione=asse_, unit=un)
        return self.__init__(P_1, P_2, P_3)

    def Cordinate_Parametriche(self, punto, base=None):
        """
        Questa funzione permette di trasformare i paramtri del piano in coordinate rispetto al piano.
        Argomenti:
            punto ==> punto da cui ottenere le cordinate parametriche.
            base:               Lista dei vettori che costituiscono la base su cui calcolare i parametri.
                                I vettori sono array; se non viene indicata nessuna base viene usata una generica base ortonormale.
        """
        # Calcolo parametri
        Cord_plan_v1 = list(self.Parametri_punto(punto, base))
        Cord_plan_v1.append(0.0)
        return np.array(Cord_plan_v1)

    def __add__(self, others):
        """
        Addizione di due rette crea una nuova retta
        """
        if type(self) == type(others):
            punto_1 = self.P1 + others.P1
            punto_2 = self.P2 + others.P2
            punto_3 = self.P3 + others.P3
            return Piano(punto_1, punto_2, punto_3)
        else:
            TypeError('operazione non supportata')

    def __sub__(self, others):
        """
        Addizione di due rette crea una nuova retta
        """
        if type(self) == type(others):
            punto_1 = self.P1 - others.P1
            punto_2 = self.P2 - others.P2
            punto_3 = self.P3 - others.P3
            return Piano(punto_1, punto_2, punto_3)
        else:
            TypeError('operazione non supportata')

    def __str__(self):
        """
        Restituisce una stringa che rappresenta il piano nello spazio quando richiamata dalla funzione print.
        """
        P = self.P1
        V = [self.Vers_1, self.Vers_2]
        Cord = ['X', 'Y', 'Z']
        stringa = "_L_'equazione parametrica del piano e':\n"
        for m in range(len(Cord)):
            stringa = stringa + ("\t \t %s = %s + u(%s) + v(%s) \n" % (Cord[m], P[m], V[0][m], V[1][m]))
        return stringa

    def __repr__(self):
        """
        REstituisce la rappresentazione della piano sottoforma di stringa.
        """
        Stringa = 'Equazione Cartesiana del piano:\n'  # Titolo della Stringa
        Stringa = Stringa + '\t %s*x + %s*y + %s*z + %s' % (self.Vers_N[0], self.Vers_N[1], self.Vers_N[2], self.d)
        return Stringa

    def __eq__(self, others):
        """
        uguaglianza ( == )
        """
        if type(self) == type(others):
            if sum(self.Vers_N == others.Vers_N) == 3:
                if self.d == others.d:
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False

    def __ne__(self, others):
        """
        disuguaglianza ( != )
        """
        if type(self) != type(others):
            if sum(self.Vers_N != others.Vers_N) == 3:
                if self.d != others.d:
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False


class Retta(object):
    """
    Questa classe crea un oggetto "retta"
    """

    def __init__(self, P_1, P_2):
        """
        _L_'inizazione di questa classe parte dal equazione parametrica della retta,
        restituendo il punto usato come rigine e il versore direzionale.
        VARIABILI
        P_1     ==>     primo
        P_2     ==>     secondo punto della retta.
        Risultati
        self.Vers       ==>         Versore direzionale della Retta.
                                        relative della retta
        """
        L_ = [P_1, P_2]
        l = Ordina_punti(L_)
        self.P1 = np.array(l[0])  # Primo punto della retta.
        self.P2 = np.array(l[1])  # Secondo punto della retta.
        # Versore direzionale della Retta.
        v = (self.P2 - self.P1)
        if v[2] < 0:
            self.Vers = Versori(v) * (0 - 1)
        else:
            self.Vers = Versori(v)
        # Punto posto sul piano z = 0 o sul piano y = 0 o sul piano x = 0 a seconda dell'inclianzione della retta
        if self.Vers[2] >= pow(10, -Vers_dec):  # Plane z = 0
            t_0 = (0 - 1) * (self.P1[2] / self.Vers[2])
        elif self.Vers[1] >= pow(10, -Vers_dec):  # Plane y = 0
            t_0 = (0 - 1) * (self.P1[1] / self.Vers[1])
        elif self.Vers[0] >= pow(10, -Vers_dec):  # Plane x = 0
            t_0 = (0 - 1) * (self.P1[0] / self.Vers[0])
        else:
            t_0 = 0.0
        P0 = self.P1 + (t_0 * self.Vers)
        self.P0 = np.round(P0, Cord_dec)

    def Eq_cart_Plane_xy(self):
        """
        Questa funzione permette di calcolare i parametri "m" e "q" per l'equazione
        della retta nel piano scritta nella forma:
                                Y = m*X + q
        VARIABILI
        P_1     ==>     Array delle cordinate del punto 1
        P_2     ==>     Array delle cordinate del punto 2
        Risultati
        m      ==>         Coeficiente angolare sul piano xy
        q      ==>         intercetta sull'asse "Y"
        """
        if (self.P2[0] - self.P1[0]) == 0.00:
            m = float('inf')
        else:
            m = round(((self.P2[1] - self.P1[1]) / (self.P2[0] - self.P1[0])), precisione)
        q = round((self.P2[1] - (m * self.P2[0])), precisione)
        return [m, q]

    def Eq_cart_(self):
        """
        Funzione per avere l'equazione cartesiana della retta nello spazio data
        dall'intersezione di due piani.

        Plane 1 ==>     (a_1 * x) + (b_1 * y) + (c_1 * z) + d_1 = 0
        Plane 2 ==>     (a_2 * x) + (b_2 * y) + (c_2 * z) + d_2 = 0

        Dove i coeficienti ( a,b,c ) rappresentono il vettore di quel piano mentre
        i coeficienti d sono i discriminanti.
        La funzione restituisce una lista contenente i due piani calcolati.
        I piani calcolati sono due degli infiniti possibili piani; i piani selezionati hanno la caratteristica di essere ortogonali tra loro.
        """
        if self.Vers[2] != 0:  # Il componente lungo l'asse z non nullo.
            v1_x = 1  # componente x del vettore del piano 1
            v1_y = 1  # componente y del vettore piano 1
            v1_z = ((self.Vers[0] * v1_x) + (self.Vers[1] * v1_y)) * (
                        (0 - 1) / self.Vers[2])  # componente z del vettore del piano 1
            v1 = np.array([v1_x, v1_y, v1_z])
        elif self.Vers[1] != 0:
            v1_x = 1  # componente x del vettore del piano 1
            v1_z = 1  # componente z del vettore del piano 1
            v1_y = ((self.Vers[0] * v1_x) + (self.Vers[2] * v1_z)) * (
                        (0 - 1) / self.Vers[1])  # componente y del vettore piano 1
            v1 = np.array([v1_x, v1_y, v1_z])
        else:
            v1_y = 1  # componente y del vettore del piano 1
            v1_z = 1  # componente z del vettore del piano 1
            v1_x = ((self.Vers[1] * v1_y) + (self.Vers[2] * v1_z)) * (
                        (0 - 1) / self.Vers[0])  # componente x del vettore piano 2
            v1 = np.array([v1_x, v1_y, v1_z])
        v2 = np.cross(self.Vers, v1)  # Vettore per il secondo piano
        P1 = self.P0 + v1  # Terso punto per creare il piano 1
        P2 = self.P0 + v2  # Terso punto per creare il piano 2
        piano_1 = Piano(P1, self.P1, self.P2)
        piano_2 = Piano(P2, self.P1, self.P2)
        return [piano_1, piano_2]

    def Punto_by_Param(self, t=None):
        """
        Funzione per calcolare le coordinate di un punto generico appartenente alla retta;
        usando l'equazione paramtrica della retta.
        VRIABILI.
        t       ==>     Indica un preciso valore del parametro per cui si vuole calcolare il punto.
        """
        # determinazione Parametro.
        if t is None:
            par_t = Rnd.randrange(-100, 100)
        else:
            par_t = t
        # Calcolo delle Coordinate
        P = self.P1 + (par_t * self.Vers)
        np.round(P, Cord_dec)
        return P

    def Punto_by_Cord(self, Cord_val, Cord=''):
        """
        Restituisce un punto partendo da una determinata coordinata.
        VARIABILI.
        Cord        ==>         Stringa indicante la coordinata scelta; se non
                                assegnata il valore predefinito sara' X
                                valori acettati: X, x , Y ,y , Z, z,
                                le altre lettere verrano tratate come paramtro generico
        Cord_val    ==>         Valore della cordinata scelta.
        """
        Cord = str(Cord).lower()
        # Evita la divisione per zero non usando le componenti nulle
        if self.Vers[0] != 0.0:
            vers = self.Vers[0]
        elif self.Vers[1] != 0.0:
            vers = self.Vers[1]
        else:
            vers = self.Vers[2]
        # Calcolo del parametro da usare nell'equazione parametrica
        if Cord == '' or Cord == 'x':  # Parametro usando la X
            t = (Cord_val - self.P1[0]) / vers
        elif Cord == 'y':  # Parametro usando la Y
            t = (Cord_val - self.P1[1]) / vers
        elif Cord == 'z':  # Parametro usando la Z.
            t = (Cord_val - self.P1[2]) / vers
        else:  # valutazione parametrica pura
            t = Cord_val
        return self.Punto_by_Param(t)

    def Param_Punto(self, punto):
        """
        Questa funzione retituisce il parametro relativo ad un punto.
        ARGOMENTI:
            punto   =>      Punto di cui determinare il parametro.
                                [ np.array() ]
        """
        T = []
        for l in range(len(self.Vers)):
            if self.Vers[l] > pow(10, -Vers_dec):
                t = round((punto[l] - self.P0[l]) / self.Vers[l], precisione)
                T.append(t)
        T.sort()
        t = np.mean(T)
        return t

    def Sposta(self, Vettore_spostamento):
        """
        Questa funzione sposta la retta secondo il vettore direzione assegnato.
        Argomenti:
            Vettore_spostamento =>      Il vettore che definisce lo spostamento lungo gli assi x,y e z.
                                        Pu? essere sia una array del modulo numpy sia una semplice lista.
        Risultato:
            Retta
        """
        p_1 = self.P1 + Vettore_spostamento
        p_2 = self.P2 + Vettore_spostamento
        return self.__init__(p_1, p_2)

    def Ruota(self, Angl, Cen_Rot=None, asse_=None, un='rad'):
        """
        Questa funzione permette di eseguire rotazioni nello spazio. La rotazione
        puo' essere eseguita attorni agli assi di una base canonica, rispetto a un
        punto e rispetto a un asse particolare.
        ARGOMENTI:
            Angl      =>      Angoli di rotazione rispetto agli assi canonici
                                                np.array() = [a,b,c]
                                    a = angolo di rotazione intorno all'asse x
                                    b = angolo di rotazione intorno all'asse y
                                    c = angolo di rotazione intorno all'asse z
                                se viene assegnato un solo valore insieme all'asse
                                di rotazione si intende come angolo di rotazione
                                intorno all'asse indicato.
        Cen_Rot    => Punto intorno a cui avviene la rotazione [ np.arra() ]
        asse_       =>      asse intorno a cui avviene la rotazione
        unit        =>      unita di misura degli angoli
                                    il valore e' impostato su radianti
                                    'dec' per decimali
                                    'cen' per centesimali
                                    'rad' per radianti
        """
        P_1 = Rotazione(self.P1, Angl, Centro_Rotazione=Cen_Rot, asse_rotazione=asse_, unit=un)
        P_2 = Rotazione(self.P2, Angl, Centro_Rotazione=Cen_Rot, asse_rotazione=asse_, unit=un)
        return self.__init__(P_1, P_2)

    def __add__(self, others):
        """
        Addizione di due rette crea una nuova retta (+)
        """
        if type(self) == type(others):
            punto_1 = self.P1 + others.P1
            punto_2 = self.P2 + others.P2
            return Retta(punto_1, punto_2)
        else:
            TypeError('operazione non supportata')

    def __sub__(self, others):
        """
        Sottrazione di due rette crea una nuova retta (-)
        """
        if type(self) == type(others):
            punto_1 = self.P1 - others.P1
            punto_2 = self.P2 - others.P2
            return Retta(punto_1, punto_2)
        else:
            TypeError('operazione non supportata')

    def __str__(self):
        """
        REstituisce la rappresentazione della reta sotto forma di stringa che si ottiene usando il comando print.
        """
        P = self.P1
        V = [self.Vers[0], self.Vers[1], self.Vers[2]]
        Cord = ['X', 'Y', 'Z']
        stringa = "_L_'equazione parametrica della retta e':\n"
        for m in range(len(Cord)):
            stringa = stringa + ("\t \t %s = %s + u(%s) \n" % (Cord[m], P[m], V[m]))
        return stringa

    def __repr__(self):
        """
        Restituisce la rappresentazione della retta sotto forma di stringa.
        """
        Piani_della_retta = self.Eq_cart_()  # Calcola l'equazione cartesiana della retta
        Stringa = 'Equazione Cartesiana della retta\n'  # Stringa che rappresenta la retta
        Piano_1 = [Piani_della_retta[0].Vers_N, Piani_della_retta[0].d]  # Versore del piano
        Piano_2 = [Piani_della_retta[1].Vers_N, Piani_della_retta[1].d]  # Versore del piano
        Stringa = Stringa + 'Plane 1 : ' + str(Piano_1[0][0]) + '*x + ' + str(Piano_1[0][1]) + '*y + ' + str(
            Piano_1[0][2]) + '*z + ' + str(Piano_1[1]) + '\n'
        Stringa = Stringa + 'Plane 2 : ' + str(Piano_2[0][0]) + '*x + ' + str(Piano_2[0][1]) + '*y + ' + str(
            Piano_2[0][2]) + '*z + ' + str(Piano_2[1])
        return Stringa

    def __eq__(self, others):
        """
        uguaglianza ( == )
        """
        if type(self) == type(others):
            if all([self.P0[n] == others.P0[n] for n in range(len(self.P0))]):
                if sum(self.Vers == others.Vers) == 3:
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False

    def __ne__(self, others):
        """
        disuguaglianza ( != )
        """
        if self.__eq__(others):
            return False
        else:
            return True


### Creazione degli assi coordinati come rette
asse_x = Retta([0.0, 0.0, 0.0], [1.0, 0.0, 0.0])
asse_y = Retta([0.0, 0.0, 0.0], [0.0, 1.0, 0.0])
asse_z = Retta([0.0, 0.0, 0.0], [0.0, 0.0, 1.0])

### Piani fondamentali
piano_xy = Piano(O, [10.0, 0.0, 0.0], [0.0, 10.0, 0.0])
piano_yz = Piano(O, [0.0, 10.0, 0.0], [0.0, 0.0, 10.0])
piano_zx = Piano(O, [10.0, 0.0, 0.0], [0.0, 0.0, 10.0])


def Punto_by_3_piani(Pia1, Pia2, Pia3):
    """
    Questa funzione restituisce il punto di intersezione tra tre piani.
    Variabili:
        Pia1        ==>     Plane numero 1  [ oggetto ]
        Pia2        ==>     Plane numero 2  [ oggetto ]
        Pia3        ==>     Plane numero 3  [ oggetto ]
    Risultati:
        P           ==>     Punto di intersezione [ oggetto ]
    """
    A = np.array([Pia1.Vers_N, Pia2.Vers_N, Pia3.Vers_N])
    B = np.array([Pia1.d, Pia2.d, Pia3.d])
    det_A = np.linalg.det(A)
    if det_A == 0.0:
        return None
    else:
        P = np.linalg.solve(A, B)
        P = np.round(P, Cord_dec)
        return P


def Punto_by_2Retta(r, s):
    """
    Questa funzione restituisce il punto in tersezione tra due rette incidenti.
    Variabili:
    s   ==>     Prima retta [ oggetto ]
    r   ==>     Seconda retta [ oggetto ]
    Risultati:
    P   ==>     Punto di intersezione [ np.array() ]
    """
    if r != s:
        V = r.P1 - s.P1  # Vettore da un punto di r a un punto di s
        rango_matrice_dei_versori = np.linalg.matrix_rank(np.array([r.Vers, s.Vers]))
        rango_matrice_completa = np.linalg.matrix_rank(np.array([r.Vers, s.Vers, V]))
        if rango_matrice_completa == 2:
            if rango_matrice_dei_versori == 2:
                R = r.Eq_cart_()  # piani per r
                S = s.Eq_cart_()  # piani per s
                # 2 piani da r e il primo da S
                P_1 = Punto_by_3_piani(R[0], R[1], S[0])
                if P_1 is None:
                    P_1 = 0
                # 2 piani da r e il secondo da S
                P_2 = Punto_by_3_piani(R[0], R[1], S[1])
                if P_2 is None:
                    P_2 = 0
                # 2 piani da s e  il primo da r
                P_3 = Punto_by_3_piani(R[0], S[1], S[0])
                if P_3 is None:
                    P_3 = 0
                # 2 piani da s e il secondo da r
                P_4 = Punto_by_3_piani(S[0], R[1], S[1])
                if P_4 is None:
                    P_4 = 0
                P = (1 / 4.0) * (P_1 + P_2 + P_3 + P_4)
                P = np.round(P, Cord_dec)
                return P
    else:
        print('Le rette sono coincidenti')
        return None


def Punto_by_plan_e_retta(retta, piano):
    """
    Questa funzione restituisce il punto di intersezione tra una retta e un piano.
    Variabili:
        retta      ==>     Retta [ oggetto ]
        piano    ==>     Plane [ oggetto ]
    Risultati:
        Punto   ==>     punnto [ np.array() ]
    """
    # Piani della retta
    Piani_r_0 = retta.Eq_cart_()[0]
    Piani_r_1 = retta.Eq_cart_()[1]
    punto = Punto_by_3_piani(piano, Piani_r_0, Piani_r_1)
    return punto


def Retta_by_P_e_Vet(P, V):
    """
    Questa funzione crea un retta partendodalsuo vettore direzionale e un punto.
    Variabili:
        P       ==>     Punto [ np.array() ]
        V       ==>     vettore direzionale [ np.arra() ]
    Risultati:
        r       ==>     Retta [ oggetto ]
    """
    Vers = Versori(V)  # Versore direzionale della retta
    P_2 = P + Vers
    return Retta(P, P_2)


def Retta_by_2Planes(piano_1, piano_2):
    """
    Questa funzione restituisce la retta generata dall'intersezione di due piani.
    Variabili:
    Pl_1    ==>     Primo piano [ oggetto ]
    Pl_2    ==>     Secondo piano [ oggetto ]
    Risultato:
    Retta   ==>     Retta inersezione tra i due piani [ oggetto ]
    """
    v_1 = piano_1.Vers_N  # Versore direzionale del primo piano
    v_2 = piano_2.Vers_N  # Versore direzionale del secondo piano
    v_r = Versori(np.cross(v_1, v_2))  # versore direzionale della retta
    A = np.array([v_1, v_2, v_r])  # Matrice dei coefficienti del sistema
    if np.linalg.det(A) != 0.00:  # se il sistema ha soluzione
        B_0 = np.array([piano_1.d, piano_2.d, 0])  # vettore termini noti per un il primo punto
        P0 = np.linalg.solve(A, B_0)
        P0 = P0.round(Cord_dec)  # Primo punto
        B_1 = np.array([piano_1.d, piano_2.d, 2])  # vettore termini noti per un il secondo punto
        P1 = np.linalg.solve(A, B_1)
        P1 = P1.round(Cord_dec)  # Secondo punto
        return Retta(P0, P1)
    else:
        return None


def Retta_by_piano_orto_e_punto(piano, punto):
    """
    Questa funzione restituisce la retta ortogonale al piano dato e passante per il
    punto assegnato.
    Variabili:
        piano     ==>         Plane [ ogetto ]
        punto   ==>         Punto [ np.array() ]
    Risultato:
        r       ==>         retta [ oggetto ]
    """
    Pia = None
    P = None
    # Ordinazione degli argomenti
    if type(punto) == type(np.zeros((3, 1))):
        if type(piano) == type(piano_xy):
            P = punto
            Pia = piano
    elif type(punto) == type(piano_xy):
        if type(piano) == type(np.zeros((3, 1))):
            P = piano
            Pia = punto
    else:
        raise TypeError()
    V_r = Pia.Vers_N  # Versore della retta.
    P_2 = P + V_r  # secondo punto retta ottenuto dall'equazione parametrica P = P_or + t*V con t = 1
    return Retta(P, P_2)


def Retta_by_retta_orto_eP(retta, punto):
    """
    Questa funzione restituisce una retta ortogonale alla retta data e passante
    per il punto scelto.
    Variabili:
        retta       ==>         Retta assegnata [ oggetto ]
        punto       ==>         Punto [ np.array() ]
    Risultati:
        s       ==>         Retta cercata [ oggetto ]
    """
    P = None
    r = None
    # Ordinazione degli argomenti
    if type(punto) == type(np.zeros((3, 1))):
        if type(retta) == type(asse_x):
            r = retta
            P = punto
    elif type(punto) == type(asse_x):
        if type(retta) == type(np.zeros((3, 1))):
            r = punto
            P = retta
    else:
        raise TypeError()
    piano = Piano(P, r.P1, r.P2)  # Plane che contiene sia la retta che il punto
    v = np.cross(piano.Vers_N, r.Vers)  # Vettore della retta
    s = Retta_by_P_e_Vet(P, v)  # Creazione della retta
    return s


def Piano_by_retta_orto_e_Punto(punto, retta):
    """
    Questa funzione crea un piano partendo da un punto appartenente al piano e
    una retta ortogonale al piano.
    Variabili:
        punto       ==>     Punto [ np.array() ]
        retta       ==>     Retta [ oggetto ]
    Risltati:
        plan    ==>     Plane [ oggetto ]
    """
    r = None
    P = None
    # Ordinazione degli argomenti
    if type(punto) == type(np.zeros((3, 1))):
        if type(retta) == type(asse_x):
            r = retta
            P = punto
    elif type(punto) == type(asse_x):
        if type(retta) == type(np.zeros((3, 1))):
            r = punto
            P = retta
    else:
        raise TypeError()
    Plan_point = []  # Punti del piano per poter ottenere il piano.
    n = 0
    if r.Vers[2] != 0:
        while n < 2:
            x = Rnd.randrange(-20, 20)
            y = Rnd.randrange(-20, 20)
            z = (0 - 1) * (((r.Vers[0] * x) + (r.Vers[1] * y)) / r.Vers[2])
            Plan_point.append(np.array([x, y, z]))
            n = n + 1
    else:
        if r.Vers_N[1] != 0:
            while n < 2:
                x = Rnd.randrange(-20, 20)
                z = Rnd.randrange(-20, 20)
                y = (0 - 1) * (((r.Vers[0] * x) + (r.Vers[2] * z)) / r.Vers_N[1])
                Plan_point.append(np.array([x, y, z]))
        else:
            while n < 2:
                y = Rnd.randrange(-20, 20)
                z = Rnd.randrange(-20, 20)
                x = (0 - 1) * (((r.Vers_N[1] * y) + (r.Vers_N[2] * z)) / r.Vers_N[0])
                Plan_point.append(np.array([x, y, z]))  # vettori ortogonali al vettore del piano (basi)
    return Piano(P, P + Plan_point[0], P + Plan_point[1])


def Piano_by_2Rette(r, s):
    """
    Questa funzione crea un piano partendo da due rette giacenti nel piano stesso.
    Variabili:
        r       ==>     Retta [ oggetto ]
        s       ==>     Retta [ oggetto ]
    Risultati:
        Plan    ==>     Plane [ oggetto ]
    """
    if all([r.Vers[n] == s.Vers[n] for n in range(len(r.Vers))]):  # Le rette sono parallele
        if all([r.P0[n] == s.P0[n] for n in range(len(r.P0))]):  # Le rette sono coincidenti
            return None  # Esistono infiniti piani quindi non possiamo sceglierne uno
        else:  # Le rette sono semplicemente parallele ma non coincidono
            P_1 = r.P1
            P_2 = r.P2
            P_3 = s.P1
            return Piano(P_1, P_2, P_3)
    elif type(Punto_by_2Retta(s, r)) == type(np.zeros((3, 1))):  # se le rette sono incidenti hanno un punto in comune.
        P_1 = r.P1
        P_2 = Punto_by_2Retta(r, s)
        P_3 = s.P2
        return Piano(P_1, P_2, P_3)
    else:  # Le rette sono sghembe
        return None


def Piano_by_P_e_Retta_gia(punto, retta):
    """
    Questa funzione crea un piano partendo da un punto e una retta giacenti nel
    piano.
    Variabili:
        punto       ==>         Punto [ np.array() ]
        retta   ==>         Retta [ oggetto ]
    Risultati:
        Plan    ==>         Plane [ oggetto ]
    """
    r = None
    P = None
    # Ordinazione degli argomenti
    if type(punto) == type(np.zeros((3, 1))):
        if type(retta) == type(asse_x):
            r = retta
            P = punto
    elif type(punto) == type(asse_x):
        if type(retta) == type(np.zeros((3, 1))):
            r = punto
            P = retta
    else:
        raise TypeError()
    P_1 = r.Punto_by_Param()
    P_2 = r.Punto_by_Param()
    return Piano(P, P_1, P_2)


def Piano_by_Vec_e_P(V, P):
    """
    Questa funzione crea un piano partendo dal suo vettore direzionale e un punto.
    Variabili:
        P       ==>     Punto [ np.array() ]
        V       ==>     vettore direzionale [ np.arra() ]
    Risultati:
        Plan    ==>     Plane [ oggetto ]
    """
    Plan_point = []  # Punti del piano per poter ottenere il piano.
    n = 0
    if V[2] != 0:
        while n < 2:
            x = Rnd.randrange(-20, 20)
            y = Rnd.randrange(-20, 20)
            z = (0 - 1) * (((V[0] * x) + (V[1] * y)) / V[2])
            Plan_point.append(np.array([x, y, z]))
            n = n + 1
    else:
        if V[1] != 0:
            while n < 2:
                x = Rnd.randrange(-20, 20)
                z = Rnd.randrange(-20, 20)
                y = (0 - 1) * (((V[0] * x) + (V[2] * z)) / V[1])
                Plan_point.append(np.array([x, y, z]))
        else:
            while n < 2:
                y = Rnd.randrange(-20, 20)
                z = Rnd.randrange(-20, 20)
                x = (0 - 1) * (((V[1] * y) + (V[2] * z)) / V[0])
                Plan_point.append(np.array([x, y, z]))
    return Piano(P, P + Plan_point[0], P + Plan_point[1])


def Dist_Punto_Retta(punto, retta):
    """
    Questa funzione restituisce la distanza tra un punto e una retta.
    Variabili:
        punto       ==>     Punto [ np.array() ]
        retta   ==>     Retta [ oggetto ]
    Risultati:
        Distanza ==>    [ numero ]
    """
    P = None
    r = None
    # Ordinazione degli argomenti
    if type(punto) == type(np.zeros((3, 1))):
        if type(retta) == type(asse_x):
            r = retta
            P = punto
    elif type(punto) == type(asse_x):
        if type(retta) == type(np.zeros((3, 1))):
            r = punto
            P = retta
    else:
        raise TypeError()
    # Metodo basato sulla retta ortogonale passante per il punto P
    Vettore_P0_P = P - r.P0  # Vettore da P0 della retta al punto P
    Vettore_P0_P1 = np.dot(Vettore_P0_P, r.Vers) * r.Vers  # Vettore da P0 alla proiezione di P sulla retta
    P_r = r.P0 + Vettore_P0_P1  # Coordinate punto proiezione
    Dis = Dist(P, P_r)  # Distanza tra il punto P e la sua proiezione sulla retta r.
    Distanza = round(Dis, precisione)
    return Distanza


def Dist_Punto_piano(punto, piano):
    """
    Funzione che restituisce la distanza tra un punto e un piano nello spazio
    Variabili:
        Punto   ==>     punto nello spazio [ array ]
        Plane   ==>     Plane nello spazio [ ogetto ]
    Risultato:
        Dist    ==>     Distanza tra il punto e la retta [ numero ]
    """
    Pia = None
    P = None
    # Ordinazione degli argomenti.
    if type(punto) == type(np.zeros((3, 1))):
        if type(piano) == type(piano_xy):
            P = punto
            Pia = piano
    elif type(punto) == type(piano_xy):
        if type(piano) == type(np.zeros((3, 1))):
            P = piano
            Pia = punto
    else:
        raise TypeError()
    # Funzione vera e propria.
    V_piano = np.round(Pia.Vers_N, Vers_dec)  # versore del piano
    r = Retta_by_P_e_Vet(P, V_piano)  # REtta ortogonale al piano e passante per P
    P1 = np.round(Punto_by_plan_e_retta(r, Pia), Cord_dec)  # Proiezione ortogonale di P sul piano
    dist_1 = round(Dist(P, P1), precisione)
    return dist_1


def Dist_tra_2_Rette(r, s):
    """
    Questa funzione calcola la distanza tra due rette parallele.
    Variabili:
        r       ==>     prima retta [ oggetto ]
        s       ==>     seconda reta [ oggetto ]
    Risultati:
        dist    ==>     distanza tra le due rette [ numero ]
    """
    p = Punto_by_2Retta(r, s)
    if p is None:
        P_1 = r.P0  # Generico punto sulla retta "r"
        dist = round(Dist_Punto_Retta(P_1, s), precisione)  # Distanza tra il punto generico su "r" e la retta "s"
    else:
        dist = 0.00
    return dist


def Dist_Retta_Piano(retta, piano):
    """
    Questa funzione calcola la distanza tra un piano e una retta ad esso parallela.
    Variabili:
        retta   ==>     Retta [ oggetto ]
        piano     ==>     Plane [ oggetto ]
    Risultati:
        Dist    ==>     Distanza [ numero ]
    """
    P_int = Punto_by_plan_e_retta(retta, piano)
    if type(P_int) == type(None):
        P_r = retta.Punto_by_Param()  # Punto su r
        s = Retta_by_piano_orto_e_punto(piano, P_r)  # Retta ortogonale a r.
        P_s = Punto_by_plan_e_retta(s, piano)  # Punto di intersezione tra s e il piano
        D = round(Dist(P_r, P_s), precisione)  # Distanza
    else:
        D = 0.0
    return D


def Dist_tra_2_piani(Pla_1, Pla_2):
    """
    Questa funzione calcola la distanza tra due piani paralleli.
    Variabili:
        Pla_1       ==>     Plane [ oggetto ]
        Pla_2       ==>     Plane [ oggetto ]
    Risultati:
        Dist        ==>     Distanza tra i piani [ numero ]
    """
    # Retta ortogonale al piano 1
    r = Retta_by_piano_orto_e_punto(Pla_1, Pla_1.P1)
    # punto di intersezione tra la retta "r" e il secondo piano.
    P_2 = Punto_by_plan_e_retta(r, Pla_2)
    Distanza = round(Dist(Pla_1.P1, P_2), precisione)  # Calcolo e arrotondamento valore distanza
    return Distanza


def r_ha_P(retta, punto):
    """
    Questa funzione ci dice se il punto appartiene alla retta indicata o meno
    Variabili:
        retta       ==>     Retta [ oggetto ]
        piano       ==>     Punto [ np.array() ]
    Risultati:
        [bolean value] True/False
    """
    P = None
    r = None
    # Ordinazione degli argomenti
    if type(punto) == type(np.zeros((3, 1))):
        if type(retta) == type(asse_x):
            r = retta
            P = punto
    elif type(punto) == type(asse_x):
        if type(retta) == type(np.zeros((3, 1))):
            r = punto
            P = retta
    else:
        r = retta
        P = punto
    #  Funzione
    D = round(Dist_Punto_Retta(P, r), precisione)  # distanza tra il punto e la retta
    if D <= 0.005:
        return True
    else:
        return False


def Retta_e_Retta(r, s, Stampa=True):
    """
    Questa funzione determina la relazione tra due rette
    Variabili:
        r       ==>     prima retta [ oggetto ]
        s       ==>     seconda reta [ oggetto ]
    Risultati:
        rette parallele:
            dist    ==>     distanza tra le due rette [ float ]
        rette incidenti:
            Punto   ==>     Punto in comune tra le rette [ np.array() ]
        rette sghembe:
            stampa le rette sono sghembe.
    """
    v_r = r.Vers  # Versore della retta r
    v_s = s.Vers  # Versore della retta s
    v_tra_rette = r.P1 - s.P1  # vettore tra le rette
    r_versori = np.linalg.matrix_rank([v_r, v_s])  # rango della matrice dei versori
    r_completa = np.linalg.matrix_rank([v_r, v_s, v_tra_rette])  # Rango della matrice completa
    if r_versori == 1:
        if r_completa == 1:
            if Stampa:
                print('Le rette sono coincidenti')
            return 0.00
        elif r_completa == 2:
            if Stampa:
                print('Le rette sono parallele')
            return Dist_tra_2_Rette(r, s)
        else:
            raise ValueError('I valori dei ranghi della matice completa sono impossibili')
    elif r_versori == 2:
        if r_completa == 2:
            if Stampa:
                print('Le rette sono incidenti')
            return Punto_by_2Retta(r, s)
        elif r_completa == 3:
            if Stampa:
                print('Le rette sono sghembe')
            return Dist_tra_2_Rette(r, s)
        else:
            raise ValueError('I valori dei ranghi della matice completa sono impossibili')


def Retta_e_Piano(retta, piano, stringa=True):
    """
    Questa funzione determina la relazione che esiste tra una retta e un piano
    nello spazio.
    Variabili:
        retta       ==>     Retta [ oggetto ]
        piano       ==>     Plane [ oggetto ]
        stringa ==>     Serve per decidere se scrivere le varie frasi che definiscono la condizione della retta. se viene chiamata all'interno di un'altra funzione e non serve la stampa di queste condizione impostare il valore su False altrimenti lasciarlo su True.
    Risultati:
        Relaz   ==>     se sono //  = distanza [ numero ]
                        se sono incidenti = Punto [np.array() ]
    """
    r = retta
    pla = piano
    D = Dist_Retta_Piano(r, pla)
    if D == 0.000:
        P = Punto_by_plan_e_retta(r, pla)
        d = Dist_Punto_piano(r.P0, pla) + Dist_Punto_piano(r.P1, pla) + Dist_Punto_piano(r.P2, pla)
        if d == 0.00:
            if stringa:
                print('La retta giace nel piano')
            return D
        else:
            if stringa:
                print("La retta e' incidente al piano")
            return P
    else:
        if stringa:
            print("La retta e' parallela al piano")
        return D


def Piano_e_Punto(punto, piano):
    """
    Questa funzione determina la relazione tra un punto e un piano.
    Se il punto appartiene al piano restituisce "True" altrimenti restituisce "False"
    Variabili:
        punto       ==>     Punto nello spazio [ oggetto ]
        piano      ==>     Plane nello spazio [ oggetto ]
    Risultati:
        True/False      ==> [bolean Value]
    """
    # Ordinazione degli argomenti
    pla = None
    P = None
    if type(punto) == type(np.zeros((3, 1))):
        if type(piano) == type(piano_xy):
            P = punto
            pla = piano
    elif type(punto) == type(piano_xy):
        if type(piano) == type(np.zeros((3, 1))):
            P = piano
            pla = punto
    else:
        raise TypeError()
    # Funzione vera e proprio
    dist = Dist_Punto_piano(P, pla)
    if dist < pow(10, -2):
        return True
    else:
        return False


def Piano_e_Piano(pla_1, pla_2, Stampa=True):
    """
    Questa funzione stabilisce i rapporti tra due piani nello spazio.
    Variabili:
        pla_1       ==>     Plane numero 1 [ oggetto ]
        pla_2       ==>     Plane numero 2 [ oggetto ]
    Risultati:
        Relaz   ==>     se sono //  = distanza [ numero ]
                        se sono incidenti = retta [ oggetto ]
    """
    Conf = Grandezza_vettore(pla_1.Vers_N - pla_2.Vers_N)
    if Conf == 0.00:
        if Stampa:
            print("I due piani sono paralleli")
        return Dist_tra_2_piani(pla_1, pla_2)
    else:
        if Stampa:
            print("I due piani sono incidenti")
        return Retta_by_2Planes(pla_1, pla_2)


def Proiezione_Punto_su_Retta(punto, retta):
    """
    Questa funzione calcola la proiezione di un punto su una generica retta.
    ARGOMENTI:
        punto   =>      Punto da proiettare [ np.array() ]
        retta   =>      Retta sulla quale proiettare il punto [ oggetto ]
    Risultato
        punto   =>      Punto derivante dalla proiezione del punto dato sulla retta.
                        [ np.array() ]
    """
    P = None
    r = None
    # Ordinazione degli argomenti
    if type(punto) == type(np.zeros((3, 1))):
        if type(retta) == type(asse_x):
            r = retta
            P = punto
    elif type(punto) == type(asse_x):
        if type(retta) == type(np.zeros((3, 1))):
            r = punto
            P = retta
    else:
        raise TypeError()
    # Funzione vera e propria
    Vettore_P0_P = P - r.P0  # Vettore da P0 della retta al punto P
    Vettore_P0_P1 = np.dot(Vettore_P0_P, r.Vers) * r.Vers  # Vettore da P0 alla proiezione di P sulla retta
    Punto_proiezione = r.P0 + Vettore_P0_P1  # Coordinate punto proiezione
    Punto_proiezione = Punto_proiezione.round(Cord_dec)
    return Punto_proiezione


def Proiezione_Punto_su_Piano(punto, piano):
    """
    questa funzione calcola la proiezione di una retta generica su un piano generico.
    ARGOMENTI:
        punto   =>      Retta da proiettare [ np.array() ]
        piano   =>      Plane sul quale eseguire la proiezione [ oggetto ]
    Risultato:
        punto   =>      Punto risultante della proiezione della punto dato sul piano.
    """
    Pia = None
    P = None
    # Ordinazione degli argomenti
    if type(punto) == type(np.zeros((3, 1))):
        if type(piano) == type(piano_xy):
            P = punto
            Pia = piano
    elif type(punto) == type(piano_xy):
        if type(piano) == type(np.zeros((3, 1))):
            P = piano
            Pia = punto
    else:
        raise TypeError()
    # Funzione vera e propria
    Retta_perp = Retta_by_piano_orto_e_punto(Pia, P)  # Creazione retta perpendicolare al piano passante per P
    Punto_proiezione = Punto_by_plan_e_retta(Retta_perp, Pia)  # intersezione tra la retta passante per P e il piano
    # La funzione restituisce il punto gia' arrotondato.
    return Punto_proiezione


def Proiezione_Retta_su_Piano(retta, piano):
    """
    questa funzione calcola la proiezione di una retta generica su un piano generico.
    ARGOMENTI:
        retta   =>      Retta da proiettare [ oggeto ]
        piano   =>      Plane sul quale eseguire la proiezione [ oggetto ]
    Risultato:
        retta   =>      Retta risultante della proiezione della retta data sul piano.
    """
    r = None
    Pla = None
    # Ordinazione degli argomenti
    if type(retta) == type(asse_x):
        if type(piano) == type(piano_xy):
            r = retta
            Pla = piano
    elif type(retta) == type(piano_xy):
        if type(piano) == type(asse_x):
            r = piano
            Pla = retta
    else:
        raise TypeError()
    # # Funzione
    Proie_P1 = Proiezione_Punto_su_Piano(r.P1, Pla)
    Proie_P2 = Proiezione_Punto_su_Piano(r.P2, Pla)
    Proiezione_retta = Retta(Proie_P1, Proie_P2)
    return Proiezione_retta


def Piano_medio(Lista_Piani):
    """
    Questa funzione calcola il piano medio tra i piani della lista.
    Il piano medio e' calcolato usando il versore della somma dei versori direzionali come versore direzionale; e come punto di partenza .
    ARGOMENTI:
        Lista_Piani     =>      Lista dei piani di cui calcolare la media.
    """
    Piani_medi = []  # Lista dove vengono messi i piani medi intermedi
    if len(Lista_Piani) != 1:
        Piani_medi = []  # Lista dove vengono messi i piani medi intermedi
        cord_x_y = np.array([0.000, 1.00])  # Coordinate da usare per calcolare i punti
        v0 = Lista_Piani[0].Vers_N  # Versore piano iniziale
        v1 = Lista_Piani[1].Vers_N  # Versore secondo piano
        if all([v0[n] == v1[n] for n in range(len(v0))]):  # Se i du piani sono paralleli
            p1 = 0.5 * (Lista_Piani[0].P1 + Lista_Piani[1].P1)
            p2 = 0.5 * (Lista_Piani[0].P2 + Lista_Piani[1].P2)
            p3 = 0.5 * (Lista_Piani[0].P3 + Lista_Piani[1].P3)
            Piani_medi.append(Piano(p1, p2, p3))  # Plane intermedio tra i due
        else:
            va = Versori(v0 + v1)  # somma dei due versori ( uno dei due piani intermedi)
            vb = Versori(v0 - v1)  # somma dei due versori ( uno dei due piani intermedi)
            r = Retta_by_2Planes(Lista_Piani[0], Lista_Piani[1])  # Retta intersezione tra i due piani
            punto = r.Punto_by_Param()  # Punto su tale retta
            piano_alfa = Piano_by_Vec_e_P(va, punto)  # uno dei due piani intermedi possibili
            piano_beta = Piano_by_Vec_e_P(vb, punto)  # uno dei due piani intermedi possibili
            z0 = Lista_Piani[0].Punto_by_Cord(cord_x_y)[2]  # Altezza punto sul primo piano
            z1 = Lista_Piani[1].Punto_by_Cord(cord_x_y)[2]  # Altezza punto sul secondo piano
            za = piano_alfa.Punto_by_Cord(cord_x_y)[2]  # Altezza punto sul piano alfa
            if min(z0, z1) <= za <= max(z0,
                                        z1):  # Valutazione tra il piano alfa e beta , quali dei due si trova tra i due piani iniziali
                Piani_medi.append(piano_alfa)
            else:
                Piani_medi.append(piano_beta)
        if len(Lista_Piani) > 2:
            for n in range(2, len(Lista_Piani)):
                v0 = Piani_medi[n - 2].Vers_N
                v1 = Lista_Piani[n].Vers_N
                if sum(v0 == v1) == 3:
                    p1 = 0.5 * (Lista_Piani[n].P1 + Piani_medi[n - 2].P1)
                    p2 = 0.5 * (Lista_Piani[n].P2 + Piani_medi[n - 2].P2)
                    p3 = 0.5 * (Lista_Piani[n].P3 + Piani_medi[n - 2].P3)
                    Piani_medi.append(Piano(p1, p2, p3))
                else:
                    va = Versori(v0 + v1)
                    vb = Versori(v0 - v1)
                    r = Retta_by_2Planes(Lista_Piani[0], Lista_Piani[1])
                    punto = r.Punto_by_Param()
                    piano_alfa = Piano_by_Vec_e_P(va, punto)
                    piano_beta = Piano_by_Vec_e_P(vb, punto)
                    z0 = Lista_Piani[0].Punto_by_Cord(cord_x_y)[2]  # Altezza punto sul primo piano
                    z1 = Lista_Piani[1].Punto_by_Cord(cord_x_y)[2]
                    za = piano_alfa.Punto_by_Cord(cord_x_y)[2]
                    if min(z0, z1) <= za <= max(z0, z1):
                        Piani_medi.append(piano_alfa)
                    else:
                        Piani_medi.append(piano_beta)
    else:
        Piani_medi.append(Lista_Piani[0])
    return Piani_medi[-1]


def Retta_media(Lista_rette):
    """
    Questa funzione di calcolare la retta che si trova a meta tra le rette presenti nella lista.

    """
    r_med = []  # Lista delle succcessive rette medie
    if len(Lista_rette) > 1:
        p1_ = Lista_rette[0].P1 + Lista_rette[1].P1  # primo punto medio
        p2_ = Lista_rette[0].P2 + Lista_rette[1].P2  # Secondo punto medio
        r_med.append(Retta(p1_, p2_))  # aggiunta della prima retta media alla lista delle rette medie
        if len(Lista_rette) > 2:
            for n in range(2, len(Lista_rette)):
                p1_ = Lista_rette[n].P1 + r_med[n - 2].P1  # primo punto medio
                p2_ = Lista_rette[n].P2 + r_med[n - 2].P2  # Secondo punto medio
                r_med.append(Retta(p1_, p2_))  # aggiunta della prima retta media alla lista delle rette medie
        return r_med[-1]
    else:
        return Lista_rette[-1]


def Elimina_punti_doppi(Lista):
    """
    Questa funzione elimina eventuali doppiuoni di punti presenti nella lista data come argomento.
    Restituisce la lista senza doppioni.
    """
    # Creazione lista dove raccogliere gli elemnti non doppi
    L_ = []
    # # Individuazione ed eliminazione doppioni
    for m in range(len(Lista)):
        n_copie = 0  # numero delle copie di quel punto
        for n in range(m + 1, len(Lista)):
            if all(Lista[m][i] == Lista[n][i] for i in range(len(Lista[m]))):
                n_copie = +1
        if n_copie == 0:
            L_.append(Lista[m])
    return L_


M = M0.Lista_punti(150, 100, 100, 100)
M = Ordina_punti(M)
pia_1 = Piano(M[0], M[5], M[2])
pia_2 = Piano(M[3], M[1], M[4])
pia_3 = Piano(M[7], M[9], M[6])
ret_1 = Retta(M[8], M[10])
ret_2 = Retta(M[10], M[132])
ret_3 = Retta(M[81], M[100])


def Crea_lista_doppioni(l, n, r=1):
    """
    Crea una lista con degli elemnti ripetuti piu' volte
    l = lunghezza finale della Lista di punti
    n = il numero di punti ripetuti
    r = il numero delle volte che i punti sono ripetuti.
    """
    L_ = []
    for doppi in range(n):
        P = Rnd.choice(M[l + 1:])
        for rip in range(r):
            L_.append(P)
    spazio_rimanente = l - len(L_)
    L_ = L_ + M[:spazio_rimanente]
    return L_


L = Crea_lista_doppioni(15, 3, 3)


def Crea_base(Lista, Tipo=None):
    """
    Questa funzione prmette di creare una base partendo da una lista di vettori.
    Argomenti:
        Lista         ==>     lista dei vettori con la quale si vuole creare la base.
        Tipo            ==>     Indica la tipologia di Base che si vule ottenere dai vettori.
                                Il parametro e' impostato su None e lascia i vettori forniti nella lista inalterati.
                                Se il parametro e' impostato su:
                                    Orto    ==> Crea na una base semplicemente ortogonale.
                                    Normale ==> Crea una una base i cui vettori componenti sono dei versori
                                    Ortonormale ==> Crea una base di vertori ortogonali.
    """
    if type(Tipo) == type(None):
        B = np.array(Lista)
        Base = B.transpose()
    elif type(Tipo) == type(''):
        if Tipo == 'Orto':
            v1 = Lista[0]
            v3 = np.cross(v1, Lista[1])
            v2 = np.cross(v1, v3)
            B = np.array([v1, v2, v3])
            Base = B.transpose()
        elif Tipo == 'Normale':
            for v in range(len(Lista)):
                vers = Versori(Lista[v])
                Lista[v] = vers
            B = np.array(Lista)
            Base = B.T
        elif Tipo == 'Ortonormale':
            for v in range(len(Lista)):
                vers = Versori(Lista[v])
                Lista[v] = vers
            v1 = Lista[0]
            v3 = np.cross(v1, Lista[1])
            v2 = np.cross(v1, v3)
            B = np.array([v1, v2, v3])
            Base = B.T
        else:
            raise Exception
    else:
        raise Exception
    return Base
