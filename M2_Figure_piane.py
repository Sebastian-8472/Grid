# -------------------------------------------------------------------------------
# Name:        modulo1
# Purpose:
#
# Author:      Sebastiano
#
# Created:     30/11/2015
# Copyright:   (c) Sebastiano 2015
# Licence:     <your licence>
# -------------------------------------------------------------------------------
from M1_Geometria_base import *  # in questo modo non devo riscrivere il nome del modulo in ogni funzione.
# Utile per mantenere i nomi dei moduli uguali e anche per copiare e incollare il codice da un modulo all'altro.


class Segmento(object):
    """
    Classe per creare un oggetto "segmento".
    Parte di una retta delimitato da due punti.
    """

    def __init__(self, P_1, P_2):
        """
        _L_'inizazione di questa classe parte dal equazione parametrica della retta su cui giace il segmento restituendo il punto usato come rigine e il versore direzionale.
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
        self.retta = Retta(P_1, P_2)
        self.Vet = self.P2 - self.P1
        self.P_medio = (self.P1 + self.P2) / 2.0
        self.Lunghezza = Dist(self.P1, self.P2)

    def __eq__(self, others):
        """
        uguaglianza ( == )
        """
        if type(self) == type(others):  # stesso Tipo
            if self.retta == others.retta:  # Stessa retta su cui giacciono i segmenti
                if self.Lunghezza == others.Lunghezza:  # Stessa lungheza
                    if sum(self.P_medio == others.P_medio) == 3.00:  # Stesso punto medio
                        return True
                    else:
                        return False  # Diverso punto medio
                else:
                    return False  # diversa lunghezza
            else:
                return False  # diversa retta
        else:
            return False  # diverso tipo

    def __ne__(self, others):
        """
        disuguaglianza ( != )
        """
        if self == others:
            return False
        else:
            return True

    def __gt__(self, others):
        """
        maggiore ( > )
        """
        if type(self) == type(others):  # stesso Tipo
            if self.retta.Vers == others.retta.Vers:  # Rette su cui giacciono i segmenti sono parallele
                if self.Lunghezza > others.Lunghezza:  # lunghezza
                    return True
                else:
                    return False  # diversa lunghezza
            else:
                return False  # diversa retta
        else:
            return False  # diverso tipo

    def __lt__(self, others):
        """
        minore ( < )
        """
        if type(self) == type(others):  # stesso Tipo
            if self.retta.Vers == others.retta.Vers:  # Rette su cui giacciono i segmenti sono parallele
                if self.Lunghezza < others.Lunghezza:  # lunghezza
                    return True
                else:
                    return False  # diversa lunghezza
            else:
                return False  # diversa retta
        else:
            return False  # diverso tipo

    def __ge__(self, others):
        """
        >=
        """
        if self.__gt__(others):  # i due segmenti sono maggiori
            return True
        else:
            if self.__eq__(others):  # i due segmenti sono identici
                return True
            else:
                return False

    def __le__(self, others):
        """
        <=
        """
        if self.__lt__(others):  # i due segmenti sono maggiori
            return True
        else:
            if self.__eq__(others):  # i due segmenti sono identici
                return True
            else:
                return False

    def __str__(self):
        """
        Restituisce la rappresentazione del segmento sotto forma di stringa che si ottiene usando il comando print.
        """
        stringa = "\nIl segmento e': \n\tLunghezza => %s \n\tPunto Medio => %s" % (self.Lunghezza, self.P_medio)
        return stringa

    def __repr__(self):
        """
        REstituisce la rappresentazione del segmento sotto forma di stringa.
        """
        Stringa = '\n\tSegmento\n\tVettore => %s \n\tLunghezza => %s' % (self.Vet, self.Lunghezza)
        return Stringa

    def Punto_by_Param(self, t=None):
        """
        Funzione per calcolare le coordinate di un punto generico appartenente alla retta del segmento;
        usando l'equazione parametrica della retta.
        VARIABILI.
        t       ==>     Indica un preciso valore del parametro per cui si vuole calcolare il punto. I valori del parametro variano tra 0 e 1
        """

        # Calcolo delle Coordinate
        if type(t) == type(None):
            t = Rnd.randrange(100) / 100.0
        P_ = self.P1 + (t * (self.P2 - self.P1))
        P_ = np.round(P_, Cord_dec)
        return P_

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
        if self.Vet[0] != 0.0:
            vers = self.Vet[0]
        elif self.Vet[1] != 0.0:
            vers = self.Vet[1]
        else:
            vers = self.Vet[2]
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

    def Proiezione(self, oggetto):
        """
        Questo metodo permette di ottenere la proiezione del segmento su piani o rette.
        La proiezione viene calcolata costruendo il segmento che congiunge le proiezione dei due vertici sull'ogggetto scelto.
        Il modulo restituisce un altro segmento, nel senso che crea un nuovo oggetto.
        ARGOMENTI:
            object      =>      Oggetto sulla quale si vuole proiettare il segmento.
                                Il modulo considera proiezioni solo su piani o rette.
        RISULTATO:
            Segmento    =>      un nuovo segmento, proiezione del segmento originale sull'oggetto scelto.
                                    [oggetto]
        """
        if type(oggetto) == type(piano_xy):
            punto_1 = Proiezione_Punto_su_Piano(self.P1, oggetto)  # Proiezione del primo punto sul piano scelto.
            punto_2 = Proiezione_Punto_su_Piano(self.P2, oggetto)  # Proiezione del secondo punto sul piano scelto.
            return Segmento(punto_1, punto_2)
        elif type(oggetto) == type(asse_x):
            punto_1 = Proiezione_Punto_su_Retta(self.P1, oggetto)  # Proiezione del primo punto sullla retta scelta.
            punto_2 = Proiezione_Punto_su_Retta(self.P2, oggetto)  # Proiezione del secondo punto sulla retta scelta.
            return Segmento(punto_1, punto_2)
        else:
            stringa = "non e' possibile eseguire la proiezione su " + str(type(oggetto))
            raise (TypeError(stringa))


seg_1 = Segmento(M[50], M[62])
seg_2 = Segmento(M[62], M[76])
seg_3 = Segmento(M[50], M[126])


def Segmento_da_Lista(Lista):
    """
    Crea una segmento partendo da una lista di due punti
    """
    s = Segmento(Lista[0], Lista[1])
    return s


def P_in_segm(Punto, segmento):
    """
    Questa funzione ci dice se un punto e' all'interno di un determinato segmento.
    ARGOMENTI:
        P           =>      Punto [np.array]
        segmento    =>      segmento [ogetto]
    Risultati:
        True        => se il punto e' all'interno del segmento
        False       => se il punto non e' all'interno del segmento
    """
    r = segmento.retta  # retta su cui giace il segmento
    if r_ha_P(r, Punto):  # il punto appartiene alla retta
        # Metodo 1
        t = round(r.Param_Punto(Punto), precisione)  # Parametro del punto scelto.
        t1 = round(r.Param_Punto(segmento.P1), precisione)  # Parametro del punto P1
        t2 = round(r.Param_Punto(segmento.P2), precisione)  # Parametro del punto P2
        if t >= min(t1, t2):  # Parametro del punto maggiore del minimo tra t1 e t2
            if t <= max(t1, t2):  # Parametro del punto minore del minimo tra t1 e t2
                metodo_1 = True
            else:
                metodo_1 = False
        else:
            metodo_1 = False
        ##                # Metodo 2
        ##        d_1 = Dist(P,segmento.P1)
        ##        d_2 = Dist(P,segmento.P2)
        ##        somma = d_1 + d_2
        ##        d = segmento.Lunghezza
        ##        if abs(somma - d) <= pow(10,-precisione):
        ##            metodo_2 = True
        ##        else:
        ##            metodo_2 = False
        return metodo_1  # ,metodo_2)
    else:
        return False  # ,False) # Il punto non e' sulla retta su cui giace il segmento


def Dist_Segme_Punto(punto, segmento):
    """
    Questa funzione permette di determinare la distanza tra un punto e un segmento.
    Tale distanza viene definita come il minimo tra la distanza tra il punto, i vertici e se esiste la proiezione del punto sul segmento.
    """
    d_1 = Dist(punto, segmento.P1)  # distanza dal primo vertice
    d_2 = Dist(punto, segmento.P2)  # Distanza dal secondo vertice
    Proiezione_punto = Proiezione_Punto_su_Retta(punto, segmento.retta)  # Punto proiezione sulla retta del segmento
    d_ = Dist(punto, Proiezione_punto)  # distanza punto - proiezione
    if P_in_segm(Proiezione_punto, segmento):
        return min(d_1, d_2, d_)
    else:
        return min(d_1, d_2)


def Segmento_e_retta(retta, segmento, stampa=True):
    """
    Questa funzione serve per capire in che rapporto stanno il segmento e la retta inseriti come argomenti.
    ARGOMENTI:
        retta       ==>     un oggetto retta.
        segmento    ==>     un oggetto segmento.
        stampa      ==>     Variabile per decidere se sampare i messaggi.
    RISULTATI:
        Dipende dal rapporto esistente tra la retta e il sgmento.
        Se il segmento giace sulla retta allora restituisce il segmento stesso.
        Se il segmento e la reta sono incindenti restituisce il punto di intersezione.
        Se il segmento e la retta non hanno niente in comune allora restituisce il valore None.
    In questo modo sara' possibile riconoscere i vari tipi di relazioni tra segmento e retta usando la funzione Type("").
    """
    # Ordinazione degli argomenti.
    if type(retta) == type(ret_1):
        r = retta
    else:
        r = segmento
    if type(segmento) == type(seg_1):
        lin = segmento
    else:
        lin = retta
    # Algoritmo vero e proprio
    r_lin = lin.retta  # linea su cui giace il segmento.
    r_e_r_lin = Retta_e_Retta(r, r_lin, stampa)  # Rapporto tra le due rette.
    if type(r_e_r_lin) == type(M[0]):  # le due rette sono incidenti.
        if P_in_segm(r_e_r_lin, lin):
            if stampa:
                print('La retta e il segmento sono incidenti')
            return r_e_r_lin
        else:
            if stampa:
                print('La rett e il segmento non sono incidenti')
            return None
    elif type(r_e_r_lin) == type(float()):  # Le due rette sono parallele o coincidenti.
        if r_e_r_lin <= pow(10, -(precisione + 1)):  # Rette coincidenti
            if stampa:
                print("Il segmento giacce sulla retta")
            return lin
        else:
            if stampa:
                print('la retta e il segmento sono paralleli')
            return None
    else:  # Le rette sono sghembe
        if stampa:
            print('la retta e il segmento sono sghembi')
        return None


def Segmento_e_Segmento(segmento_1, segmento_2, Stampa=True):
    """
    Questa funzione permette di stabilere in quale rapporto siano due segmenti.
    ARGOMENTI:
        segmento_1      ==>     Segmento numero 1
        segmento_2      ==>     Segmento numero 2
    RISULTATI:
        tuple == (intersezione trai segmenti, Tipo di rapporto tra i segmenti)
    Le stringhe che descrivono il raporto tra i segmenti sono le seguenti:
        Identici        ==>     Quando i segmenti sono uguali; sersituisce il segmento 1.
        Disgiunti       ==>     Quando i segmenti non hanno punti in comuni; resttuisce il valore None.
        Consecutivi     ==>     Quando i segmenti hanno un vertice in comune; restituisce il vertice in comune.
        Sovrapposti     ==>     Quando i segmenti hanno un tratto in comune; restituisce il segmento del tratto comune.
        Interno         ==>     Quando i segmenti sono uno dentro l'altro; restituisce il rapporto tra le lunghezze.
    In caso di semento interno l'algoritmo restituisce il rapporto tra la lunghezza del segmento 1 e la lunghezza del segmento 2.
    """
    l_vert = [segmento_1.P1, segmento_1.P2, segmento_2.P1, segmento_2.P2]  # Lista degli estremi  dei segmenti.
    if segmento_1 == segmento_2:  # Controllo se i segmenti sono uguali
        if Stampa:
            print('I segmenti sono uguali')
        return segmento_1, 'Identici'
    elif segmento_1.retta == segmento_2.retta:
        distanza_medi = Dist(segmento_1.P_medio, segmento_2.P_medio)
        if distanza_medi > 0.5 * (
                segmento_1.Lunghezza + segmento_2.Lunghezza):  # I segmenti non hanno alcuna sovraposizione.
            return None, 'Disgiunti'
        elif distanza_medi == 0.5 * (
                segmento_1.Lunghezza + segmento_2.Lunghezza):  # I segmenti hanno in comune un vertice.
            vert_seg_1 = [segmento_1.P1, segmento_1.P2]  # vertici del segmento 1.
            vert_seg_2 = [segmento_2.P1, segmento_2.P2]  # vertici del segmento 2.
            for v_1 in vert_seg_1:
                for v_2 in vert_seg_2:
                    if list(v_1) == list(v_2):
                        return v_1, 'Consecutivi'
        else:
            D_med_1 = Dist(segmento_2.P_Medio,
                           segmento_1.P1)  # Distanza tra il punto medio del segmento 2 e il vertice P1 del segmento 1.
            D_med_2 = Dist(segmento_2.P_Medio,
                           segmento_1.P2)  # Distanza tra il punto medio del segmento 2 e il vertice P2 del segmento 1.
            if D_med_1 >= (segmento_2.Lunghezza * 0.5):  # Il punto P1 del segmento 1 esterno al segmento 2
                if D_med_2 >= (segmento_2.Lunghezza * 0.5):  # Il punto P2 del segmento 1 esterno al segmento 2
                    fraz = segmento_1.Lunghezza / segmento_2.Lunghezza
                    return fraz, 'Interno'
                else:  # I due segmenti si sovrapongono parzialmente
                    if P_in_segm(segmento_2.P1, segmento_1):
                        return Segmento(segmento_2.P1, segmento_1.P1), 'Sovrapposti'
                    else:
                        return Segmento(segmento_2.P1, segmento_1.P2), 'Sovrapposti'
            else:  # Il punto P1 del segmento 1 interno al segmento 2
                if D_med_2 >= (segmento_2.Lunghezza * 0.5):  # Il punto P2 del segmento 1 esterno al segmento 2
                    if P_in_segm(segmento_2.P1, segmento_1):  # Il vertice P2 del segmento 2 si trova dentro il segmento_1
                        return Segmento(segmento_1.P1, segmento_2.P1), 'Sovrapposti'  # I segmenti si sovrappongono parzialmente.
                    else:
                        return Segmento(segmento_1.P1, segmento_2.P2)  # I segmenti si sovrappongono parzialmente.
                else:  # Il punto P2 del segmento 1 interno al segmento 2
                    fraz = segmento_1.Lunghezza / segmento_2.Lunghezza
                    return fraz, 'Interno'  # segmento 1 interno segmento 2
    else:
        Punto = Retta_e_Retta(segmento_1.retta, segmento_2.retta,
                              Stampa=False)  # relazione tra le rette su cui giaciono i segmenti
        if type(Punto) == type(M[0]):  # se le rette sono incidenti abbiamo un punto
            for vertice in l_vert:
                dist = round(Dist(vertice, Punto), precisione)
                if dist < pow(10, -precisione):
                    if Stampa:
                        print('i segmenti sono consecutivi')
                    return vertice, 'Consecutivi'
            if P_in_segm(Punto, segmento_1):  # controlla se il punto di intersezione cade sulle rette.
                if P_in_segm(Punto, segmento_2):  #
                    if Stampa:
                        print('I segmenti sono incidenti')
                    return Punto, 'Incidenti'
        else:
            if Stampa:
                print('I segmenti sono disgiunti')
            return None, 'Disgiunti'


def Segmento_e_Piano(segmento, piano, stampa=True):
    """
    Questa funzione serve per stabilire la relazione tra un segmento e un piano.
    ARGOMENTI:
        piano       ==>     un oggetto piano.
        segmento    ==>     un oggetto segmento.
        stampa      ==>     Variabile per decidere se stampare i messaggi.
    RISULTATI:
        Dipende dal rapporto esistente tra la retta e il segmento.
        Se il segmento giace nel piano allora restituisce il segmento stesso.
        Se il segmento e il piano sono incindenti restituisce il punto di intersezione.
        Se il segmento e il piano non hanno niente in comune allora restituisce il valore None se la retta ? parallela altrimenti restituisce la distanza tra la retta e il piano.
    In questo modo sara' possibile riconoscere i vari tipi di relazioni tra segmento e retta usando la funzione Type("").
    """
    # ORDINAZIONE ARGOMENTI.
    if type(segmento) == type(seg_1):
        lin = segmento
    else:
        lin = piano
    if type(piano) == type(pia_1):
        pla = piano
    else:
        pla = segmento
    # ALGORITMO VERO E PROPRIO
    R_pLA = Retta_e_Piano(lin.retta, pla, False)  # Rapporto tra il piano e la retta del segmento.
    if type(R_pLA) == type(M[0]):  # Il piano e la retta del segmento son incidenti
        if P_in_segm(R_pLA, lin):  # Il punto e' nel segmento.
            if stampa:
                print("Il piano e il segmento sono incidenti.")
            return R_pLA
        else:
            if stampa:
                print("Il segmento non intercetta il piano.")
            return None
    else:  # le rette sono esterne al piano o sono nel piano.
        if R_pLA <= pow(10, -4):  # La distanza tra piano e retta e' molto piccola.
            if stampa:
                print("Il segmento appartiene al piano.")
            return lin
        else:  # retta esterna al piano
            if stampa:
                print("Il segmento e' esterno al piano.")
            return R_pLA


class Poligono(object):
    """
    Crea un poligono usando i vertici del poligono.
    questa pu? essere considerata una classe generale dalle quale poi creremo
    """

    def __init__(self, Vertici):
        """
        Questa funzione inizializza la classe creando i dati fondamentali del poligono.
        ARGOMENTI:
            Vertici:    =>      Lista dei vertici del poligono. [ List of numpy arrays ]

        """
        self.Vertici = Vertici  # Vertici del poligono
        self.Num_lati = len(self.Vertici)  # Numero dei lati
        self.Centro = sum(self.Vertici) / 3.0  # Determina il centro tra vertici
        V = Vertici + Vertici  # Lista con i vertici ripetuti
        # Creazione dei Lati
        L_ = []  # Lista vuoto per raccogliere i lati
        for p in range(len(Vertici)):
            lato = Segmento(V[p], V[p + 1])  # lato
            L_.append(lato)
        self.Lati = L_
        # Creazione degli Angoli
        A = []  # Lista vuota per contenere gli angoli
        for p in range(len(Vertici)):
            ver_prece = Versori(V[p] - V[p - 1])  # Versore lato precedente il punto
            ver_succe = Versori(V[p] - V[p + 1])  # Versore lato successivo al punto
            cos_ang = round(np.dot(ver_prece, ver_succe), precisione)  # Coseno dell'angolo
            ang = np.arccos(cos_ang)  # Angolo interno
            if np.isnan(ang):
                ang = np.nan_to_num(ang)
            ang = round(ang, Angl_dec)  # Arrotondamento valore Angolo
            A.append(ang)
        self.Angoli = np.array(A)

    def Angoli_al_centro(self, Punto=None):
        """
        Questa funzione Calcola gli angoli al centro del poligono che sttendono i lati del poligono.
        Se viene in dicato un altro punto calcoleranno gli angoli intorno a quel punto che sottendono i lati del poligono
        """
        if type(Punto) == type(None):
            P = self.Centro
        else:
            P = Punto
        V = self.Vertici + self.Vertici
        A = []
        for n in range(self.Num_lati):
            v_1 = P - V[n]
            v_2 = P - V[n + 1]
            ang = Angoli_tra_vettori(v_1, v_2)
            A.append(ang)
        A = np.array(A)
        return (A)

    def Piano(self):
        """
        Calcola il piano su cui giace il poligono.
        """
        P_plan = Piano(self.Vertici[0], self.Centro, self.Vertici[1])
        I = []
        for n in range(2, self.Num_lati):
            I.append(Piano_e_Punto(self.Vertici[n], P_plan))
        if len(I) == sum(I):
            return (P_plan)
        else:
            print('Il piano non e\' unico.')
            return (None)

    def Perimetro(self):
        perimetro = 0.00
        for n in range(self.Num_lati):
            perimetro += self.Lati[n].Lunghezza
        perimetro = round(perimetro, precisione)
        return (perimetro)

    def Area(self):
        """
        Calcola l'area del poligono.
        """
        Vert = self.Vertici + self.Vertici  # Creazione lista ripetuta
        A = 0
        for punto in range(len(self.Vertici)):
            v_1 = self.Centro - Vert[punto]
            v_2 = self.Centro - Vert[punto + 1]
            area_parz = 0.5 * Grandezza_vettore(np.cross(v_1, v_2))  # area Parziale
            A += area_parz  # Somma delle aree parziali
        A = round(A, precisione)
        return (A)

    def Raggi(self):
        """
        Calcola i raggi del poligono ( dal centro al vertici).
        """
        R = []
        for punto in self.Vertici:
            r = Segmento(self.Centro, punto)
            R.append(r)
        return (R)

    def Diagonali(self):
        """
        Calcola tutte le diagonali del poligono. Restituisce le diagonali come una tuuple in cui sono presenti la retta e la lunghezza della diagonale.
        """
        V = self.Vertici + self.Vertici  # Accodare la lista per poter avere un circolo chiuso dei vertici
        Diag = {}  # Database vuoto dove raccogliere le diagonali
        n_d = (self.Num_lati - 1) - 2  # Numero diagonali per punto
        for e in range(self.Num_lati):
            for l in range(e + 2, e + 2 + n_d):
                d = Segmento(V[e], V[l])  # Segmento della diagonale
                # Nome della diagonale da usare come chiave del database
                if l >= self.Num_lati:
                    a = l - self.Num_lati
                else:
                    a = l
                name = str(e) + str(a)
                Diag[name] = (d)
        return (Diag)

    def Sposta(self, Spostamento, Riferimento=None):
        """
        Questa funzione permette di Spostare il poligono.
        Il punto di riferimento ? impostato sul centro del Poligono ma se si vuole si puo' indicare u altro punto generico.
        ARGOMENTI:
            Spostamento =>  Vettore di spostamento.
                            Se viene indicato un solo valore si avvra uno spostamento pari a quel valore su tutti gli assi.
                            Puo' essere  un array o un valore singolo.
            Riferimento =>  Punto di riferimento in base alla quale avviene lo spostamento. Se non viene indicato niente si usa il centro-
        """
        ### Punto di riferimento
        if type(Riferimento) == type(None):
            P_rif = self.Centro
        else:
            P_rif = Riferimento
        ### Vettore spostamento
        if type(Spostamento) != type(np.zeros((0, 3))):
            if type(Spostamento) == type(3) or type(Spostamento) == type(3.0):
                V_Spo = np.array([Spostamento, Spostamento, Spostamento])
            else:
                V_Spo = np.array(Spostamento)
        else:
            V_Spo = Spostamento
        # Calcolo dello Spostamento
        New_rif = P_rif + V_Spo  # Coordinate del punto usato come riferimento in seguito allo spostamento.
        ### Calcolo delle coordinate dei nuovi Vertici
        Nuovi_Vertici = []
        for v in range(self.Num_lati):
            pos = self.Vertici[v] - P_rif  # posizione del punto rispetto al veccio punto di riferimento
            Vert_nuovo = New_rif + pos  # Nuove coordinate del vertice
            Vert_nuovo = np.round(Vert_nuovo, Cord_dec)
            Nuovi_Vertici.append(Vert_nuovo)

    def Ruota(self, Angolo, Centro_rot=None, asse_rot=None):
        """
        Questa funzione permette di ruotare il poligono.
        _L_'angolo deve essere espresso in radici.
        """
        if type(None) == type(Centro_rot):
            Centro_rot = self.Centro
        if type(asse_rot) == type(None):
            asse_rot = self.Piano().Vers_N
        new_v = []  # Lista dei nuovi vertici
        for vertice in self.Vertici:
            nw_v = Rotazione(vertice, Angolo, Centro_rot, asse_rot)
            n_v = np.round(nw_v, Cord_dec)
            new_v.append(n_v)
        return (self.__init__(new_v))

    def Punto_Inside(self, Punto, Stampa=True):
        """
        Questa funzione stabilisce se un punto e' all'interno del Poligono
        """
        if Piano_e_Punto(Punto, self.Piano()):  # Verifica che il punto apppartenga al piano del poligono
            ## Calcolo dei raggi minimi e massimi
            # Raggio minimo
            R_min = []  # Lista dove metto i raggi relativi ai punti medi dei lati
            for lato in self.Lati:
                r_min = Dist(lato.P_medio, self.Centro)  # raggoi minimo del lato Cento--Punto medio.
                R_min.append(r_min)
            raggio_minimo = min(R_min)
            R_max = []  # Lista dove metto i raggi relativi ai vertici del poligono.
            # Raggio massimo
            for vertice in self.Vertici:
                r_max = Dist(vertice, self.Centro)  # Raggio dei vertici (tra questi si trova il massimo raggio.
                R_max.append(r_max)
            raggio_massimo = min(R_max)
            raggio_punto = Dist(self.Centro, Punto)
            if raggio_punto <= raggio_minimo:  # Il punto si trova all'interno
                return (True)
            elif raggio_punto > raggio_massimo:  # Il punto si trova di sicuro all'esterno.
                return (False)
            else:  # R_mi < R_punto < R_max
                in_vertice = False
                in_lat = False
                linee = []
                Relazioni = []
                # Verifica se il punto e' un vertice.
                for vertice in self.Vertici:
                    linee.append(Segmento(vertice, Punto))  # Segmento che unisce i? punto a un vertice.
                    if Dist(vertice, Punto) <= pow(10, -precisione):  # Il punto e' un vertice
                        in_vertice = True
                for lato in self.Lati:
                    for l in linee:
                        Relazioni.append(
                            Segmento_e_Segmento(lato, l, False))  # relazione tra i l segmento e il lato del poligono.
                    if P_in_segm(Punto, lato) == True:  # Il punto si trova nel segmento.
                        in_lat = True
                if any(vertice) == True:
                    return (True)
                elif any(in_lat) == True:
                    return (True)
                elif any(Relazioni) == 'incidenti':
                    return (True)
                else:
                    return (False)
        else:
            if Stampa == True:
                print('\nIl punto non appartiene al piano del poligono.\n')
            return (False)

    def Individua_vertice(self, punto, Stampa=True):
        """
        Questa funzione permette di individuare il numero di indice al quale si trova il vertice dalle cordinate uguale al punto fornito come argomento.
        ARGOMENTO:
            punto   =>      Coordinate del probabile vertice [ np.array() ]
        RISULTATI:
            N?      =>      Il numero d'indice del corrispondente vertice [ numero intero ]
                            Se il punto non e' un vertice del poligono il valore viene impostato sulla varibile booleana None
        """
        n = None  # Numero dell'indice
        for v in range(len(self.Vertici)):
            if Grandezza_vettore(self.Vertici[v] - punto) < pow(10, -Cord_dec):
                n = v
        if n == None:
            if Stampa == True:
                print('Il punto non e\' un vertice del Poligono')
        return (n)

    def Individua_Lato(self, segmento, Stampa=True):
        """
        Questa funzione permette di individuare il lato uguale al segmento dato, restituendo il numero d'indice alla quale corrisponde il lato nella lista dei lati del poligono.
        ARGOOMENTI:
            segmento        =>      [oggetto]
        """
        n = None  # Numero dell'indice
        for v in range(len(self.Lati)):
            if segmento == self.Lati[v]:
                n = v
        if n == None:
            if Stampa == True:
                print('Il segmento non e\' un lato del Poligono')
        return (n)

    def Individua_Raggi(self, segmento, Stampa=True):
        """
        Questa funzione permette di individuare il raggio uguale al segmento dato, restituendo il numero d'indice alla quale corrisponde il raggio nella lista dei raggi del poligono.
        ARGOOMENTI:
            segmento        =>      [oggetto]
        """
        n = None  # Numero dell'indice
        for v in range(len(self.Raggi())):
            if segmento == self.Raggi()[v]:
                n = v
        if n == None:
            if Stampa == True:
                print('Il segmento non e\' un raggio del Poligono')
        return (n)

    def Individua_Diagonale(self, segmento, Stampa=True):
        """
        Questa funzione permette di individuare la diagonale uguale al segmento dato, restituendo la chiave nel Database in cui sono contenute tutte le Diagonali del poligono.
        ARGOOMENTI:
            segmento        =>      [oggetto]
        RISULTATO:
            Chiave          =>      [stringa]
        """
        chiave = ''  # Valore iniziale della chiave
        l = [segmento.P1, segmento.P2]  # Lista dei verici dei segmenti
        if self.Num_lati > 3:  # Il Poligono ha piu' di tre lati
            if type(self.Individua_Lato(segmento, False)) == type(None):
                if type(self.Individua_Raggi(segmento, False)) == type(None):
                    n = 0  # conta il numero dei vertici
                    for punto in l:
                        indice = self.Individua_vertice(punto, False)
                        if type(indice) != type(None):
                            n += 1
                            chiave += str(indice)
            if n == 2:  # Ci sono almeno due punti
                return (chiave)
            else:
                chiave = ''
                return (chiave)
        else:  # Il poligono e' un triangolo quindi non ha diagonali
            if Stampa == True:
                print('Il segmento non e\' una diagonale del Poligono')
            return (chiave)

    def Trova_Lato(self, punto=None, vert=None, retta=None):
        """
        Questa funzione restituisce un lato del poligono che si trova vicino al punto scelto, oppure compreso tra i vertici scelti o che sia parallelo alla retta scelta. La funzione permette di assegnare anche tutti e tre i parametri di diverso tipo o di assegnare liste ai parametri.
        ARGOMENTI:
            punto       ==>     cordinate di un punto o lista di punti
            vert        ==>     vertici o lista di vertici
            retta       ==>     Retta o lista di rette
        RISULTATI:
            una segmento o una lista di segmenti.
        Quando viene assegnato un solo vertice l'algoritmo restituisce il lato precedente e quello successivo al vertice specificato, se invece viene passata una lista con piu' di due vertici l'algoritmo restituisce il lato compreso tra due vertici creando coppie di due vertici.
        Se si assegna un oggetto retta al parametro retta l'algoritmo restituisce il lato che giace su quella retta, se si inserisce una lista in cui sono contenute le rette l'algoritmo restituisce i lati paralleli a quella retta.
        """
        Risultati = []  # Lista nella quale verrano messi i lati trovati permettendo di usare le varie opzioni conteporaneamente.
        ### Punto
        if type(punto) != type(None):  # Viene inserito un punto o una lista di punti
            if type(punto) == type(M[0]):  # se il paramtro e' un solo punto
                Lista_distanze = []
                # Calcolo distanze dai lati
                for lato in self.Lati:
                    d = Dist_Segme_Punto(punto, lato)
                    Lista_distanze.append(d)
                indice = Lista_distanze.index(min(Lista_distanze))  # indice lato con minor distanza.
                Risultati.append(
                    self.Lati[indice])  # Aggiunge alla lista dei risultati il lato piu' vicino al punto indicato.
            elif type(punto) == type([]):  # Viene inserita una lista di punti.
                for point in punto:  # valuta ogni punto nella Lista dei punti.
                    Lista_distanze = []
                    # Calcolo distanze dai lati
                    for lato in self.Lati:
                        d = Dist_Segme_Punto(punto, lato)
                        Lista_distanze.append(d)
                    indice = Lista_distanze.index(min(Lista_distanze))  # indice lato con minor distanza.
                    Risultati.append(
                        self.Lati[indice])  # Aggiunge alla lista dei risultati il lato piu' vicino al punto indicato.
            else:
                raise TypeError('Assegnare un punto in formatto numpy array o una lista di punti ')
        ### Vertice
        if type(vert) != type(None):  # Vengono indicati dei vertici.
            if type(vert) == type(M[0]):  # Viene assegnat un solo vertice.
                Lista_distanze = []
                # Calcolo distanze dai lati
                for lato in self.Lati:
                    d = Dist_Segme_Punto(vert, lato)
                    Lista_distanze.append(d)
                indice = Lista_distanze.index(min(Lista_distanze))  # indice lato con minor distanza.
                Risultati.append(self.Lati[indice])  # Aggiunge alla lista dei risultati il lato precedente al vertice.
                Risultati.append(
                    self.Lati[indice + 1])  # Aggiunge alla lista dei risultati il lato sucessivo al vertice.
            elif type(vert) == type([]):  # Vengono forniti una serie di vertici
                if len(vert) == 1:  # Lista di un solo vertice ==> Errore
                    raise TypeError('La lista di vertici deve contenere almeno due vertici')
                elif len(vert) == 2:  # Lista di soli due vertici? ha senso tenerla?
                    for lato in self.Lati:  # considera i lati nella lista dei lati del poligono.
                        essere_vertice = []  # Lista per comparazione vertici del lato
                        for p in vert:  # passa attraverso tutti i punti della lista vertice
                            if Dist(p, lato.P1) <= pow(10, -4):  # Distanza dekl punto dal primo vertice del lato
                                essere_vertice.append(True)
                            elif Dist(p, lato.P2) <= pow(10, -4):  # Distanza dekl punto dal secondo vertice del lato
                                essere_vertice.append(True)
                            else:
                                essere_vertice.append(False)
                        if all(
                                essere_vertice) == True:  # Se tutti e due i vaolri sono veri il lato sitrova tra i due punti
                            Risultati.append(lato)
                else:  # Quando la lista contiene piu' di due vertici. ( si accopiano a due a due)
                    for n in range(1, len(vert)):
                        soto_vert = vert[n - 1:n + 1]  # Creazione sottolista della copia di vertici da considerare
                        for lato in self.Lati:  # considera i lati nella lista dei lati del poligono.
                            essere_vertice = []  # Lista per comparazione vertici del lato
                            for p in soto_vert:  # passa attraverso tutti i punti della lista vertice
                                if Dist(p, lato.P1) <= pow(10, -4):  # Distanza dekl punto dal primo vertice del lato
                                    essere_vertice.append(True)
                                elif Dist(p, lato.P2) <= pow(10,
                                                             -4):  # Distanza dekl punto dal secondo vertice del lato
                                    essere_vertice.append(True)
                                else:
                                    essere_vertice.append(False)
                            if all(
                                    essere_vertice) == True:  # Se tutti e due i vaolri sono veri il lato sitrova tra i due punti
                                Risultati.append(lato)
            else:
                raise TypeError('Assegnare un vertice in formatto numpy array o una lista di vertici ')
        ### Retta
        if type(retta) != type(None):  # Vengono indicate delle rette
            if type(retta) == type(asse_x):  # viene indicata una sola retta.
                for lato in self.Lati:  # considero ogni lato del poligono uno alla volta.
                    if retta == lato.retta:  # la retta e' uguale alla retta del segmento.
                        Risultati.append(lato)
            elif type(retta) == type([]):  # viene assegnata una lista di rette
                for r in retta:  # passa attraverso tutte le rette della lista
                    for lato in self.Lati:  # considero ogni lato del poligono uno alla volta.
                        v_retta = retta.Vers  # versore della retta
                        v_seg = Versori(lato.Vet)  # Versore segmento
                        if Dist(v_retta, v_seg) <= pow(10,
                                                       -4):  # La distanza tra i due versori e' inferiore a una certa precisione
                            Risultati.append(lato)
            else:
                raise TypeError('Assegnare una retta o una lista di rette')
        return (Risultati)

    def Trova_Diagonale(self, punto=None, vert=None, retta=None):
        """
        Questa funzione restituisce una diagonale del poligono che si trova vicino al punto scelto, oppure compresa tra i vertici scelti o che sia parallela alla retta scelta.
        La funzione permette di assegnare anche tutti e tre i parametri di diverso tipo o di assegnare liste ai parametri.
        ARGOMENTI:
            punto       ==>     cordinate di un punto o lista di punti
            vert        ==>     vertici o lista di vertici
            retta       ==>     Retta o lista di rette
        RISULTATI:
            una segmento o una lista di segmenti.
        Quando viene assegnato un solo vertice l'algoritmo restituisce le diagonali generate dal vertice specificato, se invece viene passata una lista con piu' di due vertici l'algoritmo restituisce la diagonale compresa tra due vertici creando coppie di due vertici.
        Se si assegna un oggetto retta al parametro retta l'algoritmo restituisce la digonale che giace su quella retta, se si inserisce una lista in cui sono contenute le rette l'algoritmo restituisce le digonali parallele a quella retta.
        """
        Risultati = []  # Lista nella quale verrano messi i diagonali trovati permettendo di usare le varie opzioni conteporaneamente.
        ### Punto
        if type(punto) != type(None):  # Viene inserito un punto o una lista di punti
            if type(punto) == type(M[0]):  # se il paramtro e' un solo punto
                Lista_distanze = []
                # Calcolo distanze dai lati
                for diag in self.Diagonali():
                    d = Dist_Segme_Punto(punto, self.Diagonali()[diag])
                    Lista_distanze.append(d, diag)
                Lista_distanze.sort()
                chiave_diagonale = Lista_distanze[0][1]
                Risultati.append(self.Diagonali()[
                                     chiave_diagonale])  # Aggiunge alla lista dei risultati la diagonale piu' vicina al punto indicato.
            elif type(punto) == type([]):  # Viene inserita una lista di punti.
                for point in punto:  # valuta ogni punto nella Lista dei punti.
                    Lista_distanze = []
                    # Calcolo distanze dai lati
                    for diag in self.Diagonali():
                        d = Dist_Segme_Punto(punto, self.Diagonali()[diag])
                        Lista_distanze.append(d, diag)
                    Lista_distanze.sort()
                    chiave_diagonale = Lista_distanze[0][1]
                    Risultati.append(self.Diagonali()[
                                         chiave_diagonale])  # Aggiunge alla lista dei risultati la diagonale piu' vicina al punto indicato.
            else:
                raise TypeError('Assegnare un punto in formatto numpy array o una lista di punti ')
        ### Vertice
        if type(vert) != type(None):  # Vengono indicati dei vertici.
            if type(vert) == type(M[0]):  # Viene assegnat un solo vertice.
                indice_vertice = None  # Variabile dove salvare l'indice a cui corrisponde il vertice indicato
                for v in range(
                        len(self.Vertici)):  # Si cerca il valore vero dell'indice del vertice indicato come argomento
                    if Dist(self.Vertici[v], vert) <= pow(10, -4):
                        indice_vertice = str(v)  # Indice del vertice pronto ad essere usato come chiave parziale
                chiavi = self.Diagonali().keys()  # lista di tutte le chiavi del database delle diagonali.
                for k in chiavi:  # considera ciascuna chiave.
                    if k.count(indice_vertice) > 0:
                        if k[0] == indice_vertice:
                            Risultati.append(self.Diagonali()[k])  # Aggiunge la diagonale alla lista dei risultati
            elif type(vert) == type([]):  # Vengono forniti una serie di vertici
                if len(vert) == 1:  # Lista di un solo vertice ==> Errore
                    raise TypeError('La lista di vertici deve contenere almeno due vertici')
                elif len(vert) == 2:  # Lista di soli due vertici? ha senso tenerla?
                    Seg = Segmento(vert[0], vert[1])  # Segmento congiungente i due vertici.
                    is_lato = False  # parametro per la verifica che il segmento non sia un lato.
                    for lato in self.Lati:
                        if Seg == lato:
                            is_lato = True
                    if is_lato == False:
                        Risultati.append(Seg)
                else:  # Quando la lista contiene piu' di due vertici. ( si accopiano a due a due)
                    for ind in range(len(vert) - 1):
                        Seg = Segmento(vert[ind], vert[ind + 1])  # Segmento congiungente i due vertici.
                        is_lato_1 = False  # parametro per la verifica che il segmento non sia un lato.
                        for lato in self.Lati:
                            if Seg == lato:
                                is_lato_1 = True
                        if is_lato_1 == False:
                            Risultati.append(Seg)
            else:
                raise TypeError('Assegnare un vertice in formatto numpy array o una lista di vertici ')
        ### Retta
        if type(retta) != type(None):  # Vengono indicate delle rette
            # cREAZIONE LISTA DEI SEGMENTI DELLE DIAGONALI.
            Lista_diagonali = []
            for item in self.Diagonali().items():
                Lista_diagonali.append(item[1])
            if type(retta) == type(asse_x):  # viene indicata una sola retta.
                for diag in Lista_diagonali:  # considero ogni lato del poligono uno alla volta.
                    if retta == diag.retta:  # la retta e' uguale alla retta del segmento.
                        Risultati.append(diag)
            elif type(retta) == type([]):  # viene assegnata una lista di rette
                for r in retta:  # passa attraverso tutte le rette della lista
                    for diag in Lista_diagonali:  # considero ogni lato del poligono uno alla volta.
                        v_retta = retta.Vers  # versore della retta
                        v_seg = Versori(diag.Vet)  # Versore segmento
                        if Dist(v_retta, v_seg) <= pow(10,
                                                       -4):  # La distanza tra i due versori e' inferiore a una certa precisione
                            Risultati.append(diag)
            else:
                raise TypeError('Assegnare una retta o una lista di rette')
        return (Risultati)

    def Trova_Raggio(self, punto=None, vert=None, retta=None):
        """
        Questa funzione restituisce un raggio del poligono che si trova vicino al punto scelto, oppure compreso tra i vertici scelti o che sia parallelo alla retta scelta. La funzione permette di assegnare anche tutti e tre i parametri di diverso tipo o di assegnare liste ai parametri.
        ARGOMENTI:
            punto       ==>     cordinate di un punto o lista di punti
            vert        ==>     vertici o lista di vertici
            retta       ==>     Retta o lista di rette
        RISULTATI:
            una segmento o una lista di segmenti.
        Se si assegna un oggetto retta al parametro retta l'algoritmo restituisce il raggio che giace su quella retta, se si inserisce una lista in cui sono contenute le rette l'algoritmo restituisce i raggi paralleli a quella retta.
        """
        Risultati = []  # Lista nella quale verrano messi i raggi trovati permettendo di usare le varie opzioni conteporaneamente.
        ### Punto
        if type(punto) != type(None):  # Viene inserito un punto o una lista di punti
            if type(punto) == type(M[0]):  # se il paramtro e' un solo punto
                Lista_distanze = []
                # Calcolo distanze dai raggi
                for rad in self.Raggi():
                    d = Dist_Segme_Punto(punto, rad)
                    Lista_distanze.append(d)
                indice = Lista_distanze.index(min(Lista_distanze))  # indice raggio con minor distanza.
                Risultati.append(
                    self.Raggi()[indice])  # Aggiunge alla lista dei risultati il raggio piu' vicino al punto indicato.
            elif type(punto) == type([]):  # Viene inserita una lista di punti.
                for point in punto:  # valuta ogni punto nella Lista dei punti.
                    Lista_distanze = []
                    # Calcolo distanze dai raggi
                    for rad in self.Raggi():
                        d = Dist_Segme_Punto(punto, rad)
                        Lista_distanze.append(d)
                    indice = Lista_distanze.index(min(Lista_distanze))  # indice raggio con minor distanza.
                    Risultati.append(
                        self.Lati[indice])  # Aggiunge alla lista dei risultati il raggio piu' vicino al punto indicato.
            else:
                raise TypeError('Assegnare un punto in formatto numpy array o una lista di punti ')
        ### Vertice
        if type(vert) != type(None):  # Vengono indicati dei vertici.
            if type(vert) == type(M[0]):  # Viene assegnat un solo vertice.
                rad = Segmento(vert, self.Centro)
                Risultati.append(rad)
            elif type(vert) == type([]):  # Vengono forniti una serie di vertici
                for v in vert:  # passa attraverso gli elementi della lista dei vertici
                    rad = Segmento(self.Centro, v)  # Segmento tra il centro del poligono e il vertice
                    Risultati.append(rad)
            else:
                raise TypeError('Assegnare un vertice in formatto numpy array o una lista di vertici ')
        ### Retta
        if type(retta) != type(None):  # Vengono indicate delle rette
            if type(retta) == type(asse_x):  # viene indicata una sola retta.
                for rad in self.Raggi():  # considero ogni raggio del poligono uno alla volta.
                    if retta == rad.retta:  # la retta e' uguale alla retta del segmento.
                        Risultati.append(rad)
            elif type(retta) == type([]):  # viene assegnata una lista di rette
                for r in retta:  # passa attraverso tutte le rette della lista
                    for rad in self.Raggi():  # considero ogni raggio del poligono uno alla volta.
                        v_retta = retta.Vers  # versore della retta
                        v_seg = Versori(rad.Vet)  # Versore segmento
                        if Dist(v_retta, v_seg) <= pow(10,
                                                       -4):  # La distanza tra i due versori e' inferiore a una certa precisione
                            Risultati.append(rad)
            else:
                raise TypeError('Assegnare una retta o una lista di rette')
        return (Risultati)

    def Proiezione(self, object):
        """
            Questo metodo permette di ottenere la proiezione del poligono su piani o rette.
            La proiezione viene calcolata costruendo il poligono che ha come vertici le proiezione dei vertici sull'ogggetto scelto.
            Il modulo restituice un altro poligono, nel senso che crea un nuovo oggetto.
            ARGOMENTI:
                object      =>      Oggetto sulla quale si vuole proiettare il segmento.
                                    Il modulo considera proiezioni solo su piani o rette.
            RISULTATO:
                Segmento    =>      un nuovo poligono, proiezione del segmento originale sull'oggetto scelto.
                                        [oggetto]
        """
        if type(object) == type(piano_xy):
            Proiezione_Vertici = []  # Lista nel quale raccogliere le proiezioni dei vertici
            for v in self.Vertici:
                punto = Proiezione_Punto_su_Piano(v, object)
                Proiezione_Vertici.append(punto)
            return (Poligono(Proiezione_Vertici))
        else:
            stringa = "non e' possibile eseguire la proiezione su " + str(type(object))
            raise (TypeError(stringa))

    def __repr__(self):
        """
        REstituisce una serie di stringhe che rappresentano il poligono.
        """
        Stringa = 'Numero lati = %s \n' % (self.Num_lati)
        for n in range(self.Num_lati):
            Stringa += 'Lato  %s: %s \n' % (n, self.Lati[n].Lunghezza)
        Stringa += 'Area: %s' % (self.Area())
        return (Stringa)

    def __str__(self):
        """
        REstituisce una serie di stringhe che rappresentano il poligono.
        """
        Stringa = 'Numero lati = %s \n' % self.Num_lati
        Stringa += 'Area: %s' % (self.Area())
        return (Stringa)

    def __eq__(self, others):
        """
        == ( I due poligono sono identici, per l'equivalenza confrontare l'area.
        """
        if type(self) == type(others):  # Stesso Tipo
            if self.Num_lati == others.Num_lati:  # Stesso numero di Lati
                stessi_punti = 0.0  # Variabile per controllare se ci sono punti diversi
                for vertice in self.Vertici:  # Controlla l'uguaglianza dei singoli vertici
                    if sum(vertice == others.Vertici) == 3.0:
                        stessi_punti += 1
                if stessi_punti == len(self.Vertici):  # Stessi vertici
                    return (True)
                else:
                    return (False)  # Diversi Vertici
            else:  # Diverso numero di lati
                return (False)
        else:  # Diverso Tipo
            return (False)

    def __ne__(self, others):
        """
        != ( I due poligono sono diversi.
        """
        if type(self) == type(others):  # Stesso Tipo
            if self.Num_lati == others.Num_lati:  # Stesso numero di Lati
                stessi_punti = 0.0  # Variabile per controllare se ci sono punti diversi
                for vertice in self.Vertici:  # Controlla l'uguaglianza dei singoli vertici
                    if sum(vertice == others.Vertici) == 3.0:
                        stessi_punti += 1
                if stessi_punti == len(self.Vertici):  # Stessi vertici
                    return (False)
                else:
                    return (True)  # Diversi Vertici
            else:  # Diverso numero di lati
                return (True)
        else:  # Diverso Tipo
            return (True)

    def __gt__(self, others):
        """
        > maggiore
        """
        if type(self) == type(others):  # Stesso Tipo
            if self.Num_lati == others.Num_lati:  # Stesso numero di Lati
                if self.Area() > others.Area():  # area maggiore
                    return (True)
                else:
                    return (False)  # area minore
            else:  # Diverso numero di lati
                return (True)
        else:  # Diverso Tipo
            return (True)

    def __lt__(self, others):
        """
        < maggiore
        """
        if type(self) == type(others):  # Stesso Tipo
            if self.Num_lati == others.Num_lati:  # Stesso numero di Lati
                if self.Area() < others.Area():  # area minore
                    return (True)
                else:
                    return (False)  # area maggiore
            else:  # Diverso numero di lati
                return (True)
        else:  # Diverso Tipo
            return (True)

    def __ge__(self, others):
        """
        >=  maggiore o uguale
        """
        if self.__eq__(others):  # Poligono uguali
            if self.__gt__(others):  # Poligoni maggiori
                return (True)
            else:
                return (False)
        else:
            return (False)

    def __le__(self, others):
        """
        < maggiore
        """
        if self.__eq__(others):  # Poligono uguali
            if self.__lt__(others):  # Poligoni minori
                return (True)
            else:
                return (False)
        else:
            return (False)


L = M[0:3]
P = Poligono(L)


class Poligono_regolare(Poligono):
    """
    Creazione Poligono regolare come classe.
    """

    def __init__(self, R=None, Lato=None, Num_lati=3, Circ='est', Centro_poligono=np.array([0.0, 0.0, 0.0])):
        self.Num_lati = Num_lati  # Numero dei lati
        self.Centro = Centro_poligono  # Determina il centro tra vertici
        alpha = (2.0 * np.pi) / float(self.Num_lati)  # Angolo al centro.
        beta = alpha / 2.0  # Semi angolo al centro.
        if Circ == 'est':  # Poligono inscritto nella cricoferenza.
            if type(Lato) == type(None):  # Si parte dal raggio.
                self.Raggio = R
                Vr = np.array([0, R, 0])  # Vettore del raggio.
            else:  # Si parte dal lato.
                R_c = Lato * (np.sin(0.5 * (np.pi - alpha)) / np.sin(alpha))  # Misura del Raggio.
                self.Raggio = R_c
                Vr = np.array([0.0, R_c, 0.0])  # Vettore del raggio.
        else:  # Poligono circoinscritto. ( La circonferenza e' dentro il poligono )
            if type(Lato) == type(None):  # Si parte dal raggio.
                R_ = 2.0 * R * (np.sin(beta) / np.sin(alpha))
                self.Raggio = R_
                Vr = np.array([0, R_, 0])  # Vettore del raggio.
            else:  # Si parte dal lato.
                gamma = (np.pi - alpha) / 2.0
                R_ = Lato * (np.sin(gamma) / np.sin(alpha))  # Misura del Raggio.
                self.Raggio = R_
                Vr = np.array([0.0, R_, 0.0])  # Vettore del raggio.
        Vertici = []  # Lista vuota nella quale metter i punti che costituiscono i vertici del poligono.
        for n in range(Num_lati):
            ang_vert = n * alpha  # Angolo per ogni punto
            rot = np.array([0.0, 0.0, ang_vert])  # array degli angoli di rotazione
            V_raggio = Rotazione(Vr, rot, self.Centro)  # Vettore raggio del vertici
            Punto = self.Centro + V_raggio  # Calcolo delle coordinate del punto
            punto = Proiezione_Punto_su_Piano(Punto,
                                              piano_xy)  # Proietta tutti i punti sul piano xy im modo che siano tutti sullo stesso piano
            punto = np.trunc(np.round(punto, Cord_dec) * pow(10, Cord_dec)) / float(
                10 ** Cord_dec)  # arrotondamento corrdinate punto
            Vertici.append(punto)  # Aggiunta del punto alla lista dei vertici.
        ## Parte di dati del poligono normale
        self.Vertici = Vertici  # Vertici del poligono
        V = Vertici + Vertici  # Lista con i vertici ripetuti
        # Creazione dei Lati
        L = []  # Lista vuoto per raccogliere i lati
        for p in range(len(Vertici)):
            lato = Segmento(V[p], V[p + 1])  # lato
            L.append(lato)
        self.Lati = L
        # Angoli ai vertici interni al poligono
        self.Angolo_al_vertice_int = round(np.pi - alpha, Angl_dec)
        self.Angolo_al_vertice_est = round(np.pi + alpha, Angl_dec)
        self.Angolo_al_centro = round(alpha, Angl_dec)

    def Piano(self):
        """
        Calcola il piano su cui giace il poligono.
        """
        P_plan = Piano(self.Vertici[0], self.Centro, self.Vertici[1])
        return (P_plan)

    def Sposta(self, Spostamento):
        """
        Questa funzione permette di Spostare il poligono.
        Il punto di riferimento ? impostato sul centro del Poligono ma se si vuole si puo' indicare u altro punto generico.
        ARGOMENTI:
            Spostamento =>  Vettore di spostamento.
                            Se viene indicato un solo valore si avvra uno spostamento pari a quel valore su tutti gli assi.
                            Puo' essere  un array o un valore singolo.
            Riferimento =>  Punto di riferimento in base alla quale avviene lo spostamento. Se non viene indicato niente si usa il centro-
        """
        ### Vettore spostamento
        if type(Spostamento) == type(np.array([])):
            V_spo = Spostamento
        elif type(Spostamento) == type([]):
            V_spo = np.array(Spostamento)
        elif type(Spostamento) == type(int()) or type(Spostamento) == type(float()):
            V_spo = np.array([Spostamento, Spostamento, Spostamento])
        else:
            raise (Exception('%s non valido come spostamento' % (type(Spostamento))))
        ### Calcolo dello spostamento.
        for v in range(len(self.Vertici)):
            V = np.round(self.Vertici[v] + V_spo, Cord_dec)
            self.Vertici[v] = np.trunc(V * pow(10, Cord_dec)) / float(10 ** Cord_dec)
        self.Centro = np.round(self.Centro + V_spo, Cord_dec)
        ## Ricalcolo dei lati
        V = self.Vertici + self.Vertici
        L = []  # Lista vuoto per raccogliere i lati
        for p in range(len(self.Vertici)):
            lato = Segmento(V[p], V[p + 1])  # lato
            L.append(lato)
        self.Lati = L

    def Ruota(self, Angolo, Centro_rot=None, asse_rotazione=None):
        """
        Questa funzione permette di ruotare il poligono.
        _L_'angolo deve essere espresso in radianti.
        """
        ## Centro di Rotazione
        if type(Centro_rot) == type(None):
            Centro_rot = self.Centro
        if type(asse_rotazione) == type(None):
            asse_rotazione = piano_xy.Vers_N
            ## Angolo di rotazione
            if type(Angolo) in Tipi_di_numeri:
                angolo = np.array([0.0, 0.0, Angolo])
            elif type(Angolo) == type(list()):
                angolo = np.array(Angolo)
            elif type(Angolo) == type(np.array([])):
                angolo = Angolo
            else:
                raise (Exception("%s non puo' essere usato come angolo" % (type(Angolo))))
        else:
            angolo = Angolo
        ## Rotazione Poligono
        for v in range(len(self.Vertici)):
            if type(self.Vertici[v]) == type(O):
                Vert = self.Vertici[v]
            elif type(self.Vertici[v]) == type([]):
                Vert = np.array(self.Vertici[v])
            else:
                raise (TypeError('Tipologia vertice non riconosciuta.\n\t%s' % (type(self.Vertici[v]))))
            V = Rotazione(Vert, angolo, Centro_rot)
            self.Vertici[v] = np.trunc(np.round(V, Cord_dec) * pow(10, Cord_dec)) / float(10 ** Cord_dec)
        self.Centro = Rotazione(self.Centro, angolo, Centro_rot)
        self.Centro.round(Cord_dec)
        ## Ricalcolo dei lati
        V = self.Vertici + self.Vertici
        L = []  # Lista vuoto per raccogliere i lati
        for p in range(len(self.Vertici)):
            lato = Segmento(V[p], V[p + 1])  # lato
            L.append(lato)
        self.Lati = L


print('3')
Triangolo = Poligono_regolare(1.0)
print('4')
Quadrilatero = Poligono_regolare(1.0, Num_lati=4)
print('5')
Pentagono = Poligono_regolare(1.0, Num_lati=5)
print('6')
Esagono = Poligono_regolare(1.0, Num_lati=6)
print('7')
Ettagono = Poligono_regolare(1.0, Num_lati=7)
print('8')
Ottagono = Poligono_regolare(1.0, Num_lati=8)
print('9')
Ennagono = Poligono_regolare(1.0, Num_lati=9)
print('10')
Decagono = Poligono_regolare(1.0, Num_lati=10)
print('11')
Endecagono = Poligono_regolare(1.0, Num_lati=11)
print('12')
Dodedecagono = Poligono_regolare(1.0, Num_lati=12)
print('100')
Cento = Poligono_regolare(10, Num_lati=100)


def Poli_e_Retta(poligono, retta, stampa=True):
    """
    Questa funzione restituisce il rapport tra un poligono e una retta.
    La retta puo' essere esterna la poligono, passare per un vertice, essere parallela a un lato oppure passare all'interno del poligono.
    ARGOMENTI:
        poligono        =>      Poligono [ oggetto ]
        retta           =>      Retta [ oggetto ]
    RISULTATI:
        Vertici     => Se la retta passa solo per un vertice o se la retta e' incidente al piano del poligono all'interno del poligono.
        Segmento    => Se la retta passa l'interno del poligono o e' parallela a un lato
        None        => Quando la retta giace nel piano ma esternamente al poligono o intercetta il piano in un punto esterno al poligono.
        Distanza    => Quando la retta si trova all'esterno del piano del poligono.
    """
    ##    poligono = Triangolo
    Plan = poligono.Piano()  # Plane del poligono
    R_e_Plan = Retta_e_Piano(retta, Plan, False)  # Relazione tra il piano e la retta
    if type(R_e_Plan) == type(M[0]):  # La retta interceta il piano in un punto.
        if poligono.Punto_Inside(R_e_Plan, False) == False:  # Il punto di incidenza rimane all'esterno del poligono.
            if stampa == True:
                print("La retta non intercetta il poligono.")
            return (None)
        else:  # Il punto di incidenza si trova dentro al poligono.
            if stampa == True:
                print(" La retta intercetta il poligono in un punto")
            return (R_e_Plan)
    else:  # Retta // al piano
        if R_e_Plan == 0.00:  # La retta giace nel piano.
            Lista_punti = []  # Lista punti di intersezione tra retta e lati poligono
            for lato in poligono.Lati:
                lato_retta = Segmento_e_retta(retta, lato, False)  # Relazione tra il lato e la retta.
                if type(lato_retta) == type(seg_1):  # La retta passa per il lato.
                    if stampa == True:
                        print("LA retta passa per un lato")
                    return (lato)
                elif type(lato_retta) == type(M[0]):  # La retta intercetta il lato.
                    if len(Lista_punti) == 0:  # Primo punto della lista.
                        Lista_punti.append(lato_retta)
                    else:
                        D = []  # Lista delle distanze.
                        for n in range(len(Lista_punti)):  # Calcolo delle distanze dai punti della lista.
                            D.append(Dist(lato_retta, Lista_punti[n]))
                        ##                        print(D)
                        if any(D) <= pow(10,
                                         -precisione):  # una della distanze ? nulla => il punto e' gia' nella lista.
                            pass
                        else:  # Il punto non e' nella lista.
                            Lista_punti.append(lato_retta)  # Aggiunta del punto
                else:
                    pass
            # Varifica dello stato tra poligono e retta
            if len(Lista_punti) == 0:
                # La retta e' esterna al punto
                if stampa == True:
                    print("La retta passa all'esterno del poligono.")
                return (None)
            elif len(Lista_punti) == 1:
                if stampa == True:
                    print("La retta passa per un vertice.")
                return (Lista_punti)
            elif len(Lista_punti) == 2:
                if stampa == True:
                    print("La retta intercetta il poligono.")
                s = Segmento(Lista_punti[0], Lista_punti[1])
                return (s)
            else:
                # Eliminazione di eventuali doppioni soppravissuti alla prima selezione
                while len(Lista_punti) > 2:
                    for e in range(len(Lista_punti)):  # prendiamo un punto
                        for n in range(e + 1, len(Lista_punti)):  # prendiamo i punti successivi a quel punto.
                            d = Dist(Lista_punti[e], Lista_punti[n])  # coalcolo distanza tra i due punti
                            if d <= pow(10, -precisione):  # distanza inferiore alla precisione.
                                Lista_punti = Lista_punti[:n] + Lista_punti[
                                                                n + 1:]  # Ricostruzione lista senza il punto.
                if stampa == True:
                    print("La retta intercetta il poligono.")
                s = Segmento(Lista_punti[0], Lista_punti[1])
                return (s)

        else:
            if stampa == True:
                print("La retta giace all'esterno del piano")
            return (R_e_Plan)


def Poligono_e_Segmento(poligono, segmento, stampa=True):
    """
    Questo algoritmo permette di calcolare l'intersezione tra un poligono e un segmento.
    ARGOMENTI:
        poligono        ==>     un oggetto poligono.
        segmento        ==>     un oggetto segmento.
    RISULTATI:
        _L_'algoritmo restituisce:
            NONE        ==>     quando il segmento e il poligono non si intersecano.
            PUNTO       ==>     quando il segmento e il poligono si intersecano in un punto.
            SEGMENTO    ==>     quando il segmento e il poligono si intersecano creando un segmento.
    """
    retta = segmento.retta  # Retta su cui giace il poligono
    ret_pol = Poli_e_Retta(poligono, retta, False)  # rapporto tra retta e poligono.
    if type(ret_pol) == type(3.98):  # poligono e retta non si intersecano.
        if stampa == True:
            print('Il poligono e il segmento non si intersecano.')
        return (None)
    elif type(ret_pol) == type(None):  # il poligono e la retta non si intersecano.
        if stampa == True:
            print('Il poligono e il segmento non si intersecano.')
        return (None)
    elif type(ret_pol) == type(M[0]):  # La retta e il poligono si intersecano in un punto.
        if P_in_segm(ret_pol) == True:  # Il punto appartiene al segmento.
            if stampa == True:
                print('Il poligono e il segmento si intersecano in un punto.')
            return (ret_pol)
        else:
            if stampa == True:
                print('Il poligono e il segmento non si intersecano.')
            return (None)
    elif type(ret_pol) == type(seg_1):  # _L_'intersezione tra la retta e il poligono e' un segmento.
        rap_seg = Segmento_e_Segmento(ret_pol, segmento, False)
        if rap_seg[1] == 'Identici':
            if stampa == True:
                print('Il poligono e il segmento si intersecano.')
            return (ret_pol)
        elif rap_seg[1] == 'Consecutivi':
            if stampa == True:
                print('Il poligono e il segmento si intersecano in un punto.')
            return (rap_seg[0])
        elif rap_seg[1] == 'Disgiunti':
            if stampa == True:
                print('Il poligono e il segmento non si intersecano.')
            return (None)
        elif rap_seg[1] == 'Interno':
            if rap_seg[0] > 1:  # segmento contenuto nel segmento di intersezione tra retta e poligono.
                if stampa == True:
                    print('Il poligono e il segmento non si intersecano.')
                return (None)
            elif rap_seg < 1:  # il segmento di intersezione tra retta e poligono e' contenuto nel segmento passato come argomento.
                if stampa == True:
                    print('Il poligono e il segmento si intersecano.')
                return (ret_pol[0])
        elif rap_seg[1] == 'Sovrapposti':
            if stampa == True:
                print('Il poligono e il segmento si intersecano.')
            return (rap_seg[0])
        else:
            if stampa == True:
                print('Il poligono e il segmento non si intersecano.')
            return (None)


def Poligono_e_Piano(poligono, piano, stampa=True):
    """
    Questo algoritmo restituisce l'interseazione ta un poligono e un piano.
    ARGOMENTI:
        poligono        ==>     oggetto poligono
        piano           ==>     oggetto piano
    Risultati:
        Dipende dal raporto tra il piano e il poligono:
            se non si intersecano il risultato sara valore None
            se si intersecano in un punto sar? un punto e cosi via
    """
    poligono = Pentagono
    pol_piano = poligono.Piano()  # Plane su cui giace il poligono.
    raporto_piani = Piano_e_Piano(pol_piano,
                                  piano)  # rapporto tra il piano del Poligono e il piano indicato come variabile.
    if type(raporto_piani) == type(asse_x):  # I piani si intercettano in una retta.
        if stampa == True:
            print('Il piano intercetta il piano del poligono.')
        return (Poli_e_Retta(poligono, raporto_piani, stampa))  # restituisce il rap
    else:
        if stampa == True:
            print('Il piano del poligono e\' paralello al piano')
        if raporto_piani != 0:  # I piani non sono coincidenti.
            if stampa == True:
                print('Il piano del poligono dista %s dal piano.' % s(raporto_piani))
            return (None)
        else:
            if stampa == True:
                print('Il piano e il piano del poligono coincidono.')
            return (poligono)


def Poligono_e_Poligono(poligono1, poligono2, stampa=True):
    """
    Funzione che stabilisce la relazione tra due poligoni, restituendo l'ogetto generato dall'intersezione tra i due poligoni.
    La funzione restituisce:
        None:       se i due poligoni non si intersecano.
        Punto:      se i due poligoni si interscano solo in un punto.
        Segmento:   se i due poligoni si intersecano in un segmento.
        Poligono:   se i due poligoni si intersecano generando una sezione di piano comune.
    """
    rela_tra_piani = Piano_e_Piano(poligono1.Piano(), poligono2.Piano())
    if type(rela_tra_piani) != type(3.20):  # I piani dei due poligoni sono incidenti.
        r_e_Poli_1 = Poli_e_Retta(poligono1, rela_tra_piani)
        r_e_Poli_2 = Poli_e_Retta(poligono2, rela_tra_piani)
        if type(r_e_Poli_1) == type(None):
            # I poligoni non si intersecano.
            return (None)
        if type(r_e_Poli_2) == type(None):
            # I poligoni non si intersecano.
            return (None)
        if type(r_e_Poli_1) == type(M[0]):
            if type(r_e_Poli_2) == type(M[0]):
                d = Dist(r_e_Poli_1, r_e_Poli_2)
                if d <= pow(10, -4):
                    # I poligoni si intersecano in un solo punto.
                    return (r_e_Poli_1)
                else:
                    # I poligoni non si intersecano.
                    return (None)
            elif type(r_e_Poli_2) == type(seg_1):
                if P_in_segm(r_e_Poli_1, r_e_Poli_2):
                    # I poligoni si intersecano in un punto.
                    return (r_e_Poli_1)
                else:
                    # I poligoni non si intersecano.
                    return (None)
            else:
                raise (TypeError('Tipologia inaspetata per la relazione tra il poligoni'))
        elif type(r_e_Poli_1) == type(seg_1):
            if type(r_e_Poli_2) == type(M[0]):
                if P_in_segm(r_e_Poli_2, r_e_Poli_1):
                    return (r_e_Poli_2)
                else:
                    return (None)
            elif type(r_e_Poli_2) == type(seg_1):
                return (Segmento_e_Segmento(r_e_Poli_1, r_e_Poli_2))
        else:
            raise (Exception())
    else:
        if rela_tra_piani != 0.0:
            # I piani dei poligoni sono parelleli.
            return (None)
        else:  # I poligoni sono complanari.
            # Ordinamento dei poligoni.
            if poligono1.Area() >= poligono2.Area():
                Poli_1 = poligono1
                Poli_2 = poligono2
            else:
                Poli_1 = poligono2
                Poli_2 = poligono1
            D_Centri = Dist(Poli_1.Centro, Poli_2.Centro)  # Distanza tra i centri.
            # Lunghezze dei raggi massimi poli_1
            L_R_1 = []
            for l in Poli_1.Raggi():
                L_R_1.append(l.Lunghezza)
            R_1 = max(L_R_1)  # Raggio massimo poligono 1
            # Lunghezze dei raggi minimi poli_1
            L_r_1 = []
            for l in Poli_1.Lati:
                L_r_1.append(Dist_Punto_Retta(Poli_1.Centro, l.retta))
            r_1 = min(L_r_1)
            # Lunghezze dei raggi massimi poli_2
            L_R_2 = []
            for l in Poli_2.Raggi():
                L_R_2.append(l.Lunghezza)
            R_2 = max(L_R_2)  # Raggio massimo poligono 1
            # Lunghezze dei raggi minimi poli_1
            L_r_2 = []
            for l in Poli_2.Lati:
                L_r_2.append(Dist_Punto_Retta(Poli_2.Centro, l.retta))
            r_2 = min(L_r_2)
            if D_Centri > (R_2 + R_1):
                # I poligono sono l'uno esterno all'altro.
                return (None)
            elif D_Centri <= abs(r_1 - r_2):
                # I poligoni sono interni.
                return (Poli_2)
            else:
                L_vert = []
                for V in Poli_2.Vertici:
                    if Poli_1.Punto_Inside(V, False):
                        L_vert.append(V)
                for V in Poli_1.Vertici:
                    if Poli_2.Punto_Inside(V, False):
                        L_vert.append(V)
                for l_1 in Poli_1.Lati:
                    for l_2 in Poli_2.Lati:
                        intersezione = Segmento_e_Segmento(l_1, l_2, False)
                        if intersezione[1] == 'Identici':
                            L_vert.append(l_1.P1)
                            L_vert.append(l_1.P2)
                        elif intersezione[1] == 'Consecutivi':
                            L_vert.append(intersezione[0])
                        elif intersezione[1] == 'Sovrapposti':
                            L_vert.append(intersezione[0].P1)
                            L_vert.append(intersezione[0].P2)
                        elif intersezione[1] == 'Interno':
                            if intersezione[0] < 1:
                                L_vert.append(l_1.P1)
                                L_vert.append(l_1.P2)
                            else:
                                L_vert.append(l_2.P1)
                                L_vert.append(l_2.P2)
                L_v = Elimina_punti_doppi(L_vert)  # Lista vertici senza doppioni.
                if len(L_v) == 0:
                    if stampa == True:
                        print('I poligoni non si intersecano.\n')
                    return (None)
                elif len(L_v) == 1:
                    if stampa == True:
                        print('I poligoni sono tangenti.\n')
                    return (L_v[0])
                elif len(L_v) == 2:
                    if stampa == True:
                        print('I poligoni si intersecano in un segmento.\n')
                    return (Segmento(L_v[0], L_v[1]))
                else:
                    Lista_vertici = Oridanazione_polare(L_v)
                    Pol_comune = Poligono(Lista_vertici)
                    if stampa == True:
                        print('I poligoni hanno una regione in comune.\n')
                    return (Pol_comune)


def Intersezione(oggetto_1, oggetto_2):
    """
    Calcola l'intersezione tra i due ogetti forniti come argomenti. Gli ogetti possono essere rette, segmenti, piani o poligoni.
    """
    if type(oggetto_1) == type(asse_x):  # il primo oggeto e' una retta
        if type(oggetto_2) == type(asse_x):  # il secondo oggetto e' una retta
            P = Punto_by_2Retta(oggetto_1, oggetto_2)
            return (P)
        elif type(oggetto_2) == type(Segmento([0, 0, 0], [1, 1, 1])):  # il secondo oggetto e' un segmento
            P = Punto_by_2Retta(oggetto_1, oggetto_2.retta)
            if P_in_segm(P, oggetto_2):
                return (P)
            else:
                print("Intersezione Nulla")
                return (None)
        elif type(oggetto_2) == type(piano_xy):  # il secondo oggetto e' un piano
            P = Punto_by_plan_e_retta(oggetto_1, oggetto_2)
            return (P)
        elif issubclass(type(oggetto_2), Poligono):  # Il secondo oggetto e' un poligono.
            return (Poli_e_Retta(oggetto_2, oggetto_1))
        else:
            Stringa = "La funzione intersezione non e' definita per " + str(type(oggetto_2))
            raise (TypeError(Stringa))
    elif type(oggetto_1) == type(seg_1):  # il primo oggetto e' un segmento
        if type(oggetto_2) == type(asse_x):  # Il secondo oggetto e' una retta
            P = Punto_by_2Retta(oggetto_1.retta, oggetto_2)
            if P_in_segm(P, oggetto_1):
                return (P)
            else:
                return (None)
        elif type(oggetto_2) == type(Segmento([0, 0, 0], [1, 1, 1])):  # Il secondo oggetto e' un segmento
            P = Punto_by_2Retta(oggetto_1.retta, oggetto_2.retta)
            if P_in_segm(P, oggetto_1):
                if P_in_segm(P, oggetto_2):
                    return (P)
                else:
                    print("Intersezione Nulla")
                    return (None)
            else:
                print("Intersezione Nulla")
                return (None)
        elif type(oggetto_2) == type(piano_xy):  # Il secondo oggetto e' un piano
            P = Punto_by_plan_e_retta(oggetto_1.retta, oggetto_2)
            if P_in_segm(P, oggetto_1):
                return (P)
            else:
                print("Intersezione Nulla")
                return (None)
        elif issubclass(type(oggetto_2), Poligono):  # Il secondo oggetto e' un poligono.
            return (Poligono_e_Segmento(oggetto_2, oggetto_1))
        else:
            Stringa = "La funzione intersezione non e' definita per " + str(type(oggetto_2))
            raise (TypeError(Stringa))
    elif type(oggetto_1) == type(piano_xy):  # il primo oggetto e' un piano
        if type(oggetto_2) == type(asse_x):  # Il secondo oggetto e' una retta
            P = Punto_by_plan_e_retta(oggetto_2, oggetto_1)
            return (P)
        elif type(oggetto_2) == type(Segmento([0, 0, 0], [1, 1, 1])):  # Il secondo oggetto e' un segmento
            P = Punto_by_plan_e_retta(oggetto_2.retta, oggetto_1)
            if P_in_segm(P, oggetto_2):
                return (P)
            else:
                print("Intersezione Nulla")
                return (None)
        elif type(oggetto_2) == type(piano_xy):  # Il secondo oggetto e' un piano
            return (Retta_by_2Planes(oggetto_1, oggetto_2))
        elif issubclass(type(oggetto_2), Poligono):  # Il secondo oggetto e' un poligono.
            return (Poligono_e_Piano(oggetto_2, oggetto_1))
        else:
            Stringa = "La funzione intersezione non e' definita per " + str(type(oggetto_2))
            raise (TypeError(Stringa))
    elif issubclass(type(oggetto_1), Poligono):  # il primo oggetto e' un poligono.
        if type(oggetto_2) == type(asse_x):  # Il secondo oggetto e' una retta
            return (Poli_e_Retta(oggetto_1, oggetto_2))
        elif type(oggetto_2) == type(seg_1):  # Il secondo oggetto e' un segmento
            return (Poligono_e_Segmento(oggetto_1, oggetto_2))
        elif type(oggetto_2) == type(piano_xy):  # Il secondo oggetto e' un piano
            return (Poligono_e_Piano(oggetto_1, oggetto_2))
        elif issubclass(type(oggetto_2), Poligono):  # Il secondo oggetto e' un poligono.
            return (Poligono_e_Poligono(oggetto_2, oggetto_1))
        else:
            Stringa = "La funzione intersezione non e' definita per " + str(type(oggetto_2))
            raise (TypeError(Stringa))
    else:
        Stringa = "La funzione intersezione non e' definita per " + str(type(oggetto_1))
        raise (TypeError(Stringa))


def Ombra(oggetto, vettore_luce, piano):
    """
    Questa funzione permette di costruire la proiezione di un oggetto secondo un determinato vettore su di un altro oggetto.
    ARGOMENTI:
        oggetto:            Oggetto di cui si vuole avere la proiezione.
        vettore_luce:       Vettore direzionale della proiezione.
        piano:              Plane su cui viene proiettata l'ombra.
    RISULTATI:
        Oggetto che descrive l'ombra proietata dall'oggetto.
    """
    if type(oggetto) == type(M[0]):  ## Punto
        retta = Retta_by_P_e_Vet(oggetto, vettore_luce)  # Retta direttrice dell'ombra
        ombra = Punto_by_plan_e_retta(retta, piano)  # Punto ombra generato dal punto a dal vertice.
        return (ombra)
    elif type(oggetto) == type(asse_x):  ## REtta
        # Creazione punti ombra:
        dire_1 = Retta_by_P_e_Vet(oggetto.P1, vettore_luce)  # direzione dell'ombra del punto 1.
        P_ombra_1 = Punto_by_plan_e_retta(dire_1, piano)  # Punto ombra 1
        dire_2 = Retta_by_P_e_Vet(oggetto.P2, vettore_luce)  # direzione dell'ombra del punto 2.
        P_ombra_2 = Punto_by_plan_e_retta(dire_2, piano)  # Punto ombra 2
        retta_ombra = Retta(P_ombra_1, P_ombra_2)  # Ombra della retta
        return (retta_ombra)
    elif type(oggetto) == type(seg_1):  ##Segmento
        # Creazione punti ombra:
        dire_1 = Retta_by_P_e_Vet(oggetto.P1, vettore_luce)  # direzione dell'ombra del punto 1.
        P_ombra_1 = Punto_by_plan_e_retta(dire_1, piano)  # Punto ombra 1
        dire_2 = Retta_by_P_e_Vet(oggetto.P2, vettore_luce)  # direzione dell'ombra del punto 2.
        P_ombra_2 = Punto_by_plan_e_retta(dire_2, piano)  # Punto ombra 2
        seg_ombra = Segmento(P_ombra_1, P_ombra_2)  # Ombra del segmento.
        return (seg_ombra)
    elif type(oggetto) == type(piano_xy):
        return (oggetto)
    elif issubclass(type(oggetto_2), Poligono):
        P_ombra = []  # Ombra dei dei vertici del poligono.
        for v in oggetto.Vertici:
            direz_ombra = Retta_by_P_e_Vet(v, vettore_luce)
            P_ombra.append(Punto_by_plan_e_retta(direz_ombra, piano))
        Pol_ombra = Poligono(P_ombra)
        return (Pol_ombra)
    else:
        stringa = 'Questa funzione non e\' supportata per %s ' % (type(oggetto))
        print(stringa)


def Poligono_by_V_and_Center(n_lati, Vertice, Centro):
    """
    Questa funzione permette di Creare un poligono partendo dal numero dei lati, da un vertice e dal centro del Poligono.
    Argomenti:
        n_lati      ==>     il numero di lati di cui ? composto il poligono.
        Vertice     ==>     Vertice del poligono.
        Centro      ==>     Centro del Poligono.
    Risultati:
        Poligono    ==> oggetto istanza della classe Poligono.
    """
    Piano_Poligono = Piano(Centro, Vertice, O)  # Plane originale in cui ? contenuto il Poligono.
    Base = Piano_Poligono.Basi_ortonormali()  # Base relativa al piano del poligono
    V1 = np.array(Base.T[0])[0]
    V2 = np.array(Base.T[1])[0]
    # occorre invertire la base ( ruotarla di 90 gradi) in modo che il vertice dato sia un vertice del Poligono
    par_vertice = np.round(Piano_Poligono.Paramtri_punto(Vertice, [V2, V1]),
                           Vers_dec)  # Coordinate del Vertice rispetto al piano.
    new_vert = O
    for c in range(len(par_vertice)):
        new_vert[c] = par_vertice[c]
    raggio = Grandezza_vettore(new_vert)  # Raggio del Poligono
    Prov = Poligono_regolare(raggio, Num_lati=n_lati)  # Poligono provisorio nel piano.
    L_Vert = []  # Lista nella quale mettere i vertici con le coordinate nel sistema iniziale.
    for V in Prov.Vertici:
        v = Piano_Poligono.Punto_by_Param(V[1], V[0])
        L_Vert.append(v)
    Poli_def = Poligono(L_Vert)
    return (Poli_def)
