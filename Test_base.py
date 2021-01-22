# -------------------------------------------------------------------------------
# Name:         T.I.M. test
#
# Purpose:     modulo per testare le funzioni presenti nel modulo TIM
#
# Author:      SEBASTIANO
#
# Created:
# -------------------------------------------------------------------------------
import Librerie_gen as B
import time
import numpy as np


def S_m_rif(Lista, riferimento):
    """
    Questo funzione calcola lo scarto medio tra gli elementi della lista
    =============================================================================
    Argomenti;
    -----------------------------------------------------------------------------
        :param Lista:           ==>     Tipo: lista
                                        Lista dei valori

        :param riferimento:     ==>     Tipo: number
                                        Valore di riferimento dalla quale si
                                        vuole calcolare lo scarto quadratico
                                        medio
    =============================================================================
        :return:    ==>     Tipo: tuple
                    Valore 1 ==> Valore dello scarto medio
                    Valore 2 ==> Lista dei valori per ogni singolo scarto
    """
    s = []
    for ele in Lista:
        if not np.isnan(ele):
            s.append(ele - riferimento)
    s_medio = np.round(np.nanmean(s), B.precisione)
    return s_medio, np.array(s)


def E_m_rif(Lista, riferimento):
    """
    Questo funzione calcola lo scarto medio tra gli elementi della lista
    =============================================================================
    Argomenti;
    -----------------------------------------------------------------------------
        :param Lista:           ==>     Tipo: lista
                                        Lista dei valori

        :param riferimento:     ==>     Tipo: number
                                        Valore di riferimento dalla quale si
                                        vuole calcolare lo scarto quadratico
                                        medio
    =============================================================================
        :return:    ==>     Tipo: tuple
                    Valore 1 ==> Valore dell' errore medio
                    Valore 2 ==> Lista dei valori per ogni singolo errore
    """
    e_ = []
    for ele in Lista:
        if not np.isnan(ele):
            val = (ele - riferimento) / riferimento
            e_.append(val)
    e_medio = np.nanmean(e_)
    return e_medio, np.array(e_)


def Varianza_rif(l):
    """
    Questa funzione calcola lo scarto quadratico medio dei valori della lista.
    =============================================================================
    Argomenti;
    -----------------------------------------------------------------------------
        :param l:               ==>     Tipo: np.array
                                        Array degli scarti di cui si vuole
                                        calcolare la varianza
    =============================================================================
        :return:    ==>     Tipo: Number
                            Valore dells varianza
    """
    quadrati = l ** 2
    return np.sqrt(np.mean(quadrati))


def Cronometro_func(funzione, *argomenti):
    """
    Questa funzione serve per cronometrare determinare il tempo impiegato
    da una determinata funzione per completare le operazioni.
    =============================================================================
    Argomenti:
    -------------------------------------------------------------------------
    :param funzione:        ==>     Tipo: oggetto funzione
                                    Funzione che si vuole cronometrare

    :param argomenti:       ==>     Tipo: Lista
                                    Lista degli argomenti per la funzione.
    =========================================================================
        :return:      (tempo medio, varianza, minimo, massimo)
    """
    i = time.time()
    a = funzione.__call__(*argomenti)
    f = time.time()
    return f-i, a


def Statistiche(array, riferimento=None):
    """
    Questa funzione permette di ottenere i basilari valori per potere
    scrivere il rapporto statistico.
    =============================================================================
    Argomenti
    -----------------------------------------------------------------------------
        :param array:       ==>         Tipo: numpy.array
                                        Array dei valori che si voglio
                                        analizzare.
        :param riferimento:     ==>     Tipo: Valore
                                        Valore di riferimento rispetto a cui
                                        effettuare l'analisi statistica.
                                        Il valore e' preimpostato su none.
    =============================================================================
    :return: Dizionario delle variabili statistiche
    """
    dizionario = {}
    dizionario.__setitem__('scarti', S_m_rif(array, np.mean(array))[1])
    dizionario.__setitem__('medio', np.mean(array))
    dizionario.__setitem__('minimo', np.min(array))
    dizionario.__setitem__('massimo', np.max(array))
    f = np.max(array)-np.min(array)
    dizionario.__setitem__('forbice', f)
    dizionario.__setitem__('SQM', np.std(array))
    dizionario.__setitem__('var', np.var(array))
    if riferimento is None:
        dizionario.__setitem__('adimensionale', None)
    else:
        adim = {}
        e_ = E_m_rif(array, riferimento)
        adim.__setitem__('errori', e_[1])
        adim.__setitem__('medio', e_[0])
        adim.__setitem__('minimo', np.min(e_[1]))
        adim.__setitem__('massimo', np.max(e_[1]))
        adim.__setitem__('forbice', np.max(e_[1]) - np.min(e_[1]))
        adim.__setitem__('SQM', np.std(e_[1]))
        adim.__setitem__('var', Varianza_rif(e_[1]))
        dizionario.__setitem__('adimensionale', adim)
    return dizionario


def Scrivi_rap(diz_val, Nome_funzione):
    """
    Questa funzione crea il rapport per l'analisi temporale della funzione
    sulla quale si e' eseguito il test.
    =============================================================================
    Argomenti
    -----------------------------------------------------------------------------
    :param diz_val:         ==>     Tipo dict
                                    dizionario con i valori statistici
    :param Nome_funzione:   ==>     Tipo: Stringa
                                    Nome della funzione che si vuole
                                    analizzare.
    =============================================================================
    :return:    inserisce il rapporto del parametro rapporti del test
    """
    s = 'Analisi statistica di %s\n' % Nome_funzione
    minimo = str(diz_val['minimo'])
    massimo = str(diz_val['massimo'])
    s += 'Valore minimo: %s\t\tValore massimo:	%s\n' %(minimo, massimo)
    forbice = str(diz_val['forbice'])
    medio = str(diz_val['forbice'])
    s += 'Forbice:	%s\t\tValore medio:	%s\n' %(forbice, medio)
    s += 'Scarto Quadratico Medio:	%s\n' % diz_val['SQM']
    s += 'Varianza: %s\n' % diz_val['var']
    if diz_val['adimensionale'] is not None:
        a = diz_val['adimensionale']
        s += '\nAdimensionale|\n'
        forbice = str(a['medio'])
        medio = str(a['forbice'])
        s += 'Errore medio:	%s\tForbice:	%s\n' % (medio, forbice)
        minimo = str(a['minimo'])
        massimo = str(a['massimo'])
        s += 'Errore minimo:%s\tErrore massimo:	%s\n' % (minimo, massimo)
        s += 'Varianza: % s\n' % a['var']
    return s


'''
Rapporto accuratezza su %s
minimo:		%s	
forbice:	%s
media:		%s'''


class TestBase(object):
    """ Classe per creare un forma di test basilare in modo da poter avere
    una certa uniformita in tutti i test che eseguiremo sulle funzioni che
    creiamo
    """
    def __init__(self):
        """
        Funzione che crea il nucleo del test nel quale verrano delle variabili
        statistiche che ci permettono di stabilire se l'algoritmo che stiamo
        costruendo corrisponde alle nostre esigenze..
        """
        self.tempo = {}
        self.risultati = {}
        self.rapporti = {}
        self.lista_tempi = np.zeros(1)

    def Tempo_test(self, lista_dei_tempi):
        """
        Questa funzione permette di eseguire l'analisi statistica dei tempi
        di esecuzione della funzione testata.
        =========================================================================
        Argomenti
        -------------------------------------------------------------------------
        :param lista_dei_tempi:     ==>     Tipo: np.array
                                            Array dei tempi di esecuzione
        =========================================================================
        :return:    Scrive i risultai dell'analisi nella variabile self.tempo
        """
        if not isinstance(lista_dei_tempi, type(np.array)):
            np.array(lista_dei_tempi)
        self.tempo = Statistiche(lista_dei_tempi)
        lista_t = list(self.lista_tempi)
        for ele in lista_dei_tempi:
            lista_t.append(ele)
        self.lista_tempi = np.array(lista_t)

    def Analisi_accuratezza(self, valori, riferimento):
        """
        Questo metodo permette di valutare l'accuratezza di una determinata
        funzione.
        =========================================================================
        Argomenti
        -------------------------------------------------------------------------
            :param valori:          ==>     Tipo: Numpy.array
                                            Array dei valori che si sono
                                            ricavati dalla funzione analizzata.
            :param riferimento:     ==>     Tipo: Vario
                                            Valore che la funzione dovrebbe
                                            restituire in condizioni ideali.
        =========================================================================
        :return:        dizionario con i valori dell'analisi statistica che
        descrivono l'accuratezza della funzione.
        """
        # trasformazione da array valori in array percentuale.
        v_percentuali = valori/riferimento
        diz_stat = Statistiche(valori, riferimento)
        self.rapporti.__setitem__('accuratezza', diz_stat)


if __name__ == '__main__':
    v = np.random.random(25000)
    T = TestBase()
    T.Tempo_test(v)
    d = Statistiche(v, np.mean(v))
    print(Scrivi_rap(d, 'Prova'))
