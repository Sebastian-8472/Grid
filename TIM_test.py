# -------------------------------------------------------------------------------
# Name:         T.I.M. test
#
# Purpose:     modulo per testare le funzioni presenti nel modulo TIM
#
# Author:      SEBASTIANO
#
# Created:
# -------------------------------------------------------------------------------
import Test_base
import numpy as np
import TIM


class TestMediaMatrici(Test_base.TestBase):
    """
    Questa classe serve per generare i test sulle funzione che fanno la media
    tra due elementi consecutivi in una matrice.
    """
    def Test_shape(self, ripetizioni=100):
        """
        Questa funzione fa il test per controllare che la forma delle matrice
        che viene restituita da questa funzione e' quella aspettata.
        media su righe riduce il numero delle righe di uno
        media sulle colonne riduce il numero delle colonne di uno
        medie sulle diagonali riduce sia il numero delle righe e delle colonne di uno.
        =========================================================================
        Arogomenti:
        -------------------------------------------------------------------------
            :param ripetizioni:     ==>     Tipo: numero intero
                                            Numero che si vuole ripetere la
                                            struttura base del test.
                                            Il valore predefinito e'
                                            impostato su 100
        =========================================================================
        :return:
        """
        tempi = []
        risultati = []
        for n in range(ripetizioni):
            r = np.random.randint(1, 20)
            c = np.random.randint(1, 20)
            matrice = np.ones((r, c))
            tempi.append(Test_base.Cronometro_func(TIM.Media_su_riga,
                                                   matrice)[0])

        self.Tempo_test(tempi)

T_ = TestMediaMatrici()
T_.Test_shape()
