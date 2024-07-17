#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#

'''
author: Bill Kayser
date: 02/2021
descritpion: functions class of Parabolic Equation method
'''

# =========================================================================
# Importation des packages nécessaires
# =========================================================================

import numpy as np
import math as mt
import random as rnd
import scipy.sparse as sparse
import scipy.sparse.linalg as splinalg
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline

# =========================================================================
# Création de la calss de fonctions nécessaires au modèle Eq. Parabolique
# =========================================================================

class fonctions():

# =========================================================================
# Définition des paramètres globaux du modèle
# =========================================================================

    def __init__(self,freq,z_source,z_mic,x_dim,z_dim,discrx,discrz,T):
        '''
        Initialisation des paramètres communs à plusieurs fonctions de la classe.
        Les entrées du modèle sont définies dans le main et peuvent être modifiées par l'utilisateur.
        Les constantes et paramètres numériques ci dessous ne doivent pas être modifiés.
        '''
        ######################################
        ## Paramètres variables du modèle
        ######################################

        #----- La Source -----#
        self.freq = freq                                              # fréquence acoustique de la source (Hz)
        self.z_source = z_source                                    # hauteur de la source (m)
        self.z_mic = z_mic        #utilisé que pour plot et pour une autre variable, pas la peine                                      # hauteur du microphone (m)

        #-----  Le Domaine -----#
        self.x_dim = x_dim                                           # longueur du domaine suivant x (m)
        self.z_dim = z_dim                                            # hauteur du domaine suivant z (m)

        #----- Constantes -----#
        self.Patm = 101300                                          # pression atmosphérique (Pa)
        self.Rg = 286.6896552                               # constante des gaz parfait
        self.gam = 1.41                                          # rapport des chaleurs spécifiques
        self.temp0 = T + 273.15       #use T instead of self.T              # température (K)
        self.cel = mt.sqrt(self.gam*self.Rg*self.temp0)             # célérité adiabatique dans le milieu (m/s)
        self.k0 = 2*mt.pi*self.freq/self.cel        # nombre d'onde adiabatique

        # pression de saturation de la vapeur d'eau en fonction de la température
        #self.Psat = np.array([[-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40], [103,133.4,150,158.4,172.4,190,203.7,221.2,240,260.3,260,300,330,370,386,416,460,484,522,560,611,651,706,752,813,866,935,995,1073,1141,1228,1312,1402,1497,1598,1705,1818,1937,2063,2197,2338,2487,2643,2809,2983,3167,3360,3564,3780,4005,4243,4492,4755,5030,5319,5623,5941,6275,6625,6992,7375]])
        #self.idx = (np.abs(self.Psat-self.T)).argmin()          # cherche l'indice dans le tableau Psat correspondant à la bonne pression de saturation en fonction de T
        #self.Psat = self.Psat[1,self.idx]                                 # Valeur de la pression de staturation de la vapeur d'eau (Pa)
        #self.rho = ((1-(0.3783*self.hr*self.Psat/self.Patm))/(self.Rg*self.temp0)*self.Patm)  # densité de l'air (kg.m-3)

        Lambda = self.cel/self.freq                # longueur d'onde acoustique
        self.delta_x = Lambda/discrx       # pas de discretisation horizontale
        self.delta_z = Lambda/discrz       # pas de discretisation vertical

        self.nb_x = int(mt.floor(self.x_dim/self.delta_x)+1)   # nombre de point de discrétisation selon x
        self.nb_z = int(mt.floor(self.z_dim/self.delta_z)+1)   # nombre de point de discrétisation selon z (incluant la zone d'amortissement)

        self.z = np.arange(0,self.z_dim,self.delta_z) # vecteur hauteur z (m)

        #----- Ecran acoustique -----#
        # self.xe =  ??                                   # position de l'écran acoustique (m)
        # self.he = 0                                                     # hauteur de l'écran (m)
        # self.Je = int(mt.floor(self.xe/self.delta_x))               #indice de position des écrans
        # self.Jh = int(mt.floor(self.he/self.delta_z))               # indice de hauteur des écrans

    # end def __init__

# =========================================================================
# Calcul de l'admittance effective du sol b_eff = b_rug + b
# =========================================================================
    def impedance(self,x_imp,cgs1,cgs2,lc1,lc2,sigmah1,sigmah2):
        '''
        Calcul des propriétés acoustiques de sol à l'aide du modèle d'impédance de Miki (miki 1990)
        &  prise en compte de l'effet de la rugosité du sol avec le calcul de l'admittance effective (Bourlier, O.Faure 2014)
        Les hauteurs de rugosité du sol suivent un spectre Gaussien

        cgs : vecteur 1D contenant les valeurs de resistivité de chaque coté de la discontinuité de sol
        lc : vecteur 1D contenant les valeurs de longueur de corrélation de chaque coté de la discontinuité de sol
        sigmah : vecteur 1D contenant les valeurs des écarts-types de hauteur de rugosité de chaque coté de la discontinuité de sol

        Z : vecteur 1D contenant les valeurs d'impédances de sol de part et d'autre de la discontinuité
        B_rug : vecteur 1D contenant les valeurs d'admittance rugueuse de chaque coté de la discontinuité de sol
        B_eff : vecteur 1D contenant les valeurs d'admittance effective de chaque coté de la discontinuité de sol
        Admittance : vecteur 1D de taille nb_x qui contient les valeurs d'admittance effective en chaque point du sol
        kappa : scalaire, angle d'incidence / à la normale

        Les autres variables sont des variables de stockage ou de calcul mathématique (ex:intégrale)
        '''

        # initialisation des paramètres
        cgs = np.array([cgs1, cgs2])
        lc = np.array([lc1, lc2])
        sigmah = np.array([sigmah1, sigmah2])

        Z = np.zeros(np.shape(cgs),dtype=complex)
        B_rug = np.zeros(np.shape(cgs),dtype=complex)
        Beta_eff = np.zeros(np.shape(cgs),dtype=complex)
        self.Admittance = np.zeros(self.nb_x,dtype=complex)

        # angle d'incidence de l'onde sur le sol, pi/2 : incidence rasante
        incidence = mt.pi/2
        # angle par rapport à la normale
        kappa = self.k0*mt.sin(incidence)

        # paramètres numériques pour le calcul de B_rug
        s1 = 1
        s2 = -1

        # boucle sur la longueur du vecteur cgs
        for kk in range(len(cgs)):

            if lc[kk] == 0: #test si sol rugueux ou non

                B_rug[kk] = 0

            else: # si sol rugueux

                # Correction pour prise en compte de la rugosité (thèse Olivier Faure Annexe A)
                d_int_alpha = mt.sqrt(self.k0) / 100
                u_alpha = np.arange(0,mt.sqrt(self.k0)+d_int_alpha,d_int_alpha)

                integrande_alpha1 = np.real((1/(self.k0*np.sqrt(-u_alpha**2+2*self.k0)))*((self.k0**2+s1*kappa*(self.k0-u_alpha**2))**2)*(((sigmah[kk]**2)*(lc[kk]/(2*np.sqrt(mt.pi))))*np.exp(-(((kappa+s1*(self.k0-u_alpha**2))*lc[kk])**2)/4)))
                integrande_alpha2 = np.real((1/(self.k0*np.sqrt(-u_alpha**2+2*self.k0)))*((self.k0**2+s2*kappa*(self.k0-u_alpha**2))**2)*(((sigmah[kk]**2)*(lc[kk]/(2*np.sqrt(mt.pi))))*np.exp(-(((kappa+s2*(self.k0-u_alpha**2))*lc[kk])**2)/4)))

                alpha1 = np.trapz(integrande_alpha1,u_alpha)
                alpha2 = np.trapz(integrande_alpha2,u_alpha)
                alpha = alpha1 + alpha2

                d_int_beta =(6/lc[kk])/100
                u_beta = np.arange(0,6/lc[kk]+d_int_beta,d_int_beta)

                integrande_beta1 = np.real((1/(self.k0*np.sqrt(self.k0**2+u_beta**2)))*((self.k0**2+s1*kappa*np.sqrt(self.k0**2+u_beta**2))**2)*(((sigmah[kk]**2)*(lc[kk]/(2*np.sqrt(mt.pi))))*np.exp(-(((kappa+s1*np.sqrt(self.k0**2+u_beta**2))*lc[kk])**2)/4)))
                integrande_beta2 = np.real((1/(self.k0*np.sqrt(self.k0**2+u_beta**2)))*((self.k0**2+s2*kappa*np.sqrt(self.k0**2+u_beta**2))**2)*(((sigmah[kk]**2)*(lc[kk]/(2*np.sqrt(mt.pi))))*np.exp(-(((kappa+s2*np.sqrt(self.k0**2+u_beta**2))*lc[kk])**2)/4)))

                beta1 = -np.trapz(integrande_beta1,u_beta)
                beta2 = -np.trapz(integrande_beta2,u_beta)
                beta = beta1 + beta2

                B_rug[kk] = alpha + 1.j*beta


            if cgs[kk] >= 100000: # cas parfaitement refléchissant

                Beta_eff[kk] = 0 + B_rug[kk]

            else: # cas avec absorption

                Z[kk] = 1+5.50*(self.freq/cgs[kk])**(-0.632)+1.j*8.43*(self.freq/cgs[kk])**(-0.632) # modèle d'impédance de Miki
                Beta_eff[kk] = 1/Z[kk] + B_rug[kk]


            self.Rp = (1-Beta_eff[0])/(1+Beta_eff[0])           # coefficient de reflexion sol 1
            self.coeff_absorption = 1-abs(self.Rp)**2                     #coefficient d'absorption sol 1
            kimp = int(round(x_imp/self.delta_x))           # indice suivant x du changement d'impedance, doit être un int car c'est un indice

            self.Admittance[0:kimp] = Beta_eff[0]                    # admittance effective du 1er sol (jusqu'à kimp)
            self.Admittance[kimp+1:self.nb_x] = Beta_eff[1]  # admittance effective du 2e sol (après kimp)

    # end def impedance

# =========================================================================
# Calcul du starter WAPE
# =========================================================================
    def starterWAPE(self):
        '''
        # Starter pour la Wide-angle PE
        # voir livre Salomons computational atmospheric acoustics p.179

        A0,A2,B : scalaires qui sont des coefficients pour le calcul de U0
        U0 : vecteur 1D sur la dimension z qui correspond au champ de pression initial en x = 0
        '''
        A0 = 1.3717
        A2 = -0.3701
        B = 3

        self.U0 = mt.sqrt(self.k0)*(A0+A2*self.k0**2*(self.z - self.z_source)**2)*np.exp(-(self.k0**2 *(self.z-self.z_source)**2)/B) + self.Rp*mt.sqrt(self.k0)*(A0+A2*self.k0**2 *(self.z + self.z_source)**2)*np.exp(-(self.k0**2*(self.z+self.z_source)**2)/B)

        # end def starterWAPE

# =========================================================================
# Calcul des profils météorologiques moyens verticaux (sur z) forme logarithmique
# =========================================================================
    def profils_meteo(self,au,aT,theta,hv):
        '''
        Fonction permettant de calculer les profils verticaux de vent, de température et de célérité effective.
        L'indice de réfraction acoustique du milieu est retourné, pour calcul du terme epsilon lors de la résolution matricielle
        (voir fonction calc)

        au : scalaire, coefficient de profil log (m/s)
        aT : scalaire, coefficient de profil log (K/m)
        theta : scalaire, angle de propagation (°)
        hv : scalaire, hauteur de végétation (m)
        refraction : scalaire, indice acoustique du milieu
        ceff : vecteur de taille nb_z, profil vertical de célérité effective
        u_z : vecteur de taille nb_z, profil vertical de vent
        t_z : vecteur de taille nb_z, profil vertical de température
        '''

        d = 0.66*hv                # hauteur de déplacement des profils météo (m)
        z0 = 0.13*hv + 0.01      # hauteur de rugosité des profils météo (m)
        
        if au == 0 and aT == 0: # cas sans profil micro-meteo (homogène)

            self.refraction = 1                     # indice de refraction constant = 1
            self.ceff = np.ones(np.shape(self.z))*self.cel                    # pas de profil de célérité, ceff = cel

        else: # cas hétérogène avec profils log

            indice_d = int(np.argmin(abs(self.z-d)))        # trouve l'indice de z correspondant à la hauteur de déplacement, doit être un int

            u_z = au*(np.log((self.z-d)/z0 +1))                        # profil vertical de vent
            t_z = self.temp0 + aT*(np.log((self.z-d)/z0 +1))  # profil vertical de température

            self.ceff = np.sqrt(self.gam*self.Rg*(t_z))+u_z*mt.cos(mt.pi*theta/180)    # profil vertical de célérité effective
            self.ceff[0:indice_d] = self.cel                        # ceff = cel jusqu'à la hauteur de déplacement d (lié à l'influence de la végétation)
            self.refraction = (self.cel/self.ceff)              # indice de refraction acoustique du milieu

    # end def profils
    
# =========================================================================
# Fonction qui permet d'initialiser les paramètres nécessaires au calcul
# du champ thermique turbulent
# =========================================================================

    def init_turb_therm(self,nmode,mu2,bigl,littlel):
        '''
        initialisation des parametres servant au calcul du champ turbulent thermique
        avec un spectre de Von Karman (cf thèse P.Chevret 1994, ECN Lyon)
        
        nmode : scalaire, nombre de modes de Fourier aléatoire
        mu2   : scalaire, intensité de la turbulence
        bigl  : scalaire, échelle externe de la turbulence
        littlel : scalaire, échelle interne de la turbulence
        kmin    : scalaire, petit nombre d'onde turbulent
        kmax    : scalaire, grand nombre d'onde turbulent
        
        Les paramètres de sortie sont :
        kct, kst, tmcf, tmsf  : paramètres de la turbulence thermique
        '''
        
        # initialisation des paramètres
        # Coordonnées sphériques :
        theta = np.zeros(nmode)# 1er angle de propagation du mode, aléatoire sur [0;2pi] pour isotropie de la turbulence
        phi = np.zeros(nmode) # 2e agle de propagation du mode, aléatoire sur [0;2pi] pour assurer invariance par translation des propriétés du milieu
        
        kmin = 0.1/bigl # petit nombre d'onde turbulente
        kmax = 6/bigl # grand nombre d'onde turbulente
        Km = 5.92/littlel # constante voir eq1.5 p13 thèse Chevret

        # discrétisation logarithmique sur Nmode, entre kmin et kmax, du spectre d'énergie G(k) cf fig. 4.5 p48 thèse Chevret
        dkl = (mt.log(kmax)-mt.log(kmin))/(nmode-1) # utilisation d'une discrétisation log car les nombres d'onde turbulents bas sont plus influents sur l'acoustique (cf article Karweit 1991)
        
        dk = kmin*(mt.exp(dkl)- 1)*np.exp(np.arange(0,nmode,1)*dkl) # distance entre chaque mode de Fourier aléatoire , dans l'espace des nombres d'onde K

        k = kmin*np.exp(np.arange(0,nmode,1)*dkl) # nombres d'onde turbulent
        
        r = k/Km # variables pour calcul

        #tm = np.sqrt(np.exp(-r*r)*dk*k*(1+k*k/(bigl**2))**(-11/6)) # voir eq. 4.11 - 4.13 p47 thèse Chevret. Soucis 1 + K2/LO2 ??
        tm = np.sqrt(np.exp(-r*r)*dk*k*(k*k+(1/(bigl**2)))**(-11/6)) # voir eq. 4.11 - 4.13 p47 thèse Chevret. Soucis 1 + K2/LO2 ??

        for imode in range(0, nmode):
            
            #rnd.seed(imode) ############
            #coordonnées sphériques dans l'espace des nombres d'onde k
            theta[imode] = rnd.random()*2*mt.pi # 1er angle de propagation du mode cf fig4.1 p 43 these Chevret
            
            #rnd.seed(imode) ############
            phi[imode] = rnd.random()*2*mt.pi # 2e angle de propagation du mode cf fig4.1 p 43 these Chevret

        # end for
        
        kct = k*np.cos(theta)
        kst = k*np.sin(theta)
        tmcf = tm*np.cos(phi)
        tmsf = tm*np.sin(phi)
        somme = 0
        
        for i in range(0,1000):
            
            #rnd.seed(i) ############
            x = 100*rnd.random()
            
            #rnd.seed(i) ############
            z = 100*rnd.random()
            
            # calcul la temperature de la turbulence en un point x,z
            phase = kct*x + kst*z
            cp = np.cos(phase)
            sp = np.sin(phase)
            ttemp = np.sum(tmcf*cp-tmsf*sp)

            somme = somme+ttemp**2
        # end for
        
        inorm = 1/np.sqrt(somme/1000)
        tm = tm*inorm*mt.sqrt(mu2)
        
        self.kct = k*np.cos(theta)
        self.kst = k*np.sin(theta)
        self.tmcf = tm*np.cos(phi)
        self.tmsf = tm*np.sin(phi)

    # end def init_turb_therm

# =========================================================================
# Mode de Fourier Aléatoires pour le calcul de la partie stochastique de l'indice de réfraction epsilon
# =========================================================================    

    def MFA(self,x):
        '''
        Calcul de la temperature turbulente au point x,z pour la hauteur complète du domaine
        (cf thèse P.Chevret 1994, ECN Lyon)
        
        z  : vecteur, axe vertical du domaine en cours
        x  : scalaire, pas d'avancement suivant l'axe x
        kct,kst,tmcf,tmsf : vecteurs, paramètres de calcul de la turbulence (cf init_turb_therm)

        '''
        dim = len(self.z)
        phase = np.zeros((dim,len(self.kct)))

        for i in range(0,dim):
            phase[i,:] = self.kct*x+self.kst*self.z[i]
        #end
        cp = np.cos(phase)
        sp = np.sin(phase)
        tmcf_mat = np.ones((dim,1))*self.tmcf
        tmsf_mat = np.ones((dim,1))*self.tmsf   
        temp_turb = np.sum((tmcf_mat*cp-tmsf_mat*sp),axis=1)

        return temp_turb

    # end def MFA
# =========================================================================
# Calcul de l'amortissement en haut du domaine (z)
# =========================================================================
    def amortissement(self):
        '''
        Permet de prendre en compte l'amortissement en haut du domaine pour éviter les reflexions parasites

        haut_a : scalaire 0<haut_a<1
        coeff_a : scalaire, coeff d'amortissement
        za : scalaire, indice sur z à partir duquel l'amortissement commence
        idex : vecteur 1D sur z, contient les indices où l'amortissement a lieu
        amort : vecteur 1D sur z, coefficient d'amortissement sur une"tranche" verticale du domaine, il faut multiplier chaque "tranche" p_ij par ce vecteur
        '''

        haut_a = 0.6                  # hauteur du domaine à laquelle on commence l'amortissement (0.8 = 80%)
        coeff_a = 10                 # coefficient d'amortissement

        za = haut_a*self.nb_z          # hauteur à partir de laquelle on applique l'amortissement
        idex = np.arange(mt.ceil(za)-1,self.nb_z-1) # renvoie les indices du vecteur z où on applique l'amortissement (> za)
        idex = idex.astype(int)         # tableau qui doit contenir des entiers car ce sont des indices

        self.amort = np.ones(self.nb_z)

        self.amort[idex] = self.amort[idex]* np.exp( -((idex-za)/(coeff_a*(self.nb_z-idex)))**2) # calcul de l'amortissement pour chaque z
        self.amort[self.nb_z-1] = 0        # dernier élement nul

    # end def amortissement

# =========================================================================
# Coeur de calcul suivant la méthode de Crank Nicholson et développement Padé(1,1)
# =========================================================================
    def calc_pade_one_one(self):
        '''
        Résolution matricielle utilisant le schéma numérique de Crank-Nicholson
        et l'approximation du terme Q par un développement Padé(1,1)
        (blairon_phd2002, Eqs.(2.15), p.44
        salomons computational atmospheric acoustic)

        p_ij : matrice 2D de calcul, de taille (nb_z,nb_x)
        pp : matrice 2D,champ de pression 3D après correction /sqrt(x)
        epsilon : vecteur, indice acoustique du milieu suivant z

        p1, q1 : scalaire, coefficient Padé(1,1)
        a,b,c,d,e,f : scalaires complexes, termes des diagonales des matrices Aij Bij
        Aij & Bij : matrice 2D de la méthode Crank Nicholson, taille (nb_z,nb_z)
        nb_x, nb_z : scalaires, nombres de points de discrétisations spatiaux
        '''

        self.pp =  np.zeros((self.nb_z,self.nb_x),dtype = complex)     # initialisation de la variable globale qui correspond au champ de pression final
        p_ij = np.zeros((self.nb_z,self.nb_x),dtype = complex)     # initialisation de la variable champ de pression pour calcul dans les boucles
        self.pp[:,0] = self.U0  # à x = 0 la pression correspond au starter U0 sur tout z
        p_ij[:,0] = self.U0  # idem

        epsilon = (self.cel/self.ceff)**2 -1                        # indice de refraction acoustique du milieu

        # Padé(1,1)
        sig = 1j * self.k0 * self.delta_x
        p1 = (1. + sig) / 4.
        q1 = (1. - sig) / 4.

        ## Initialisation des coefficients des matrices Aij & Bij
        # matrice Aij
        aj = q1 * (1. / (self.k0 * self.delta_z) ** 2)                          # below diagonal term
        bj = 1. + q1 * (epsilon -2. / (self.k0 * self.delta_z) ** 2)    # diagonal term
        cj = aj                                                                                         # above diagonal term
        # matrice Bij
        dj = p1 * (1. / (self.k0 * self.delta_z) ** 2)                          # below diagonal term
        ej = 1. + p1 * (epsilon -2. / (self.k0 * self.delta_z) ** 2)    # diagonal term
        fj = dj                                                                                         # above diagonal term

        # remplissage des matrices avec méthode de Cranck-Nicholson       sparse.csr_matrix
        Aij = np.zeros((self.nb_z, self.nb_z), dtype = complex)     # ??? utilisation de matrices sparses pour éviter surcharge mémoire
        Bij = np.zeros((self.nb_z, self.nb_z), dtype = complex)
        rng = np.arange(self.nb_z)

        Aij[rng[:-1], rng[:-1] + 1] = aj
        Aij[rng, rng] = bj
        Aij[rng[1:], rng[1:] - 1] = cj

        Bij[rng[:-1], rng[:-1] + 1] = dj
        Bij[rng, rng] = ej
        Bij[rng[1:], rng[1:] - 1] = fj

        # condition aux limites (sol)
        bj_0 = bj[0]+2.*1j*self.k0*self.Admittance[0]*self.delta_z * q1 * (1. / (self.k0 * self.delta_z) ** 2)
        cj_0 = 2. * q1 * (1. / (self.k0 * self.delta_z) ** 2)

        Aij[0, 0] = bj_0
        Aij[0, 1] = cj_0

        ej_0 = ej[0]+2.*1j*self.k0*self.Admittance[0]*self.delta_z * p1 * (1. / (self.k0 * self.delta_z) ** 2)
        fj_0 = 2. * p1 * (1. / (self.k0 * self.delta_z) ** 2)
        Bij[0, 0] = ej_0
        Bij[0, 1] = fj_0

        Aij = sparse.csr_matrix(Aij) #csr
        Bij = sparse.csr_matrix(Bij)

        # résolution matriciel à chaque pas ix, pour tout z
        for ix in range(0, self.nb_x-1): # boucle sur les indices x, indice commence à 0, donc de 0 à nb_x -1

            p_ij[:,ix+1] = Bij.dot(p_ij[:,ix])
            p_ij[:,ix+1] = splinalg.spsolve(Aij,p_ij[:,ix+1])

            p_ij[:,ix+1] = p_ij[:,ix+1] * self.amort   # prise en compte de l'amortissement (pour eviter reflexion parasite en haut du domaine) commence à l'indice ix+1 (voir matrice pp qui n'est pas bien amortie en haut du domaine si ix)

            #Passage pression 3D
            self.pp[:,ix] = np.abs(p_ij[:,ix+1])/mt.sqrt((ix+1)*self.delta_x) #p3D = p2/sqrt(x) hypothèse d'axisymétrie suivant y
            # /!\ si on écrit p_ij[:,ix] il y de mauvais résultats à la colone 0, & un décalage  d'une colone // modèle Matlab...??

        # calcul de la dernière colone de la matrice pp, elle  n'est pas calculée dans la boucle précédente
        self.pp[:,self.nb_x-1] = np.abs(p_ij[:,self.nb_x-1])/mt.sqrt((self.nb_x)*self.delta_x)


        return self.pp
    # end def calc_pade_one_one
    
# =========================================================================
# Coeur de calcul suivant la méthode de Crank Nicholson et développement Padé(1,1)
# avec prise en compte de la turbulence thermique à l'aide des modes de Fourrier aléatoires
# =========================================================================
    def calc_pade_one_one_turbulence(self, nb_ind,nmode,mu2,bigl,littlel):
        '''
        Résolution matricielle utilisant le schéma numérique de Crank-Nicholson
        et l'approximation du terme Q par un développement Padé(1,1)
        (blairon_phd2002, Eqs.(2.15), p.44
        salomons computational atmospheric acoustic)
        
        Les modes de Fourrier aléatoires (thèse chevret 1994) permettent de prendre
        en compte la turbulence thermiques. L'indice du milieu epsilon doit être 
        recalculé à chaque pas d'avancement suivant x, ce qui implique un temps de calcul long car
        les matrices Aij et Bij doivent être re-remplies à chaque pas dx

        nb_ind : scalaire, nombre d'itération de la turbulence sur lequel on moyenne le champ de pression pp
        
        p_ij : matrice 2D de calcul, de taille (nb_z,nb_x)
        pp : matrice 2D,champ de pression 3D après correction /sqrt(x)
        epsilon : vecteur, indice acoustique du milieu (partie moyenne + aléatoire MFA) suivant z

        p1, q1 : scalaire, coefficient Padé(1,1)
        a,b,c,d,e,f : scalaires complexes, termes des diagonales des matrices Aij Bij
        Aij & Bij : matrice 2D de la méthode Crank Nicholson, taille (nb_z,nb_z)
        nb_x, nb_z : scalaires, nombres de points de discrétisations spatiaux
        '''

        self.pp =  np.zeros((self.nb_z,self.nb_x),dtype = complex)     # initialisation de la variable globale qui correspond au champ de pression final
        p_ij = np.zeros((self.nb_z,self.nb_x),dtype = complex)     # initialisation de la variable champ de pression pour calcul dans les boucles
        self.pp[:,0] = self.U0  # à x = 0 la pression correspond au starter U0 sur tout z
        p_ij[:,0] = self.U0  # idem
        pp_turb = np.zeros((self.nb_z,self.nb_x),dtype = complex) # permet de sommer les champs de pressions lros de l'itération de la turbulence
                
        # Padé(1,1)
        sig = 1j * self.k0 * self.delta_x
        p1 = (1. + sig) / 4.
        q1 = (1. - sig) / 4.

        ## Initialisation des coefficients des matrices Aij & Bij
        # matrice Aij
        aj = q1 * (1. / (self.k0 * self.delta_z) ** 2)                          # below diagonal term
        # b dans la boucle iturb
        cj = aj                                                                                         # above diagonal term
        # matrice Bij
        dj = p1 * (1. / (self.k0 * self.delta_z) ** 2)                          # below diagonal term
        # e dans la boucle iturb
        fj = dj                                                                                         # above diagonal term

        # remplissage des matrices avec méthode de Cranck-Nicholson       sparse.csr_matrix
        Aij = np.zeros((self.nb_z, self.nb_z), dtype = complex)     # ??? utilisation de matrices sparses pour éviter surcharge mémoire
        Bij = np.zeros((self.nb_z, self.nb_z), dtype = complex)
        rng = np.arange(self.nb_z)

        Aij[rng[:-1], rng[:-1] + 1] = aj
        # b dans la boucle iturb
        Aij[rng[1:], rng[1:] - 1] = cj

        Bij[rng[:-1], rng[:-1] + 1] = dj
        # e dans la boucle iturb
        Bij[rng[1:], rng[1:] - 1] = fj

        # condition aux limites (sol) (la suite dans la boucle sur iturb car epsilon change à chaque itération)
        cj_0 = 2. * q1 * (1. / (self.k0 * self.delta_z) ** 2)
        fj_0 = 2. * p1 * (1. / (self.k0 * self.delta_z) ** 2)
        
        Aij = sparse.csr_matrix(Aij) #csr
        Bij = sparse.csr_matrix(Bij)
        
        # boucle sur le nombre d'itérations de la turbulence
        for iturb in range(nb_ind):

            #initialisation du champ thermique turbulent avec les MFA
            self.init_turb_therm(nmode,mu2,bigl,littlel)     

            # résolution matriciel à chaque pas ix, pour tout z
            for ix in range(0, self.nb_x-1): # boucle sur les indices x, indice commence à 0, donc de 0 à nb_x -1

                # ------------------------------------------------------------------
                # La turbulence nécessite de calculer l'indice acoustique epsilon à chaque pas x
                # et donc de re-remplir les termes b, e, et les matrices A, B
            
                temp_turb = self.MFA(ix*self.delta_x)
                epsilon = ((self.cel/self.ceff)-temp_turb)**2 -1                        # indice de refraction acoustique du milieu

                bj = 1. + q1 * (epsilon -2. / (self.k0 * self.delta_z) ** 2)    # diagonal term
                bj_0 = bj[0]+2.*1j*self.k0*self.Admittance[0]*self.delta_z * q1 * (1. / (self.k0 * self.delta_z) ** 2)
                
                Aij[rng, rng] = bj
                Aij[0, 0] = bj_0
                Aij[0, 1] = cj_0

                ej = 1. + p1 * (epsilon -2. / (self.k0 * self.delta_z) ** 2)    # diagonal term
                ej_0 = ej[0]+2.*1j*self.k0*self.Admittance[0]*self.delta_z * p1 * (1. / (self.k0 * self.delta_z) ** 2)
                
                Bij[rng, rng] = ej
                Bij[0, 0] = ej_0
                Bij[0, 1] = fj_0
                # ------------------------------------------------------------------

                p_ij[:,ix+1] = Bij.dot(p_ij[:,ix])
                p_ij[:,ix+1] = splinalg.spsolve(Aij,p_ij[:,ix+1])

                p_ij[:,ix+1] = p_ij[:,ix+1] * self.amort   # prise en compte de l'amortissement (pour eviter reflexion parasite en haut du domaine) commence à l'indice ix+1 (voir matrice pp qui n'est pas bien amortie en haut du domaine si ix)

                #Passage pression 3D
                self.pp[:,ix] = np.abs(p_ij[:,ix+1])/mt.sqrt((ix+1)*self.delta_x) #p3D = p2/sqrt(x) hypothèse d'axisymétrie suivant y
                # /!\ si on écrit p_ij[:,ix] il y de mauvais résultats à la colone 0, & un décalage  d'une colone // modèle Matlab...??

            # calcul de la dernière colone de la matrice pp, elle  n'est pas calculée dans la boucle précédente
            self.pp[:,self.nb_x-1] = np.abs(p_ij[:,self.nb_x-1])/mt.sqrt((self.nb_x)*self.delta_x)
            
            # addition des champs de pression pour chaque itération de turbulence
            pp_turb = np.add(pp_turb,self.pp)

        self.pp = pp_turb/nb_ind # moyenne par le nombre d'itération pour avoir un champ de pression moyen
        
    # end def calc_pade_one_one

# =========================================================================
# Calcul de la divergence géométrique en tout point du domaine
# =========================================================================
    def divergence_geometrique(self,stock_x,stock_z):
        '''
        fonction qui calcule la correction divergence géométrique
        en tout point (x,z) du domaine

        x1 : vecteur 1D (m) de 0 à x_dim par pas delta_x
        z1 : vecteur 1D (m) de 0 à z_dim par pas delta_z

        Zi, Xi : matrices 2D (m) permettant le calcul matriciel de la divergence géométrique Ri
        Ri : matrice 2D (z,x) (m) de distance à la source en tout point du domaine

        '''
        x = np.arange(0,self.x_dim,stock_x)
        z = np.arange(0,self.z_dim,stock_z)

        Zi,Xi = np.meshgrid(z,x)
        Zi = Zi.transpose()
        Xi = Xi.transpose()

        self.Ri = np.sqrt((Xi)*(Xi)+(Zi-self.z_source)*(Zi-self.z_source))

    # end def divergence_geometrique

# =========================================================================
# Calcul de l'absorption atmosphérique
# =========================================================================
    def abs_atmos(self,hr):
        '''
        fonction permettant de prendre en compte l'absorption atmosphérique
        d'après la norme iso 9613.
        Dépend de la fréquence, température, hygrométrie et pression.

        pr : scalaire en kPa
        pa : scalaire en kPa
        T0 : scalaire en K
        C : scalaire
        hr : scalaire 0<hr<1
        h : scalaire

        fro : scalaire
        frn : scalaire
        alpha : scalaire

        self.Ri : matrice 2D qui est la divergence géométrique en tout point du domaine
        
        atmos: matrice 2D (nb_z,nb_x), coefficient d'absorption en tout point du domaine (dépend de Ri et alpha)
        '''

        pr = 101.325 # pression atmosphérique de référence (kPa)
        pa = self.Patm/1000 # pression atmosphérique (kPa)
        Tr = 293.15   # température de référence (20°) en K

        T0 = 273.15 # 0 degré (K)

        C = -6.8346*(T0/self.temp0)**1.261+4.6151
        h = hr*10**C*pr/pa      # fraction molaire de vapeur d'eau

        # Eq (3) (4) et (5) p3 de la norme ISO 9613-1

        fro = pa/pr*(24+4.04*1e4*h*(0.02+h)/(0.391+h)) # Frequence de relaxation O2

        frn = pa/pr*(self.temp0/Tr)**(-0.5)*(9+280*h*mt.exp(-4.170*((self.temp0/Tr)**(-1/3)-1))) # Frequence de relaxation N

        # Attenuation en dB/m
        alpha = 8.686*self.freq**2*((1.84*1e-11*(pr/pa)*(self.temp0/Tr)**(0.5))+(self.temp0/Tr)**(-5/2)*(0.01275*(mt.exp(-2239.1/self.temp0))*(fro+(self.freq**2/fro))**(-1) + 0.1068*(mt.exp(-3352.0/self.temp0))*(frn+(self.freq**2/frn))**(-1)))
        
        # calcul de l'absorption en tout point du domaine z,x (Ri est la divergence géométrique)      
        self.atmos = self.Ri*alpha 

    # end def abs_atmos
    
# =========================================================================
# Permet d'interpoler les points de la carte de bruit tous les 0.5m et non en lambda/20
# =========================================================================
    
    def interpolation(self,stock_x,stock_z):

        '''
        Interpolation linéaire equivalente à interp2(z, xi, yi,'linear') en MATLAB
        oInterpP : objet Interpolation Pression
        x : vecteur de coordonnées des données en x
        z : vecteur de coordonnées des données en z   
        xi : vecteur de coordonnées où l'interpolation est souhaitée en x
        zi : vecteur de coordonnées où l'interpolation est souhaitée en z
        stock_x : pas de stockage en x entre chaque point interpolé (m)
        stock_z : pas de stockage en z entre chaque point interpolé (m)
        self.pp_interp : matrice de pression interpolée
        '''
        
        x = np.arange(0,self.x_dim,self.delta_x)
        z = np.arange(0,self.z_dim,self.delta_z)
        
        oInterpP = RectBivariateSpline(z, x, self.pp)
                
        xi =  np.arange(0,self.x_dim,stock_x)
        zi =  np.arange(0,self.z_dim,stock_z)
        
        self.pp_interp = oInterpP(zi,xi)
        
  #end def interpolation      
        
# =========================================================================
# Calcul des grandeurs en dB
# =========================================================================
    def to_dB(self):
        '''
        Fonction permettant de passer le champ de pression en dB
        Prends en compte l'absorption atmosphérique self.atmos en tout point du domaine (z,x)

        DL = normalisation // p_source - div_geo
        '''
        
        # Passage en dB du champ de pression
        self.pp_interp = 20. * np.ma.log10(np.abs(self.pp_interp)/2e-5) # passage en dB, .ma = masked array (éviter division /0 dans le log)

        # Calcul de la divergence géométrique en dB        
        div_geo = 20. * np.ma.log10(np.abs(self.Ri))        

        p_max = np.amax(20. * np.ma.log10(np.abs(self.pp)/2e-5))                 # trouve le maximum de la carte de bruit
        
        self.pp_interp = self.pp_interp - p_max             # normalisation par rapport au niveau max
        
        self.pp_interp = self.pp_interp - self.atmos          # prise en compte de l'absorption atmosphérique
        
        # Atténuation // champ libre
        self.DeltaL = self.pp_interp + div_geo # on compense la perte liée à la divergence géométrique pour calcul d'atténuation / champ libre

    # end def to_dB

# =========================================================================
# Lancement du modèle en faisant appel à chaque fonction définie précédemment
# =========================================================================

    def launch(self,x_imp,cgs1,cgs2,lc1,lc2,sigmah1,sigmah2,au,aT,theta,hv,nmode,mu2,bigl,littlel,nb_ind,stock_x,stock_z,hr,turbulence):
        
        'propriétés de sol'
        self.impedance(x_imp,cgs1,cgs2,lc1,lc2,sigmah1,sigmah2) 

        'champ de pression initial en x = 0'
        self.starterWAPE()

        'calcul des profils météorologiques et de célérité effective'
        self.profils_meteo(au,aT,theta,hv)

        'calcul de la zone d amortissement en haut du domaine pour éviter les réflexions parasites'
        self.amortissement()

        'coeur de calcul de la méthode PE : résolution matricielle'
        # avec turbulence (True) ou sans turbulence (False). Sans turbulence bien plus rapide.
        if turbulence:
            self.calc_pade_one_one_turbulence(nb_ind,nmode,mu2,bigl,littlel)
        else:
            self.calc_pade_one_one()

        'calcul de la divergence géométrique Ri'
        self.divergence_geometrique(stock_x,stock_z)

        'calcul de l absorption atmosphérique'
        self.abs_atmos(hr)

        'interpolation linéaires des données (en Pa) tous les "stock_x, stock_z" (et non lambda/n)'
        self.interpolation(stock_x,stock_z)

        'permet de passer le champ de pression pp en (dB)'
        self.to_dB()
        
        return self.Ri
        
    # end def launch
    
# =========================================================================
# Affichage des résultats (cartes de bruit, courbe d'atténuation)
# =========================================================================
    def plot_map(self,stock_x,stock_z):

        distance = np.arange(0,self.x_dim,stock_x) # vecteur distance en m
        idx_mic = int(mt.floor(self.z_mic/stock_z) +1)
        
        #distance = np.arange(0,self.x_dim,self.x_dim/len(self.pp[0])) # vecteur distance en m
        #idx_mic = int(mt.floor(self.z_mic/self.delta_z) +1)

        # affichage de l'atténuation // champ libre
        # plt.xlabel('distance (m)')
        # plt.ylabel('heigh (m)')
        # plt.title(r'DeltaL, $f = $' + str(self.freq) +' (Hz)', fontsize = 16)
        # plt.imshow(self.DeltaL,vmin = -5, vmax = 10, extent=[0,self.x_dim,0,self.z_dim],aspect='auto',cmap='jet',origin='lower')
        # plt.colorbar()
        # #plt.contour(self.p_ij, levels = 4, colors='white', extent=[0,self.x_dim,0,self.z_dim])
        # plt.show()

        #affichage attenuation // champ libre
        # plt.plot(distance,self.DeltaL[idx_mic,:],label = "z_mic = " + str(self.z_mic) + " (m)")
        # plt.ylim(-5,10)
        # plt.xlabel('distance (m)')
        # plt.ylabel('gain (dB)')
        # plt.title(r'Attenuation / free field, $f = $' + str(self.freq) +' (Hz)', fontsize = 16)
        # plt.legend()
        # plt.show()


        # affichage de la carte de bruit
        plt.figure()
        plt.xlabel('distance (m)')
        plt.ylabel('heigh (m)')
        plt.title(r'2D noise map, $f = $' + str(self.freq) +' (Hz)', fontsize = 16)
        plt.imshow(self.pp_interp,vmin = -120, vmax = 0, extent=[0,self.x_dim,0,self.z_dim],aspect='auto',cmap='jet',origin='lower')
        plt.colorbar()
        #plt.contour(self.p_ij, levels = 4, colors='white', extent=[0,self.x_dim,0,self.z_dim])
        plt.show()

        #affichage d'une coupe de pression
        plt.figure()
        plt.plot(distance,self.pp_interp[idx_mic,:],label = "z_mic = " + str(self.z_mic) + " (m)")
        plt.ylim(-120,0)
        plt.xlabel('distance (m)')
        plt.ylabel('gain (dB)')
        plt.title(r'Attenuation / source, $f = $' + str(self.freq) +' (Hz)', fontsize = 16)
        plt.legend()
        plt.show()

     # end def plot_map

#=============================================================================
#Export en .txt de la carte de bruit
#=============================================================================

    def save(self,titre):
           np.savetxt( str(titre), self.pp_interp)
           #np.savetxt(r'WAPE_DeltaL_' + str(self.freq) +'Hz_hs' + str(self.z_source) + 'm', self.DeltaL)

    #end def save

# end class
