 #! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#

'''
author: Bill Kayser
contributors: Benoit Gauvreau, Pierre Chobeau, Bertrand Lihoreau
date: 02/2021
description: main script of Parabolic Equation method that uses functions defined in PE.py
                      this code is adapted from a Matlab model previously developped at UMRAE, Univ.G.Eiffel
'''

import time
tps1 = time.time()

'============================================================================='
''' Importation du fichier.py contenant la class de fonctions du modèle PE'''
'============================================================================='
import class_PE

'============================================================================='
''' Définition des paramètres d'entrée du modèle'''
'============================================================================='
#--- La source ---#
freq = 500                # fréquence de la source (Hz)
z_source = 10          # hauteur de la source (m)
z_mic = 1.5      # hauteur du microphone (m)

#--- Le domaine ---#
x_dim = 300            # longueur du domaine suivant x (m)
z_dim = 50             # hauteur du domaine suivant z (m)
discrx = 10             # discretisation spatiale de calcul suivant x = lambda/discrx
discrz = 10             # discretisation spatiale de calcul suivant z = lambda/discrz
stock_x = 0.5          # discrétisation spatiale de stockage suivant x (m)
stock_z = 0.5          # discrétisation spatiale de stockage suivant z (m)

#--- Le milieu ---#
T = 10                       # température atmosphérique de surface (°C)
au = 0                      # 0 homogène, coefficient du profil log de vent (m/s) 0<au<1.7
aT = 0                     # 0 homogène, coefficient du profil log de température (K/m) -0.5<aT<0.25
theta = 0               # angle du vent (degré) par rapport à la propagation (0 : dans l'axe, 180 : opposé)
hr = 0.8                  # hygrométrie de l'air (%)

#--- Le sol ---#
hv = 0                       # hauteur de végétation (m)
lc1 = 0                     # longueur de corrélation de la rugosité (m), si 0 > pas de rugosité, 0.05<lc<1.5
sigmah1 = 0                # écart-type des hauteur de rugosité (m), 0.01<sigma<0.05
cgs1 = 100000             # resistivité du sol (kNsm-4) traduit l'effet d'absorption (cgs > 100000 => beta = 1)
x_imp = x_dim           # position de la discontinuité d'impédance
lc2 = lc1                  # propriété du 2e type de sol                 
sigmah2 = sigmah1       # propriété du 2e type de sol           
cgs2 = cgs1            # propriété du 2e type de sol

#--- La turbulence ---#
turbulence = False  # Booléens pour lancer la turbulence (True) ou non (False). Beaucoup plus rapide sans turbulence (False)
nb_ind = 2             # nb d'itérations (calculs) de la turbulence
nmode = 100             # nombre de modes de Fourier pour le calcul de la turbulence
mu2 = 8e-06           # mu2 classiquement dans la biblio : 2e-006 < vart < 8e-006 (contestable ?)
bigl = 1.1              # échelle externe de la turbulence (m) spectre von Karman
littlel = 0.001     # échelle interne de la turbulence (m) spectre von Karman

#--- Sauvegarde ---#
titre = (r'turb_' + str(freq) +'Hz_hs' + str(z_source) + 'm')

'============================================================================='
'''Instantiation / Initialisation de l'objet Calcul'''
'============================================================================='
oCalc = class_PE.fonctions(freq,z_source,z_mic,x_dim,z_dim,discrx,discrz,T)

'============================================================================='
'''Lancement du modèle : différentes fonctions à la suite (cf class PE.py)'''
'============================================================================='
oCalc.launch(x_imp,cgs1,cgs2,lc1,lc2,sigmah1,sigmah2,au,aT,theta,hv,nmode,mu2,bigl,littlel,nb_ind,stock_x,stock_z,hr,turbulence)

tps2 = time.time()
print(tps2 - tps1)
'============================================================================='
''' Affichage de la carte de bruit & de la courbe d'atténuation'''
'============================================================================='
oCalc.plot_map(stock_x,stock_z)

'============================================================================='
'''Export en .txt de la carte de bruit'''
'============================================================================='
#oCalc.save(titre)
