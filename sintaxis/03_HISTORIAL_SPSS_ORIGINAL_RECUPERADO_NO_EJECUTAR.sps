* TRANSCRIPCION DEL HISTORIAL RECUPERADO DE PERDIZROJA.spv.
* CONTENIDO DOCUMENTAL; NO EJECUTAR SIN ADAPTAR RUTAS.
* LA SELECCION ORIGINAL USO FSTEP(COND), 68 VARIABLES CANDIDATAS.
* LA F_AMB se calculo con PRE_1 y la razon 892/4443.
* Los comandos de seleccion y los de exportacion fueron copiados del SPV.

* ===== HISTORIAL SPV - BLOQUE 2 =====.
LOGISTIC REGRESSION VARIABLES PerdizRoja
  /METHOD=BSTEP(COND) La Lo La2 Lo2 LaLo
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20) CUT(.5).

* ===== HISTORIAL SPV - BLOQUE 4 =====.
COMPUTE Y_ESPACIAL=-22.934+1.325 * La-0.436 * Lo-0.020 * La2-0.010 * Lo2+0.010 * LaLo.
EXECUTE.
LOGISTIC REGRESSION VARIABLES PerdizRoja
  /METHOD=FSTEP(COND) DP10Spr DP1Win DP30Spr DP30Sum DTN20 DTN20Aut DTX25Win SISSpr ETP_rec ETR_rec
    Alt_rec RA_rec DAut_rec Dens_pob_rec U500_rec U100_rec CTI Ene_sab Pina_abe_ab Pina_abe_dens
    Haya_ab Haya_den Cast_ab Cast_den Chopo_ab Chopo_den Fres_ab Fres_den Enc_alq_ab Enc_alq_den
    Planif_ab Planif_den Plani_con_ab Plani_con_den Euca Pal Lauri_motver Riberas_arb Deforest Mat_ab
    Mat_den Herb_ralos Herb_altos Riberas_desarb Masas_agua Cul_herb Cul_reg Cult_inund Vid Frutales
    Mosaico_agri Ciudad Urbanizacion Otros_urba Carreteras Area_degradadas Roquedos Arenales
    otros_habitat Shannon C_forest_den C_forest_ab C_forest_cadu C_forest_mixta C_agric_total
    C_agric_arb C_amb_acuat C_amb_urbani
  /SAVE=PRED
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20) CUT(.5).

* ===== HISTORIAL SPV - BLOQUE 6 =====.
COMPUTE F_AMB=(((PRE_1)/(1-PRE_1))/((892/4443)+(PRE_1/(1-PRE_1)))).
EXECUTE.
LOGISTIC REGRESSION VARIABLES PerdizRoja
  /METHOD=FSTEP(COND) DP10Spr DP1Win DP30Spr DP30Sum DTN20 DTN20Aut DTX25Win SISSpr ETP_rec ETR_rec
    Alt_rec RA_rec DAut_rec Dens_pob_rec U500_rec U100_rec CTI Ene_sab Pina_abe_ab Pina_abe_dens
    Haya_ab Haya_den Cast_ab Cast_den Chopo_ab Chopo_den Fres_ab Fres_den Enc_alq_ab Enc_alq_den
    Planif_ab Planif_den Plani_con_ab Plani_con_den Euca Pal Lauri_motver Riberas_arb Deforest Mat_ab
    Mat_den Herb_ralos Herb_altos Riberas_desarb Masas_agua Cul_herb Cul_reg Cult_inund Vid Frutales
    Mosaico_agri Ciudad Urbanizacion Otros_urba Carreteras Area_degradadas Roquedos Arenales
    otros_habitat Shannon C_forest_den C_forest_ab C_forest_cadu C_forest_mixta C_agric_total
    C_agric_arb C_amb_acuat C_amb_urbani Y_ESPACIAL
  /SAVE=PRED
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20) CUT(.5).

* ===== HISTORIAL SPV - BLOQUE 8 =====.
COMPUTE F_AMBESP=(((PRE_2)/(1-PRE_2))/((892/4443)+(PRE_2/(1-PRE_2)))).
EXECUTE.

SAVE OUTFILE='C:\Users\Lenovo\Desktop\ARTEMISAN\modelofavorabilidad\PERDIZROJA.sav'
  /COMPRESSED.

SAVE TRANSLATE OUTFILE='C:\Users\Lenovo\Desktop\ARTEMISAN\modelofavorabilidad\PERDIZROJA.csv'
  /TYPE=CSV
  /ENCODING='UTF8'
  /MAP
  /REPLACE
  /FIELDNAMES
  /CELLS=VALUES.



* ===== HISTORIAL SPV - BLOQUE 9 =====.
GET
  FILE='C:\Users\Lenovo\Desktop\ARTEMISAN\FAVORABILIDAD\Perdiz Roja\PERDIZROJA.sav'.
DATASET NAME ConjuntoDatos1 WINDOW=FRONT.
ROC F_AMB BY PerdizRoja (1)
  /PLOT=CURVE
  /CRITERIA=CUTOFF(INCLUDE) TESTPOS(LARGE) DISTRIBUTION(FREE) CI(95)
  /MISSING=EXCLUDE.