* PERDIZ ROJA - VERIFICACION SIN REAJUSTE DE SPSS (5 NOV 2026).
* Fuente: comandos originales recuperados de PERDIZROJA.spv.
* AVISO: Cambiar la ruta siguiente a la del archivo .sav en tu ordenador.
* NO cambia ni guarda el archivo original.
GET FILE='PERDIZROJA.sav'.
DATASET NAME PerdizVerificacion WINDOW=FRONT.

* Constante + 20 coeficientes del modelo AMBIENTAL real, precision completa.
COMPUTE LOGIT_RECON = (-0.11595901957801) +
  (-0.29263518923873 * DP10Spr) +
  (+0.28518108569134 * DP30Spr) +
  (-0.45565477097057 * DP30Sum) +
  (-0.27224706067233 * DTX25Win) +
  (+0.0016555718012446 * ETP_rec) +
  (-0.0035823827047674 * U100_rec) +
  (-0.18864575371958 * CTI) +
  (+0.027608008183601 * Pina_abe_ab) +
  (-0.6863202466027 * Cast_den) +
  (+0.01312801918463 * Enc_alq_ab) +
  (+0.015942739355276 * Enc_alq_den) +
  (+0.057219490192534 * Planif_den) +
  (+0.034166699125999 * Euca) +
  (+0.018311472009083 * Mat_den) +
  (+0.020615955313468 * Herb_altos) +
  (-0.079562690835534 * Cult_inund) +
  (-0.090972741165053 * Urbanizacion) +
  (+0.45773828153739 * Shannon) +
  (-0.017529881398651 * C_forest_cadu) +
  (+0.028569794627774 * C_agric_total).
COMPUTE PRE_RECON = 1/(1+EXP(-LOGIT_RECON)).
COMPUTE FAV_RECON = PRE_RECON/(PRE_RECON+(892/4443)*(1-PRE_RECON)).
COMPUTE DELTA_PRE = ABS(PRE_RECON-PRE_1).
COMPUTE DELTA_FAV = ABS(FAV_RECON-F_AMB).
EXECUTE.

DESCRIPTIVES VARIABLES=DELTA_PRE DELTA_FAV
  /STATISTICS=MEAN STDDEV MIN MAX.

* Matriz de confusion con favorabilidad 0.5.
COMPUTE CLASE_RECON = (FAV_RECON >= .5).
IF (MISSING(FAV_RECON)) CLASE_RECON=$SYSMIS.
EXECUTE.
CROSSTABS /TABLES=PerdizRoja BY CLASE_RECON
 /CELLS=COUNT ROW COLUMN.

* No se ejecuta SAVE para preservar el archivo original.
