# MODELO DE FAVORABILIDAD DE PERDIZ ROJA: REPRODUCCION DEL ARCHIVO SPSS
# Fecha: 5 noviembre 2026. Solo necesita R base (ningun paquete adicional).
# Ejecutar en RStudio: Source de este archivo. NO sobreescribe datos originales.
# La sintaxis original de SPSS se recupero del registro de PERDIZROJA.spv.
# Esta reconstruccion permite comprobar el modelo, no demuestra que las
# 4.443 celdas codificadas 0 sean ausencias observadas en el campo.

localiza_carpeta <- function() {
  rutas <- lapply(sys.frames(), function(x) x$ofile)
  rutas <- Filter(function(x) !is.null(x) && length(x)==1 && file.exists(x), rutas)
  if(length(rutas)) return(dirname(normalizePath(tail(rutas,1)[[1]])))
  return(getwd())
}
RAIZ <- dirname(localiza_carpeta())
archivo <- file.path(RAIZ,'datos','PERDIZROJA_extraido.csv')
if(!file.exists(archivo)) stop('No encuentro datos/PERDIZROJA_extraido.csv.')
d <- read.csv(archivo, check.names=FALSE, stringsAsFactors=FALSE)

# Las 68 candidatas corresponden a METHOD=FSTEP(COND) del historial de SPSS.
candidatas <- c(
    'DP10Spr', 'DP1Win', 'DP30Spr', 'DP30Sum', 'DTN20', 'DTN20Aut',
  'DTX25Win', 'SISSpr', 'ETP_rec', 'ETR_rec', 'Alt_rec', 'RA_rec',
  'DAut_rec', 'Dens_pob_rec', 'U500_rec', 'U100_rec', 'CTI', 'Ene_sab',
  'Pina_abe_ab', 'Pina_abe_dens', 'Haya_ab', 'Haya_den', 'Cast_ab', 'Cast_den',
  'Chopo_ab', 'Chopo_den', 'Fres_ab', 'Fres_den', 'Enc_alq_ab', 'Enc_alq_den',
  'Planif_ab', 'Planif_den', 'Plani_con_ab', 'Plani_con_den', 'Euca', 'Pal',
  'Lauri_motver', 'Riberas_arb', 'Deforest', 'Mat_ab', 'Mat_den', 'Herb_ralos',
  'Herb_altos', 'Riberas_desarb', 'Masas_agua', 'Cul_herb', 'Cul_reg', 'Cult_inund',
  'Vid', 'Frutales', 'Mosaico_agri', 'Ciudad', 'Urbanizacion', 'Otros_urba',
  'Carreteras', 'Area_degradadas', 'Roquedos', 'Arenales', 'otros_habitat', 'Shannon',
  'C_forest_den', 'C_forest_ab', 'C_forest_cadu', 'C_forest_mixta', 'C_agric_total', 'C_agric_arb',
  'C_amb_acuat', 'C_amb_urbani'
)

# Estas son las 20 variables que contribuyen exactamente a PRE_1 (modelo ambiental).
seleccionadas <- c(
    'DP10Spr', 'DP30Spr', 'DP30Sum', 'DTX25Win', 'ETP_rec', 'U100_rec',
  'CTI', 'Pina_abe_ab', 'Cast_den', 'Enc_alq_ab', 'Enc_alq_den', 'Planif_den',
  'Euca', 'Mat_den', 'Herb_altos', 'Cult_inund', 'Urbanizacion', 'Shannon',
  'C_forest_cadu', 'C_agric_total'
)
stopifnot(all(c('PerdizRoja','PRE_1','F_AMB') %in% names(d)))
stopifnot(all(candidatas %in% names(d)), all(seleccionadas %in% names(d)))
stopifnot(nrow(d)==5335, sum(d$PerdizRoja==1,na.rm=TRUE)==892,
          sum(d$PerdizRoja==0,na.rm=TRUE)==4443)

# SPSS ajusto el FSTEP sobre los 5293 casos sin faltantes en todas las
# 68 candidatas (892 celdas positivas + 4401 codificadas 0).
calibracion <- d[complete.cases(d[,c('PerdizRoja',candidatas)]),]
stopifnot(nrow(calibracion)==5293, sum(calibracion$PerdizRoja==1)==892)

# Reajuste reproducible del MISMO modelo final, sin repetir la seleccion stepwise.
formula_original <- reformulate(seleccionadas, response='PerdizRoja')
modelo <- glm(formula_original, family=binomial(link='logit'),
              data=calibracion, control=glm.control(maxit=150))

# Predecir sobre TODAS las celdas con las 20 variables finales disponibles.
# Asi se obtienen 5306 predicciones, como en el archivo SPSS.
pred <- predict(modelo, newdata=d, type='response')
prevalencia <- 892/4443  # Se usa exactamente la razon presente en SPSS.
favorabilidad <- pred / (pred + prevalencia*(1-pred))

validas <- is.finite(pred) & is.finite(d$PRE_1) & is.finite(d$F_AMB)
stopifnot(sum(validas)==5306)
error_p <- max(abs(pred[validas]-d$PRE_1[validas]))
error_f <- max(abs(favorabilidad[validas]-d$F_AMB[validas]))
cat('\n--- REPRODUCCION DEL MODELO AMBIENTAL ---\n')
cat('Celdas totales:',nrow(d),' | Positivos:',sum(d$PerdizRoja==1),
    ' | Ceros codificados:',sum(d$PerdizRoja==0),'\n')
cat('Muestra ajuste:',nrow(calibracion),' | Predicciones:',sum(validas),'\n')
cat('N variables finales:',length(seleccionadas),'\n')
cat('Diferencia maxima vs PRE_1 archivado:',signif(error_p,5),'\n')
cat('Diferencia maxima vs F_AMB archivado:',signif(error_f,5),'\n')
if(error_p > 1e-7 || error_f > 1e-7) warning('Discrepancia numerica: revisar version de R y entrada.')

# Tabla 1 COMPLETA: incluye las 20 variables y el intercepto.
b <- summary(modelo)$coefficients
tabla1 <- data.frame(Variable=rownames(b), B=b[,1], SE=b[,2],
                     Wald=(b[,1]/b[,2])^2, p=b[,4], row.names=NULL)
write.csv(tabla1, file.path(RAIZ,'resultados','TABLA_1_REGENERADA_R.csv'),row.names=FALSE)

# Tabla 2 para favorabilidad >=0.5, mismo modelo guardado.
obs <- d$PerdizRoja[validas]
est <- favorabilidad[validas]>=0.5
TP <- sum(est & obs==1); FN <- sum(!est & obs==1)
TN <- sum(!est & obs==0); FP <- sum(est & obs==0)
N <- TP+TN+FP+FN
sens <- TP/(TP+FN); esp <- TN/(TN+FP); acc <- (TP+TN)/N
p_obs <- ((TP+FN)/N)*((TP+FP)/N) + ((TN+FP)/N)*((TN+FN)/N)
kappa <- (acc-p_obs)/(1-p_obs)
cat('Matriz: TP=',TP,' FN=',FN,' TN=',TN,' FP=',FP,'\n',sep='')
cat('Sensibilidad:',round(sens,4),' Especificidad:',round(esp,4),
    ' CCR:',round(acc,4),' Kappa:',round(kappa,4),'\n')
write.csv(data.frame(TP,FN,TN,FP,Sensibilidad=sens,Especificidad=esp,
                     CCR=acc,Kappa=kappa),
          file.path(RAIZ,'resultados','TABLA_2_REGENERADA_R.csv'),row.names=FALSE)

# Guardar las predicciones reconstruidas para trazabilidad, sin tocar CSV original.
write.csv(data.frame(UTMCODE=d$UTMCODE, PerdizRoja=d$PerdizRoja,
                     PRE_1_original=d$PRE_1, PRE_1_reconstruido=pred,
                     F_AMB_original=d$F_AMB,F_AMB_reconstruido=favorabilidad),
          file.path(RAIZ,'resultados','PREDICCIONES_RECONSTRUIDAS_R.csv'),row.names=FALSE)
cat('OK: archivos nuevos guardados junto al script.\n')
cat('ATENCION: falta verificar si las celdas 0 tienen evidencia de muestreo.\n')
