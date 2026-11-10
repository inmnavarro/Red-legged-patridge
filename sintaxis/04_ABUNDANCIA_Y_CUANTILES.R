# REVISIÓN RÁPIDA — PERDIZ ROJA, NOVIEMBRE 2026
# Ejecutar de principio a fin en RStudio (Source) o con Rscript.
# NO edita los datos originales. Sólo dependencia: quantreg.
# Los datos necesarios están incluidos en la carpeta 'datos'.
# ADVERTENCIA: índices sin corrección de detectabilidad; no representan IC.

localiza_raiz <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep('^--file=', args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub('^--file=', '', file_arg[1]), winslash='/', mustWork=TRUE)))
  if (requireNamespace('rstudioapi', quietly=TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error=function(e) '')
    if (nzchar(p) && file.exists(p)) return(dirname(normalizePath(p, winslash='/', mustWork=TRUE)))
  }
  getwd()
}
RAIZ <- dirname(localiza_raiz())
if (!dir.exists(file.path(RAIZ,'datos'))) stop('No encuentro la carpeta datos junto al script. Abre 01_REVISION_PERDIZ_EJECUTAR.R en RStudio y usa Source.')
if (!requireNamespace('quantreg', quietly = TRUE)) stop("Instala 'quantreg' una sola vez: install.packages('quantreg'); luego ejecuta otra vez.")
OUT <- file.path(RAIZ, 'resultados')
FIG <- file.path(RAIZ, 'figuras')
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
dir.create(FIG, showWarnings=FALSE, recursive=TRUE)
lee <- function(f) read.csv(file.path(RAIZ,'datos',f), sep=';', fileEncoding='UTF-8-BOM', stringsAsFactors=FALSE, check.names=FALSE, na.strings=c('', 'NA','NaN'))
norm <- function(x) gsub('[[:space:]]+','_',tolower(trimws(as.character(x))))

a <- lee('transectos_abundancia_original.csv')
t <- lee('transectos_tendencias_incluye_ceros.csv')
fav <- lee('favorabilidad_espana.csv')
orig <- lee('celdas_modelo_original.csv')

# 1. Reconstruir el criterio aplicado a las observaciones positivas.
stopifnot(nrow(a)==2050, nrow(t)==2634)
rep_a <- norm(a$Max_Habitat) == norm(a$X_majority)
rep_t <- norm(t$Max_Habitat_t) == norm(t$X_majority)
rep_a[is.na(rep_a)] <- FALSE
rep_t[is.na(rep_t)] <- FALSE
stopifnot(all(rep_a == (a$CHECK == 'SI')))
positivos <- a[rep_a, ,drop=FALSE]
ceros <- t[!is.na(t$perdiz) & t$perdiz==0 & rep_t, ,drop=FALSE]
stopifnot(nrow(positivos)==1387, nrow(ceros)==399)
stopifnot(sum(rep_t & !is.na(t$perdiz) & t$perdiz>0)==1387)

# 2. Agregar la densidad observada por celda UTM.
positivos$Densidad <- as.numeric(positivos$Densidad)
stopifnot(all(is.finite(positivos$Densidad)))
dens <- aggregate(Densidad ~ Grid, data=positivos, FUN=mean)
np <- aggregate(Densidad ~ Grid, data=positivos, FUN=length)
names(np)[2] <- 'n_positivos'
nz <- aggregate(rep(1L,nrow(ceros)), by=list(Grid=ceros$Grid), FUN=sum)
names(nz)[2] <- 'n_ceros'
cel <- merge(merge(dens,np,by='Grid',all=TRUE),nz,by='Grid',all=TRUE)
cel$n_positivos[is.na(cel$n_positivos)] <- 0L
cel$n_ceros[is.na(cel$n_ceros)] <- 0L
names(cel)[names(cel)=='Densidad'] <- 'densidad_original'
cel$densidad_con_ceros <- ifelse(cel$n_positivos+cel$n_ceros>0,
  ifelse(is.na(cel$densidad_original),0,cel$densidad_original*cel$n_positivos) / (cel$n_positivos+cel$n_ceros), NA_real_)
cel$Total_original <- 100*cel$densidad_original
cel$Total_con_ceros <- 100*cel$densidad_con_ceros
cel <- merge(cel, fav[c('Grid','F_AMB')], by='Grid', all.x=TRUE)
stopifnot(nrow(cel)==694)
comp <- merge(orig[c('Grid','Densidad')], cel[c('Grid','densidad_original')],by='Grid')
stopifnot(max(abs(comp$Densidad-comp$densidad_original),na.rm=TRUE)<1e-5)
write.csv(cel,file.path(OUT,'01_auditoria_por_celda.csv'),row.names=FALSE,fileEncoding='UTF-8')

# 3. Ajustar los cuantiles: se conserva el umbral del script ORIGINAL,
#    Total < 4600 individuos por celda, para una comparación estricta.
qs <- c(.25,.50,.70,.75,.90,.95)
proy <- fav$F_AMB[is.finite(fav$F_AMB)]
esc <- list(
  original = cel[cel$n_positivos>0,c('Grid','F_AMB','Total_original')],
  ceros_celdas_existentes = cel[cel$n_positivos>0,c('Grid','F_AMB','Total_con_ceros')],
  ceros_y_celdas_nuevas = cel[,c('Grid','F_AMB','Total_con_ceros')]
)
mods <- list(); resultados <- list(); usados <- list()
for (nombre in names(esc)) {
  d <- esc[[nombre]]
  names(d)[3] <- 'Total'
  d <- d[is.finite(d$Total) & is.finite(d$F_AMB) & d$Total<4600,,drop=FALSE]
  usados[[nombre]] <- d
  for (q in qs) {
    fit <- quantreg::rq(Total ~ F_AMB,tau=q,data=d,method='br')
    cf <- as.numeric(coef(fit))
    pred <- cf[1] + cf[2]*proy
    clave <- sprintf('%s_%.2f',nombre,q)
    mods[[clave]] <- fit
    resultados[[length(resultados)+1L]] <- data.frame(
      escenario=nombre, cuantil=q, n_celdas=nrow(d), intercepto=cf[1],
      pendiente=cf[2], total_nacional=sum(pred), n_pred_negativas=sum(pred<0),
      n_celdas_proyectadas=length(pred))
  }
}
res <- do.call(rbind,resultados)
write.csv(res,file.path(OUT,'02_cuantiles_sensibilidad.csv'),row.names=FALSE,fileEncoding='UTF-8')
write.csv(usados[['ceros_y_celdas_nuevas']], file.path(OUT,'03_celdas_utilizadas_analisis_revisado.csv'),row.names=FALSE)
rr <- subset(res, cuantil %in% c(.5,.7))
write.csv(rr,file.path(OUT,'04_resultados_principales.csv'),row.names=FALSE)

# Control de reproducibilidad: las proyecciones originales publicadas.
a50 <- subset(res, escenario=='original' & cuantil==.5)$total_nacional
a70 <- subset(res, escenario=='original' & cuantil==.7)$total_nacional
if (abs(a50-1775949)>100 || abs(a70-2830213)>100) stop('No se reproducen los totales originales. Revisar versión R/quantreg y datos de entrada.')

# 4. Figura 3 revisada. Eje Y: individuos por 100 km2 (NO individuos/km2).
d <- usados[['ceros_y_celdas_nuevas']]
f50 <- mods[['ceros_y_celdas_nuevas_0.50']]
f70 <- mods[['ceros_y_celdas_nuevas_0.70']]
png(file.path(FIG,'Figura_3_cuantiles_revisada.png'),width=1800,height=1200,res=200)
par(mar=c(5,5,2,1))
plot(d$F_AMB,d$Total,pch=16,cex=.5,col=adjustcolor('gray20',alpha.f=.35),
 xlab='Environmental favourability',ylab='Observed index (individuals per 100 km2 grid cell)')
x <- seq(0,1,length.out=101)
lines(x,coef(f50)[1]+coef(f50)[2]*x,lwd=3,col='#185A8B')
lines(x,coef(f70)[1]+coef(f70)[2]*x,lwd=3,col='#AF5636',lty=2)
legend('topleft',c('Q50','Q70'),lty=c(1,2),col=c('#185A8B','#AF5636'),lwd=3,bty='n')
dev.off()

# 5. Figura 5 revisada: sensibilidad, NO intervalo de confianza.
p <- subset(res,escenario=='ceros_y_celdas_nuevas')
p <- p[order(p$cuantil),]
png(file.path(FIG,'Figura_5_sensibilidad_revisada.png'),width=1800,height=1100,res=200)
par(mar=c(5,5,2,1))
plot(p$cuantil,p$total_nacional/1e6,type='b',pch=19,col='#185A8B',lwd=2,
 xlab='Quantile (tau)',ylab='Sum of projected cell values (millions of individuals)',
 xlim=c(.2,1),ylim=range(c(0,p$total_nacional/1e6)))
points(p$cuantil[p$cuantil %in% c(.5,.7)],p$total_nacional[p$cuantil %in% c(.5,.7)]/1e6,
  pch=19,cex=1.6,col='#AF5636')
dev.off()

# 6. Validación espacial exploratoria con 5 particiones por bloques de 100 km;
#    comparar la pérdida pinball frente a predecir un cuantil constante.
pinball <- function(y,yhat,tau) mean(ifelse(y-yhat>=0,tau*(y-yhat),(tau-1)*(y-yhat)))
cv <- list()
set.seed(20261102)
for (nm in c('original','ceros_y_celdas_nuevas')) {
 d <- usados[[nm]]
 bloques <- unique(substr(d$Grid,1,5))
 fold_bloques <- setNames(sample(rep(1:5,length.out=length(bloques))),bloques)
 fold <- fold_bloques[substr(d$Grid,1,5)]
 for (q in c(.5,.7)) {
  perdidas_modelo <- c(); perdidas_constante <- c(); pesos <- c()
  for (k in 1:5) {
    train <- d[fold!=k,,drop=FALSE]; test <- d[fold==k,,drop=FALSE]
    if (!nrow(test) || nrow(train)<10) next
    mod <- quantreg::rq(Total~F_AMB,tau=q,data=train,method='br')
    pm <- as.numeric(predict(mod,newdata=test))
    pc <- rep(as.numeric(quantile(train$Total,probs=q,type=1)),nrow(test))
    perdidas_modelo <- c(perdidas_modelo,pinball(test$Total,pm,q))
    perdidas_constante <- c(perdidas_constante,pinball(test$Total,pc,q))
    pesos <- c(pesos,nrow(test))
  }
  Lm <- weighted.mean(perdidas_modelo,pesos)
  Lc <- weighted.mean(perdidas_constante,pesos)
  cv[[length(cv)+1L]] <- data.frame(escenario=nm,cuantil=q,n=nrow(d),
    perdida_modelo=Lm,perdida_constante=Lc,mejora_porcentaje=100*(Lc-Lm)/Lc)
 }
}
write.csv(do.call(rbind,cv),file.path(OUT,'05_validacion_espacial_exploratoria.csv'),row.names=FALSE)

cat('\n----------------------- RESUMEN -----------------------\n')
cat('Transectos representativos originales:',nrow(positivos),'\n')
cat('Ceros representativos recuperados:',nrow(ceros),'\n')
cat('Celdas con ceros y/o positivos:',nrow(cel),'\n')
print(rr,row.names=FALSE)
cat('\nResultados y figuras escritos en:',RAIZ,'\n')
cat('NOTA: el ajuste de favorabilidad y su matriz de confusión NO se recalculan aquí.\n')
cat('NOTA: la extrapolación NO corrige detectabilidad, sesgo de voluntariado ni sueltas.\n')
cat('NOTA: Q50-Q70 son ESCENARIOS, no un intervalo estadístico para la abundancia real.\n')
