# FIGURA 5 CON IC95 - PERDIZ ROJA
# Archivo listo: pulsar Source en RStudio. No modifica archivos originales.
# Reproduce las 1.000 muestras espaciales calculadas para el manuscrito.
# Necesita quantreg, ya utilizado en la primera revisión.

localiza_raiz <- function() {
  frames <- sys.frames()
  paths <- vapply(frames, function(f) {
    z <- f$ofile
    if (is.null(z) || length(z) != 1L) "" else as.character(z)
  }, character(1))
  paths <- paths[nzchar(paths)]
  candidatos <- c(if (length(paths)) dirname(normalizePath(tail(paths, 1), winslash="/", mustWork=FALSE)), getwd())
  ok <- candidatos[file.exists(file.path(dirname(candidatos), "datos", "datos_celdas_686.csv"))]
  if (!length(ok)) stop("No encuentro datos_celdas_686.csv. Extrae TODO el ZIP y abre este script desde la carpeta extraida.")
  dirname(ok[1])
}
RAIZ <- localiza_raiz()
if (!requireNamespace("quantreg", quietly=TRUE)) stop("Falta el paquete quantreg. Ejecuta install.packages('quantreg') y repite.")

D <- read.csv(file.path(RAIZ, "datos", "datos_celdas_686.csv"), stringsAsFactors=FALSE)
F <- read.csv(file.path(RAIZ, "datos", "favorabilidad_espana.csv"), sep=";", fileEncoding="UTF-8-BOM", stringsAsFactors=FALSE)
B <- read.csv(file.path(RAIZ, "datos", "bootstrap_bloques_rep.csv"), check.names=FALSE)
Qs <- c(.25, .50, .70, .75, .90, .95)
F$F_AMB <- suppressWarnings(as.numeric(F$F_AMB))
fav <- F$F_AMB[is.finite(F$F_AMB)]
stopifnot(nrow(D)==686L, length(fav)==5306L, nrow(B)==1000L)

# El identificador de bloque son los 5 primeros caracteres del codigo UTM.
D$bloque <- substr(D$Grid, 1, 5)
filas <- split(seq_len(nrow(D)), D$bloque)
bloques <- sort(names(filas))
stopifnot(identical(bloques, names(B)))

ajustar <- function(datos, tau) {
  fit <- quantreg::rq(Total ~ F_AMB, tau=tau, data=datos, method="br")
  beta <- unname(coef(fit))
  # Mismo calculo algebraico que el estudio original: sin recortar negativos.
  sum(beta[1] + beta[2]*fav)
}

original <- vapply(Qs, function(q) ajustar(D, q), numeric(1))
stopifnot(abs(original[Qs==.50] - 1129259)<5,
          abs(original[Qs==.70] - 2202340)<5)

sims <- matrix(NA_real_, nrow=nrow(B), ncol=length(Qs))
for (i in seq_len(nrow(B))) {
  indices <- unlist(Map(function(ix, k) rep(ix, k), filas, as.integer(unlist(B[i, bloques, drop=FALSE], use.names=FALSE))), use.names=FALSE)
  dboot <- D[indices, c("Total","F_AMB"), drop=FALSE]
  sims[i,] <- vapply(Qs, function(q) ajustar(dboot, q), numeric(1))
  if (i %% 100 == 0) cat("Bootstrap espacial:", i, "de", nrow(B), "replicas\n")
}

ci <- apply(sims, 2, function(x) as.numeric(quantile(x, c(.025,.975), na.rm=FALSE)))
out <- data.frame(cuantil=paste0("Q",round(100*Qs)), tau=Qs,
                  proyeccion=round(original), ic95_inferior=round(ci[1,]),
                  ic95_superior=round(ci[2,]), n_bootstrap=nrow(B), n_bloques=length(bloques))
write.csv(out, file.path(RAIZ,"resultados","FIGURA_5_IC95_RESULTADOS_R.csv"), row.names=FALSE)

png(file.path(RAIZ,"figuras","FIGURA_5_IC95_REPRODUCIDA_R.png"), width=2400, height=1600, res=300)
par(mar=c(5,5,2,2))
x <- Qs*100; y <- original/1e6; low <- ci[1,]/1e6; high <- ci[2,]/1e6
plot(x,y,type="b",pch=19,lwd=2,col="#155c84",ylim=range(c(0,low,high))*c(1,1.05),
     xaxt="n",xlab="Conditional quantile",ylab="Sum of projected cell indices (million individuals)")
axis(1,at=x,labels=paste0("Q",x))
arrows(x,low,x,high,angle=90,code=3,length=.04,col="#155c84")
points(x[2:3],y[2:3],pch=19,col="#b34c36",cex=1.5)
mtext("95% percentile intervals: spatial block bootstrap (1,000 replicates; 74 blocks)", side=3, line=0.25,cex=.8)
dev.off()
cat("\nTERMINADO. Revisa FIGURA_5_IC95_RESULTADOS_R.csv y FIGURA_5_IC95_REPRODUCIDA_R.png\n")
print(out, row.names=FALSE)
cat("NOTA: los intervalos son condicionales al modelo de favorabilidad; NO son IC de la poblacion real.\n")
cat("NOTA: Q25 genera predicciones negativas en algunas celdas; no es una estimacion biologica.\n")
