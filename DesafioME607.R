##Função AR_MM 
AR_MM <- function(serie_simulada, n, p){
  gamma0 <- var(serie_simulada)
  dimensao <- p
  soma_ro_phi <- c()
  intervalos_confianca <- c()
  limite_superior <- c()
  limite_inferior <- c()
  M1 <- matrix(NA, nrow = dimensao, ncol = dimensao)
  acf_vals <- acf(serie_simulada, lag.max = dimensao, plot = FALSE)$acf
  M2 <- matrix(acf_vals[-1], ncol = 1)
  
  #Matriz M1
  for(i in 1:dimensao){
    for(j in 1:dimensao){
      M1[i,j] <- acf_vals[abs(i-j)+1]
    }
  }
  #Encontrar phi_hat
  phi_hat <- solve(M1) %*% M2 
  
  #Encontrar o valor de sigma2_e_hat (usando rho e phi_hat)
  for (i in 1: nrow(phi_hat)){
    soma_ro_phi <- append(soma_ro_phi, phi_hat[i] * acf_vals[i+1])
  }
  sigma2_e_hat <- gamma0 * (1 - sum(soma_ro_phi))
  
  #Variância estimada de phi estimado
  m_Gamma <- gamma0*M1
  ACOV <- n^{-1} * sigma2_e_hat * solve(m_Gamma)
  
  #Intervalo de confiança 
  for (i in 1:nrow(ACOV)){
    limite_inferior <- append(limite_inferior, phi_hat[i] - 1.96*sqrt(ACOV[i,i]))
    limite_superior <- append(limite_superior, phi_hat[i] + 1.96*sqrt(ACOV[i,i]))
    ic <- c(phi_hat[i] - 1.96*sqrt(ACOV[i,i]), phi_hat[i] + 1.96*sqrt(ACOV[i,i]))
    intervalos_confianca <- append(intervalos_confianca, paste0("phi_",i, ":",
                                                                "[", round(ic[1],4), ", ", round(ic[2],4), "]"))
  }
  
  #Resultados
  return(list(phi_hat = phi_hat,
              limite_inferior = limite_inferior,
              limite_superior = limite_superior,
              intervalos_confianca = intervalos_confianca))
}


# Exemplo de como utilizar
serie_ar3 <- arima.sim(n = 500, list(ar = c(0.5, -0.3, 0.2)), sd = 1)
modelo <- AR_MM(serie_simulada = serie_ar3, n = 500, p =3)
modelo

# Resultados da simulação
set.seed(607)
n <- 500
phi1_sup <- c()
phi2_sup <- c()
phi3_sup <- c()
phi1_inf <- c()
phi2_inf <- c()
phi3_inf <- c()
phi1_estimado <- c()
phi2_estimado <- c()
phi3_estimado <- c()
for (i in 1:1000){
  serie_ar2 <- arima.sim(n = n, list(ar = c(0.5, -0.3, 0.2)), sd = 1)
  modelo <- AR_MM(serie_ar2, n, 3)
  phi1_inf <- append(phi1_inf, modelo$limite_inferior[1])
  phi2_inf <- append(phi2_inf, modelo$limite_inferior[2])
  phi3_inf <- append(phi3_inf, modelo$limite_inferior[3])
  phi1_sup <- append(phi1_sup, modelo$limite_superior[1])
  phi2_sup <- append(phi2_sup, modelo$limite_superior[2])
  phi3_sup <- append(phi3_sup, modelo$limite_superior[3])
  phi1_estimado <- append(phi1_estimado , modelo$phi_hat[1])
  phi2_estimado <- append(phi2_estimado , modelo$phi_hat[2])
  phi3_estimado <- append(phi3_estimado , modelo$phi_hat[3])
}
phi_1 <- data.frame(phi1_inf, phi1_sup)
phi_2 <- data.frame(phi2_inf, phi2_sup)
phi_3 <- data.frame(phi3_inf, phi3_sup)
phi_1$avaliar <- ifelse(phi_1$phi1_inf < 0.5 & phi_1$phi1_sup > 0.5, NA, "*")
phi_2$avaliar <- ifelse(phi_2$phi2_inf < -0.3 & phi_2$phi2_sup > -0.3, NA, "*")
phi_3$avaliar <- ifelse(phi_3$phi3_inf < 0.2 & phi_3$phi3_sup > 0.2, NA, "*")

#Vies
vies_phi1 <- mean(phi1_estimado) - 0.5
vies_phi2 <- mean(phi2_estimado) - (-0.3)
vies_phi3 <- mean(phi3_estimado) - 0.2

#RMSE
rmse_phi1 <- sqrt(mean((phi1_estimado - 0.5)^2))
rmse_phi2 <- sqrt(mean((phi2_estimado - (-0.3))^2))
rmse_phi3 <- sqrt(mean((phi3_estimado - 0.2)^2))

#Quantidade de elementos fora do IC
sum(!is.na(phi_1$avaliar))
sum(!is.na(phi_2$avaliar))
sum(!is.na(phi_3$avaliar))
