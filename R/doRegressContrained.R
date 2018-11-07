#' obtain residulas afer monitonically constrained regression
#'
#' @param single_level_bc text
#' @param fcb_df text
#' @param Loc text
#' @param Weight text
#' @param trans text
#' @param val3 text
#' @param constrained_flag text
#' @param columns text
#' @param monodir text
#' @param cofactor text
#' @return some text
#'
#' @seealso text
#' @examples text

doRegressContrained <- function(single_level_bc, fcb_df = NULL,
                                Loc, weight, trans,
                                fsc_ssc = c(fsc = 'FSC-A', ssc = 'SSC-A'),
                                val3 = NULL, constrained_flag = 1,
                                columns = NULL, monodir = NULL, cofactor = NULL) {
  #print(cofactor)

  print('drc')
  single_level_bc <- as.data.frame(single_level_bc)
  fcb_df <- as.data.frame(fcb_df)

  fsc <- fsc_ssc[1]
  ssc <- fsc_ssc[2]


  lo <- 1 #log offset

  fsc_limits <- c(min(single_level_bc[fsc]), max(single_level_bc[fsc]))
  ssc_limits <- c(min(single_level_bc[ssc]), max(single_level_bc[ssc]))
  data.corr.ls <- vector("list", length = length(columns))

  for (n in seq_along(columns)){
    ind <- append(c(fsc, ssc), columns[n])
    D <- single_level_bc[,ind]
    if (is.null(fcb_df)) {
      D2 <- D
    } else {
      D2 <- fcb_df[,ind]
    }

    rm <- regression_model_matrix()

    switch (trans,
            none = {
              D <- D
              D2 <- D2
            },
            logF = {
              D <- cbind(D[1:2], log(D[3] + lo))
              D2 <- cbind(D2[1:2], log(D2[3] + lo))
            },
            arcsinh = {
              D <- cbind(D[1:2], asinh(D[3]/cofactor))
              D2 <- cbind(D2[1:2], asinh(D2[3]/cofactor))
            }
    )

    #generate regressors
    regressors_1 <- generate_regressors(D, rm)
    regressors_2 <- generate_regressors(D2, rm)

    D <- regressors_1[['D']]
    D2 <- regressors_2[['D']]

    Y <- regressors_1[['Y']]
    Y2 <- regressors_2[['Y']]

    X <- regressors_1[['X']]
    X2 <- regressors_2[['X']]

    #intial unconstrained regression
    B <- lm(Y ~ X)
    OFFSET <- B$coefficients[1]
    Bx <- B$coefficients[2:9]
    mod.resid <- as.numeric(B$residuals)

    if (constrained_flag == 1) {
      con_reg <- constrained_regression(X, Y, fsc_limits, ssc_limits, val3, D,
                                        OFFSET, Bx, rm, monodir)
      Bc <- con_reg[['Bc']]
      NOFFSET <- con_reg[['NOFFSET']]
      NB <- con_reg[['NB']]
      Nresidual <- Y - (X%*%Bc[1:8] + NOFFSET)
      Nresidual2 <- Y2 - (X2%*%Bc[1:8] + NOFFSET)
    }
    else {
      NB <- Bx
      NOFFSET <- OFFSET
      Nresidual <- mod.resid
    }

    XO <- matrix(nrow= length(weight), ncol = 4)
    XO[,1] <- Loc[,1]
    XO[,2] <- Loc[,2]
    XO[,3] <- Loc[,1]*Loc[,2]
    XO[,4] <- NOFFSET

    for (i in 1:ncol(rm)){
      XO[,4] <- XO[,4] + NB[i]*XO[,rm[1,i]]^rm[2,i]
    }

    R2 <- as.numeric(Nresidual2) + as.vector(t(weight)%*%XO[,4]/sum(weight))

    if (trans == 'logF') {
      data.corr <- exp(R2) - lo
    } else if (trans == 'none') {
      data.corr <- R2
    } else if (trans == 'arcsinh') {
      data.corr <- sinh(R2)*cofactor
    }


    print(paste(columns[n], "complete"))
    data.corr.ls[[n]] <- data.corr
  }
  data.corr.df <- as.data.frame(data.corr.ls)
  colnames(data.corr.df) <- columns

  return(list(df = data.corr.df, coefs = Bc))
}
