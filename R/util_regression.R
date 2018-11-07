
#Performing the Constrained regression
#----------------------------------------

#vars (1= FSC, 2= SSC, 3= FSC*SSC)
#powers (indicates power to which each predictor is raised)
regression_model_matrix <- function(nrow = 2, ncol = 8,
                                    vars = c(1, 1, 1, 2, 2, 2, 3, 3),
                                    powers = c(1, 0.5, 2, 1, 0.5, 2, 1, 0.5)){

    rm <- matrix(nrow = nrow, ncol = ncol)
    rm[1,] <- vars
    rm[2,] <- powers

    return(rm)
}

generate_regressors <- function(D, rm){

    Y <- D[,3]
    D[,3] <- D[,1]*D[,2]
    X <- matrix(0, nrow = length(Y), ncol = ncol(rm))

    for (m in 1:ncol(rm)){
        X[,m] = D[,rm[1,m]]^rm[2,m]
    }

    return(list('Y' = Y, 'D' = D, 'X' = X))
}

constrained_regression <- function(X, Y, fsc_limits, ssc_limits, val3, D,
                                   OFFSET, Bx, rm, monodir){
    S = 5
    area = 0.5
    xlimr <- c(floor(min(fsc_limits)), ceiling(max(fsc_limits)))
    ylimr <- c(floor(min(ssc_limits)), ceiling(max(ssc_limits)))
    P <- MASS::kde2d(D[,1], D[,2], n=2^S, lims= c(xlimr, ylimr))
    XX <- matrix(rep(P[[1]],2^S), ncol = 2^S, nrow = 2^S, byrow = TRUE)
    YY <- matrix(rep(P[[2]],2^S), ncol = 2^S, nrow = 2^S, byrow = FALSE)
    PP <- P[[3]]
    PP <- PP/sum(PP)

    ft <- function(t) {abs(sum(PP[PP>t])- area)}
    tmin <- optimize(ft, c(0,max(PP)))$minimum

    ZZ= matrix(nrow=sum(PP>tmin), ncol=4)
    ZZ[,1] <- as.numeric(XX[PP>tmin])
    ZZ[,2] <- as.numeric(YY[PP>tmin])
    ZZ[,3] <- ZZ[,1]*ZZ[,2]
    ZZ[,4] <- OFFSET

    for (i in 1:ncol(rm)){
        ZZ[,4] <- ZZ[,4] + Bx[i]*ZZ[,rm[1,i]]^rm[2,i]
    }

    dir1 <- lm(ZZ[,4] ~ ZZ[,1])
    dir2 <- lm(ZZ[,4] ~ ZZ[,2])

    if(is.null(monodir)){
        monodir <- +c(dir1$coefficients[2]>=0, dir2$coefficients[2]>=0)
    } else {
        monodir <- monodir
    }

    #2. Select points to evaluate on monotonicity
    S= 3
    P <- MASS::kde2d(D[,1], D[,2], n=2^S, lims= c(xlimr, ylimr))
    XX <- matrix(rep(P[[1]],2^S), ncol = 2^S, nrow = 2^S, byrow = TRUE)
    YY <- matrix(rep(P[[2]],2^S), ncol = 2^S, nrow = 2^S, byrow = FALSE)
    PP <- P[[3]]

    noc <- 2*(2^S-1)*2^S
    Q <- matrix(nrow=noc, ncol = 4)
    p <- 0
    for (i in 1:2^S) {
        for (j in 1:(2^S-1)) {
            p <- p + 1
            if (monodir[1] == 1) {
                Q[p,] <- c(XX[1,j], YY[i,1], XX[1,j+1], YY[i,1])
            } else {
                Q[p,] <- c(XX[1, j+1], YY[i, 1], XX[1, j], YY[i,1])
            }
        }
    }

    for (i in 1:2^S) {
        for (j in 1:(2^S-1)) {
            p <- p + 1
            if (monodir[2] == 1) {
                Q[p,] <- c(XX[1,i], YY[j,1], XX[1,i], YY[j+1,1])
            } else {
                Q[p,] <- c(XX[1, i], YY[j+1, 1], XX[1,i], YY[j,1])
            }
        }
    }

    XX1 <- matrix(nrow = noc, ncol = 3)
    XX1[,1] <- as.numeric(Q[,1])
    XX1[,2] <- as.numeric(Q[,2])
    XX1[,3] <- XX1[,1]*XX1[,2]
    A1 <- matrix(nrow=noc, ncol = ncol(rm))
    for (m in 1:ncol(rm)) {
        A1[,m] <- XX1[,rm[1,m]]^rm[2,m]
    }
    XX2 <- matrix(nrow = noc, ncol = 3)
    XX2[,1] <- as.numeric(Q[,3])
    XX2[,2] <- as.numeric(Q[,4])
    XX2[,3] <- XX2[,1]*XX2[,2]
    A2 <- matrix(nrow=noc, ncol = ncol(rm))
    for (m in 1:ncol(rm)) {
        A2[,m] <- XX2[,rm[1,m]]^rm[2,m]
    }


    if (!is.null(val3)){
        A <- matrix(nrow = noc + 2, ncol = ncol(rm))
        A[1:noc,] <- A1-A2

        XBP <- matrix(nrow = 2, ncol=3)
        if (monodir[1]==1 & monodir[2]==1){
            XBP[,1]<- XX[1, c(1, 2^S)]
            XBP[,2]<- YY[c(1, 2^S), 1]
            "FSC up, SSC up"
        } else if (monodir[1]==0 & monodir[2]==1){
            XBP[,1]<- XX[1, c(2^S,1)]
            XBP[,2]<- YY[c(1, 2^S), 1]
            "FSC down, SSC up"
        } else if (monodir[1]==1 & monodir[2]==0){
            XBP[,1]<- XX[1, c(1, 2^S)]
            XBP[,2]<- YY[c(2^S, 1), 1]
            "FSC up, SSC down"
        } else if (monodir[1]==0 & monodir[2]==0){
            XBP[,1]<- XX[1, c(2^S, 1)]
            XBP[,2]<- YY[c(2^S, 1), 1]
            "FSC down, SSC down"
        }

        XBP[,3] <- XBP[,1]*XBP[,2]
        AX <- matrix(nrow = 2, ncol= ncol(rm))
        for (m in 1:ncol(rm)){
            AX[,m] <- XBP[,rm[1,m]]^rm[2,m]
        }

        A[noc+1, ] <- -AX[1,]
        A[noc+2, ] <- AX[2,]

        b <- matrix(0, nrow=noc+2, 1)
        b(noc+1) <- val3[1]
        b(noc+2) <- val3[2]

    } else {
        A <- A1-A2
        b <- matrix(0, nrow=noc, ncol=1)
    }

    # #constrained regression analysis
    # if(FALSE) {
    #   ## This seems silly, but prcma is dependent on quadprog,
    #   ## but shiny wasn't recognizing this dependency and won't install
    #   ## the package on the shinyapps.io server otherwise.
    #   ## This block of code shouldn't run, but in case it does,
    #   ## this is just the example code from ?quadprog::solve.QP
    #
    # }





    NB <- numeric(length = length(Bx))
    Bc <- pracma::lsqlincon(cbind(X, 1),Y,cbind(A,0),b)
    NOFFSET <- Bc[9]
    NB <- Bc[1:8]


    return(list('NB' = NB, 'Bc' = Bc, 'NOFFSET' = NOFFSET))
}

