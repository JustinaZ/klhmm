# Example with Gaussian emissions: Extracting latent variables that satisfy min segmental length constraints


#=======Simulate some data===========================================================
library("KLHMM")

m <- c(-1, 0, 1) # Gaussian means
s <- rep(0.1,length=length(m)) # standard deviation
prm <- list(mean = m,sd = s)
# latent variable set-up:
z <-c(rep(2,length=500), rep(1, length=5),rep(3, length=5),
rep(2,length=500), rep(3, length=10),rep(1, length=10),
rep(2,length=500), rep(3, length=20),rep(1, length=20),
rep(2,length=500), rep(1, length=50),rep(3, length=50),
rep(2,length=500), rep(1, length=100),rep(3, length=100))

nu <- c(0,1,0) # initial probabilities

S <- length(m)
kl <- 0 * diag(S) # matrix with Kullback-Liebler distances
for (i in 1:S){
    for (j in 1:S){
        kl[i,j] <- KL(c(m[i], m[j]), c(s[i], s[j]))
    }
}
N <- length(z)
TT <- list()
L_seg <- c(5, 10, 20, 50, 100)
epsilon <- 0.001
alpha <- 0.03
for (it in 1:length(L_seg)){
    ST <- c()
    L <- L_seg[it]
    repeat{

        x <- c()
        for (n in 1:N){
            x <- c(x, rnorm(1, m[z[n]], s[z[n]]) ) # simulate data
        }
        if (L==5){
            lam <- matrix(runif(9), ncol=3) # starting values for transition matrix
        }else{
            lam <- -A5  # starting value
        }
        x1 <- GD(kl, L, lam, epsilon, alpha, plt=FALSE)
        A0 <- -x1 # keep it in log space for Viterbi
        
        llk <- c()
        for (i in 1:N){
          llk <- rbind(llk, dnorm(x[i], prm$mean, prm$sd, log=TRUE))
        }
        
        z1 <- viterbi(llk, A0, nu)
        ST <- rbind(ST, z1)
        if (L==5){
            A5 <- A0
        }
        if (nrow(ST)==100){
            break
        }
    }
    TT[[it]] <- apply(ST, 2, Mode)
}

# plot results
TT1 <- matrix(unlist(TT), nrow=5, byrow=T)

pdf("GausImg.pdf")
image(t(TT1), bty='n', yaxt="n", xaxt = "n", col=gray.colors(2, start = 0.9, end = 0.3, gamma = 2.2, alpha = NULL), main = "")
dev.off()


pdf("L5_100.pdf")
par(mfrow = c(3,2))

plot(z, type="l", main="Latent Variable (True)",lwd=2, col= "black", cex=1.7, xlab="Index", ylab="States",  cex.lab=1.5, cex.axis=1.3, font.lab=2, bty='n')

plot(TT[[1]], type="l", main="Latent Variable (L=5)",lwd=2, col= "grey50", cex=1.7, xlab="", ylab="",  cex.lab=1.5, cex.axis=1.3, font.lab=2, bty='n')

plot(TT[[2]], type="l", main="Latent Variable (L=10)",lwd=2, col= "grey50", cex=1.7, xlab="", ylab="",  cex.lab=1.5, cex.axis=1.3, font.lab=2, bty='n')

plot(TT[[3]], type="l", main="Latent Variable (L=20)",lwd=2, col= "grey50", cex=1.7, xlab="", ylab="",  cex.lab=1.5, cex.axis=1.3, font.lab=2, bty='n')

plot(TT[[4]], type="l", main="Latent Variable (L=50)",lwd=2, col= "grey50", cex=1.7, xlab="", ylab="",  cex.lab=1.5, cex.axis=1.3, font.lab=2, bty='n')

plot(TT[[5]], type="l", main="Latent Variable (L=100)",lwd=2, col= "grey50", cex=1.7, xlab="", ylab="",  cex.lab=1.5, cex.axis=1.3, font.lab=2, bty='n')
dev.off()

clr <- c("grey55", "grey85", "grey20")
pdf("Gaussian_data.pdf")
plot(x,  type="p", pch=21, cex=1.3, main="", bg=clr[z], col = "black", xlab="Location", ylab="Data",  cex.lab=1.5, cex.axis=1.3, font.lab=2, bty='n')
dev.off()



