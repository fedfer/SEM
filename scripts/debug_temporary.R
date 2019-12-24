n <- 2
k <- 3
l <- 4
m <- 5


eta <- matrix(rnorm(n*k),n,k)     # Initial values
eta.T <- t(eta)

Z <- matrix(rnorm(n * l), n, l)

Deltas <- array(data = rnorm(m * k * l), c(m, k, l))



tmp <- cbind(eta, Z)
inter_chem_cov <- t(apply(tmp, c(1), function(tmp_i){
  apply(Deltas, c(1), etai_Delta_zi, etai = tmp_i[1:ncol(eta)], zi = tmp_i[(ncol(eta) + 1):ncol(tmp)])
}))

