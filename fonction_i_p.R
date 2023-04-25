i_p <- function(p){
  exp(- qnorm(1-p)^2 / 2) / (p * sqrt(2*pi) )
}

i_p(1)
i_p(0.5)
i_p(0.07)
i_p(0.26)
i_p(0.3)
i_p(0.6)




i <- sapply(seq(0.01,1,0.01) , FUN = i_p)

plot(x = seq(0.01,1,0.01) , y = i)


