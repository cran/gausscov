#' Generates sin(pi*j*(1:n)/n) and cos(pi*j*(1:n)/n) for j=1,...,m for a given sample size n
#'
#' @param n Sample size.
#' @param m Maximum order of sine and cosine functions. 
#' @return x  The functions sin(pi*j*(1:n)/n) (odd) and cos(pi*j*(1:n)/n) (even) for j=1,...,m.
#' @examples 
#'trig<-fgentrig(36,36)
fgentrig<-function(n,m){
        tmp<-.Fortran(
                "triggen",
                as.integer(n),
                as.integer(m),
                double(2*n*m)
                )
        x<-tmp[[3]]
        dim(x)<-c(n,2*m)
        list(x)
}
