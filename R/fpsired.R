#' Calculation of Hampel's redescending psi function
#'
#' @param x The point at which the psi function is evaluated.
#' @param cnr The parameters of Hampel's redescending psi function}
#' @return rpsi The value of the psi function
#' @examples 
#' fpsired(1,c(1,3,5))
fpsired<-function(x,cnr){
	xx<-abs(x)
      	if(xx<=cnr[1]){rpsi<-xx}
      	if((xx>cnr[1])&(xx<=cnr[2])){rpsi=cnr[1]}
      	if((xx>cnr[2])&(xx<=cnr[3])){rpsi<-cnr[1]*(cnr[3]-xx)/(cnr[3]-cnr[2])}
      	if(xx>cnr[3]){rpsi<-0}
	if(x<0){rpsi<- -rpsi}	
      return(rpsi)
}
 