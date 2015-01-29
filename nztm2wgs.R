nztm2wgs <-
function(ce,cn){
#
# R Function to convert NZTM Easting and Northing co-ordinates to 
# WGS84 latlon (in decimal degrees)
#
# original function in Matlab, Adapted from LINZ C algorithms by Ben Knight 2007
# see following website for details of method:
# http://www.linz.govt.nz/core/surveysystem/geodeticinfo/geodeticdatums/index.html
#
# coded in R by WeiminJ, 2011
#
# usage: latlon<-nztm2wgs(ce,cn)
#
if (length(ce)!=length(cn))  stop('Inputs must be the same size')

# Constants, used to define the parameters for the International Ellipsoid
# used for the NZGD2000 datum (and hence for NZTM) 
 a <-6378137 
 rf <- 298.257222101
 
 cm  <-   173.0*pi/180  #in radians
 lto  <-  0.0*pi/180   #in radians
 sf   <-  0.9996 
 fe  <-   1600000.0 
 fn   <-  10000000.0
 utom <- 1
 
## ellipsoid calcs
 if(rf!= 0.0)   f <- 1.0/rf 
 else   f <- 0.0
e2 <- 2.0*f - f*f
ep2 <- e2/( 1.0 - e2 )
om <- meridian_arc(lto,e2,a)

 
## main code  
cn1  <-  (cn - fn)*utom/sf + om
fphi <- foot_point_lat(cn1,f,a)
slt <- sin(fphi)
clt <- cos(fphi)

eslt <- (1.0-e2*slt*slt)
eta <- a/sqrt(eslt)
rho <- eta * (1.0-e2) / eslt
psi <- eta/rho

E <- (ce-fe)*utom
x <- E/(eta*sf)
x2 <- x*x

t <- slt/clt
t2 <- t*t
t4 <- t2*t2

trm1 <- 1.0/2.0

trm2 <- ((-4.0*psi 
             +9.0*(1-t2))*psi 
             +12.0*t2)/24.0

trm3 <- ((((8.0*(11.0-24.0*t2)*psi 
              -12.0*(21.0-71.0*t2))*psi 
              +15.0*((15.0*t2-98.0)*t2+15))*psi 
              +180.0*((-3.0*t2+5.0)*t2))*psi + 360.0*t4)/720.0

trm4 <- (((1575.0*t2+4095.0)*t2+3633.0)*t2+1385.0)/40320.0

lat <- fphi+(t*x*E/(sf*rho))*(((trm4*x2-trm3)*x2+trm2)*x2-trm1)

trm1 <- 1.0

trm2 <- (psi+2.0*t2)/6.0

trm3 <- (((-4.0*(1.0-6.0*t2)*psi 
           +(9.0-68.0*t2))*psi 
           +72.0*t2)*psi 
           +24.0*t4)/120.0

trm4 <- (((720.0*t2+1320.0)*t2+662.0)*t2+61.0)/5040.0

lon <- cm - (x/clt)*(((trm4*x2-trm3)*x2+trm2)*x2-trm1)

# convert to degrees
lat<-lat*180/pi
lon<-lon*180/pi 
return(cbind(lat,lon))
}



meridian_arc<-function(lt,e2,a) {   
  # /***************************************************************************/
  # /*                                                                         */
  # /*  meridian_arc                                                           */
  # /*                                                                         */
  # /*  Returns the length of meridional arc (Helmert formula)                 */
  # /*  Method based on Redfearn's formulation as expressed in GDA technical   */
  # /*  manual at http://www.anzlic.org.au/icsm/gdatm/index.html               */
  # /*                                                                         */
  # /*  Parameters are                                                         */
  # /*    projection                                                           */
  # /*    latitude (radians)                                                   */
  # /*                                                                         */
  # /*  Return value is the arc length in metres                               */
  # /*                                                                         */
  # /***************************************************************************/    
  e4 <- e2*e2
  e6 <- e4*e2
  
  A0 <- 1 - (e2/4.0) - (3.0*e4/64.0) - (5.0*e6/256.0)
  A2 <- (3.0/8.0) * (e2+e4/4.0+15.0*e6/128.0)
  A4 <- (15.0/256.0) * (e4 + 3.0*e6/4.0)
  A6 <- 35.0*e6/3072.0
  
  out<-a*(A0*lt-A2*sin(2*lt)+A4*sin(4*lt)-A6*sin(6*lt))
  return(out)
}
foot_point_lat<-function(m,f,a){
  #  /*************************************************************************/
  # /*                                                                       */
  # /*   foot_point_lat                                                      */
  # /*                                                                       */
  # /*   Calculates the foot point latitude from the meridional arc          */
  # /*   Method based on Redfearn's formulation as expressed in GDA technical*/
  # /*   manual at http://www.anzlic.org.au/icsm/gdatm/index.html            */
  # /*                                                                       */
  # /*   Takes parameters                                                    */
  # /*      tm definition (for scale factor)                                 */
  # /*      meridional arc (metres)                                          */
  # /*                                                                       */
  # /*   Returns the foot point latitude (radians)                           */                                                                        */
  # /*************************************************************************/
  
  n  <- f/(2.0-f)
  n2 <- n*n
  n3 <- n2*n
  n4 <- n2*n2
  
  g <- a*(1.0-n)*(1.0-n2)*(1+9.0*n2/4.0+225.0*n4/64.0)
  sig <- m/g
  
  phio <- sig + (3.0*n/2.0 - 27.0*n3/32.0)*sin(2.0*sig) 
  + (21.0*n2/16.0 - 55.0*n4/32.0)*sin(4.0*sig) 
  + (151.0*n3/96.0) * sin(6.0*sig) 
  + (1097.0*n4/512.0) * sin(8.0*sig)   
  return(phio)   
}