wgs2nztm <-
function(lon,lat){

# Constants, used to define the parameters for the International Ellipsoid
# used for the NZGD2000 datum (and hence for NZTM) 
#
# translated to R by WeiminJ, based on matlabe function WGS2NZTM by BenK 2011


if(length(lon)!=length(lat)) stop("Inputs must be the same size")
 a      <- 6378137  
 rf     <- 298.257222101 
 
 cm     <- 173.0*pi/180   #in radians
 lto    <- 0.0*pi/180  
 sf     <- 0.9996  
 fe     <- 1600000.0  
 fn     <- 10000000.0 
 utom   <- 1 
 
# ellipsoid calcs
 if (rf != 0.0){
     f = 1.0/rf 
 } else{
     f = 0.0
}
e2 = 2.0*f - f*f
ep2 = e2/( 1.0 - e2 ) 
om = meridian_arc(lto,e2,a) 

# convert to radians
lat=lat*pi/180
lon=lon*pi/180

# main code  
dlon  =  lon - cm 
while ( dlon > pi ) dlon = dlon-2*pi    
while (dlon < -pi ) dlon = dlon+2*pi    
m = meridian_arc(lat,e2,a)
slt = sin(lat)
eslt = (1.0-e2*slt*slt)
eta = a/sqrt(eslt)
rho = eta * (1.0-e2)/ eslt
psi = eta/rho

rm(eslt, rho)

clt = cos(lat)
rm(lon, lat)

wc = clt*dlon
wc2 = wc*wc
rm(wc)

t = slt/clt
#note commented out by bk to save memory, but slower!
# t2 = t.*t 
# t4 = t2.*t2 
# t6 = t2.*t4 

trm1 = (psi-t^2)/6.0

trm2 = (((4.0*(1.0-6.0*t^2)*psi 
              + (1.0+8.0*t^2))*psi 
              - 2.0*t^2)*psi+t^4)/120.0

trm3 = (61 - 479.0*t^2 + 179.0*t^4 - t^6)/5040.0

gce = (sf*eta*dlon*clt)*(((trm3*wc2+trm2)*wc2+trm1)*wc2+1.0)
rm( clt)
ce = gce/utom+fe
rm( gce, fe)

trm1 = 1.0/2.0

trm2 = ((4.0*psi+1)*psi-t^2)/24.0

trm3 = ((((8.0*(11.0-24.0*t^2)*psi 
            -28.0*(1.0-6.0*t^2))*psi 
            +(1.0-32.0*t^2))*psi  
            -2.0*t^2)*psi 
            +t^4)/720.0

trm4 = (1385.0-3111.0*t^2+543.0*t^4-t^6)/40320.0

gcn = (eta*t)*((((trm4*wc2+trm3)*wc2+trm2)*wc2+trm1)*wc2)
rm(trm1, trm2, trm3, trm4)
cn = (gcn+m-om)*sf/utom+fn
return(cbind(ce,cn))
}


meridian_arc<-function(lt,e2,a) {   
  # /.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*/
  # /.*                                                                         .*/
  # /.*  meridian_arc                                                           .*/
  # /.*                                                                         .*/
  # /.*  Returns the length of meridional arc (Helmert formula)                 .*/
  # /.*  Method based on Redfearn's formulation as expressed in GDA technical   .*/
  # /.*  manual at http://www.anzlic.org.au/icsm/gdatm/index.html               .*/
  # /.*                                                                         .*/
  # /.*  Parameters are                                                         .*/
  # /.*    projection                                                           .*/
  # /.*    latitude (radians)                                                   .*/
  # /.*                                                                         .*/
  # /.*  Return value is the arc length in metres                               .*/
  # /.*                                                                         .*/
  # /.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*/    
  e4 = e2*e2 
  e6 = e4*e2 
  
  A0 = 1 - (e2/4.0) - (3.0*e4/64.0) - (5.0*e6/256.0) 
  A2 = (3.0/8.0) * (e2+e4/4.0+15.0*e6/128.0) 
  A4 = (15.0/256.0) * (e4 + 3.0*e6/4.0) 
  A6 = 35.0*e6/3072.0 
  
  out=a*(A0*lt-A2*sin(2*lt)+A4*sin(4*lt)-A6*sin(6*lt)) 
  return(out)
}
foot_point_lat<-function(m,f,a){
  #  /.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*/
  # /.*                                                                       .*/
  # /.*   foot_point_lat                                                      .*/
  # /.*                                                                       .*/
  # /.*   Calculates the foot point latitude from the meridional arc          .*/
  # /.*   Method based on Redfearn's formulation as expressed in GDA technical.*/
  # /.*   manual at http://www.anzlic.org.au/icsm/gdatm/index.html            .*/
  # /.*                                                                       .*/
  # /.*   Takes parameters                                                    .*/
  # /.*      tm definition (for scale factor)                                 .*/
  # /.*      meridional arc (metres)                                          .*/
  # /.*                                                                       .*/
  # /.*   Returns the foot point latitude (radians)                           .*/                                                                        .*/
  # /.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*/
  
  n  = f/(2.0-f) 
  n2 = n*n 
  n3 = n2*n 
  n4 = n2*n2 
  
  g = a*(1.0-n)*(1.0-n2)*(1+9.0*n2/4.0+225.0*n4/64.0) 
  sig = m/g 
  
  phio = sig + (3.0*n/2.0 - 27.0*n3/32.0)*sin(2.0*sig) 
  + (21.0*n2/16.0 - 55.0*n4/32.0)*sin(4.0*sig) 
  + (151.0*n3/96.0) * sin(6.0*sig)
  + (1097.0*n4/512.0) * sin(8.0*sig)   
  return(phio)
}