DEGREE= 0.0174533;
INCONSISTY= "Inconsisty in the ";

func ysil_uniform_disc(x, rads, centrs, angs=)
{
  // Set-up
  y= transpose(x);
  a= rads(1, );
  b= rads(2, );
  n= dimsof(rads)(0);
  if (!is_array(centrs) && n== 1) centrs= [0, 0];
  if (is_array(angs)) angs*= DEGREE;
  else angs= array(double, n);
  dims= dimsof(x);
  discs= array(double, dims);
  disc= xc= yc= xr= yr= rho= array(double, dims)(, , -:1:n);

  for (i= 1; i<= n; ++i) {
    // Centres of discs
    xc(, , i)= x - centrs(, i)(1);
    yc(, , i)= y - centrs(, i)(2);
    // Rotate discs
    xr(, , i)=  cos(angs(i))*xc(, , i) + sin(angs(i))*yc(, , i);
    yr(, , i)= -sin(angs(i))*xc(, , i) + cos(angs(i))*yc(, , i);
    // Build discs
    rho(, , i)= sqrt(xr(, , i)*xr(, , i) / (a(i)*a(i)) + yr(, , i)*yr(, , i) / (b(i)*b(i)));
    idx= where(rho(, , i) <= 1.0);
    dsc= disc(, , i);
    dsc(idx)= 1.0;
    disc(, , i)= dsc;
    discs += disc(, , i);
  }

  // Normalise
  discs/= sum(discs);

  return discs;
}
/** ----------------------------------------- oOo ----------------------------------------- **/



func ysil_gaussian_disc(..)
{
  
}
/** ----------------------------------------- oOo ----------------------------------------- **/



func ysil_limb_darkened_disc(..)
{
}
/** ----------------------------------------- oOo ----------------------------------------- **/



func ysil_gaussian_ring(..)
{
  
}
/** ----------------------------------------- oOo ----------------------------------------- **/
