DEGREE= 0.0174533;

/* ================
 *  Error messages
 * ================ */
local BAD_DATA_TYPE, INCONSISTENCY, ORIGIN;
BAD_DATA_TYPE= "BAD DATA TYPE for ";
INCONSISTENCY= "INCONSISTENCY between ";
ORIGIN= " MUST CONTAIN THE ORIGIN.";
/** ---------------- **/


func ysil_gaussian_disc(x, sigs, centrs, angs=, amps=, norm=)
/* DOCUMENT discs= ysil_gaussian_disc(x, sigs, centrs, angs=, amps=, norm=);
 *
 * Returns a square array containing Gaussian elliptical discs, eventually rotated by some
 * angles. The array is built upon a (X, transpose(X)) mesh.
 * CENTRS defaults to [0, 0] if there is only one disc and no CENTRS array is input.
 * Otherwise, the CENTRS array shall have dimensions [2, 2, n], where n is the number of
 * discs. CENTRS(1, ) and CENTRS(2, ) correspond respectively to the X- and Y-coordinates of
 * the centres of all discs.
 * The keyword AMPS is facultative. If it is ommited, the amplitude of Gaussian i is equal to
 *
 *                                   1.0
 *                       -----------------------------.
 *                          -------------------------
 *                        \/ 2 * π * σ_x(i) * σ_y(i)
 *
 * Otherwise, each Gaussian i is equal to
 *
 *                            1/2*(xr^2/σx^2 + yr^2/σy^2)
 *                amps(i) * e^                           ,
 *
 *
 * where 'xr' and 'yr' are the X- and Y-coordinates rotated by angle ANGS and displaced by
 * CENTRS (expressed in X-units).
 *
 * VARIABLES
 * x= double, square array of X-coordinates.
 * sigs= double, the HWHM of the discs (semi-minor and semi-major axes), in units of X.
 * centrs= double, the centres of the discs, in units of X.
 *
 * KEYWORDS
 * angs= Double, angles of rotation of the ellipses, counted anticlockwise from the X-axis,
 *       in degrees. Default: angs= array(double, number_of_discs).
 * amps= Double, amplitudes of the Gaussian discs.
 * norm= int, if TRUE (1, 1.0 or 1n) the final array is normalised to the total flux.
 *       Default: norm= 1n.
 *
 * EXAMPLES
 * 1. Uniform disc, FOV of 10mas with 500 pixels per side of the square array, X-axis of 4mas,
 *    Y-axis of 2mas, centred at the origin of the coordinate system, rotated 10 degrees.
 *
 *    x= span(-5, 5, 500)(, -:1:500);
 *    disc= ysil_gaussian_disc(x, [4., 2.], angs= 10);
 *
 * 2. Three uniform discs, FOV of 20mas with 300 pixels per side of the square array, X-axes
 *    of 2, 3.5 and 5.6mas, Y-axes of 4, 3.5 and 1.3mas, centred at [-6, 0]mas, [5, 5.5]mas
 *    and [4, -7]mas, rotated by 30, 0 and -160degrees.
 *
 *    x= span(-10, 10, 300)(, -:1:300);
 *    discs= ysil_gaussian_disc(x, [[2, 4], [3.5, 3.5], [5.6, 1.3]], [[-6, 0], [5, 5.5],
 *                              [4, -7]], angs= [30, 0, -160]);
 *
 * HISTORY
 * Revision 0.1 2016/12/17 19:30:53 Nuno Gomes
 *  - First release of the function.
 *
 * SEE ALSO: ysil_uniform_disc.
 ** --------------------------------------------------------------------------------------- **/
{
  // Set-up
  norm= (is_void(norm) ? 1n : int(norm));
  local y, sigx, sigy, n;
  y= transpose(x);
  sigx= sigs(1, );
  sigy= sigs(2, );
  if (dimsof(sigs)(1)== 1) sigs= sigs(, -);
  n= dimsof(sigs)(0);
  if (!is_array(centrs) && n== 1) centrs= [0, 0];
  if (is_array(angs)) angs*= DEGREE;
  else angs= array(double, n);
  // Checks
  if (!is_numerical(x)     ) error, BAD_DATA_TYPE+"'x'.";
  if (max(x(, 1)) * min(x(, 1))> 0.0) error, "'x'"+ORIGIN;
  if (!is_numerical(sigx)  ) error, BAD_DATA_TYPE+"'sigs'.";
  if (!is_numerical(sigy)  ) error, BAD_DATA_TYPE+"'sigs'.";
  if (dimsof(centrs)(1)!= 1 && dimsof(centrs)(0)!= n)
    error, INCONSISTENCY+"'sigs' and 'centrs'.";
  if (numberof(angs)!= n   ) error, INCONSISTENCY+"'sigs' and 'angs'.";
  if (is_array(amps)) {
    if (!is_numerical(amps)) error, BAD_DATA_TYPE+"'amps'.";
    if (numberof(amps)!= n ) error, INCONSISTENCY+"'sigs' and 'amps'.";
  }
  
  // Build Gaussians
  local dims, discs, disc, xc, yc, xr, yr, rho;
  dims= dimsof(x);
  discs= array(double, dims);
  disc= xc= yc= xr= yr= rho= array(double, dims)(, , -:1:n);

  for (i= 1; i<= n; ++i) {
    // centres of discs
    xc(, , i)= x - centrs(, i)(1);
    yc(, , i)= y - centrs(, i)(2);
    // rotate discs
    xr(, , i)=  cos(angs(i))*xc(, , i) + sin(angs(i))*yc(, , i);
    yr(, , i)= -sin(angs(i))*xc(, , i) + cos(angs(i))*yc(, , i);
    // build Gaussian discs
    disc(, , i)= exp(-0.5*(xr(, , i)*xr(, , i)/(sigx(i)*sigx(i)) +
                           yr(, , i)*yr(, , i)/(sigy(i)*sigy(i))));
    if (is_array(amps)) {
      disc(, , i)*= amps(i);
    } else {
      local amp;
      amp= 1.0 / sqrt(2.0 * PI * sigx(i) * sigy(i));
      disc(, , i)*= amp;
    }

    discs += disc(, , i);
  }

  // Normalise
  if (norm) discs/= sum(discs);

  return discs;
}
/** ----------------------------------------- oOo ----------------------------------------- **/



func ysil_limb_darkened_disc(x, radi, centrs, l0=, u=, norm=)
/* DOCUMENT stars= ysil_limb_darkened_disc(x, radi, centrs, l0=, u=, norm=);
 *
 * Returns square array containing limb-darkened discs, that can be used to simulate stellar
 * photospheres.
 * The equation used for limb darkening was adapted from "Physics Topics" by J. B. Tatum and
 * "Astrophysical Quantities" by Arthur N. Cox, and it is written as
 *
 *   I(r) = I(0) * [1 - u * (1 - sqrt[(rad² - l²) / rad²])].
 *
 * where 'r' is the radial distance from the centre of the disc, 'u' is the limb darkening
 * coefficient, 'rad' is the radius of the star, 'l' is the darkening radial distance from the
 * centre of the star, and 'I' is the brightness distribution at the surface.
 * 'l0' indicates the radius of the surface where the limb darkening is observed (behind it
 * there is an uniform disc of radius RAD).
 *
 * The same value of U is used for all stars if only on value is input.
 *
 * VARIABLES
 * x= double, square array of X-coordinates.
 * radi= double, the radii of the discs, in units of X.
 * centrs= double, the centres of the discs, in units of X.
 *
 * KEYWORDS
 * l0= double, darkening radial distance from the centre of the star, in units of X.
 *     Default: l0= radi.
 * u= double, limb darkening coefficients.
 * norm= int, if TRUE (1, 1.0 or 1n) the final array is normalised to the total flux.
 *       Default: norm= 1n.
 *
 * EXAMPLES
 * 1. Single stellar photosphere with radius of 3mas, centred at the FOV, with u= 0.2.
 *
 *    x= span(-10, 10, 300)(, -:1:300);
 *    star= ysil_limb_darkened_disc(x, 3, [0, 0], u= 0.2, norm= 1n);
 *
 * 2. Three stellar photospheres with radii [2.8, 4.3, 1.5]mas, at positions [[-5, -5],
 *    [3.5, 2], [-6.5, 7], with u= [0.1, 0.3, 0.7].
 *
 *    x= span(-5, 5, 400)(, -:1:400);
 *    stars= ysil_limb_darkened_disc(x, [2.8, 4.3, 1.5], [[-5, -5], [3.5, 2], [-6.5, 7]],
 *                                   u= [0.1, 0.3, 0.7], norm= TRUE);
 * 
 * HISTORY
 * Revision 0.1 2016/12/18 09:01:46 Nuno Gomes
 *  - First release of the function.
 *
 * SEE ALSO: ysil_uniform_disc.
 ** --------------------------------------------------------------------------------------- **/ 
{
  // Set-up
  if (is_void(l0)) l0= radi;
  norm= (is_void(norm) ? 1n : int(norm));
  local y, dims, n, rads;
  y= transpose(x);
  dims= dimsof(x);
  n= numberof(radi);
  rads= transpose(radi(, -:1:2));
  // Checks
  if (!is_numerical(x)        ) error, BAD_DATA_TYPE+"'x'.";
  if (max(x(, 1)) * min(x(, 1))> 0.0) error, "'x'"+ORIGIN;
  if (!is_numerical(rads(1, ))) error, BAD_DATA_TYPE+"'radi'.";
  if (!is_numerical(rads(2, ))) error, BAD_DATA_TYPE+"'radi'.";
  if (dimsof(centrs)(1)!= 1 && dimsof(centrs)(0)!= n)
    error, INCONSISTENCY+"'radi' and 'centrs'.";
  if (!is_numerical(l0)       ) error, BAD_DATA_TYPE+"'l0'.";
  if (numberof(u)< n) u= u(, -:1:n);
  if (numberof(u)!= n         ) error, INCONSISTENCY+"'radi' and 'u'.";
  
  // Build stars
  local star, stars, mask, cx, cy, xpos, ypos, l2, idx, mu;
  star=  array(double, dims)(, , -:1:n);
  stars= array(double, dims);
  
  for (i= 1; i<= n; ++i) {
    star(, , i)= ysil_uniform_disc(x, rads(, i), centrs(, i), norm= 0n);
    mask= double(star(, , i) > 0.0); // mask of ones where stars exists
    // displace darkened region
    cx= centrs(, i)(1);
    cy= centrs(, i)(2);
    xpos= x - cx;
    ypos= y - cy;
    // apply equation
    l2= xpos*xpos + ypos*ypos;
    idx= where(l2 > l0(i)*l0(i));
    l2(idx)= 0.0;
    mu= sqrt((radi(i)*radi(i) - l2) / (radi(i)*radi(i)));
    star(, , i)*= (1.0 - u(i)*(1.0 - mu));
    stars+= star(, , i);
  }

  // Normalise
  if (norm) stars/= sum(stars);
  
  return stars;
}
/** ----------------------------------------- oOo ----------------------------------------- **/



func ysil_uniform_disc(x, radi, centrs, angs=, norm=)
/* DOCUMENT discs= ysil_uniform_disc(x, radi, centrs, angs=, norm=);
 *
 * Returns a square array containing uniform elliptical discs, eventually rotated by some
 * angle. The array is built upon a (X, transpose(X)) mesh.
 * CENTRS defaults to [0, 0] if there is only one disc and no CENTRS array is input.
 * Otherwise, the CENTRS array shall have dimensions [2, 2, n], where n is the number of
 * discs. CENTRS(1, ) and CENTRS(2, ) correspond respectively to the X- and Y-coordinates of
 * the centres of all discs.
 *
 * VARIABLES
 * x= double, square array of X-coordinates.
 * radi= double, the radii of the discs (semi-minor and semi-major axes), in units of X.
 * centrs= double, the centres of the discs, in units of X.
 *
 * KEYWORDS
 * angs= double, angles of rotation of the ellipses, counted anticlockwise from the X-axis,
 *       in degrees. Default: angs= array(double, number_of_discs).
 * norm= int, if TRUE (1, 1.0 or 1n) the final array is normalised to the total flux.
 *       Default: norm= 1n.
 *
 * EXAMPLES
 * 1. Uniform disc, FOV of 10mas with 500 pixels per side of the square array, X-axis of 4mas,
 *    Y-axis of 2mas, centred at the origin of the coordinate system, rotated 10 degrees.
 *
 *    x= span(-5, 5, 500)(, -:1:500);
 *    disc= ysil_uniform_disc(x, [4., 2.], angs= 10);
 *
 * 2. Three uniform discs, FOV of 20mas with 300 pixels per side of the square array, X-axes
 *    of 2, 3.5 and 5.6mas, Y-axes of 4, 3.5 and 1.3mas, centred at [-6, 0]mas, [5, 5.5]mas
 *    and [4, -7]mas, rotated by 30, 0 and -160degrees.
 *
 *    x= span(-10, 10, 300)(, -:1:300);
 *    discs= ysil_uniform_disc(x, [[2, 4], [3.5, 3.5], [5.6, 1.3]], [[-6, 0], [5, 5.5],
 *                             [4, -7]], angs= [30, 0, -160]);
 *
 * HISTORY
 * Revision 0.1 2016/12/16 19:05:28 Nuno Gomes
 *  - First release of the function.
 *
 * SEE ALSO: ysil_gaussian_disc.
 ** --------------------------------------------------------------------------------------- **/
{
  // Set-up
  norm= (is_void(norm) ? 1n : int(norm));
  local y, a, b, n;
  y= transpose(x);
  a= radi(1, );
  b= radi(2, );
  if (dimsof(radi)(1)== 1) radi= radi(, -);
  n= dimsof(radi)(0);
  if (!is_array(centrs) && n== 1) centrs= [0, 0];
  if (is_array(angs)) angs*= DEGREE;
  else angs= array(double, n);
  // Checks
  if (!is_numerical(x)     ) error, BAD_DATA_TYPE+"'x'.";
  if (max(x(, 1)) * min(x(, 1))> 0.0) error, "'x'"+ORIGIN;
  if (!is_numerical(a)     ) error, BAD_DATA_TYPE+"'radi'.";
  if (!is_numerical(b)     ) error, BAD_DATA_TYPE+"'radi'.";
  if (dimsof(centrs)(1)!= 1 && dimsof(centrs)(0)!= n)
    error, INCONSISTENCY+"'radi' and 'centrs'.";
  if (numberof(angs)!= n   ) error, INCONSISTENCY+"'radi' and 'angs'.";

  // Build discs
  local dims, discs, disc, xc, yc, xr, yr, rho;
  dims= dimsof(x);
  discs= array(double, dims);
  disc= xc= yc= xr= yr= rho= array(double, dims)(, , -:1:n);

  for (i= 1; i<= n; ++i) {
    // centres of discs
    xc(, , i)= x - centrs(, i)(1);
    yc(, , i)= y - centrs(, i)(2);
    // rotate discs
    xr(, , i)=  cos(angs(i))*xc(, , i) + sin(angs(i))*yc(, , i);
    yr(, , i)= -sin(angs(i))*xc(, , i) + cos(angs(i))*yc(, , i);
    // build discs
    local idx;
    rho(, , i)= sqrt(xr(, , i)*xr(, , i) / (a(i)*a(i)) + yr(, , i)*yr(, , i) / (b(i)*b(i)));
    idx= where(rho(, , i)<= 1.0);
    dsc= disc(, , i);
    dsc(idx)= 1.0;
    disc(, , i)= dsc;
    discs += disc(, , i);
  }

  // Normalise
  if (norm) discs/= sum(discs);

  return discs;
}
/** ----------------------------------------- oOo ----------------------------------------- **/
