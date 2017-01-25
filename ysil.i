/*********************************************************************************************
 *********************************************************************************************
 ****                                                                                     ****
 **                                          YSIL                                           **
 **                             Yorick Synthetic Image Library                              **
 **                                                                                         **
 ****                                                                                     ****
 *********************************************************************************************
 *********************************************************************************************
 ****                                                                                     ****
 **                                                                                         **
 ** DOCUMENT YSIL v0.3                                                                      **
 **                                                                                         **
 ** A Yorick Synthetic Image Library for the interferometric image reconstruction cookbook. **
 ** The package contains basic functions that can be used to build astronomical objects,    **
 ** such as stellar clusters, stellar photospheres or young stellar objects.                **
 **                                                                                         **
 ** Total number of functions in the package: 5.                                            **
 **                                                                                         **
 ** LIST OF FUNCTIONS                                                                       **
 ** =================                                                                       **
 ** ysil_convolve..................... Computes convolution between arrays.                 **
 ** ysil_gaussian_disc................ Creates discs with Gaussian profile.                 **
 ** ysil_gaussian_ring................ Creates rings with Gaussian profile.                 **
 ** ysil_limb_darkened_disc........... Creates discs with limb darkened profile.            **
 ** ysil_uniform_disc................. Creates uniform discs.                               **
 **                                                                                         **
 ** HISTORY                                                                                 **
 ** Revision 0.3 2017/01/24 22:39:21 Nuno Gomes                                             **
 **  - ADDED function 'ysil_gaussian_ring'.                                                 **
 ** Revision 0.2 2017/01/14 19:35:17 Nuno Gomes                                             **
 **  - ADDED function 'ysil_convolve'.                                                      **
 ** Revision 0.1 2016/12/23 13:20:20 Nuno Gomes                                             **
 **  - FIRST RELEASE of 'ysil.i'.                                                           **
 **                                                                                         **
 ****                                                                                     ****
 *********************************************************************************************
 *********************************************************************************************/


DEGREE= 0.0174533; // conversion from degrees to radians


/* ================
 *  Error messages
 * ================ */
local BAD_DATA_TYPE, INCONSISTENCY, ORIGIN;
BAD_DATA_TYPE= "BAD DATA TYPE for ";
INCONSISTENCY= "INCONSISTENCY between ";
ORIGIN= " MUST CONTAIN THE ORIGIN.";
/** ---------------- **/




/* ===========
 *  FUNCTIONS
 * =========== */
func ysil_convolve(img, filt, norm, rmnd=)
/* DOCUMENT new_img= ysil_convolve(img, filt, norm, rmnd=);
 *
 * Computes the discrete convolution between array IMG and kernel FILT by means of a FFT.
 * Note that both arrays must have the same dimensions.
 * If NORM is not nil or zero, the convolution is normalised.
 * 'rmnd' controls how much is added to the sum of the convolution in case it is zero and
 * NORM is equal to 1. It defaults to 1e-8.
 *
 * VARIABLES
 * img= array of doubles, to be convolved with FILT.
 * filt= array of doubles, to be convolved with IMG.
 * norm= integer, works as a switch. If on (non nil), the convolution is normalised.
 *       Default: norm= 0.
 *
 * KEYWORDS
 * rmnd= double, number to be added to the sum of the convolution array in case it is zero
 *       and NORM is 1. Default: rmnd= 1e-8.
 *
 * HISTORY
 * Revision 0.1 2017/01/14 19:02:41 Nuno Gomes
 *  - FIRST RELEASE of the function.
 *
 * SEE ALSO: ysil_gaussian_ring, fft, roll.
 ** --------------------------------------------------------------------------------------- **/
{
  local real, dims, wspc, conv, sumconv;
  norm= (!norm ? 0 : 1);
  rmnd= (!rmnd ? 1e-8 : double(rmnd));
  real= (!is_complex(img) && !is_complex(filt));
  dims= dimsof(img);
  wspc= fft_setup(dims);
  conv= fft( fft(img, setup= wspc) * fft(roll(filt), setup= wspc), -1, setup= wspc );

  if (real) conv= double(conv);

  conv= 1.0/numberof(conv) * conv;
  sumconv= sum(conv);

  if (norm) conv/= (sumconv== 0 ? sumconv + rmnd : sumconv);

  return conv;
}
/** ----------------------------------------- oOo ----------------------------------------- **/



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
 *  - FIRST RELEASE of the function.
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



func ysil_gaussian_ring(x, out_radi, in_radi, centrs_out, centrs_in, sigs, angs=, norm=)
/* DOCUMENT rings= ysil_gaussian_ring(x, out_radi, in_radi, centrs_out, centrs_in, sigs,
 *                                    angs=, norm=);
 *
 * Returns a square array containing Gaussian elliptical rings, eventually rotated by some
 * angle. The array is built upon a (X, transpose(X)) mesh.
 * The widths of the rings are defined by variables OUT_RADI and IN_RADI, ie., the ring is
 * built by subtracting a uniform elliptical disc of semi- major and minor axes IN_RADI to
 * a uniform elliptical disc of semi- major and minor axes OUT_RADI.
 * CENTRS_OUT and CENTRS_IN arrays shall have dimensions [2, 2, n], where n is the number of
 * rings. CENTRS_OUT/IN(1, ) and CENTRS_OUT/IN(2, ) correspond respectively to the X- and
 * Y-coordinates of the centres of all outer/inner rings.
 * The ring is created by convolving an uniform ring with a Gaussian point spread function
 * (PSF).
 *
 * VARIABLES
 * x= double, square array of X-coordinates.
 * out_radi= double, the radii of the outer discs (semi-minor and semi-major axes), in units
 *           of X.
 * in_radi= double, the radii of the inner discs (semi-minor and semi-major axes), in units
 *          of X.
 * centrs_out= double, the centres of the outer discs, in units of X.
 * centrs_in= double, the centres of the inner discs, in units of X.
 * sigs= double, HWHM of the Gaussian used as PSF.
 *
 * KEYWORDS
 * angs= double, angles of rotation of the ellipses, counted anticlockwise from the X-axis,
 *       in degrees. Default: angs= array(double, number_of_discs).
 * norm= int, if TRUE (1, 1.0 or 1n) the final array is normalised to the total flux.
 *       Default: norm= 1n.
 *
 * EXAMPLES
 * 1. Gaussian ring, HWHM of 1mas, FOV of 10mas with 300 pixels per side of the square array,
 *    outer semi-X and -Y axes of 4.5mas and 3.5mas, respectively, inner semi-X and -Y axes of
 *    3mas and 2mas, centred at the origin of the coordinate system, inner void vertically
 *    displaced -0.5mas, all structure rotated 15 degrees.
 *
 *    x= span(-5, 5, 400)(, -:1:400);
 *    ring= ysil_gaussian_ring(x, [4.5, 3.5], [3., 2.], [0, 0], [0, -0.5], 1.0, angs= 15);
 *
 * 2. Two normalised Gaussian rings, PSFs HWHMs of 0.5 and 0.4mas, FOV of 20mas with 500 pixels
 *    per side of the square array, outer semi-X axes of 6.1 and 3mas, outer semi-Y axes of 3.8
 *    and 1.7mas, inner semi-X axes of 5 and 2mas, inner semi-Y axes of 2.9 and 0.9mas, centred
 *    at [-0.3, -3.5]mas and [1, 6]mas (outer discs), and [-0.2, -3.3]mas and [1, 6.45]mas
 *    (inner discs), rings rotated by -40 and 0 degrees.
 *
 *    x= span(-10, 10, 500)(, -:1:500);
 *    rings= ysil_gaussian_ring(x, [[6.1, 3.8], [3.0, 1.7]], [[5.0, 2.9], [2.0, 0.9]],
 *                                 [[-0.3, -3.5], [1.0, 6.0] ], [[-0.2, -3.3], [1.0, 6.45]],
 *                                 [0.5, 0.4], angs= [-40.0, 0.0], norm= 1);
 *
 * HISTORY
 * Revision 0.2 2017/01/24 21:44:19 Nuno Gomes
 *  - UPDATED code and documentation.
 * Revision 0.1 2016/12/26 17:48:39 Nuno Gomes
 *  - FIRST RELEASE of the function.
 *
 * SEE ALSO: ysil_gaussian_disc.
 ** --------------------------------------------------------------------------------------- **/
{
  // Set-up
  norm= (is_void(norm) ? 1n : int(norm));
  local y, a, b, n;
  y= transpose(x);
  a_out= out_radi(1, );
  b_out= out_radi(2, );
  a_in= in_radi(1, );
  b_in= in_radi(2, );
  if (dimsof(out_radi)(1)== 1) out_radi= out_radi(, -);
  if (dimsof(in_radi)(1)== 1) in_radi= in_radi(, -);
  n= dimsof(out_radi)(0);
  if (dimsof(centrs_out)(1)== 1) centrs_out= centrs_out(, -);
  if (dimsof(centrs_in)(1)== 1) centrs_in= centrs_in(, -);
  if (is_array(angs)) angs*= DEGREE;
  else angs= array(double, n);
  // Checks
  if (!is_numerical(x)) error, BAD_DATA_TYPE+"'x'.";
  if (max(x(, 1)) * min(x(, 1))> 0.0) error, "'x'"+ORIGIN;
  if (!is_numerical(a_out)) error, BAD_DATA_TYPE+"'out_radi'.";
  if (!is_numerical(b_out)) error, BAD_DATA_TYPE+"'out_radi'.";
  if (!is_numerical(a_in)) error, BAD_DATA_TYPE+"'in_radi'.";
  if (!is_numerical(b_in)) error, BAD_DATA_TYPE+"'in_radi'.";
  if (dimsof(centrs_out)(1)!= 1 && dimsof(centrs_out)(0)!= n)
    error, INCONSISTENCY+"'radi' and 'centrs'.";
  if (!is_numerical(sigs)) error, BAD_DATA_TYPE+"'sigs'.";
  if (numberof(sigs)!= n) error, INCONSISTENCY+"'radi' and 'sigs'.";
  if (numberof(angs)!= n) error, INCONSISTENCY+"'radi' and 'angs'.";

  // Build discs
  local dims, rings, disc_out, disc_in, ring, psf;
  local xc_out, yc_out, xc_in, yc_in, xr_out, yr_out, xr_in, yr_in, rho_out, rho_in;
  dims= dimsof(x);
  rings= array(double, dims);
  disc_out= disc_in= ring= psf= array(double, dims)(, , -:1:n);
  xc_out= yc_out= xc_in= yc_in= xr_out= yr_out= xr_in= yr_in= rho_out= rho_in=
    array(double, dims)(, , -:1:n);

  for (i= 1; i<= n; ++i) {
    // centres of rings
    xc_out(, , i)= x - centrs_out(, i)(1);
    yc_out(, , i)= y - centrs_out(, i)(2);
    xc_in(, , i)= x - centrs_in(, i)(1);
    yc_in(, , i)= y - centrs_in(, i)(2);
    // rotate rings
    xr_out(, , i)=  cos(angs(i))*xc_out(, , i) + sin(angs(i))*yc_out(, , i);
    yr_out(, , i)= -sin(angs(i))*xc_out(, , i) + cos(angs(i))*yc_out(, , i);
    xr_in(, , i)=  cos(angs(i))*xc_in(, , i) + sin(angs(i))*yc_in(, , i);
    yr_in(, , i)= -sin(angs(i))*xc_in(, , i) + cos(angs(i))*yc_in(, , i);
    // build rings
    local idx_out, idx_in, dsc_out, dsc_in;
    rho_out(, , i)= sqrt(xr_out(, , i)*xr_out(, , i) / (a_out(i)*a_out(i)) +
                         yr_out(, , i)*yr_out(, , i) / (b_out(i)*b_out(i)));
    rho_in(, , i)= sqrt(xr_in(, , i)*xr_in(, , i) / (a_in(i)*a_in(i)) +
                         yr_in(, , i)*yr_in(, , i) / (b_in(i)*b_in(i)));
    idx_out= where(rho_out(, , i)<= 1.0);
    idx_in= where(rho_in(, , i)<= 1.0);
    dsc_out= disc_out(, , i);
    dsc_in= disc_in(, , i);
    dsc_out(idx_out)= 1.0;
    dsc_in(idx_in)= 1.0;
    disc_out(, , i)= dsc_out;
    disc_in(, , i)= dsc_in;
    ring(, , i)= disc_out(, , i) - disc_in(, , i);
    psf(, , i)= ysil_gaussian_disc(x, [sigs(i), sigs(i)]);
    ring(, , i)= ysil_convolve(ring(, , i), psf(, , i));
    ring(, , i)+= abs(min(ring(, , i)));
    rings += ring(, , i);
  }

  // Normalise
  if (norm) rings/= sum(rings);

  return rings;
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
 *  - FIRST RELEASE of the function.
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
 *  - FIRST RELEASE of the function.
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
