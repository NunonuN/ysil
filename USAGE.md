# Installation of YSIL

## Prerequisites
You must have [Yorick](http://dhmunro.github.io/yorick-doc) installed in order to use YSIL.

## Installation
Select a [version](https://github.com/NunonuN/ysil/releases) of YSIL, download it to any desired folder and unpack it.
For instance:

        wet https://github.com/NunonuN/ysil/archive/v${VERSION}.tar.gz
        tar xzvf ysil-${VERSION}.tar.gz

To use YSIL inside Yorick, it is sufficient to write in the corresponding environment:

        include, "ysil.i";

Make sure you indicate the full path to the `"ysil.i"` file.


# Quick instructions on how to use YSIL

## Using YSIL inside the Yorick interpreter
Launch Yorick and load `"ysil.i"`:

        include, "ysil.i";
        
### Examples
* Normalised uniform disc, FOV of 10mas with 500 pixels per side of the square array, X-axis of 4mas, Y-axis of 2mas, centred at the origin of the coordinate system, rotated 10 degrees:

        x= span(-5, 5, 500)(, -:1:500);
        disc= ysil_uniform_disc(x, [4., 2.], angs= 10);

<img src="/figures/uniform_disc.png" width="500" />

* Same as previous example, but with a Gaussian disc, X- and Y- axes respectively of 3 and 4mas, rotated -20 degrees:

        x= span(-5, 5, 500)(, -:1:500);
        disc= ysil_gaussian_disc(x, [3., 2.], angs= -20);

<img src="/figures/gaussian_disc.png" width="500" />

* Three normalised Gaussian rings, PSFs HWHMs of 0.5, 0.3 and 0.4mas, FOV of 20mas with 500 pixels per side of the square array, outer semi-X axes of 6.1, 3.5 and 3mas, outer semi-Y axes of 3.8, 1.5 and 1.7mas, inner semi-X axes of 5, 1.9 and 2mas, inner semi-Y axes of 2.9, 0.6 and 0.9mas, centred at [4.0, 0.0]mas, [-5.5, -5.5] and [-5.0, 6.0]mas (outer discs), and [4.2, 0.2]mas, [-5.5, -5.9]mas, [-5.0, 6.45]mas (inner discs), rings rotated by -60, 30 and 0 degrees:

        x= span(-10, 10, 500)(, -:1:500);
        rings= ysil_gaussian_ring(x,
                            [[6.1, 3.8], [3.5, 1.5], [3.0, 1.7]], // out_radi
                            [[5.0, 2.9], [1.9, 0.6], [2.0, 0.9]], // in_radi
                            [[4.0, 0.0], [-5.5, -5.5], [-5.0, 6.0] ], // centrs_out
                            [[4.2, 0.2], [-5.5, -5.9], [-5.0, 6.45]], // centrs_in
                            [0.5, 0.3, 0.4], // sigs
                            angs= [-60.0, 30.0, 0.0]
                           );

<img src="/figures/gaussian_rings.png" width="500" />

* Single stellar photosphere with radius of 3mas, centred at the FOV, with u= 0.2:

        x= span(-5, 5, 300)(, -:1:300);
        star= ysil_limb_darkened_disc(x, 3, [0, 0], u= 0.2);

<img src="/figures/star.png" width="500" />
