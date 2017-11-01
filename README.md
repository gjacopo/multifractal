[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.322420.svg)](https://doi.org/10.5281/zenodo.322420)
multifractal
============

Tools (C/Matlab) for multifractal analysis of 1D (time-series) and 2D (images) signals
---

Low-level image processing filters used in multifractal experiments/analyses over signals and images.

**About**

<table align="center">
    <tr> <td align="left"><i>version</i></td> <td align="left">1.0 <i>(non-active development)</i> </td> </tr> 
    <tr> <td align="left"><i>since</i></td> <td align="left">2002</td> </tr> 
    <tr> <td align="left"><i>license</i></td> <td align="left"><a href="https://joinup.ec.europa.eu/sites/default/files/eupl1.1.-licence-en_0.pdfEUPL">EUPL</a>  <i>(cite the source code or any of the references below!)</i> </td> </tr> 
</table>

**Description**

The code source files will enable you to reproduce some of the experiments mentioned in the references [below](References), namely for:
* multifractal analysis of signals/images using the so-called micro-wavelet,
* calculation of the multifractal singularity exponents and the multifractal spectrum,
* multifracal decomposition of signals/images,
* extraction of the so-called Most Singular Manifold,
* characterisation of the so-called chromatically reduced signal/image,
* multifractal reconstruction of signals/images.

See also [`singularity` project](https://github.com/gjacopo/singularity) for multifractal singularity analysis of 1D (time-series) signals.

**Usage**

To launch the programs:
* in `Matlab`, see also [`Contents.m`](matlab/Contents.m) for the whole list of functions; call [`run_multifractal.m`](matlab/run_multifractal.m) for an example of application; 
* in `C`, use the `Makefile` to build from source. 

**<a name="References"></a>References**

* Yahia H., Turiel A., Chrysoulakis N., Grazzini J., Prastacos P., and Herlin I. (2008): [**Application of the multifractal microcanonical formalism to the detection of fire plumes in NOAA-AVHRR data**](http://www.tandfonline.com/doi/abs/10.1080/01431160701840174), _International Journal of Remote Sensing_, 29(14):4189-4205, doi:[10.1080/01431160701840174](http://dx.doi.org/10.1080/01431160701840174).
* Grazzini J., Turiel A., Yahia H., and Herlin I. (2007): [**A multifractal approach for extracting relevant textural areas in satellite meteorological images**](http://www.sciencedirect.com/science/article/pii/S1364815205001970), _Environmental Modelling & Software_, 22(3):323-334, doi: [10.1016/j.envsoft.2005.07.032](http://dx.doi.org/10.1016/j.envsoft.2005.07.032).
* Turiel A., Pérez-Vicente C.J., and Grazzini J. (2006): [**Numerical methods for the estimation of singularity spectra on sampled data: a comparative study**](http://www.sciencedirect.com/science/article/pii/S0021999105005565), _Journal of Computational Physics_, 216(1):362-390, doi: [10.1016/j.jcp.2005.12.004](http://dx.doi.org/10.1016/j.jcp.2005.12.004).
* Grazzini J., Turiel A., and Yahia H. (2006): **Multifractal formalism for remote sensing: A contribution to the description and the understanding of meteorological phenomena in satellite images**, _Complexus Mundi. Emergent Patterns in Nature_, pp.247-256, _World Scientific_, doi: [10.1142/9789812774217_0022](http://dx.doi.org/10.1142/9789812774217_0022).
* Grazzini J. and Chrysoulakis N. (2005): [**Extraction of surface properties from a high accuracy DEM using multiscale remote sensing techniques**](http://enviroinfo.eu/sites/default/files/pdfs/vol111/0352.pdf), in _Proc. EnviroInfo_, pp.352-356.
* Turiel A., Grazzini J., and Yahia H. (2005): **Multiscale techniques for the detection of precipitation using thermal IR satellite images**, _IEEE Geoscience and Remote Sensing Letters_, 2(4):447-450, doi: [10.1109/LGRS.2005.852712](https://doi.org/10.1109/LGRS.2005.852712).
* Grazzini J., Turiel A., and Yahia H. (2005): **Presegmentation of high-resolution satellite images with a multifractal reconstruction scheme based on an entropy criterium**, in _Proc. IEEE ICPR_, pp.I-649-52, doi: [10.1109/ICIP.2005.1529834](https://doi.org/10.1109/ICIP.2005.1529834).
* Grazzini J., Turiel A., Yahia H., and Herlin I. (2004): [**Edge-preserving smoothing of high-resolution images with a partial multifractal reconstruction scheme**](http://www.isprs.org/proceedings/XXXV/congress/comm3/papers/435.pdf), in _Proc. ISPRS_, pp.1125-1129.
* Grazzini J. (2003): [**Analyses multiéchelle et multifractale d'images météorologiques : application à la détection de zones précipitantes**](https://tel.archives-ouvertes.fr/tel-00005940/file/tel-00006264.pdf) (in French), _Université de Marne-la-Vallée_, PhD thesis.
* Grazzini J., and Turiel A., and Yahia H. (2002): [**Entropy estimation and multiscale processing in meteorological satellite images**](http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1048103), in _Proc. IEEE ICPR_, vol. 3, pp.764-768, doi: [10.1109/ICPR.2002.1048103](https://doi.org/10.1109/ICPR.2002.1048103).
* Turiel A. and Parga, N. (2000): [**Multifractal wavelet filter of natural images**](https://arxiv.org/pdf/cond-mat/0107335.pdf), _Physical Review Letters_, 85(15):3325-8, doi:[10.1103/PhysRevLett.85.3325](https:/doi.org/10.1103/PhysRevLett.85.3325).
* Turiel A. and Parga, N. (2000): [**The multifractal structure of contrast changes in natural images: From sharp edges to textures**](http://www.lps.ens.fr/~risc/rescomp/Antonio/PAPERS/paper4.pdf), _Neural Computation_, 12(4):763-793, doi:[10.1162/089976600300015583](https://doi.org/10.1162/089976600300015583).
