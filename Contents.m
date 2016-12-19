% Multifractal analysis 
% 
% license: European Union Public License
%
% Main functions :
%  reduced_msm        - Compute the MSM and the chromatically reduced image of a greyscale image.
%  unitary            - Compute the MSM and the chromatically reduced image of an image with a naive 
%                       unitary gradient instead of the true one.
%  source             - Compute the source of an image, starting from the estimation of the MSM and  
%                       of the derivatives of the reconstructed image.
%  reconstruction     - Compute the reconstruction of an image.
%
% Analysis and reconstruction with the micro-wavelet:
%  fourier_cross      - Cross Fourier Transform or Inverse Cross  Fourier Transform with a micro-wavelet.
%  init_cross         - Initialise the micro-wavelet for analysis-reconstruction. 
%  convolution_cross  - Convolution with the micro-wavelet.
%  propagation        - Reconstruction of a signal suing the essential density and the universal propagator. 
%  propagation_cross  - Reconstruction of a local patch with the micro-wavelet.
%  derive_spectral    - Compute derivatives of an image by considering the frequency domain and by use of
%                       the micro-wavelet.
%  derive_msm         - Transform the Hausdorff measure of the MSM in the frequency space.
%  derive_msm_unitary - Compute the reduced_from_msm unitary essential vectorial field.
%  mask_gradient      - Apply the mask of MSM on the gradient  images.
%  mask_gradient_msm  - Compute an orientated MSM (or not, depending on the orientation of the gradient
%                       over MSM) and an orientated gradient vectors.
%  filter_spectral    - Transformation in the frequency domain: computes the value of a matrice after 
%                       having processed a scalar product with the norm of the frequency vector in the 
%                       frequency domain.
%                                  
% Estimation of the singularity exponents and the MSM:
%  msm                - Compute the MSM.
%  upm                - Compute the MSM and the density of the MSM in the image according to some 
%                       adjusting parameters.
%  distribution_UPM   - Compute the values and the distribution of singularity exponents.
%  quantile_threshold - Compute the adaptative threshold giving the value of the singularity exponents of
%                       the pixels of the MSM, according to the desired density.
%  msm_from_sing      - Compute the MSM and its density according to the distribution of the singularity 
%                       exponents.
%
% Computation of unitary and essential gradients:                                 
%  dummy              - Compute the reconstruction of the image starting from the reduced_from_msm unitary 
%                       essential vectorial field.
%  reduced_unitary    - Compute the chromatically reduced_from_msm image reconstructed using an unitary 
%                       gradient instead of the original one.
%  gradient_essential - Compute the essential gradient field 
%                                given the chromatically reduced_from_msm image.
%
% Extraction and representation of the source:
%  represent_source   - Represent the vectorial field of source by its norm and its orientation.
%  source_angle       - Compute the orientation of a vectorial field. Used for the representation of 
%      		        the source.
%  source_norm        - Compute the norm of a vectorial field (eventually log-representation).
%                                 
% Routines :
%  duplicate          - Duplicate the image edges so that it is periodic in both directions.
%  bits               - Return the nearest upper powers of 2.
%  clique             - Extract a patch from an image.
%  psnr               - Compute the PSNR ratio between an image and its approximation.
%  vec_normalise      - Normalise an image by its module (pixel by pixel).
%  vec_divide         - Compute the complex division of two vectorial fields in the complex domain.
%  vec_multiply       - Compute the complex multiplication of two vectorial fields in the complex domain.
%
% Exemple of use : launch :
%   1)    reduced_msm or unitary  
%   2)    source 
%   3)    reconstruction
