function [yFWHM, PWb, OnAxisIntensity, CumulativeFWHM, Overall_Dithered_xz_OTF2, DitheredxzPSFCrossSection, ...
    PlotSLMPhase, MaskIntensity, SampleIntensityAtFocus, DitheredIntensity, OverallCrossSection, SheetCrossSection, DetPSFCrossSection, xzSIMOTF5] = ...
    Calc_and_Plot_3D_LLS_PSFs_and_OTFs_v9e_cleaned(algorithm, lightsheet_descrip, xyPol, PW, NAlattice, NAsigma, NAannulus, crop_factor, NAdet, ...
    xz_det_PSF, xz_det_OTF, gamma, yplot_range, ystepsize_in, outputdir, subfolder)
%
%CALC_AND_PLOT_3D_LLS_PSFs_and_OTFs_v9e_cleaned
%Calculates the PSFs, OTFs, and their propagation characterisitics for of any 
%    desired light sheet. Light sheets come in one of four types, based on
%    the value of the input string algorithm.
%The program first requires 'plane wave sets for common lattices.mat' to be
%   loaded into the workspace
%IF algorithm = 'original' 
%   The lattice is calculated using the original algorithm given in 
%   Plot_Bound_2DLattice_Images_and_OTFs_2021 where the ideal 2D lattice is
%   bound in z, and the binarized version of the real part of the bound
%   lattice E field is written to the SLM.  This creates an E field
%   impinging on the mask that has equal kz bounding for all k vectors.
%   This field is then filtered by an annular mask of NAmax ans NAmin
%   determined by the input values of these variables, and approximately
%   equal to the kz bounding when fill_factor = 1.
%ELSE IF algorithm = 'multi-Bessel':
%   The theoretical pupil field before the mask is determined by replacing each  
%   point of the ideal lattice with a gaussian stripe in kz of sigma = NAsigma
%   large compared to NAmax of the annulus, and then computationally cropped
%   by the annulus to produce a series of stripes of relatively uniform
%   intensity.  The real part of the inverse FT of this then normalized, cropped to set
%   values below crop_factor to zero, and binarized to determine the pattern
%   applied to the SLM.  This is the same procedure followed for multi-Bessel
%   light sheets in Chen et al Science 2014.
%ELSE IF algorithm = "axially confined"
%   The theoretical pupil field before the mask is determined by replacing each  
%   point of the ideal lattice with a gaussian stripe in kz of sigma = NAsigma
%   generally smaller than the width NAmax-NAmin of the annulus, and then 
%   computationally cropped (generally weakly) by the annulus. The real part of
%   the inverse FT of this then normalized, cropped to set values below 
%   crop_factor to zero, and binarized to determine the pattern applied to 
%   the SLM.  This is the same procedure followed for axially confined light
%   sheets in Chen et al Science 2014.
%ELSE IF algorithm = 'cosine-sinc' 
%   The theoretical SLM pattern is calculated in exactly the same way as in
%   the multi-Bessel case.  However, to enforce the conditions of a
%   cosine-sinc light sheet, the pupil field after the annular mask is
%   simulated by replacing each point of the ideal lattice with a uniform stripe
%   in kz large compared to NAmax of the annulus, which is then cropped by
%   the annulus.  Becuase the manipulation is done post-SLM, this case is
%   only used for simulations, and not actual experiments.
%ELSE IF algorithm = 'harmonic balanced'
%   The theoretical pupil field before the mask is determined by first
%   finding the range deltaNAb (Eqs. 41a,b of the paper) of each beamlet b,
%   assuming that NAsigmab for each beamlet is identical and given by the
%   input variable NAsimga.  The actual NAsigmab for the beamlet at the
%   position of highest kz in the pupil (defined as the reference beamlet)
%   is set equal to the input NAsigma, and for all other beamlets NAsigmab
%   = (deltaNAref/deltaNAb)*NAsigma.  This theoretical pupil field is then 
%   computationally cropped (generally weakly) by the annulus. The real part
%   of the inverse FT of this then normalized, cropped to set values below
%   crop_factor to zero, rescaled to the full bit depth of the SLM, and then
%   applied to the SLM. 
%ELSE
%   set algorithm = 'grayscale'
%   The symmetry of the light sheet is assumed to be defined by PW, the 
%   pupil beamlet properties are defined by lightsheet_descrip and NAsigma, 
%   and the real part of the inverse FT of the pupil field is normalized, 
%   cropped to set values below crop_factor to zero, rescaled to the full 
%   bit depth of the SLM, and then applied to the SLM.
%END
%
%After determining the post mask E field, the program calculates
%the overall xz PSFs and OTFs at various positions along the propagation
%direction y to see how the PSFs and the relative amplitudes of the 
%excitation spatial frequencies vary form their initital amplitudes on focus
%
%
%valid choices for lightsheet_descrip are:
%Any desired lattice type (e.g. Hex or Sq45)
%FST: field synthesis theorem.  All point wavevectors of the ideal lattice
%     are replaced with uniform stripes along kz spanning the pupil.  This 
%     pattern then repalces the E field impinging on the annular mask. 
%Ideal Gaussian: DC Gaussian in the pupil creates a Gaussian beam at the focal plane
%ExpGaussian: a laterally offset pair of Gaussian lines in the pupil creates a 
%             Gaussian bound lateral standing wave at the focal plane
%IdealSincz: a DC rect ("top-hat") pupil beam to create a sinc(z) light sheet 
%            cross section
%ExpSincz: a DC rect ("top-hat") pupil beam to create a sinc(z) light sheet
%            cross section
%ExpSincy: same as IdealSLMSincy, but with a pair of latterally offset
%          beams that create a lateral standing wave at the focal plane 
%
%IF lightsheet_descrip = "Gaussian"
%   a.  PW is overwritten with PW = [1 0 0 0 0 0] to let the program know
%       the only k vector is one along the axis of the objective
%   b.  A circular mask of NA = 0.6 is used, with no DC block (since the Gaussian beam is DC) 
%   c.  Enforces an E field impinging on the mask given by a Gaussian
%       function
%END
%IF lightsheet_descrip = "Sinc" (sinc at sample, rect = top-hat in pupil)
%   a.  PW is overwritten with PW = [1 0 0 0 0 0] to let the program know
%       the only k vector is one along the axis of the objective
%   b.  A circular mask of NAmax = NAsigma is used, to block unwanted
%       diffraction orders outside the desired pupil band of uniform intensity
%   c   Enforces an E field impinging on the mask given by a 1D beam in kz 
%       with uniform intensity between +/-NAsinc = NAmax
%END
%
%   PW is an N x 6 matrix containing the directions of the N plane
%      waves that constitute the lattice.  Each row refers to a single 
%      plane wave.  The properties of each plane wave are entered in the
%      columns of PW as follows:
%      PW(n, 1:3) = complex electric field vector.  
%         IMPORTANT:  FOR THIS PROGRAM ONLY, PW(n,1) contains the relative
%         E field amplitude that should be used for each k vector, and 
%         PW(n,2:3) are ignored, since the final E field compoments
%         PW(n,1:3) are calculated within the program itself
%      PW(n, 4:6) = direction ek
%   NAlattice is the NA of the ideal lattice form which the LLS is created
%   NAsigma is the 1/e^2 width of the Gaussian INTENSITY bounding envelope 
%       applied to each plane wave of the LLS before any adjustment 
%       for location within the annulus, meaasured in terms of the NA span 
%       sigma within the rear pupil.
%   NAannulus is a vector giving [NAmax NAmin] of the annular mask
%   crop_factor gives the value of the normalized real(E) at the SLM below
%      which the phase on the SLM should be overwritten with a grating pattern
%      diffract the largely DC signal beyond this level outside the pupil 
%      so that the DC beam block at the annular mask is not overloaded
%
%create the folder where the results will be stored:
% 
% 
%Input Variables:
%    algorithm: see above
%    PW:    an N x 6 matrix containing the directions of the N plane
%           waves that constitute the ideal lattice.  Each row refers to a single 
%           plane wave.  The properties of each plane wave are entered in the
%           columns of PW as follows:
%           PW(n, 1:3) = complex electric field vector.  THIS PART IS IGNORED
%           PW(n, 4:6) = real vector, which in unit form gives propagation
%             direction ek
%           E and k vector dot product must be zero (true for all plane waves)
%           examples for different lattice symmetries are given in 'PW_OTF_PSF.mat'
%    xyPol: a two element row vector giving the complex input electric
%           field, which must be in the xy plane, since the input illumination
%           is along the z axis
%    NAlattice: numerical aperture of the ideal optical lattice upon which
%        the LLS is based, i.e., the numerical aperture at the center of each
%        beamlet of the desired LLS
%    NAsigma is the 1/e^2 width of the Gaussian INTENSITY bounding envelope 
%        applied to each plane wave of the LLS before any adjustment 
%        for location within the annulus, meaasured in terms of the NA span 
%        sigma within the rear pupil.
%    NAannulus is a vector giving [NAmax NAmin] of the annular mask
%    lightsheet_descrip: a string giving the type of light sheet to be
%        calculated:
%        IF lightsheet_descrip = 'IdealGaussian' calculates the
%            properties of the Gaussian light sheetcreated by a pupil field
%            consisting of a 1D Gaussian stripe centered in the pupil and
%            having sigma = NAsigma
%        ELSE IF lightsheet_descrip = 'ExpGaussian' calculates the properties
%            of a Gaussian bound lateral standing wave light sheet created
%            by a pair of laterally offset 1D Gaussian stripes of 
%            sigma = NAsigma. The gap between the stripes allows for a DC
%            beam block to exclude undiffracted light from the SLM.
%        ELSE IF lightsheet_descrip = 'IdealSinc' calculates the
%            properties of the sinc light sheet created by a pupil field
%            consisting of a 1D stripe of uniform amplitude in the
%            pupil with endpoints at kz = +/-NAsigma
%        ELSE IF lightsheet_descrip = 'ExpSinc' calculates the properties
%            of a sinc bound lateral standing wave light sheet created
%            by a pair of laterally offset 1D uniform stripes with endpoints
%            at kz = +/-NAsigma. The gap between the stripes allows for a DC
%            beam block to exclude undiffracted light from the SLM.
%        ELSE IF lightsheet_descrip = 'AxialSW' calculates the properties
%            properties of an axially confined axial standing wave.  If
%            NAlattice + NAsigma > NAmin, the pupil field consists of a 
%            pair of 1D Gaussian stripes of sigma = NAsigma vertically
%            offset by +/- NAlattice.  Otherwise, the pupil field consists
%            of two such pairs of vertically offset Gaussian stripes, with
%            each pair additionally offset laterally by NA = 0.22, and
%            NAmin = 0.2 for the annulus.  The lateral offset allows for a DC
%            beam block to exclude undiffracted light from the SLM.
%        ELSE
%            calculates the properties of a lattice light sheet defined by
%            the wavevectors in PW according to the procedure of the input
%            variable algorithm.  The value input to lightsheet_descrip
%            should reflect the light sheet symmetry (e.g., axial SW,
%            square, hex, or hexrect).
%        END
%    crop_factor determines the level below which the normalized real part 
%        of the calculated E field at the SLM is set to zero.  crop_factor 
%        is >0 and generally <0.2.  A high crop factor leads to a light sheet 
%        more tightly bound in z, at the expense of a shorter propagation length  
%    NAdet is the numerical aperture of the detection objective used in the
%        simulation
%    xz_det_PSF is the detection PSF in the xz plane as given by the variable
%        xz_PSF_RW_510nm_NA1p0 in 'PW_OTF_PSF.mat'
%    xz_det_OTF is the detection OTF in the xz plane as given by the variable
%        xz_OTF_RW_510nm_NA1p0 in 'PW_OTF_PSF.mat'
%    gamma is the factor Iplot = I.^gamma used in plotting the OTF curves
%        to better visualize weak regions within the support
%    yplot_range: wave number plotted in ky
%    ystepsize_in: input wave number interval for the plot 
%    outputdir is the path to the directory to where the simulation results are stored 
%    subfolder: user defined subfolders within outputdir to save the simulation results
%
%Output Variables:
%   yFWHM is an estimate of the +/- extent of the light sheet in the y 
%      direction, given by the 50% intensity points of the sheet
%   OnAxisIntensity is a vector giving the intensity of the dithered lattice 
%      at the z = 0 central plane every 3 lambda along the propagation vector y
%   CumulativeFWHM is a vector giving the full width at half minimum of the 
%      of the cumulative intensity in the dithered LLS with increasing distance
%      from z = 0, calculated every 3 lambda along the propagation vector y
%   PWb is an N x 7 matrix where the first six columns are PW from above,
%      and the 7th column is the fill factor ratios between the different k
%      vectors
%   Overall_Dithered_xz_OTF2: overall dithered xz OTF
%   DitheredxzPSFCrossSection: overall dithered xz PSF cross section
%   PlotSLMPhase: SLM Phase
%   MaskIntensity: Annular mask intensity
%   SampleIntensityAtFocus: the LLS xz excitation PSF at the focal plane of the excitation objective
%   DitheredIntensity: intensity of yz slice at x = 0 through the dithered exc PSF
%   OverallCrossSection: Overall PSF axial linecuts at focus
%   SheetCrossSection: Excitation PSF axial linecuts at focus
%   DetPSFCrossSection: Detection PSF axial linecuts at focus
%   xzSIMOTF5: the xz lattice SIM OTF at at multiple points along the propagation axis
%
% Author: Eric Betzig
% Reorganized by Xiongtao Ruan

if nargin < 16 || isempty(subfolder)
    subfolder = ['NAlattice', num2str(NAlattice, '%1.2f'), filesep, lightsheet_descrip, filesep, 'NAAnnulusMax', num2str(NAannulus(1), '%1.2f'), filesep, 'NAsigma', num2str(NAsigma, '%1.2f'), filesep];
end
% subfolder = ['NAlattice', num2str(NAlattice, '%1.2f'), filesep, lightsheet_descrip, filesep, 'NAAnnulusMax', num2str(NAannulus(1), '%1.2f'), filesep, 'NAsigma', num2str(NAsigma, '%1.2f'), filesep, 'Crop', num2str(crop_factor, '%1.2f')];
mkdir(outputdir, subfolder);
imgpath2 = [outputdir, subfolder];
%
NAmax = NAannulus(1);
NAmin = NAannulus(2);
index = 1.33;
%normalize the NA values to the refractive index:
namax = NAmax./index;
naideal = NAlattice./index;
namin = NAmin./index;
nasigma = NAsigma./index;  %the initial E field bounding is given by %E = Eo*exp(-(kz/(nasigma))^2)
%
%find the plot range for the PSFs and OTFs, and pupil field:
pixsize = 0.1;
plot_range = 100;  %calculate over +/-100 media wavelength FOV in pixel 
                   %steps of 0.1 media wavelengths
halfpix = round(plot_range./pixsize);    
nasigma_pix = halfpix.*nasigma./5;  %Gaussian beamlet envelope in pixels
A = [1 1].*(2.*halfpix + 1);  %image size for PSFs and OTFs 
%
%
%now find the range dNA of illumination associated with the +/-kz edges
%   of each bound k vector:
B = size(PW);
dNA = zeros(1,B(1));
namax_beamlet = dNA;
namin_beamlet = dNA;
for p = 1:B(1)
    namax_beamlet(p) = sqrt(naideal.^2 + 2.*naideal.*nasigma.*abs(PW(p,5)) + nasigma.^2);
    if naideal.*abs(PW(p,5)) >= nasigma  
        %the current illumination stripe does not pass the equator
        namin_beamlet(p) = sqrt(naideal.^2 - 2.*naideal.*nasigma.*abs(PW(p,5)) + nasigma.^2);
    else
        %the current illumination stripe does intersect the equator
        namin_beamlet(p) = naideal;
    end
    dNA(p) = namax_beamlet(p) - namin_beamlet(p);
end
%
%find the minimum dNA among the k vectors, corresponding to the beam with
%   the least bounding and hence largest beam waist at the focus:
[dNAmin, kindex] = min(dNA);
kymax = 2.*pi.*sqrt(1-(max(namin_beamlet(kindex),namin)).^2);
kymin = 2.*pi.*sqrt(1-(min(namax_beamlet(kindex),namax)).^2);
kydiff = kymax - kymin;

%first find the THEORETICAL DESIRED pupil electric field based on the 
%    algorithm and lattice description:
DesiredPupilEField = zeros(A(1),A(2));  %array for the desired pupil E field
B = size(PW);  %number of wavevectors (beamlets) comprising the lattice
switch lightsheet_descrip
    case 'IdealGaussian'
        PW = [1 0 0 0 0 0];  %enforce a single k vector at the center of the objective for Gaussian or sinc
        NAmax = 0.6;
        NAlattice = 0;
        naideal = NAlattice./index;
        NAmin = 0;   %enforces a centrally located 1D Gaussian beam impinigng on 
                     %the mask having a 1/e^2 half-width = NAsigma. This creates
                     %a Gaussian profile at the focal plane.
        crop_factor = 0.01;    %minimal cropping for grayscale SLM
        algorithm = 'grayscale';
        centerx = halfpix + 1;
        centerz = halfpix + 1;
        EAmp = PW(1,1);  %manual adjustment of k vector stored in first index of PW
        EFieldkz = EAmp.*exp(-(((1:1:A(2))-centerz)./nasigma_pix).^2)';
        DesiredPupilEField(:,centerx) = DesiredPupilEField(:,centerx) + EFieldkz;
    case 'IdealSinc'
        PW = [1 0 0 0 0 0];  %enforce a single k vector at the center of the objective for Gaussian or sinc
        NAmax = 0.6;
        NAlattice = 0;
        naideal = NAlattice./index;
        NAmin = 0;   %enforces a centrally located 1D uniform beam of length
                     %+/-NAsigma. This creates a sinc profile at the focal plane.
        crop_factor = 0.01;    %minimal cropping for grayscale SLM
        algorithm = 'grayscale';
        centerx = halfpix + 1;
        centerz = halfpix + 1;
        EAmp = PW(1,1);  %manual adjustment of k vector stored in first index of PW
        EFieldz = zeros(1,A(2));
        EFieldz((centerz-round(nasigma_pix)):(centerz+round(nasigma_pix))) = ones(1,(2.*round(nasigma_pix)+1));
        EFieldz = EAmp.*EFieldz';
        DesiredPupilEField(:,halfpix+1) = DesiredPupilEField(:,halfpix+1) + EFieldz;        
    case 'ExpGaussian'        
        %to avoid the undiffracted light near DC, the experimental Gaussian at the
        %focal plane is created by two Gaussians offset 
        %by +/-NAlattice > NAmin in the pupil plane:
        PW = [1 0 0 1 0 0; 1 0 0 -1 0 0]; 
        NAmin = max(NAmin,0.2);
        naideal = NAlattice./index;
        NAmin = 0.2;  %blocks the DC beam
        crop_factor = 0.01;  %minimal cropping for grayscale SLM
        algorithm = 'grayscale';
        %use a pair of 1D Gaussian pupil beams to create a lateral standing
        %wave at the focal plane axially bound by a Gaussian:
        %first add the left beam:
        centerx = halfpix + 1 + round(halfpix.*naideal.*PW(1,4)./5);
        centerz = halfpix + 1 + round(halfpix.*naideal.*PW(1,5)./5);
        EAmp = PW(1,1);  %manual adjustment of k vector stored in first index of PW
        EFieldz = EAmp.*exp(-(((1:A(1))-centerz)./nasigma_pix).^2)';        
        DesiredPupilEField(:,centerx) = DesiredPupilEField(:,centerx) + EFieldz; 
        %now add the right beam:
        centerx = halfpix + 1 + round(halfpix.*naideal.*PW(2,4)./5);
        centerz = halfpix + 1 + round(halfpix.*naideal.*PW(2,5)./5);
        EAmp = PW(2,1);  %manual adjustment of k vector stored in first index of PW
        EFieldz = EAmp.*exp(-(((1:A(1))-centerz)./nasigma_pix).^2)';        
        DesiredPupilEField(:,centerx) = DesiredPupilEField(:,centerx) + EFieldz;
    case 'ExpSinc'
        %to avoid the undiffracted light near DC, the experimental sinc beam at the
        %focal plane is created by two rect beams offset by +/-NAlattice > NAmin
        %in the pupil plane:
        PW = [1 0 0 1 0 0; 1 0 0 -1 0 0]; 
        NAmax = 0.6;  %max OD of annulus to permit Gaussian tails to be trnamitted
                      %as much as possible
        NAmin = max(NAmin,0.2);
        naideal = NAlattice./index;
        NAmin = 0.2;  %blocks the DC beam
        crop_factor = 0.01;  %minimal cropping for grayscale SLM
        algorithm = 'grayscale';
        %add the first rect beam:
        centerx = halfpix + 1 + round(halfpix.*naideal.*PW(1,4)./5);
        centerz = halfpix + 1 + round(halfpix.*naideal.*PW(1,5)./5);
        EAmp = PW(1,1);  %manual adjustment of k vector stored in first index of PW
        EFieldz = zeros(1,A(2));
        EFieldz((centerz-round(nasigma_pix)):(centerz+round(nasigma_pix))) = ones(1,(2.*round(nasigma_pix)+1));
        EFieldz = EAmp.*EFieldz';
        DesiredPupilEField(:,centerx) = DesiredPupilEField(:,centerx) + EFieldz;        
        %now add the right beam:
        centerx = halfpix + 1 + round(halfpix.*naideal.*PW(2,4)./5);
        centerz = halfpix + 1 + round(halfpix.*naideal.*PW(2,5)./5);
        EAmp = PW(2,1);  %manual adjustment of k vector stored in first index of PW
        EFieldz = zeros(1,A(2));
        EFieldz((centerz-round(nasigma_pix)):(centerz+round(nasigma_pix))) = ones(1,(2.*round(nasigma_pix)+1));
        EFieldz = EAmp.*EFieldz';
        DesiredPupilEField(:,centerx) = DesiredPupilEField(:,centerx) + EFieldz;                 
    case 'ExpAxialSW'
        % crop_factor = 0.01;  %minimal cropping for grayscale SLM
        algorithm = 'original';
        if (NAlattice - NAsigma) > NAmin
            %the beamlets are not clipped by the DC part of the annulus, so
            %they can be placed at the lateral center of the pupil:
            PW = [1 0 0 0 1 0; 1 0 0 0 -1 0];
        else
           %use two pairs of laterally and axially offset beamlets to avoid the
           %DC beam block:
           NAmin = max(NAmin,0.2);
           NALateralOffset = 0.22; 
           xna = NALateralOffset;
           yna = NAlattice;
           PW = [1 0 0 xna yna 0; 1 0 0 -xna yna 0; 1 0 0 xna -yna 0; 1 0 0 -xna -yna 0]; 
           PW(:,4:5) = PW(:,4:5)./sqrt(xna.^2 + yna.^2);
           % crop_factor = 0.01;
           lightsheet_descrip = 'ExpAxialSW';
        end
        naideal = NAlattice./index;
    otherwise
end

switch algorithm
    case 'multi-Bessel'
        %add an infinite stripe in kz at the kx location for each beamlet:
        for p = 1:B(1)  %loop thru all beamlets
          %find the kx location of  vector:
          centerx = halfpix + 1 + round(halfpix.*naideal.*PW(p,4)./5);
          EAmp = PW(p,1);  %manual adjustment of k vector stored in first index of PW
          EFieldkz = EAmp.*ones(A(1),1);
          % xruan update DesiredPupilEField
          % DesiredPupilEField(:,centerx) = DesiredPupilEField(:,centerx) + EFieldkz;
          DesiredPupilEField(:,centerx) = EFieldkz;
        end
    case 'cosine-sinc'
        %add an infinite stripe in kz at the kx location for each beamlet:
        for p = 1:B(1)  %loop thru all beamlets
          %find the kx location of  vector:
          centerx = halfpix + 1 + round(halfpix.*naideal.*PW(p,4)./5);
          EAmp = PW(p,1);  %manual adjustment of k vector stored in first index of PW
          EFieldkz = EAmp.*ones(A(1),1);
          DesiredPupilEField(:,centerx) = EFieldkz;
        end                
    case 'harmonic balanced'
        %find the range dNA of illumination associated with the +/-kz edges
        %of each bound k vector:
        dNA = zeros(1,B(1));
        namax_beamlet = dNA;
        namin_beamlet = dNA;
        for p = 1:B(1)
            namax_beamlet(p) = sqrt(naideal.^2 + 2.*naideal.*nasigma.*abs(PW(p,5)) + nasigma.^2);
            if naideal.*abs(PW(p,5)) >= nasigma  
                %the current illumination stripe does not pass the equator
                namin_beamlet(p) = sqrt(naideal.^2 - 2.*naideal.*nasigma.*abs(PW(p,5)) + nasigma.^2);
            else
                %the current illumination stripe does intersect the equator
                namin_beamlet(p) = naideal;
            end
            dNA(p) = namax_beamlet(p) - namin_beamlet(p);
        end
        %
        %find the minimum dNA among the k vectors, corresponding to the beam with
        %the least bounding and hence largest beam waist at the focus:
        [dNAmin, kindex] = min(dNA);
        %
        %adjust NAsigma for every k vector so that every beam has this same minimum
        %dNA and hence the same beam waist at the focus:
        nasigma_adjusted = nasigma.*dNAmin./dNA;
        %
        %calc the rear pupil E field for this condition where the bunding of the 
        %k vectors have been adjusted so that they all have the same dNA and
        %hence the same beam waist size at the focus:
        %express kzsigma_adjusted in terms of pixels:
        nasigma_adjusted_pix = halfpix.*nasigma_adjusted./5;
        B = size(PW);
        for p = 1:B(1)
          %find the ideal center location of each k vector:
          centerx = halfpix + 1 + round(halfpix.*naideal.*PW(p,4)./5);
          centerz = halfpix + 1 + round(halfpix.*naideal.*PW(p,5)./5);
          EAmp = PW(p,1);  %manual adjustment of k vector stored in first index of PW
          EFieldkz = EAmp.*exp(-(((1:A(1))-centerz)./nasigma_adjusted_pix(p)).^2)';
          %adjust the Efield amp so that all beams have the same Efield amp at the
          %focus:
          EFieldkz = (dNA(p)./dNAmin).*EFieldkz;
          DesiredPupilEField(:,centerx) = DesiredPupilEField(:,centerx) + EFieldkz;
        end
    case 'axially confined'
        %algorithm = 'axially confined'
        %load at Gaussian stripe in kz at the center position of each beamlet:
        for p = 1:B(1)  %loop thru all beamlets
          %find the ideal center location of each k vector:
          centerx = halfpix + 1 + round(halfpix.*naideal.*PW(p,4)./5);
          centerz = halfpix + 1 + round(halfpix.*naideal.*PW(p,5)./5);
          EAmp = PW(p,1);  %manual adjustment of k vector stored in first index of PW
          EFieldkz = EAmp.*exp(-(((1:1:A(2))-centerz)./nasigma_pix).^2)';
          DesiredPupilEField(:,centerx) = DesiredPupilEField(:,centerx) + EFieldkz;
        end
    case 'original'
        %
        %insure that the k vectors lie in a plane (i.e., cone angle = 90 deg):
        PW(:,4:5) = PW(:,4:5) ./ sqrt(1 - PW(1,6)^2);
        PW(:,6) = 0;
        %
        %now modify each k vector to reflect the cone angle upon which they lie:
        PW(:,4:5) = naideal .* PW(:,4:5);
        PW(:,6) = sqrt(1-naideal.^2);
        %
        %normalize the input polarization and make it a 1 x 3 vector:
        InputPol = [xyPol(1) xyPol(2) 0];
        InputPol = InputPol ./ norm(InputPol);
        %
        %find the electric field for each lattice wavevector after passing through
        %the objective:
        B = size(PW);
        for n = 1:B(1)
            %find the orthonormal vectors defining the coordinate system for the
            %nth beam when entering the rear pupil:
            phivec = cross(squeeze(PW(n,4:6)), [0 0 1]);
            phivec = phivec ./ norm(phivec);  %azimuthal unit vector
            radvec = cross([0 0 1], phivec);  %radial unit vector
            %
            %the azimuthal component of the electric field is unaffected when passing
            %through the objective:
            ephi = dot(phivec, InputPol);
            %
            %the radial component is tilted by refraction when passing through the
            %objective to point in the theta direction as defined by a spherical
            %coordinate system centered at the focal point:
            thetavec = cross(squeeze(PW(n,4:6)), phivec);
            etheta = dot(radvec, InputPol);
            %
            %save the desired electric field amplitude stored in PW(n,1) before 
            %loading the complexelectric field into PW(n,1:3)
            EAmp = PW(n,1);
            %the sum of the azimuthal and theta components gives the total electric
            %field for the nth plane wave of the lattice:
            PW(n,1:3) = ephi .* phivec + etheta .* thetavec;
            %
            %confirm that the electric field is of unit strength:
            PW(n,1:3) = PW(n,1:3) ./ norm(PW(n,1:3));
            %adjust the eletric field to the desired amplitude:
            PW(n,1:3) = EAmp.*PW(n,1:3);
        end
    %
        %calculate the complete electric field of the ideal 2D lattice:
        x = 0:pixsize:plot_range;
        y = 0:pixsize:plot_range;
        [X Y] = ndgrid(x, y);
        %now calculate the E field at each mesh point:
        A = size(X);
        Ex = zeros(A(1), A(2));
        Ey = zeros(A(1), A(2));
        Ez = zeros(A(1), A(2));
        for q = 1:B(1)   %loop thru all plane waves
            phase = exp(2 .* pi .* i .* (PW(q, 4) .* X + PW(q, 5) .* Y));
            Ex = Ex + PW(q, 1) .* phase;
            Ey = Ey + PW(q, 2) .* phase;
            Ez = Ez + PW(q, 3) .* phase;
        end
        %expand through all quadrants:
        Extemp = zeros(2 .* A(1), 2 .* A(2));
        Eytemp = zeros(2 .* A(1), 2 .* A(2));
        Eztemp = zeros(2 .* A(1), 2 .* A(2));
        %load the original data into the first quadrant:
        Extemp((A(1) + 1):end, (A(2) + 1):end) = Ex;
        Eytemp((A(1) + 1):end, (A(2) + 1):end) = Ey;
        Eztemp((A(1) + 1):end, (A(2) + 1):end) = Ez;
        %now mirror along each dimension and use parities to fill other quadrants:
        %simply mirror the data since parity is always even for magnitudes:
        Extemp(1:A(1), (A(2) + 1):end) = flipdim(Ex, 1);
        Eytemp(1:A(1), (A(2) + 1):end) = flipdim(Ey, 1);
        Eztemp(1:A(1), (A(2) + 1):end) = flipdim(Ez, 1);
        Extemp(:, 1:A(2)) = flipdim(Extemp(:, (A(2) + 1):end), 2);
        Eytemp(:, 1:A(2)) = flipdim(Eytemp(:, (A(2) + 1):end), 2);
        Eztemp(:, 1:A(2)) = flipdim(Eztemp(:, (A(2) + 1):end), 2);
        %delete the extra vector from mirroring in each dimension:
        Extemp(A(1), :) = [];
        Eytemp(A(1), :) = [];
        Eztemp(A(1), :) = [];
        Extemp(:, A(2)) = [];
        Eytemp(:, A(2)) = [];
        Eztemp(:, A(2)) = [];
        Ex = Extemp;
        Ey = Eytemp;
        Ez = Eztemp;
        clear Extemp Eytemp Eztemp;
        A = size(Ex);
        %
        %find the ideal 2D lattice intensity:  %%%%%%
        EComp = zeros(3,A(1),A(2));
        EComp(1, :, :) = Ex;
        EComp(2, :, :) = Ey;
        EComp(3, :, :) = Ez;
        ESq = EComp .* conj(EComp);
        ESqTot = squeeze(sum(ESq, 1));
        maxval = max(max(ESqTot));
        ESqTot = ESqTot ./ maxval;
        %
        %calc and plot the ideal 2D lattice of the component of the real electric field 
        %projected onto the state of the input electric field polarization:
        RealE = real(conj(Ex) .* InputPol(1) + conj(Ey) .* InputPol(2));
        %                                  
        %calc and plot the bound 2D electric field lattice at the SLM:
        z = -plot_range:pixsize:plot_range;
        kxdiff = 2.*pi*nasigma;
        lattice_full_width = pi./kxdiff;  %approximate half width of the function limiting
            %the extent of the bound lattice, in media wavelengths
        sigma = lattice_full_width ./ sqrt(2.*log(2));
        envelope = exp(-2.*(z./sigma).^2)' * ones(1,A(1));
        RealE = RealE .* envelope';
        RealE = RealE';     
end        
%
%calc the rear pupil E field for this nominal condition where all k vectors 
%   have the same bounding:
NominalPupilEfield = zeros(2001,2001); %match the pupil field of other 
%   programs that use pixsize = 0.1 media wavelengths and +/-100 media 
%   wavelgnths field of view. In this case, the pupil field calculation 
%   runs from +/-2.5 times the max spatial frequency defined by two
%   counterpropagating plane waves
%calc and plot the function for the filtering provided by the annular mask:
halfpix = round((A(1)-1)./2);
x = -halfpix:halfpix;
x = 5.*x./halfpix;   %match the spectral range of the electric field at the annular mask
y = x;
[X Y] = ndgrid(x, y);
R = sqrt(X.*X + Y.*Y);
MaxRad = namax;
MinRad = namin;
AnnularFilter = (R <= MaxRad) & (R >= MinRad);
switch lightsheet_descrip
    case {'ExpGaussian', 'ExpSinc', 'ExpAxialSW'}
    %   %block the spurious center stripe that results from diffraction from
    %   %the SLM:
       AnnularFilter((halfpix-10):(halfpix+10),:) = 0;
end
figure  %create a new figure window for the plots
set(gcf, 'Position', [50 100 600 600]);
axes_h = axes('Position', [0.13, 0.12, 0.8, 0.8]);
%crop the result to the +/- the max spatial freq defined by
%conterpropagating plane waves:
AnnularFilter2 = AnnularFilter((0.6.*halfpix):(1.4.*halfpix),(0.6.*halfpix):(1.4.*halfpix));
image(256 .* AnnularFilter2');
colormap hot(256);
B = size(AnnularFilter2);
axis([1 B(1) 1 B(2)]);
axis square;
set(gca, 'XTick', [1:(B(1)-1)./4:B(1)]);
set(gca, 'XTickLabel', [-1:0.5:1], 'FontSize', 12);
xlabel(['k_{x} / (4\pi/\lambda)'], 'FontSize', 14);
set(gca, 'YTick', [1:(B(2)-1)./4:B(2)]);
set(gca, 'YTickLabel', [-1:0.5:1], 'FontSize', 12);
ylabel(['k_{z} / (4\pi/\lambda)'], 'FontSize', 14);
text(0.2 .*B(1), -0.04 .* B(2), ['Annular Mask Transmission, NA ', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')], 'FontSize', 12);
%
%
% apply the mask to the current DesiredPupilEField to find its final form:
switch algorithm
    case {'multi-Bessel', 'harmonic balanced'}
        DesiredPupilEField = DesiredPupilEField.*(AnnularFilter');
end
% calc and plot the DesiredPupilIntensity:
DesiredPupilIntensity = DesiredPupilEField.*conj(DesiredPupilEField);
DesiredPupilIntensity = DesiredPupilIntensity./max(max(DesiredPupilIntensity));
figure  %create a new figure window for the plots
set(gcf, 'Position', [10 -200 1000 1000]);
axes_h = axes('Position', [0.13, 0.12, 0.8, 0.8]);
%crop the result to the +/- the max spatial freq defined by conterpropagating 
%plane waves:
DesiredPupilIntensity2 = DesiredPupilIntensity((0.6.*halfpix):(1.4.*halfpix),(0.6.*halfpix):(1.4.*halfpix));
image(256 .* (DesiredPupilIntensity2.^gamma));
colormap hot(256);
C = size(DesiredPupilIntensity2);
axis([1 C(1) 1 C(2)]);
axis square;
set(gca, 'XTick', [1:(C(1)-1)./4:C(1)]);
set(gca, 'XTickLabel', [-1:0.5:1], 'FontSize', 16);
xlabel(['k_{x} / (4\pi/\lambda)'], 'FontSize', 16);
set(gca, 'YTick', [1:(C(2)-1)./4:C(2)]);
set(gca, 'YTickLabel', [-1:0.5:1], 'FontSize', 16);
ylabel(['k_{z} / (4\pi/\lambda)'], 'FontSize', 16);
text(0.1 .*C(1), -0.06 .* C(2), ['Desired Pupil Intensity, ', algorithm, ' ', lightsheet_descrip, ', gamma = ', num2str(gamma, '%1.2f')], 'FontSize', 15);
lightsheet_descrip2 = ['lattice NA = ', num2str(NAlattice, '%1.2f'), ', annulus NA', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')];
if (strcmp(algorithm, 'grayscale') || strcmp(algorithm, 'axially confined')) || strcmp(algorithm, 'harmonic balanced')
    lightsheet_descrip2 = ['lattice NA = ', num2str(NAlattice, '%1.2f'), ', NAsigma = ', num2str(NAsigma, '%1.2f'), ', annulus NA = ', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')];
end
text(0.2 .*C(1), -0.02 .* C(2), lightsheet_descrip2, 'FontSize', 15);
%text(0.22 .*C(1), -0.06 .* C(2), ['Desired Pupil Intensity, annulus NA', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')], 'FontSize', 16);
%text(0.08 .*C(1), -0.02 .* C(2), [lightsheet_descrip, ', lattice NA = ', num2str(NAlattice, '%1.2f'), ', NAsigma = ', num2str(NAsigma, '%1.2f'), ', gamma = ', num2str(gamma, '%1.2f')], 'FontSize', 16);
saveas(gcf ,[imgpath2 'Desired Pupil Intensity.png']);    %
%
%
%project the desired pupil E field aback to the SLM:
if ~strcmp(algorithm, 'original')
    ESLM = fftshift(ifft2(ifftshift(DesiredPupilEField)));  %complex electric field impinging on the annular mask    
    RealE = real(ESLM);
end
maxval = max(max(RealE));
minval = min(min(RealE));
PlotRealE = (RealE - minval) ./ (maxval - minval);
figure  %create a new figure window for the plots
set(gcf, 'Position', [10 -200 1000 1000]);
axes_h = axes('Position', [0.13, 0.1, 0.8, 0.8]);
image(256 .* PlotRealE);
colormap hot(256);
A = size(PlotRealE);
axis([1 A(1) 1 A(2)]);
axis square;
set(gca, 'XTick', [1:(A(1)-1)./4:A(1)]);
set(gca, 'XTickLabel', plot_range .* [-1:0.5:1], 'FontSize', 16);
xlabel(['x / \lambda'], 'FontSize', 16);
set(gca, 'YTick', [1:(A(1)-1)./4:A(1)]);
set(gca, 'YTickLabel', plot_range .* [-1:0.5:1], 'FontSize', 16);
ylabel(['z / \lambda'], 'FontSize', 16);
text(0.10 .*A(1), -0.06 .* A(2), ['Bounded E Field Amplitude at SLM, annulus NA', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')], 'FontSize', 16);
text(0.25 .*A(1), -0.02 .* A(2), [lightsheet_descrip, ', lattice NA = ', num2str(NAlattice, '%1.2f'), ', NAsigma = ', num2str(NAsigma, '%1.2f')], 'FontSize', 16);
saveas(gcf ,[imgpath2 'Bounded E Field Amplitude at SLM.png']);    %
%
%
%calc and plot the SLM pattern:
if strcmp(algorithm, 'multi-Bessel') || strcmp(algorithm, 'axially confined') || strcmp(algorithm, 'original') || strcmp(lightsheet_descrip, 'AxialSW') || strcmp(lightsheet_descrip, 'ExpAxialSW')
    %SLM is binary
    %calc and plot the binary phase pattern at SLM 
    %truncate any values below the level of crop_factor:
    RealE = RealE ./ (max(max(RealE)));
    RealE = RealE .* (abs(RealE) > crop_factor);
    %now create the binary phase function (0 or pi) for the SLM:
    SLMPhase = (sign(RealE' + eps).*pi./2 + pi./2);
    %SLMPhase = mod(RealE', 2.*pi) - pi;  %%%
    PWb = PW;
else 
    %SLM is grayscale
    RealE = RealE ./ (max(max(abs(RealE))));
    RealE = RealE .* (abs(RealE) > crop_factor); 
    SLMPhase = RealE'.*pi;  
    %SLMPhase = mod(RealE', 2.*pi) - pi;  %%%
    if strcmp(lightsheet_descrip, 'ExpGaussian') || strcmp(lightsheet_descrip, 'ExpSincy')      
       SLMPhase = mod(RealE', 2.*pi) - pi;
    end    
    PWb = PW;     
end
%plot the SLM pattern:
maxval = max(max(SLMPhase));
minval = min(min(SLMPhase));
PlotSLMPhase = (SLMPhase' - minval) ./ (maxval - minval);
figure  %create a new figure window for the plots
set(gcf, 'Position', [10 -200 1000 1000]);
%set(gcf, 'Position', [850 100 600 600]);
axes_h = axes('Position', [0.13, 0.1, 0.8, 0.8]);
image(256 .* PlotSLMPhase);
colormap hot(256);
axis([1 A(1) 1 A(2)]);
axis square;
set(gca, 'XTick', [1:(A(1)-1)./4:A(1)], 'FontSize', 16);
set(gca, 'XTickLabel', plot_range .* [-1:0.5:1]);
xlabel(['x / \lambda'], 'FontSize', 16);
set(gca, 'YTick', [1:(A(1)-1)./4:A(1)], 'FontSize', 16);
set(gca, 'YTickLabel', plot_range .* [-1:0.5:1]);
ylabel(['z / \lambda'], 'FontSize', 16);
text(0.22 .*A(1), -0.06 .* A(2), ['Phase at SLM, annulus NA', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')], 'FontSize', 16);
text(0.20 .*A(1), -0.02 .* A(2), [lightsheet_descrip, ', lattice NA = ', num2str(NAlattice, '%1.2f'), ', NAsigma = ', num2str(NAsigma, '%1.2f')], 'FontSize', 16);
saveas(gcf ,[imgpath2 'Phase at SLM.png']); 
%
%plot the intensity impinging on the annular mask:
if strcmp(algorithm, 'cosine-sinc')
    %for cosine-sinc, we override the pattern from the SLM to enforce a 
    %pupil field where all ideal point wavevectors are replaced with
    %lines spanning the pupil:
    B = size(PW);
    EImpingingOnMask = zeros(A(1),A(2));
    for p = 1:B(1)
      %find the ideal center location of each k vector:
      centerx = halfpix + 1 + round(halfpix.*naideal.*PW(p,4)./5);
      EAmp = PW(p,1);  %manual adjustment of k vector stored in first index of PW
      EFieldkz = ones(1,A(1));
      EImpingingOnMask(centerx,:) = EFieldkz;
    end 
else
    ESLM = exp(1i .* SLMPhase);
    EImpingingOnMask = fftshift(fft2(ifftshift(ESLM)));  %complex electric 
                                      %field impinging on the annular mask
    %there is a singularity in EImpiningOnMask at center pixel 501,501 which we remove by
    %   replacing it with the value at the adjacent pixel 501,500:
    EImpingingOnMask(halfpix+1,halfpix+1) = EImpingingOnMask(halfpix+1,halfpix);
end
IntensityOnMask = EImpingingOnMask .* conj(EImpingingOnMask);   %intensity impinging on the annular mask
maxval = max(max(IntensityOnMask));
IntensityOnMask = IntensityOnMask ./ maxval;
%crop the result to the +/- the max spatial freq defined by
%conterpropagating plane waves:
IntensityOnMask = IntensityOnMask((0.6.*halfpix):(1.4.*halfpix),(0.6.*halfpix):(1.4.*halfpix));
figure  %create a new figure window for the plots
%set(gcf, 'Position', [50 100 600 600]);
set(gcf, 'Position', [10 -200 1000 1000]);
axes_h = axes('Position', [0.13, 0.12, 0.8, 0.8]);
image(256 .* (IntensityOnMask'.^gamma));
colormap hot(256);
B = size(IntensityOnMask);
axis([1 B(1) 1 B(2)]);
axis square;
set(gca, 'XTick', [1:(B(1)-1)./4:B(1)]);
set(gca, 'XTickLabel', [-1:0.5:1], 'FontSize', 16);
xlabel(['k_{x} / (4\pi/\lambda)'], 'FontSize', 16);
set(gca, 'YTick', [1:(B(2)-1)./4:B(2)]);
set(gca, 'YTickLabel', [-1:0.5:1], 'FontSize', 16);
ylabel(['k_{z} / (4\pi/\lambda)'], 'FontSize', 16);
text(0.13 .*B(1), -0.06 .* B(2), ['Intensity Impinging on Mask, annulus NA', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')], 'FontSize', 16);
text(0.10 .*B(1), -0.02 .* B(2), [lightsheet_descrip, ', lattice NA = ', num2str(NAlattice, '%1.2f'), ', NAsigma = ', num2str(NAsigma, '%1.2f'), ', gamma = ', num2str(gamma, '%1.2f')], 'FontSize', 16);
saveas(gcf ,[imgpath2 'Intensity Impinging on Mask.png']);
%
%
%the original and new algorithms for lattice generation merge here, except
%    for the annulus filtering
maxval = max(max(SLMPhase));
minval = min(min(SLMPhase));
PlotSLMPhase = (SLMPhase' - minval) ./ (maxval - minval);
figure  %create a new figure window for the plots
set(gcf, 'Position', [10 -200 1000 1000]);
%set(gcf, 'Position', [850 100 600 600]);
axes_h = axes('Position', [0.13, 0.1, 0.8, 0.8]);
image(256 .* PlotSLMPhase);
colormap hot(256);
axis([1 A(1) 1 A(2)]);
axis square;
set(gca, 'XTick', [1:(A(1)-1)./4:A(1)], 'FontSize', 16);
set(gca, 'XTickLabel', plot_range .* [-1:0.5:1]);
xlabel(['x / \lambda'], 'FontSize', 16);
set(gca, 'YTick', [1:(A(1)-1)./4:A(1)], 'FontSize', 16);
set(gca, 'YTickLabel', plot_range .* [-1:0.5:1]);
ylabel(['z / \lambda'], 'FontSize', 16);
text(0.22 .*A(1), -0.06 .* A(2), ['Phase at SLM, annulus NA', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')], 'FontSize', 16);
text(0.20 .*A(1), -0.02 .* A(2), [lightsheet_descrip, ', lattice NA = ', num2str(NAlattice, '%1.2f'), ', NAsigma = ', num2str(NAsigma, '%1.2f')], 'FontSize', 16);
saveas(gcf ,[imgpath2 'Phase at SLM.png']); 
%
%now find the pupil field at the mask created by this SLM pattern:
if strcmp(lightsheet_descrip, 'IdealGaussian')
    nasigma_pix = halfpix.*nasigma./5;
    EImpingingOnMask = zeros(A(1),A(2));
    EImpingingOnMask(halfpix+1,:) = exp(-((-halfpix:1:halfpix)./nasigma_pix).^2);
elseif strcmp(lightsheet_descrip, 'FST') 
    %for field synthesis theorem, replace all ideal point wavevectors with
    %lines spanning the pupil:
    B = size(PW);
    EImpingingOnMask = zeros(A(1),A(2));
    for p = 1:B(1)
      %find the ideal center location of each k vector:
      centerx = halfpix + 1 + round(halfpix.*naideal.*PW(p,4)./5);
      EAmp = PW(p,1);  %manual adjustment of k vector stored in first index of PW
      EFieldkz = ones(1,A(1));
      EImpingingOnMask(centerx,:) = EFieldkz;
    end 
else
    ESLM = exp(1i .* SLMPhase);
    EImpingingOnMask = fftshift(fft2(ifftshift(ESLM)));  %complex electric field impinging on the annular mask
    %there is a singularity in EImpiningOnMask at center pixel 501,501 which we remove by
    %   replacing it with the value at the adjacent pixel 501,500:
    EImpingingOnMask(halfpix+1,halfpix+1) = EImpingingOnMask(halfpix+1,halfpix);
end
%plot the intensity impinging on the annular mask:
IntensityOnMask = EImpingingOnMask .* conj(EImpingingOnMask);   %intensity impinging on the annular mask
maxval = max(max(IntensityOnMask));
IntensityOnMask = IntensityOnMask ./ maxval;
%crop the result to the +/- the max spatial freq defined by
%conterpropagating plane waves:
IntensityOnMask = IntensityOnMask((0.6.*halfpix):(1.4.*halfpix),(0.6.*halfpix):(1.4.*halfpix));
figure  %create a new figure window for the plots
%set(gcf, 'Position', [50 100 600 600]);
set(gcf, 'Position', [10 -200 1000 1000]);
axes_h = axes('Position', [0.13, 0.12, 0.8, 0.8]);
image(256 .* (IntensityOnMask'.^gamma));
colormap hot(256);
B = size(IntensityOnMask);
axis([1 B(1) 1 B(2)]);
axis square;
set(gca, 'XTick', [1:(B(1)-1)./4:B(1)]);
set(gca, 'XTickLabel', [-1:0.5:1], 'FontSize', 16);
xlabel(['k_{x} / (4\pi/\lambda)'], 'FontSize', 16);
set(gca, 'YTick', [1:(B(2)-1)./4:B(2)]);
set(gca, 'YTickLabel', [-1:0.5:1], 'FontSize', 16);
ylabel(['k_{z} / (4\pi/\lambda)'], 'FontSize', 16);
text(0.13 .*B(1), -0.06 .* B(2), ['Intensity Impinging on Mask, annulus NA', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')], 'FontSize', 16);
text(0.10 .*B(1), -0.02 .* B(2), [lightsheet_descrip, ', lattice NA = ', num2str(NAlattice, '%1.2f'), ', NAsigma = ', num2str(NAsigma, '%1.2f'), ', gamma = ', num2str(gamma, '%1.2f')], 'FontSize', 16);
saveas(gcf ,[imgpath2 'Intensity Impinging on Mask.png']);    %
%
%calc the E field immediately after transmission through the annulus:
EAfterMask = EImpingingOnMask .* AnnularFilter;
%shift EAfterMask by one pixel to correctly center it within the pupil:
EAfterMask(:,1:(A(2)-1)) = EAfterMask(:,2:A(2));
%plot the intensity immediately after the annular mask:
MaskIntensity = EAfterMask .* conj(EAfterMask);
MaskIntensity = MaskIntensity ./ max(max(MaskIntensity));
MaskIntensity = MaskIntensity((0.6.*halfpix):(1.4.*halfpix),(0.6.*halfpix):(1.4.*halfpix));
figure  %create a new figure window for the plots
set(gcf, 'Position', [10 -200 1000 1000]);
%set(gcf, 'Position', [50 100 600 600]);
axes_h = axes('Position', [0.13, 0.12, 0.8, 0.8]);
image(256 .* (MaskIntensity'.^gamma));
colormap hot(256);
B = size(MaskIntensity);
axis([1 B(1) 1 B(2)]);
axis square;
set(gca, 'XTick', [1:(B(1)-1)./4:B(1)]);
set(gca, 'XTickLabel', [-1:0.5:1], 'FontSize', 16);
xlabel(['k_{x} / (4\pi/\lambda)'], 'FontSize', 16);
set(gca, 'YTick', [1:(B(2)-1)./4:B(2)]);
set(gca, 'YTickLabel', [-1:0.5:1], 'FontSize', 16);
ylabel(['k_{z} / (4\pi/\lambda)'], 'FontSize', 16);
text(0.10 .*B(1), -0.06 .* B(2), ['Post-Mask Intensity at Rear Pupil, annulus NA', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')], 'FontSize', 16);
text(0.10 .*B(1), -0.02 .* B(2), [lightsheet_descrip, ', lattice NA = ', num2str(NAlattice, '%1.2f'), ', NAsigma = ', num2str(NAsigma, '%1.2f'), ', gamma = ', num2str(gamma, '%1.2f')], 'FontSize', 16);
saveas(gcf ,[imgpath2 'Post-Mask Intensity at Rear Pupil.png']);    %
%
%calc and plot the LLS xz excitation PSF at the focal plane of the excitation objective:
ESample = fftshift(ifft2(ifftshift(EAfterMask)));
SampleIntensity = ESample .* conj(ESample);
SampleIntensity = SampleIntensity ./ max(max(SampleIntensity));
SampleIntensityAtFocus = SampleIntensity;
figure  %create a new figure window for the plots
set(gcf, 'Position', [850 100 600 600]);
axes_h = axes('Position', [0.13, 0.1, 0.8, 0.8]);
image(256 .* SampleIntensity');
colormap hot(256);
axis([1 A(1) 1 A(2)]);
axis square;
set(gca, 'XTick', [1:(A(1)-1)./4:A(1)]);
set(gca, 'XTickLabel', plot_range .* [-1:0.5:1]);
xlabel(['x / \lambda'], 'FontSize', 14);
set(gca, 'YTick', [1:(A(1)-1)./4:A(1)]);
set(gca, 'YTickLabel', plot_range .* [-1:0.5:1]);
ylabel(['z / \lambda'], 'FontSize', 14);
text(0.08 .*A(1), -0.06 .* A(2), ['Stationary LLS xz PSF at Focus, annulus NA', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')], 'FontSize', 12);
text(0.10 .*A(1), -0.02 .* A(2), [lightsheet_descrip, ', lattice NA = ', num2str(NAlattice, '%1.2f'), ', NAsigma = ', num2str(NAsigma, '%1.2f'), ', gamma = ', num2str(gamma, '%1.2f')], 'FontSize', 12);
saveas(gcf ,[imgpath2 'Stationary LLS xz PSF at Focus.png']);    %
%
if strcmp(lightsheet_descrip,'ideal')
    yFWHM = [];
    OnAxisIntensity = [];
    CumulativeFWHM = [];
    OverallDitheredPSF1 = [];
    Overall_Dithered_xz_OTF2 = [];
    DitheredxzPSFCrossSection = [];
    DitheredIntensity = [];
    OverallCrossSection = [];
    SheetCrossSection = [];
    DetPSFCrossSection = [];
    return
end
%now calculate the LLS xz excitation pattern and overall OTF at points y > 0
%   along the propagation direction
%to speed the 2D integration over all points in the rear pupil, create a
%   list of all points in the pupil at which the amplitude of the
%   Efield is above some minimum threshold:
EAfterMask = EAfterMask';
EAmpAfterMask = abs(EAfterMask);
EAmpAfterMask = EAmpAfterMask./max(max(EAmpAfterMask));
EAmpThreshhold = 0.05;
% 
%loop through all rows and columns of EAmpAfterMask within the square
%   bounding the annulus:
halfrange = ceil(namax.*halfpix./2);
NumPointsWithinBoundingSquare = (2.*halfrange + 1).^2;
%define matrices to store the info about the illuminated points in the pupil:
NAx = zeros(1,NumPointsWithinBoundingSquare);
NAy = NAx;
PupilEfield = NAx;
numpoints = 1;
for p = (halfpix-halfrange):(halfpix+halfrange)
    if max(EAmpAfterMask(p,:)) > EAmpThreshhold
        for q = (halfpix-halfrange):(halfpix+halfrange)
            if EAmpAfterMask(p,q) > EAmpThreshhold
                NAx(numpoints) = p;
                NAy(numpoints) = q;
                PupilEfield(numpoints) = EAfterMask(p,q);
                numpoints = numpoints + 1;
            end
        end
    end
end
numpoints = numpoints - 1;
NAx = 2.5.*namax.*(NAx(1:numpoints) - halfpix)./halfrange;
NAy = 2.5.*namax.*(NAy(1:numpoints) - halfpix)./halfrange;
PupilEfield = PupilEfield(1:numpoints);
%
%find the k vector direction for all illuminated points in the pupil:
sinth = sqrt(NAx.^2 + NAy.^2);
costh = sqrt(1 - sinth.^2);  %angle of k vector to y axis
sinphi = NAy./sinth;
cosphi = NAx./sinth;
kx = cosphi.*sinth;
kz = sinphi.*sinth;
ky = costh;
%
%
%plot a yz slice at x = 0 through the dithered exc PSF;
yplot_range_xzPSF = 300;
ycalc_range = 400;
ypixsize = 0.25;
yp = 0:ypixsize:ycalc_range;
zp = -plot_range:pixsize:plot_range;
[YP, ZP] = ndgrid(yp,zp);
Q = size(YP);
DitheredIntensity = zeros(Q(1),Q(2));
%In the dithered mode, only the intensity contributed by electric fields in
%the same z column remain after dither, since they have kx = 0 in
%intensity:
vidfile = VideoWriter([imgpath2, 'Components of Dithered xz PSF'],'Uncompressed AVI');
vidfile.FrameRate = 1;
open(vidfile);
NumberComponentsDitheredLightSheet = 0;
for p = (halfpix-halfrange):(halfpix+halfrange)
    if max(EAmpAfterMask(:,p)) > EAmpThreshhold
        %the current x column has E field above threshold
        numpoints = 1;
        NAx = zeros(1,NumPointsWithinBoundingSquare);
        NAy = NAx;
        PupilEfield = NAx;
        for q = (halfpix-halfrange):(halfpix+halfrange)
            if EAmpAfterMask(q,p) > EAmpThreshhold
                NAx(numpoints) = p;
                NAy(numpoints) = q;
                PupilEfield(numpoints) = EAfterMask(q,p);
                numpoints = numpoints + 1;
            end
        end
        numpoints = numpoints - 1;
        NAx = 2.5.*namax.*(NAx(1:numpoints) - halfpix)./halfrange;
        NAy = 2.5.*namax.*(NAy(1:numpoints) - halfpix)./halfrange;
        PupilEfield = PupilEfield(1:numpoints);
        sinth = sqrt(NAx.^2 + NAy.^2);
        costh = sqrt(1 - sinth.^2);  %angle of k vector to y axis
        sinphi = NAy./sinth;
        cosphi = NAx./sinth;
        kx = cosphi.*sinth;
        kz = sinphi.*sinth;
        ky = costh;
        EField = zeros(Q(1),Q(2));
        for q = 1:numpoints
            EField = EField + PupilEfield(q).*exp(2 .* pi .* 1i .* (ky(q).*YP + kz(q).*ZP));
        end
        PartialDitheredIntensity = EField.*conj(EField);
        NumberComponentsDitheredLightSheet = NumberComponentsDitheredLightSheet + 1;
        DitheredIntensity = DitheredIntensity + PartialDitheredIntensity;
        PartialDitheredIntensity = PartialDitheredIntensity'./max(max(PartialDitheredIntensity));        
        PartialDitheredIntensity2 = PartialDitheredIntensity(round(0.5.*halfpix):round(1.5.*halfpix),1:round(yplot_range_xzPSF./ypixsize+1));
        A = size(PartialDitheredIntensity2);
        %find the y HWHM propagation length of the current contribution to the total 
        %   dithered intensity from the intensity along the y axis (i.e., z = 0):
        midpix = (A(1)+1)./2;
        OnAxisIntensity = PartialDitheredIntensity2(midpix,:);
        OnAxisIntensity = OnAxisIntensity./max(OnAxisIntensity);
        YThreshold = (OnAxisIntensity <= 0.5.*OnAxisIntensity(1));
        [threshval, threshindex] = max(YThreshold);
        yFWHM_PartialIntensity = 2.*threshindex.*ypixsize;
        %now plot the yz view of the current contribution to the dithered light sheet:
        figure  %create a new figure window for the plots
        set(gcf, 'Position', [150 200 1300 1300.*plot_range./yplot_range_xzPSF]);
        if ~ispc
            set(gcf, 'visible', 'off');
        end

        axes_h = axes('Position', [0.08, 0.12, 0.85, 0.75]);
        image(256 .* PartialDitheredIntensity2);
        colormap hot(256);
        axis([1 A(2) 1 A(1)]);
        set(gca, 'XTick', [1:(A(2)-1)./6:A(2)]);
        set(gca, 'XTickLabel', yplot_range_xzPSF .* [0:(1./6):1]);
        xlabel(['y / \lambda'], 'FontSize', 14);
        set(gca, 'YTick', [1:(A(1)-1)./4:A(1)]);
        set(gca, 'YTickLabel', plot_range .* [-1:0.5:1]./2);
        ylabel(['z / \lambda'], 'FontSize', 14);
        text(0.25 .*A(1), -0.1 .* A(2), ['Propagation of Component ', num2str(NumberComponentsDitheredLightSheet), ' of the Dithered Light Sheet, annulus NA', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')], 'FontSize', 12);
        text(0.43 .*A(1), -0.05 .* A(2), [lightsheet_descrip, ', lattice NA = ', num2str(NAlattice, '%1.2f'), ', NAsigma = ', num2str(NAsigma, '%1.2f')], 'FontSize', 12);
        text(0.9.*A(1), 0.06.*A(2), ['y FWHM = ', num2str(yFWHM_PartialIntensity, '%1.2f'), '\lambda'], 'Color', 'white', 'FontSize', 14);
%        saveas(gcf ,[imgpath2 'Propagation of Dithered Component ', num2str(NumberComponentsDitheredLightSheet),'.png']);    %
        F(p) = getframe(gcf); 
        writeVideo(vidfile,F(p));
    end
    if mod(p,20) == 0
        close all
    end
end
DitheredIntensity = DitheredIntensity'./max(max(DitheredIntensity));
DitheredIntensity2 = DitheredIntensity(round(0.5.*halfpix):round(1.5.*halfpix),1:round(yplot_range_xzPSF./ypixsize+1));
A = size(DitheredIntensity2);
%
%
%refine the estimate of the y HWHM propagation length by evalauting the
%intensity of the dithered light sheet along the y axis (i.e., z = 0):
midpix = (A(1)+1)./2;
OnAxisIntensity = DitheredIntensity2(midpix,:);
OnAxisIntensity = OnAxisIntensity./max(OnAxisIntensity);
YThreshold = (OnAxisIntensity <= 0.5.*OnAxisIntensity(1));
[threshval, threshindex] = max(YThreshold);
yHWHM = threshindex.*ypixsize;
yFWHM = 2.* yHWHM;
%
%save the OnAxisIntensity:
save([imgpath2 'OnAxisIntensity.mat'], 'OnAxisIntensity');
%
%
%now plot the yz view of the dithered light sheet:
figure  %create a new figure window for the plots
set(gcf, 'Position', [150 200 1300 1300.*plot_range./yplot_range_xzPSF]);
set(gcf, 'visible', 'off');
axes_h = axes('Position', [0.08, 0.12, 0.85, 0.75]);
image(256 .* DitheredIntensity2);
colormap hot(256);
axis([1 A(2) 1 A(1)]);
set(gca, 'XTick', [1:(A(2)-1)./6:A(2)]);
set(gca, 'XTickLabel', yplot_range_xzPSF .* [0:(1./6):1]);
xlabel(['y / \lambda'], 'FontSize', 14);
set(gca, 'YTick', [1:(A(1)-1)./4:A(1)]);
set(gca, 'YTickLabel', plot_range .* [-1:0.5:1]./2);
ylabel(['z / \lambda'], 'FontSize', 14);
text(0.35 .*A(1), -0.1 .* A(2), ['Dithered LLS Propagation in x = 0 Plane, annulus NA', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')], 'FontSize', 12);
text(0.43 .*A(1), -0.05 .* A(2), [lightsheet_descrip, ', lattice NA = ', num2str(NAlattice, '%1.2f'), ', NAsigma = ', num2str(NAsigma, '%1.2f')], 'FontSize', 12);
text(0.9.*A(1), 0.06.*A(2), ['y FWHM = ', num2str(yFWHM, '%1.2f'), '\lambda'], 'Color', 'white', 'FontSize', 14);
%saveas(gcf ,[imgpath2 'Dithered LLS Propagation.png']);    %
f = getframe(gcf); 
writeVideo(vidfile,f);
close(vidfile);
%
%
%plot the intensity along the optical axis:
DitheredIntensity2 = DitheredIntensity(round(0.5.*halfpix):round(1.5.*halfpix),1:round(yplot_range./ypixsize+1));
G = size(DitheredIntensity2);
Symmmetric_yz_profile = zeros(G(1), (2.*G(2)-1));
Symmmetric_yz_profile(:,1:G(2)) = flip(DitheredIntensity2,2);
Symmmetric_yz_profile(:,(G(2):(2.*G(2)-1))) = DitheredIntensity2;
H = size(Symmmetric_yz_profile);
halfpix2 = (H+1)./2;
yaxis_Intensity = Symmmetric_yz_profile(halfpix2(1),:);
figure  %create a new figure window for the plots
set(gcf, 'Position', [450 100 600 600]);
axes_h = axes('Position', [0.12, 0.1, 0.8, 0.8]);
axis square;
plot(yaxis_Intensity, 'r', 'LineWidth', 2);
%now calculate the ideal sinc^2 intensity function corresponding to a
%    rect (top-hat) ky function for the light sheet.  If NAlattice = 
%    max(NAexc) = NAsinc and thetasinc = asin(NAsinc/index), then the
%    rect function is positive for ky = 1 to k*cos(thetasinc), and the
%    ideal function for the intenisty along the optical axis is:
%    I(y) = sinc^2(k*(1-cos(thetasinc))*y)
costhetasinc = sqrt(1 - nasigma.^2);
%the location of the first zero occurs when k*cos(thetasinc)*ynode = pi,
%   Thus ynode = pi/(k*(1-cos(thetasinc))) = lambda_n/(2.*(1-cos(thetasinc))), where
%   lambda_n is the wavelength in the media.  Thus, I(y) = sinc^2(pi*y/ynode)
y = -yplot_range:ypixsize:yplot_range;
sincsq = (sin(pi.*y.*(1-costhetasinc)+eps)./(pi.*y.*(1-costhetasinc)+eps)).^2;
hold on
plot(sincsq, 'b', 'LineWidth', 1, 'LineStyle', '--');
set(gca, 'XTick', [1:((H(2)-1)./10):H(2)]);
set(gca, 'XTickLabel', yplot_range .* [-1:0.2:1]);
xlabel(['z / \lambda'], 'FontSize', 14);
set(gca, 'YTick', [0:0.25:1]);
set(gca, 'YTickLabel', [0:0.25:1]);
ylabel(['Intensity'], 'FontSize', 14);
grid on;
text(0.23 .*H(2), 1.06, ['y Axis Intensity, annulus NA', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')], 'FontSize', 12);
text(0.20 .*H(2), 1.02, [lightsheet_descrip, ', lattice NA = ', num2str(NAlattice, '%1.2f'), ', NAsigma = ', num2str(NAsigma, '%1.2f')], 'FontSize', 12);
saveas(gcf ,[imgpath2 'y Axis Intensity.png']); 
%
%calc and plot the Fourier transform of the y axis intensity:
fftyaxisIntensity = abs(fftshift(fft(ifftshift(yaxis_Intensity))));
fftyaxisIntensity = fftyaxisIntensity./max(fftyaxisIntensity);
fftsincsq = abs(fftshift(fft(ifftshift(sincsq))));
fftsincsq = fftsincsq./max(fftsincsq);
%extract the part from +/-0.1 kzmax:
G = size(fftyaxisIntensity);
halfpix2 = (G(2)-1)./2;
fftyaxisIntensity2 = fftyaxisIntensity(0.9.*(halfpix2):(1.1.*halfpix2)+1);
fftsincsq2 = fftsincsq(0.9.*(halfpix2):(1.1.*halfpix2)+1);
G = size(fftyaxisIntensity2);
figure  %create a new figure window for the plots
set(gcf, 'Position', [450 100 600 600]);
axes_h = axes('Position', [0.12, 0.1, 0.8, 0.8]);
axis square;
plot(fftyaxisIntensity2, 'r', 'LineWidth', 2)
hold on;
plot(fftsincsq2, 'b', 'LineWidth', 1, 'LineStyle', '--');
set(gca, 'YTick', 0:0.2:1);
set(gca, 'YTickLabel', [0:0.2:1]);
xlabel(['k_{z} / (4\pi/\lambda)'], 'FontSize', 14);
set(gca, 'XTick', [1:(G(2)-1)./4:G(2)]);
set(gca, 'XTickLabel', [-0.1:0.05:0.1], 'FontSize', 12);
ylabel(['Amplitude'], 'FontSize', 14);
grid on;
text(0.23 .*G(2), 1.06, ['y Axis OTF, annulus NA', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')], 'FontSize', 12);
text(0.20 .*G(2), 1.02, [lightsheet_descrip, ', lattice NA = ', num2str(NAlattice, '%1.2f'), ', NAsigma = ', num2str(NAsigma, '%1.2f')], 'FontSize', 12);
saveas(gcf ,[imgpath2 'y Axis OTF.png']); 
%
%calc and plot the Fourier transform of the yz view:
figure  %create a new figure window for the plots
set(gcf, 'Position', [450 100 600 600]);
axes_h = axes('Position', [0.13, 0.12, 0.8, 0.8]);
axis square;
fftyzview = abs(fftshift(fft2(ifftshift(Symmmetric_yz_profile))));
fftyzview = fftyzview./max(max(fftyzview));
halfpix2 = (H-1)./2;
fftyzview = fftyzview(((0.6.*halfpix2):(1.4.*halfpix2)),:);
image(256 .* (fftyzview.^0.5));
colormap hot(256);
G = size(fftyzview);
set(gca, 'XTick', [1:(G(2)-1)./4:G(2)]);
set(gca, 'XTickLabel', [-1:0.5:1], 'FontSize', 12);
xlabel(['k_{y} / (4\pi/\lambda)'], 'FontSize', 14);
set(gca, 'YTick', [1:(G(1)-1)./4:G(1)]);
set(gca, 'YTickLabel', [-1:0.5:1], 'FontSize', 12);
ylabel(['k_{z} / (4\pi/\lambda)'], 'FontSize', 14);
text(0.3 .*G(2), -0.06 .* G(1), ['yz OTF, annulus NA', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')], 'FontSize', 12);
text(0.23.*G(2), -0.02 .* G(1), [lightsheet_descrip, ', lattice NA = ', num2str(NAlattice, '%1.2f'), ', NAsigma = ', num2str(NAsigma, '%1.2f')], 'FontSize', 12);
saveas(gcf ,[imgpath2 'yz OTF.png']);    %
%
%plot the xz excitation PSF at at multiple points along the propagation axis,
%   with a step size determined by the yFWHM:
numslices = 16;
ystepsize = ceil(1.5.*yFWHM/(numslices-1)/6).*ystepsize_in;  %ystepsize is a multiple of 3 lambda
yp = ystepsize.*(0:1:(numslices-1));
%to calculate the stationary xz PSF vs propagation distance y, first find
%   all points of significant E field within the rear pupil:
%loop through all rows and columns of EAmpAfterMask within the square
%   bounding the annulus:
halfrange = ceil(namax.*halfpix./2);
NumPointsWithinBoundingSquare = (2.*halfrange + 1).^2;
%define matrices to store the info about the illuminated points in the pupil:
NAx = zeros(1,NumPointsWithinBoundingSquare);
NAy = NAx;
PupilEfield = NAx;
numpoints = 1;
EAmpThreshhold = 0.05;
for p = (halfpix-halfrange):(halfpix+halfrange)
    if max(EAmpAfterMask(p,:)) > EAmpThreshhold
        for q = (halfpix-halfrange):(halfpix+halfrange)
            if EAmpAfterMask(p,q) > EAmpThreshhold
                NAx(numpoints) = p;
                NAy(numpoints) = q;
                PupilEfield(numpoints) = EAfterMask(p,q);
                numpoints = numpoints + 1;
            end
        end
    end
end
numpoints = numpoints - 1;
NAx = 2.5.*namax.*(NAx(1:numpoints) - halfpix)./halfrange;
NAy = 2.5.*namax.*(NAy(1:numpoints) - halfpix)./halfrange;
PupilEfield = PupilEfield(1:numpoints);
%find the k vector direction for all illuminated points in the pupil:
sinth = sqrt(NAx.^2 + NAy.^2);
costh = sqrt(1 - sinth.^2);  %angle of k vector to y axis
sinphi = NAy./sinth;
cosphi = NAx./sinth;
kx = cosphi.*sinth;
kz = sinphi.*sinth;
ky = costh;
%define the plot parameters for the xz excitation PSF at each y position
xp = -plot_range:pixsize:plot_range;
zp = xp;
%the above values are in media wavelengths
[XP, ZP] = ndgrid(xp, zp);
Q = size(XP);
xzExcPSF = zeros(numslices,Q(1),Q(2));  %stores the xz intensity at the 
                   %calculation points yp along the propagation direction
for p = 1:numslices
    EField = zeros(Q(1),Q(2));  %initialize the field before integration 
                                %over all illuminated points in the pupil
    for q = 1:numpoints
        EField = EField + PupilEfield(q).*exp(2 .* pi .* 1i .* (kx(q).*XP + kz(q).*ZP + ky(q).*yp(p)));
    end
    xzExcPSF(p,:,:) = EField.*conj(EField);    
end
%normalize all the xz images along y to the maximum value, presumably at y = 0:
xzExcPSF = xzExcPSF./max(max(max(xzExcPSF)));
%loop thru and plot the xz excitation PSF at all calculated postions
%   yp along the propagation direction:
J = size(xzExcPSF);
midzpix = (J(3)+1)./2;
%imgpath3 = [imgpath2, 'Stationary LLS xz PSF\'];
%mkdir(imgpath2, 'Stationary LLS xz PSF\');
vidfile = VideoWriter([imgpath2, 'Stationary LLS xz PSF'],'Uncompressed AVI');
vidfile.FrameRate = 2;
open(vidfile);
a = size(OnAxisIntensity);
for p = 1:J(1)
    SampleIntensity = squeeze(xzExcPSF(p,:,:));
    figure  %create a new figure window for the plots
    set(gcf, 'Position', [50 100 600 600]);
    if ~ispc
        set(gcf, 'visible', 'off');
    end
    axes_h = axes('Position', [0.13, 0.1, 0.8, 0.8]);
    image(256 .* SampleIntensity);
    colormap hot(256);
    axis([1 J(3) 1 J(2)]);
    axis square
    set(gca, 'XTick', [1:(J(2)-1)./4:J(2)]);
    set(gca, 'XTickLabel', max(xp) .* [-1:0.5:1]);
    xlabel(['x / \lambda'], 'FontSize', 14);
    set(gca, 'YTick', [1:(J(3)-1)./4:J(3)]);
    set(gca, 'YTickLabel', max(zp) .* [-1:0.5:1]);
    ylabel(['z / \lambda'], 'FontSize', 14);
    text(0.0 .*J(2), -0.06 .* J(3), ['Stationary LLS xz PSF ', num2str(yp(p), '%1.2f'), '\lambda From Focus, annulus NA', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')], 'FontSize', 12);
    text(0.23 .*J(2), -0.02 .* J(3), [lightsheet_descrip, ', lattice NA = ', num2str(NAlattice, '%1.2f'), ', NAsigma = ', num2str(NAsigma, '%1.2f')], 'FontSize', 12);
if round(yp(p)./ypixsize + 1) < floor(yplot_range./ypixsize)
    text(0.6.*J(2), 0.05.*J(3), ['On Axis Intensity = ', num2str(OnAxisIntensity(round(yp(p)./ypixsize + 1)), '%1.3f')], 'Color', 'white',  'FontSize', 12);
end
    F(p) = getframe(gcf); 
    writeVideo(vidfile,F(p));
end
close(vidfile);
close all
%
%
%plot the xz excitation OTF at at multiple points along the propagation axis,
%   with a step size determined by the yFWHM:
vidfile = VideoWriter([imgpath2, 'Stationary LLS xz OTF'],'Uncompressed AVI');
vidfile.FrameRate = 2;
open(vidfile);
xzExcOTF = zeros(J(1),J(2),J(3));
K = (J(2)-1).*0.4 + 1;
xzExcOTF2 = zeros(J(1),K,K);  %this version extracts the portion of xzExcOTF
                              % between +/-kmax
halfpix = (J(2)-1)./2;
C = [K K];
for p = 1:J(1)
    %calc and plot the xz excitation OTF of the stationary lattice light sheet
    xzExcOTF(p,:,:) = abs(fftshift(fft2(ifftshift(squeeze(xzExcPSF(p,:,:))))));
    xzExcOTF(p,:,:) = xzExcOTF(p,:,:)./max(max(xzExcOTF(p,:,:)));
    %extract the portion of xzExcOTF from +/-kmax:
    xzExcOTF2(p,:,:) = xzExcOTF(p,(0.6.*halfpix):(1.4.*halfpix),(0.6.*halfpix):(1.4.*halfpix));
    figure  %create a new figure window for the plots
    set(gcf, 'Position', [450 100 600 600]);
    if ~ispc
        set(gcf, 'visible', 'off');
    end
    
    axes_h = axes('Position', [0.13, 0.12, 0.8, 0.8]);
    image(256 .* squeeze(xzExcOTF2(p,:,:)).^gamma);
    colormap hot(256);
%    midpt = (C(1)+1)./2;
%    FattestColumn = round(midpt + C(1).*NAdet./1.33./4);
%    line([1 C(1)], [midpt midpt], 'LineWidth', 1, 'Color', [0 0 1]);
%    line([midpt midpt], [1 C(1)], 'LineWidth', 1, 'Color', [1 0 0]);
%    line([FattestColumn FattestColumn], [1 C(1)], 'LineWidth', 1, 'Color', [0 1 0]);
    axis([1 C(1) 1 C(2)]);
    axis square;
    set(gca, 'XTick', [1:(C(1)-1)./4:C(1)]);
    set(gca, 'XTickLabel', [-1:0.5:1], 'FontSize', 12);
    xlabel(['k_{x} / (4\pi/\lambda)'], 'FontSize', 14);
    set(gca, 'YTick', [1:(C(2)-1)./4:C(2)]);
    set(gca, 'YTickLabel', [-1:0.5:1], 'FontSize', 12);
    ylabel(['k_{z} / (4\pi/\lambda)'], 'FontSize', 14);
    text(0.0 .*C(1), -0.06 .* C(2), ['Stationary LLS xz OTF ', num2str(yp(p), '%1.2f'), '\lambda From Focus, annulus NA', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')], 'FontSize', 12);
    text(0.23.*C(1), -0.02 .* C(2), [lightsheet_descrip, ', lattice NA = ', num2str(NAlattice, '%1.2f'), ', NAsigma = ', num2str(NAsigma, '%1.2f')], 'FontSize', 12);
    F(p) = getframe(gcf); 
    writeVideo(vidfile,F(p));
end
close(vidfile);
%
%
%plot the xz lattice SIM OTF at at multiple points along the propagation axis,
%   with a step size determined by the yFWHM:
vidfile = VideoWriter([imgpath2, 'Lattice SIM xz OTF'],'Uncompressed AVI');
vidfile.FrameRate = 2;
open(vidfile);
%xzSIMOTF = zeros(J(1),J(2),J(3));
%The number of pixels in the widefield OTF must be doubled for the OTF pixel size
%%used here (based on extending the plot_range from +/-50 to +/-100 lambda)
C = size(xz_det_OTF);
[xq,zq] = meshgrid(1:0.5:C(1),1:0.5:C(2));
xz_det_OTF2 = interp2(xz_det_OTF,xq,zq);
%extract the portion of xz_det_OTF2 from +/-kmax:
C = size(xz_det_OTF2);
halfpix = (C(1)-1)./2;
xz_det_OTF2 = xz_det_OTF2((0.6.*halfpix):(1.4.*halfpix),(0.6.*halfpix):(1.4.*halfpix));
D = size(xz_det_OTF2);
xzSIMOTF3 = zeros(J(1),D(1),D(2));
xzSIMOTF5 = zeros(J(1), D(1) * 2 - 1, D(2) * 2 - 1);
for p = 1:J(1)
    %calc and plot the xz excitation OTF of the stationary lattice light sheet
    xzExcOTF3 = squeeze(xzExcOTF2(p,:,:));
    xzSIMOTF = abs(conv2(xzExcOTF3, xz_det_OTF2));
    xzSIMOTF = xzSIMOTF./max(max(xzSIMOTF));
    %extract the data withtin +/-kmax:    
    C = size(xzSIMOTF);
    halfpix = (C(1)-1)./2;
    xzSIMOTF2 = xzSIMOTF((0.5.*halfpix):(1.5.*halfpix),(0.5.*halfpix):(1.5.*halfpix));
    xzSIMOTF3(p,:,:) = xzSIMOTF2;
    xzSIMOTF5(p,:,:) = xzSIMOTF;
    C = size(xzSIMOTF2);
    figure  %create a new figure window for the plots
    set(gcf, 'Position', [450 100 600 600]);
    if ~ispc
        set(gcf, 'visible', 'off');
    end

    axes_h = axes('Position', [0.13, 0.12, 0.8, 0.8]);
    image(256 .* xzSIMOTF2.^gamma);
    colormap hot(256);
    midpt = (C(1)+1)./2;
    FattestColumn = round(midpt + C(1).*NAdet./1.33./4);    
    line([1 C(1)], [midpt midpt], 'LineWidth', 1, 'Color', [0 0 1]);
    line([midpt midpt], [1 C(1)], 'LineWidth', 1, 'Color', [1 0 0]);
    line([FattestColumn FattestColumn], [1 C(1)], 'LineWidth', 1, 'Color', [0 1 0]);
    axis([1 C(1) 1 C(2)]);
    axis square;
    set(gca, 'XTick', [1:(C(1)-1)./4:C(1)]);
    set(gca, 'XTickLabel', [-1:0.5:1], 'FontSize', 12);
    xlabel(['k_{x} / (4\pi/\lambda)'], 'FontSize', 14);
    set(gca, 'YTick', [1:(C(2)-1)./4:C(2)]);
    set(gca, 'YTickLabel', [-1:0.5:1], 'FontSize', 12);
    ylabel(['k_{z} / (4\pi/\lambda)'], 'FontSize', 14);
    text(0.0 .*C(1), -0.06 .* C(2), ['Lattice SIM xz OTF ', num2str(yp(p), '%1.2f'), '\lambda From Focus, annulus NA', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')], 'FontSize', 12);
    text(0.23.*C(1), -0.02 .* C(2), [lightsheet_descrip, ', lattice NA = ', num2str(NAlattice, '%1.2f'), ', NAsigma = ', num2str(NAsigma, '%1.2f')], 'FontSize', 12);
    F(p) = getframe(gcf); 
    writeVideo(vidfile,F(p));
end
close(vidfile);
close all    
%
%extract linecuts through the widefield detection OTF for comparison:
D = size(xz_det_OTF);
halfpix2 = (D(1)-1)./2;
xz_det_OTF3 = xz_det_OTF((0.6.*halfpix2):(1.4.*halfpix2),(0.6.*halfpix2):(1.4.*halfpix2));
%xz_det_OTF = xz_det_OTF3./max(max(xz_det_OTF3));
G = size(xz_det_OTF2);
midpt = (G(1)+1)./2;
FattestColumn = round(midpt + G(1).*NAdet./1.33./4);
LateralWFCrossSection = squeeze(xz_det_OTF2(midpt,:));
AxialWFCrossSection = squeeze(xz_det_OTF2(:,midpt));
BowtieWFCrossSection = squeeze(xz_det_OTF2(:,FattestColumn));
%
%plot orthogonal linecuts through the lattice overall OTF, as well as  
%    an axial linecut at the fattest part of the xz OTF:
vidfile = VideoWriter([imgpath2, 'Lattice SIM OTF Linecuts'],'Uncompressed AVI');
vidfile.FrameRate = 2;
open(vidfile);
for p = 1:J(1)
    %find the location of the peak of the OTF through which to make the axial
    %   and lateral linecuts:
    xzSIMOTF4 = squeeze(xzSIMOTF3(p,:,:));
    D = size(xzSIMOTF4);
    [OTFMax,I] = max(xzSIMOTF4(:));
    I_row = floor(I./D(1));
    I_col = I - D(1).*I_row;
    %now extract the linecuts through the OTF:
    LateralOTFCrossSection = squeeze(xzSIMOTF4(I_col,:));
    LateralOTFCrossSection = LateralOTFCrossSection ./ OTFMax;
    AxialOTFCrossSection = squeeze(xzSIMOTF4(:,I_row))';
    AxialOTFCrossSection = AxialOTFCrossSection ./ OTFMax;
    OffsetAxialLinecut = squeeze(xzSIMOTF4(:,FattestColumn))';
    OffsetAxialLinecut = OffsetAxialLinecut./OTFMax;
    figure  %create a new figure window for the plots
    set(gcf, 'Position', [650 100 600 600]);
    axes_h = axes('Position', [0.13, 0.1, 0.8, 0.77]);
    if ~ispc
        set(gcf, 'visible', 'off');
    end
    D = size(AxialOTFCrossSection);
    hold on
    plot(log10(AxialOTFCrossSection), 'r', 'LineWidth', 2);
    plot(log10(AxialWFCrossSection), 'r', 'LineStyle', '--', 'LineWidth', 2);    hold on
    plot(log10(LateralOTFCrossSection), 'b', 'LineWidth', 2);
    plot(log10(LateralWFCrossSection), 'b', 'LineStyle', '--', 'LineWidth', 2);
    plot(log10(OffsetAxialLinecut), 'g', 'LineWidth', 2);
    plot(log10(BowtieWFCrossSection), 'g', 'LineStyle', '--', 'LineWidth', 2);
    axis([1 D(2) -3 0]);
    axis square;
    grid on;
    set(gca, 'XTick', [1:(D(2)-1)./10:D(2)]);
    set(gca, 'XTickLabel', [-1:0.2:1]);
    xlabel(['k / (4\pi/\lambda)'], 'FontSize', 14);
    set(gca, 'YTick', [-3:1:0]);
    set(gca, 'YTickLabel', 10.^[-3:1:0]);
    ylabel(['OTF Strength'], 'FontSize', 14);
    text(-0.1 .*D(2), 0.23, ['Lattice SIM OTF linecuts ', num2str(yp(p), '%1.2f'), '\lambda From Focus, annulus NA', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')], 'FontSize', 12);
    text(0.23 .*D(2), 0.1, [lightsheet_descrip, ', lattice NA = ', num2str(NAlattice, '%1.2f'), ', NAsigma = ', num2str(NAsigma, '%1.2f')], 'FontSize', 12);
    text(0.02.*D(2), -0.1, 'OTF(k_{x}) at k_{z} = 0', 'Color', 'blue', 'FontSize', 12);
    text(0.7.*D(2), -0.1, '----', 'Color', 'blue', 'FontSize', 12);
    text(0.02.*D(2), -0.24, 'OTF(k_{z}) at k_{x} = 0', 'Color', 'red', 'FontSize', 12);
    text(0.7.*D(2), -0.24, '----', 'Color', 'red', 'FontSize', 12);
    text(0.76.*D(2), -0.19, 'widefield', 'Color', 'black', 'FontSize', 12);    
    text(0.76.*D(2), -0.29, 'versions', 'Color', 'black', 'FontSize', 12);    
    text(0.02.*D(2), -0.38, 'OTF(k_{z}) at k_{x} = 0.5k_{max}', 'Color', 'green', 'FontSize', 12);
    text(0.7.*D(2), -0.38, '----', 'Color', 'green', 'FontSize', 12);
%    saveas(gcf ,[imgpath3 'Overall OTF Linecuts' num2str(p,'%d') '.png'])
    F(p) = getframe(gcf); 
    writeVideo(vidfile,F(p));
end
close(vidfile);
%
%
%now calculate and plot the xz dithered light sheet cross section and cumulative 
%   intensity at every 3 lambda along the propagation direction:
%calculate every 3 lambda
% xruan change 3 lambda with user input 
% yplot_range = 300;  %%%%%
yp = 0:ystepsize_in:(ceil(yplot_range./ystepsize_in))*ystepsize_in;
K = size(yp);
%create a yz matrix of the slice-be-slice normalized dithered cross section 
%   at various longitudinal positions y:
DitheredxzPSFCrossSection = zeros(K(2),J(3));
for p = 1:K(2)
    DitheredxzPSFCrossSection(p,:) = DitheredIntensity(:,round(yp(p)./ypixsize + 1))';
    %normalize the xzPSF cross section for plotting:
    DitheredxzPSFCrossSection(p,:) = DitheredxzPSFCrossSection(p,:)./max(DitheredxzPSFCrossSection(p,:));
end
%
%now find the cumulative intensity with increasing distance from the z = 0
%center of the light sheet:
midpt = (J(3)+1)./2;
CumulativeIntensity = zeros(K(2),J(3));
CumulativeIntensity(:,midpt) = DitheredxzPSFCrossSection(:,midpt);
for zp = (midpt+1):J(3)
    CumulativeIntensity(:,zp) = CumulativeIntensity(:,zp-1) + DitheredxzPSFCrossSection(:,zp);
end
CumulativeIntensity(:,1:(midpt-1)) = flipdim(CumulativeIntensity(:,(midpt+1):J(3)),2);
%normalize to the cumulative intensity from z = 0 to the edge of the z FOV 
[MaxCumulativeIntensity, maxindex] = max(CumulativeIntensity, [], 2);
CumulativeFWHM = zeros(1,K(2));
for p = 1:K(2)
    CumulativeIntensity(p,:) = CumulativeIntensity(p,:) ./ MaxCumulativeIntensity(p);
    %find the full width at half minimum of the cumulative intensity:
    YThreshold = (CumulativeIntensity(p,midpt:J(3)) >= 0.5);
    [threshval, threshindex] = max(YThreshold);
    CumulativeFWHM(p) = 2.*pixsize.*threshindex;
end
vidfile = VideoWriter([imgpath2, 'Cross Sectional and Cumulative Intensity'],'Uncompressed AVI');
vidfile.FrameRate = 5;
open(vidfile);
for p = 1:K(2)
    figure  %create a new figure window for the plots
    set(gcf, 'Position', [250 100 1200 600]);
    if ~ispc
        set(gcf, 'visible', 'off');
    end
    
    axes_h = axes('Position', [0.13, 0.1, 0.8, 0.8]);
    plot(squeeze(DitheredxzPSFCrossSection(p,:)), 'r', 'LineWidth', 1.5);
    hold on
    plot((squeeze(CumulativeIntensity(p,:))), 'b', 'LineWidth', 1.5);
    axis([1 J(3) 0 1]);
    grid on;
    set(gca, 'XTick', [1:(J(3)-1)./4:J(3)]);
    set(gca, 'XTickLabel', max(xp) .* [-1:0.5:1]);
    xlabel(['z / \lambda'], 'FontSize', 14);
    set(gca, 'YTick', [0:0.25:1]);
    set(gca, 'YTickLabel', [0:0.25:1]);
    ylabel(['Light Sheet Cross Section'], 'FontSize', 14);
    text(0.8.*J(2), 0.92, 'Intensity(z)', 'Color', 'red', 'FontSize', 14);
    text(0.8.*J(2), 0.87, 'Cumulative(z)', 'Color', 'blue', 'FontSize', 14);
if round(yp(p)./ypixsize + 1) < floor(yplot_range./ypixsize)
    text(0.06.*J(2), 0.92, ['On Axis Intensity = ', num2str(OnAxisIntensity(round(yp(p)./ypixsize + 1)), '%1.3f')], 'Color', 'black',  'FontSize', 14);
end
    text(0.06.*J(2), 0.87, ['Cumulative FWHM = ', num2str(CumulativeFWHM(p)), '\lambda'], 'Color', 'black',  'FontSize', 14);
    text(0.2.*J(2), 1.08, ['Dithered LLS Cross Sectional Intensity and Cumulative Intensity, ', num2str(yp(p), '%1.2f'), '\lambda From Focus'], 'FontSize', 12);
    text(0.35 .*J(2), 1.04, [lightsheet_descrip, ', lattice NA = ', num2str(NAlattice, '%1.2f'), ', NAsigma = ', num2str(NAsigma, '%1.2f')], 'FontSize', 12);
%    saveas(gcf ,[imgpath3 'Cross Sectional and Cumulative Intensity' num2str(p,'%d') '.png'])
    F(p) = getframe(gcf); 
    writeVideo(vidfile,F(p));
    if mod(p,20) == 0
        close all
    end
end
close(vidfile);
%
%
%now plot the axial excitation OTF of the dithered lattice light sheet along
%   the kx = 0 line at multiple positions y along the propagation direction:
vidfile = VideoWriter([imgpath2, 'Axial Excitation OTF'],'Uncompressed AVI');
vidfile.FrameRate = 5;
open(vidfile);
halfpix = (J(2)-1)./2;
AxialDitheredExcOTF2 = zeros(K(2),(0.4*halfpix+1));
for p = 1:K(2)
    %calc and plot the xz excitation OTF of the dithered lattice light sheet,
    %   using only the cross section data withtin z = +/-50 lambda:
    AxialDitheredExcOTF = abs(fftshift(fft(ifftshift(squeeze(DitheredxzPSFCrossSection(p,(0.5.*halfpix):(1.5.*halfpix)))))));
    AxialDitheredExcOTF = AxialDitheredExcOTF((0.6.*halfpix./2):(1.4.*halfpix./2));
    AxialDitheredExcOTF = AxialDitheredExcOTF ./ max(AxialDitheredExcOTF);
    AxialDitheredExcOTF2(p,:) = AxialDitheredExcOTF;
    figure  %create a new figure window for the plots
    set(gcf, 'Position', [850 100 600 600]);
    if ~ispc
        set(gcf, 'visible', 'off');
    end
    
    axes_h = axes('Position', [0.13, 0.1, 0.8, 0.77]);
    B = size(AxialDitheredExcOTF);
    plot(log10(AxialDitheredExcOTF), 'r', 'LineWidth', 2);
    axis([1 B(2) -3 0]);
    axis square;
    grid on;
    set(gca, 'XTick', [1:(B(2)-1)./10:B(2)]);
    set(gca, 'XTickLabel', [-1:0.2:1]);
    xlabel(['k_{x} / (4\pi/\lambda)'], 'FontSize', 14);
    set(gca, 'YTick', [-3:1:0]);
    set(gca, 'YTickLabel', 10.^[-3:1:0]);
    ylabel(['OTF Strength'], 'FontSize', 14);
    text(-0.05 .*B(2), 0.3, ['Dithered LLS axial excitation OTF ', num2str(yp(p), '%1.2f'), '\lambda From Focus, annulus NA', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')], 'FontSize', 12);
    text(0.23 .*B(2), 0.15, [lightsheet_descrip, ', lattice NA = ', num2str(NAlattice, '%1.2f'), ', NAsigma = ', num2str(NAsigma, '%1.2f')], 'FontSize', 12);
%    saveas(gcf ,[imgpath3 'Axial Excitation OTF' num2str(p,'%d') '.png'])
    F(p) = getframe(gcf); 
    writeVideo(vidfile,F(p));
    if mod(p,20) == 0
        close all
    end
end  
close(vidfile);
%
%
%now calc and plot the xz overall dithered light sheet OTF at multiple positions 
%    y along the propagation direction:
vidfile = VideoWriter([imgpath2, 'Overall OTF'],'Uncompressed AVI');
vidfile.FrameRate = 5;
open(vidfile);
B = size(AxialDitheredExcOTF2);
C = size(xz_det_OTF);
%extract the portion of xz_det_OTF from +/-kmax:
halfpix2 = (C(1)-1)./2;
xz_det_OTF3 = xz_det_OTF((0.6.*halfpix2):(1.4.*halfpix2),(0.6.*halfpix2):(1.4.*halfpix2));
C = size(xz_det_OTF3);
%%The number of pixels in the OTF must be doubled for the OTF pixel size
%%used here (based on extending the plot_range from +/-50 to +/-100 lambda)
%[xq,zq] = meshgrid(1:0.5:C(1),1:0.5:C(2));
halfpix3 = (C(1)-1)./2;
Overall_Dithered_xz_OTF2 = zeros(K(2),C(1),C(2));
for p = 1:K(2)
    %convolve the excitation and detection OTFs to find the overall OTF by
    %summing the detection OTF repeatedly across z, weighted by the axial
    %excitation OTF:
    Overall_Dithered_xz_OTF = zeros(C(1),C(2));
    for q = 1:B(2)
        detrowi = max((halfpix3 - q + 1), 1);
        detrowf = min((C(1) + halfpix3 - q), C(1));
        overallrowi = max((q - halfpix3), 1);
        overallrowf = overallrowi + detrowf - detrowi;
        Overall_Dithered_xz_OTF(overallrowi:overallrowf,:) = Overall_Dithered_xz_OTF(overallrowi:overallrowf,:) + AxialDitheredExcOTF2(p,q).*xz_det_OTF3(detrowi:detrowf,:);
    end
    Overall_Dithered_xz_OTF = Overall_Dithered_xz_OTF./max(max(Overall_Dithered_xz_OTF));
    Overall_Dithered_xz_OTF2(p,:,:) = Overall_Dithered_xz_OTF;
    C = size(Overall_Dithered_xz_OTF);
    figure  %create a new figure window for the plots
    set(gcf, 'Position', [450 100 600 600]);
    if ~ispc
        set(gcf, 'visible', 'off');
    end
    
    axes_h = axes('Position', [0.13, 0.12, 0.8, 0.8]);
    image(256 .* Overall_Dithered_xz_OTF.^gamma);
    colormap hot(256);
    midpt = (C(1)+1)./2;
    FattestColumn = round(midpt + C(1).*NAdet./1.33./4);
    line([1 C(1)], [midpt midpt], 'LineWidth', 1, 'Color', [0 0 1]);
    line([midpt midpt], [1 C(1)], 'LineWidth', 1, 'Color', [1 0 0]);
    line([FattestColumn FattestColumn], [1 C(1)], 'LineWidth', 1, 'Color', [0 1 0]);
    axis([1 C(1) 1 C(2)]);
    axis square;
    set(gca, 'XTick', [1:(C(1)-1)./4:C(1)]);
    set(gca, 'XTickLabel', [-1:0.5:1], 'FontSize', 12);
    xlabel(['k_{x} / (4\pi/\lambda)'], 'FontSize', 14);
    set(gca, 'YTick', [1:(C(2)-1)./4:C(2)]);
    set(gca, 'YTickLabel', [-1:0.5:1], 'FontSize', 12);
    ylabel(['k_{z} / (4\pi/\lambda)'], 'FontSize', 14);
    text(0.0 .*C(1), -0.06 .* C(2), ['Dithered LLS Overall OTF ', num2str(yp(p), '%1.2f'), '\lambda From Focus, annulus NA', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')], 'FontSize', 12);
    text(0.23.*C(1), -0.02 .* C(2), [lightsheet_descrip, ', lattice NA = ', num2str(NAlattice, '%1.2f'), ', NAsigma = ', num2str(NAsigma, '%1.2f')], 'FontSize', 12);
%    saveas(gcf ,[imgpath3 'Overall OTF' num2str(p,'%d') '.png'])
    F(p) = getframe(gcf); 
    writeVideo(vidfile,F(p));
    if mod(p,20) == 0
        close all
    end
end
close(vidfile)
%
%
%extract linecuts through the widefield detection OTF for comparison:
D = size(xz_det_OTF);
halfpix2 = (D(1)-1)./2;
xz_det_OTF3 = xz_det_OTF((0.6.*halfpix2):(1.4.*halfpix2),(0.6.*halfpix2):(1.4.*halfpix2));
xz_det_OTF3 = xz_det_OTF3./max(max(xz_det_OTF3));
LateralWFCrossSection = squeeze(xz_det_OTF3(midpt,:));
AxialWFCrossSection = squeeze(xz_det_OTF3(:,midpt));
BowtieWFCrossSection = squeeze(xz_det_OTF3(:,FattestColumn));
%
%plot orthogonal linecuts through the lattice overall OTF, as well as  
%    an axial linecut at the fattest part of the xz OTF:
vidfile = VideoWriter([imgpath2, 'Overall OTF Linecuts'],'Uncompressed AVI');
vidfile.FrameRate = 5;
open(vidfile);
for p = 1:K(2)
    %find the location of the peak of the OTF through which to make the axial
    %   and lateral linecuts:
    Overall_Dithered_xz_OTF3 = squeeze(Overall_Dithered_xz_OTF2(p,:,:));
    D = size(Overall_Dithered_xz_OTF3);
    [OTFMax,I] = max(Overall_Dithered_xz_OTF3(:));
    I_row = floor(I./D(1));
    I_col = I - D(1).*I_row;
    %now extract the linecuts through the OTF:
    LateralOTFCrossSection = squeeze(Overall_Dithered_xz_OTF3(I_col,:));
    LateralOTFCrossSection = LateralOTFCrossSection ./ OTFMax;
    AxialOTFCrossSection = squeeze(Overall_Dithered_xz_OTF3(:,I_row))';
    AxialOTFCrossSection = AxialOTFCrossSection ./ OTFMax;
    OffsetAxialLinecut = squeeze(Overall_Dithered_xz_OTF3(:,FattestColumn))';
    OffsetAxialLinecut = OffsetAxialLinecut./OTFMax;
    figure  %create a new figure window for the plots
    set(gcf, 'Position', [650 100 600 600]);
    if ~ispc
        set(gcf, 'visible', 'off');
    end    
    axes_h = axes('Position', [0.13, 0.1, 0.8, 0.77]);
    D = size(AxialOTFCrossSection);
    hold on
    plot(log10(AxialOTFCrossSection), 'r', 'LineWidth', 2);
    plot(log10(AxialWFCrossSection), 'r', 'LineStyle', '--', 'LineWidth', 2);    hold on
    plot(log10(LateralOTFCrossSection), 'b', 'LineWidth', 2);
    plot(log10(LateralWFCrossSection), 'b', 'LineStyle', '--', 'LineWidth', 2);
    plot(log10(OffsetAxialLinecut), 'g', 'LineWidth', 2);
    plot(log10(BowtieWFCrossSection), 'g', 'LineStyle', '--', 'LineWidth', 2);
    axis([1 D(2) -3 0]);
    axis square;
    grid on;
    set(gca, 'XTick', [1:(D(2)-1)./10:D(2)]);
    set(gca, 'XTickLabel', [-1:0.2:1]);
    xlabel(['k / (4\pi/\lambda)'], 'FontSize', 14);
    set(gca, 'YTick', [-3:1:0]);
    set(gca, 'YTickLabel', 10.^[-3:1:0]);
    ylabel(['OTF Strength'], 'FontSize', 14);
    text(-0.1 .*D(2), 0.23, ['Dithered LLS Overall OTF linecuts ', num2str(yp(p), '%1.2f'), '\lambda From Focus, annulus NA', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')], 'FontSize', 12);
    text(0.23 .*D(2), 0.1, [lightsheet_descrip, ', lattice NA = ', num2str(NAlattice, '%1.2f'), ', NAsigma = ', num2str(NAsigma, '%1.2f')], 'FontSize', 12);
    text(0.02.*D(2), -0.1, 'OTF(k_{x}) at k_{z} = 0', 'Color', 'blue', 'FontSize', 12);
    text(0.7.*D(2), -0.1, '----', 'Color', 'blue', 'FontSize', 12);
    text(0.02.*D(2), -0.24, 'OTF(k_{z}) at k_{x} = 0', 'Color', 'red', 'FontSize', 12);
    text(0.7.*D(2), -0.24, '----', 'Color', 'red', 'FontSize', 12);
    text(0.76.*D(2), -0.19, 'widefield', 'Color', 'black', 'FontSize', 12);    
    text(0.76.*D(2), -0.29, 'versions', 'Color', 'black', 'FontSize', 12);    
    text(0.02.*D(2), -0.38, 'OTF(k_{z}) at k_{x} = 0.5k_{max}', 'Color', 'green', 'FontSize', 12);
    text(0.7.*D(2), -0.38, '----', 'Color', 'green', 'FontSize', 12);
%    saveas(gcf ,[imgpath3 'Overall OTF Linecuts' num2str(p,'%d') '.png'])
    F(p) = getframe(gcf); 
    writeVideo(vidfile,F(p));
    if mod(p,20) == 0
        close all
    end
end
close(vidfile);
%
% xruan 
% calc xz overall PSF at multiple positions y along propagation direction
% with size 1001 X 10001
PSF_out_range = 50;
PSFOuthalfpix = round(PSF_out_range./pixsize);
PSFOutnumpix = 2.*PSFOuthalfpix + 1;
DitheredPSF1 = zeros(K(2),PSFOutnumpix,PSFOutnumpix);

xz_det_PSF1 = circshift(xz_det_PSF, [1, 1]);
OverallDitheredPSF1 = zeros(K(2), 1001, 1001);
for p = 1 : K(2)
    %calculate the dithered PSF over the whole range from the dithered
    %    xz PSF cross section:
    DitheredxzPSF = ones(J(3),1)*squeeze(DitheredxzPSFCrossSection(p,:));    
    %find the peak of the current dithered PSF along the z axis near z = 0:
    zoffset = 1;    
    DitheredPSF1(p,:,:) = DitheredxzPSF((halfpix - PSFOuthalfpix):(halfpix + PSFOuthalfpix),(halfpix - PSFOuthalfpix + zoffset):(halfpix + PSFOuthalfpix + zoffset))';
    OverallDitheredPSF1(p, :, :) = OnAxisIntensity(round(yp(p)./ypixsize + 1)).*squeeze(DitheredPSF1(p,:,:)) .* xz_det_PSF1;    
end
OverallDitheredPSF1 = OverallDitheredPSF1 .* (OverallDitheredPSF1 > 0);
%
%now calc and plot the xz overall PSF at multiple positions y along the 
%   propagation direction:
PSF_plot_range = 10;
PSFhalfpix = round(PSF_plot_range./pixsize);
PSFnumpix = 2.*PSFhalfpix + 1;
DitheredPSF2 = zeros(K(2),PSFnumpix,PSFnumpix);
OverallDitheredPSF = DitheredPSF2;

%imgpath3 = [imgpath2, 'Overall PSF\'];
%mkdir(imgpath2, 'Overall PSF\');
vidfile = VideoWriter([imgpath2, 'Overall PSF'],'Uncompressed AVI');
vidfile.FrameRate = 5;
open(vidfile);
for p = 1:K(2)
    %calculate the dithered PSF over the PSF plot range from the dithered
    %    xz PSF cross section:
    DitheredxzPSF(:,:) = ones(J(3),1)*squeeze(DitheredxzPSFCrossSection(p,:));
    %find the peak of the current dithered PSF along the z axis near z = 0:
    zoffset = 1;
    DitheredPSF2(p,:,:) = DitheredxzPSF((halfpix - PSFhalfpix):(halfpix + PSFhalfpix),(halfpix - PSFhalfpix + zoffset):(halfpix + PSFhalfpix + zoffset))';
    OverallPSF = OnAxisIntensity(round(yp(p)./ypixsize + 1)).*squeeze(DitheredPSF2(p,:,:)) .* xz_det_PSF((500 - PSFhalfpix):(500 + PSFhalfpix),(500 - PSFhalfpix):(500 + PSFhalfpix));
    C = size(OverallPSF);
    OverallDitheredPSF(p,:,:) = OverallPSF;
    figure  %create a new figure window for the plots
    set(gcf, 'Position', [850 100 600 600]);
    if ~ispc
        set(gcf, 'visible', 'off');
    end    
    axes_h = axes('Position', [0.13, 0.1, 0.8, 0.8]);
    image(256 .* OverallPSF);
    colormap hot(256);
    axis([1 C(1) 1 C(2)]);
    axis square;
    set(gca, 'XTick', [1:(C(1)-1)./4:C(1)]);
    set(gca, 'XTickLabel', PSF_plot_range .* [-1:0.5:1]);
    xlabel(['x / \lambda'], 'FontSize', 14);
    set(gca, 'YTick', [1:(C(1)-1)./4:C(1)]);
    set(gca, 'YTickLabel', PSF_plot_range .* [-1:0.5:1]);
    ylabel(['z / \lambda'], 'FontSize', 14);
    text(0.02 .*C(1), -0.06 .* C(2), ['Dithered LLS Overall PSF ', num2str(yp(p), '%1.2f'), '\lambda From Focus, annulus NA', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')], 'FontSize', 12);
    text(0.23.*C(1), -0.02 .* C(2), [lightsheet_descrip, ', lattice NA = ', num2str(NAlattice, '%1.2f'), ', NAsigma = ', num2str(NAsigma, '%1.2f')], 'FontSize', 12);
    text(0.6.*C(1), 0.06.*C(2), ['Max Intensity = ', num2str(OnAxisIntensity(round(yp(p)./ypixsize + 1)), '%1.2f')], 'Color', 'white', 'FontSize', 14);
%    text(0.05.*C(1), 0.06.*C(2), ['On Axis Intensity = ', num2str(OnAxisIntensity(round(yp(p)./ypixsize + 1)), '%1.3f')], 'Color', 'white',  'FontSize', 14);
    text(0.05.*C(1), 0.06.*C(2), ['Cumulative FWHM = ', num2str(CumulativeFWHM(p)), '\lambda'], 'Color', 'white',  'FontSize', 14);
%    saveas(gcf ,[imgpath3 'Overall PSF' num2str(p,'%d') '.png'])
    F = getframe(gcf); 
    writeVideo(vidfile,F);
    if mod(p,20) == 0
        close all
    end
end
close(vidfile)
%mkdir(imgpath2, 'Overall PSF\');
vidfile = VideoWriter([imgpath2, 'Overall PSF gamma0p5'],'Uncompressed AVI');
vidfile.FrameRate = 5;
open(vidfile);
for p = 1:K(2)
    OverallPSF = squeeze(OverallDitheredPSF(p,:,:));
    figure  %create a new figure window for the plots
    set(gcf, 'Position', [850 100 600 600]);
    if ~ispc
        set(gcf, 'visible', 'off');
    end    
    axes_h = axes('Position', [0.13, 0.1, 0.8, 0.8]);
    image(256 .* OverallPSF.^gamma);
    colormap hot(256);
    axis([1 C(1) 1 C(2)]);
    axis square;
    set(gca, 'XTick', [1:(C(1)-1)./4:C(1)]);
    set(gca, 'XTickLabel', PSF_plot_range .* [-1:0.5:1]);
    xlabel(['x / \lambda'], 'FontSize', 14);
    set(gca, 'YTick', [1:(C(1)-1)./4:C(1)]);
    set(gca, 'YTickLabel', PSF_plot_range .* [-1:0.5:1]);
    ylabel(['z / \lambda'], 'FontSize', 14);
    text(0.02 .*C(1), -0.06 .* C(2), ['Dithered LLS Overall PSF ', num2str(yp(p), '%1.2f'), '\lambda From Focus, annulus NA', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')], 'FontSize', 12);
    text(0.23.*C(1), -0.02 .* C(2), [lightsheet_descrip, ', lattice NA = ', num2str(NAlattice, '%1.2f'), ', NAsigma = ', num2str(NAsigma, '%1.2f')], 'FontSize', 12);
    text(0.6.*C(1), 0.06.*C(2), ['Max Intensity = ', num2str(OnAxisIntensity(round(yp(p)./ypixsize + 1)), '%1.2f')], 'Color', 'white', 'FontSize', 14);
%    text(0.05.*C(1), 0.06.*C(2), ['On Axis Intensity = ', num2str(OnAxisIntensity(round(yp(p)./ypixsize + 1)), '%1.3f')], 'Color', 'white',  'FontSize', 14);
    text(0.05.*C(1), 0.06.*C(2), ['Cumulative FWHM = ', num2str(CumulativeFWHM(p)), '\lambda'], 'Color', 'white',  'FontSize', 14);
%    saveas(gcf ,[imgpath3 'Overall PSF 0.5Gamma' num2str(p,'%d') '.png'])
    G = getframe(gcf); 
    writeVideo(vidfile,G);
    if mod(p,20) == 0
        close all
    end
end
close(vidfile);
%
%plot a cross-section through the swept sheet excitation PSF, the 
%    axial detection PSF, and the axial overall PSF at multiple positions 
%    along the propagation direction:
%imgpath3 = [imgpath2, 'PSF Axial Linecuts\'];
%mkdir(imgpath2, 'PSF Axial Linecuts\');
vidfile = VideoWriter([imgpath2, 'PSF Axial Linecuts'],'Uncompressed AVI');
vidfile.FrameRate = 2;
open(vidfile);
for p = 1:K(2)
    SheetCrossSection = squeeze(DitheredPSF2(p,:,(PSFhalfpix+1)));
    SheetCrossSection = SheetCrossSection ./ max(SheetCrossSection);
    OverallCrossSection = squeeze(OverallDitheredPSF(p,:,(PSFhalfpix+1)));
    OverallCrossSection = OverallCrossSection./max(OverallCrossSection);
    DetPSFCrossSection = xz_det_PSF(400:600,501);
    DetPSFCrossSection = DetPSFCrossSection ./ max(DetPSFCrossSection);
    figure  %create a new figure window for the plots
    set(gcf, 'Position', [850 100 600 600]);
    if ~ispc
        set(gcf, 'visible', 'off');
    end    
    axes_h = axes('Position', [0.13, 0.1, 0.8, 0.77]);
    A = size(OverallCrossSection);
    plot(OverallCrossSection, 'r', 'LineWidth', 2);
    hold on
    plot(DetPSFCrossSection, 'b', 'LineWidth', 2);
    plot(SheetCrossSection, 'g', 'LineWidth', 2);
    axis([1 A(2) 0 1]);
    axis square;
    grid on;
    set(gca, 'XTick', [1:(A(2)-1)./10:A(2)]);
    set(gca, 'XTickLabel', 10 .* [-1:0.2:1]);
    xlabel(['z / \lambda'], 'FontSize', 14);
    set(gca, 'YTick', [0:0.25:1]);
    set(gca, 'YTickLabel', [0:0.25:1]);
    ylabel(['Intensity'], 'FontSize', 14);
    text(0.05 .*A(2), 1.06, ['PSF axial linecuts ', num2str(yp(p), '%1.2f'), '\lambda From Focus, annulus NA', num2str(NAmax, '%1.2f'), '/', num2str(NAmin, '%1.2f')], 'FontSize', 12);
    text(0.23 .*A(2), 1.02, [lightsheet_descrip, ', lattice NA = ', num2str(NAlattice, '%1.2f'), ', NAsigma = ', num2str(NAsigma, '%1.2f')], 'FontSize', 12);
    text(0.60  *A(2), 0.95, ['Max Intensity = ', num2str(max(max(OverallDitheredPSF(p,:,:))), '%1.2f')], 'Color', 'black', 'FontSize', 14);
%    saveas(gcf ,[imgpath3 'PSF Axial Linecuts' num2str(p,'%d') '.png'])
    F(p) = getframe(gcf); 
    writeVideo(vidfile,F(p));
    if mod(p,20) == 0
        close all
    end
end
close(vidfile);
close all;

% output line cuts at p = 1
p = 1;
SheetCrossSection = squeeze(DitheredPSF2(p,:,(PSFhalfpix+1)));
SheetCrossSection = SheetCrossSection ./ max(SheetCrossSection);
OverallCrossSection = squeeze(OverallDitheredPSF(p,:,(PSFhalfpix+1)));
OverallCrossSection = OverallCrossSection./max(OverallCrossSection);
DetPSFCrossSection = xz_det_PSF(400:600,501);
DetPSFCrossSection = DetPSFCrossSection ./ max(DetPSFCrossSection);


