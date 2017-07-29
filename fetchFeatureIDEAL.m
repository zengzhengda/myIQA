%% ========================================================================
% IDEAL Implementation, Version 1.0
% Copyright(c) 2016 Dohyoung Lee
% All Rights Reserved.
%
% This is an implementation of the quality-aware feature (QAF) extraction routine of 
% Invariance DEscriptor-based ALgorithm (IDEAL) 
%
% For detailed explanation of the algorithm, see reference:
%       D. Lee and K. N. Plataniotis, "Towards a No-Reference Image Quality Assessment
%       Using Statistics of Perceptual Color Descriptors", in IEEE Transaction on Image Processing, 2016.
%
% Input :   (1) img : The test RGB image to be evaluated
%
% Output:   (1) qaf : The set of quality-aware features (QAFs) extracted from input image
%
% Basic Usage :
%   Given an input RGB image 'distorted.jpg'
%   >>
%   >> img = imread('distorted.jpg');
%   >> qaf = extract_feature_IDEAL(img);
%   
% Code Dependency:
%   This code makes use of the 'BRISQUE Software Release' implementation by Anish Mittal. 
%   (http://live.ece.utexas.edu/research/quality/BRISQUE_release.zip)
%
%========================================================================

%function [l_feat, a_feat, h_feat, o_feat, s_feat] = extract_feature_IDEAL(img)
function qaf = extract_feature_IDEAL(img)
      
    img = uint8(img);
       
    window = fspecial('gaussian', 7, 1); 
    window = window/sum(sum(window));
   
    l_feat = [];
    h_feat = [];
    o_feat = [];
    s_feat = [];
    a_feat = [];
        
    shifts        = [ 0 1; 1 0; 1 1; -1 1];
    
    n_scale = 2;    % scale parameter for quality-aware feature extraction
    for i_scale = 1:n_scale
                
        %% Achromatic QAFs - BRISQUE features
        img_Y = double(rgb2gray(uint8(img)));    
        
        mu            = filter2(window, img_Y, 'same');
        mu_sq         = mu.*mu;
        sigma         = sqrt(abs(filter2(window, img_Y.*img_Y, 'same') - mu_sq));
        structdis     = (img_Y-mu)./(sigma+1);

        for itr_shift = 1:4
            shifted_structdis           = circshift(structdis,shifts(itr_shift,:));
            pair                        = structdis(:).*shifted_structdis(:);
            [alpha, leftstd, rightstd]  = estimateaggdparam(pair);
            const                       = (sqrt(gamma(1/alpha))/sqrt(gamma(3/alpha)));
            meanparam                   = (rightstd-leftstd)*(gamma(2/alpha)/gamma(1/alpha))*const;
            l_feat                      = [l_feat alpha meanparam leftstd^2 rightstd^2];               
        end
      
        %% Chromatic QAFs are only extracted from the original dimension       
        if(i_scale == 1) 
            
            df_sigma = 1;                        
            [im_s, im_h, im_o, im_a] = get_Invariance_Descriptors(double(img),df_sigma);
            
            % Consider only 2 primary orientations
            n_start = 1;    n_end = 2;  
                            
            for i = n_start:n_end

                %% Relative saturation features  
                shifted_s = circshift(im_s, shifts(i,:));
                diff_s = im_s - shifted_s; 
                diff_s = diff_s(2:end-1,2:end-1); diff_s = diff_s(:);
                [paramEsts(1), paramEsts(2)] = estimateggdparam(diff_s); 
                s_feat = [s_feat paramEsts];
                clear diff_s paramEsts                   
                
                %% Relative hue features  
                shifted_h = circshift(im_h, shifts(i,:));
                diff_h = get_circ_diff(im_h, shifted_h); 
                diff_h = diff_h(2:end-1,2:end-1); diff_h = diff_h(:);
                [paramEsts(1), paramEsts(2)] = circ_wcpar_estimate(diff_h);
                cku = get_circ_kurtosis(diff_h);
                h_feat = [h_feat paramEsts cku];
                clear diff_h paramEsts                    

                %% Relative opponent angle features               
                shifted_o = circshift(im_o, shifts(i,:));
                diff_o = get_circ_diff(im_o, shifted_o); 
                diff_o = diff_o(2:end-1,2:end-1); diff_o = diff_o(:);
                [paramEsts(1), paramEsts(2)] = circ_wcpar_estimate(diff_o);
                cku = get_circ_kurtosis(diff_o);
                o_feat = [o_feat paramEsts cku];
                clear diff_o paramEsts 

                %% Relative spherical angle features               
                shifted_a = circshift(im_a, shifts(i,:));
                diff_a = get_circ_diff(im_a, shifted_a); 
                diff_a = diff_a(2:end-1,2:end-1); diff_a = diff_a(:);
                [paramEsts(1), paramEsts(2)] = circ_wcpar_estimate(diff_a);
                cku = get_circ_kurtosis(diff_a);
                a_feat = [a_feat paramEsts cku];
                clear diff_a paramEsts 

            end

        end
        
        %% For the extraction of achromatic QAFs in lower dimension
        img = imresize(img,0.5);
        
    end
     
    qaf = [l_feat, s_feat, h_feat, o_feat, a_feat];
end

%% ================================================================
%
% get_circ_diff() evaluates difference in angle. Output lies within [-pi pi)
%
%============================================================================
function diffA = get_circ_diff(A1, A2)
    diffA = A1 - A2;
    diffA(diffA >= pi) = diffA(diffA >= pi) - (2*pi); 
    diffA(diffA < -pi) = diffA(diffA < -pi) + (2*pi);
end

%% ================================================================
%
% get_circ_kurtosis() evaluates circular kurtosis of input angular samples
%
%============================================================================
function c_kurt = get_circ_kurtosis(theta)
    cosI1 = sum(cos(theta(:))); 
    sinI1 = sum(sin(theta(:)));
    mu    = atan2(sinI1,cosI1);
    c_kurt = sum(cos(2*get_circ_diff(theta,mu))) / numel(theta);

end

%% ================================================================
%
% get_Invariance_Descriptors() evaluates color invariance descriptor maps 
% from provided input image
%
%============================================================================
function [im_s, im_h, im_ox, im_ax] = get_Invariance_Descriptors(in,sig)

    if(nargin == 1)
        sig = 1;
    end

    R=double(in(:,:,1));
    G=double(in(:,:,2));
    B=double(in(:,:,3));

    im_s = 1 - ((3*(min(min(R,G),B)) ./ (R + G + B + eps)));
    
    im_h = atan2((R-G)/sqrt(2),(R+G-(2*B))/sqrt(6));
    im_h(im_h < 0) = im_h(im_h < 0) + (2*pi);
    im_h(isnan(im_h)) = 0;  

    sigma = sig;
    Rx=gDer(R,sigma,1,0);   %R=gDer(R,sigma,0,0);
    Gx=gDer(G,sigma,1,0);   %G=gDer(G,sigma,0,0);
    Bx=gDer(B,sigma,1,0);   %B=gDer(B,sigma,0,0);    
    
    f_O1_x = (Rx-Gx)/sqrt(2);
    f_O2_x = (Rx+Gx-2*Bx)/sqrt(6);
    im_ox = atan2(f_O1_x,f_O2_x);
    im_ox(im_ox < 0) = im_ox(im_ox < 0) + (2*pi);
    im_ox(isnan(im_ox)) = 0;      
    
    t_01_x = (Gx.*R - Rx.*G) ./ sqrt(R.*R + G.*G + eps);
    t_02_x = (Rx.*R.*B + Gx.*G.*B - Bx.*R.*R - Bx.*G.*G) ./ sqrt((R.*R + G.*G + eps) .* (R.*R + G.*G + B.*B + eps));
    im_ax = atan2(t_01_x,t_02_x);
    im_ax(im_ax < 0) = im_ax(im_ax < 0) + (2*pi);
    im_ax(isnan(im_ax)) = 0;      
    
end

%% ================================================================
%
% circ_wcpar_estimate() computes the wrapped Cauchy distribution parameters
% from provided samples
%
%============================================================================
function [a, rho] = circ_wcpar_estimate(alpha)

    alpha = alpha(:);

    % Initialize two variables
    u1 = 0.3; u2 = 0.3;

    thres = 0.000001;
    k = 0;
    while (k < 500)  
        mu1 = u1;
        mu2 = u2;
        w = 1./(1 - mu1*cos(alpha) - mu2*sin(alpha) + eps);
        num1 = w.*cos(alpha);
        num2 = w.*sin(alpha);    
        u1 = sum(num1(:))/sum(w(:));
        u2 = sum(num2(:))/sum(w(:));
        k = k + 1;
        if ((abs(u1-mu1)< thres) && (abs(u2-mu2)< thres))
            break;
        end
    end
    if k == 500
        error('data do not converge');
    end

    a = atan2(u2,u1); 
    
    rho = (1 - sqrt(1 - u1^2 - u2^2))/sqrt(u1^2 + u2^2);
end

%% ================================================================
%
% estimateaggdparam() computes the assymetric generalized Gaussian distribution 
% parameters from provided samples
%
%   This code makes use of the 'BRISQUE Software Release' implementation by Anish Mittal. 
%   (http://live.ece.utexas.edu/research/quality/BRISQUE_release.zip)
%
%============================================================================
function [alpha, leftstd, rightstd] = estimateaggdparam(vec)
    gam   = 0.2:0.001:10;
    r_gam = ((gamma(2./gam)).^2)./(gamma(1./gam).*gamma(3./gam));

    leftstd            = sqrt(mean((vec(vec<0)).^2));
    rightstd           = sqrt(mean((vec(vec>0)).^2));
    gammahat           = leftstd/rightstd;
    rhat               = (mean(abs(vec)))^2/mean((vec).^2);
    rhatnorm           = (rhat*(gammahat^3 +1)*(gammahat+1))/((gammahat^2 +1)^2);
    [min_difference, array_position] = min((r_gam - rhatnorm).^2);
    alpha              = gam(array_position);
end

%% ================================================================
%
% estimateggdparam() computes generalized Gaussian distribution parameters
% from provided samples
%
%   This code makes use of the 'BRISQUE Software Release' implementation by Anish Mittal. 
%   (http://live.ece.utexas.edu/research/quality/BRISQUE_release.zip)
%
%============================================================================
function [gamparam, sigma] = estimateggdparam(vec)
    gam                              = 0.2:0.001:10;
    r_gam                            = (gamma(1./gam).*gamma(3./gam))./((gamma(2./gam)).^2);

    sigma_sq                         = mean((vec).^2);
    sigma                            = sqrt(sigma_sq);
    E                                = mean(abs(vec));
    rho                              = sigma_sq/E^2;
    [min_difference, array_position] = min(abs(rho - r_gam));
    gamparam                         = gam(array_position);  
end

%% ================================================================
%
% gDer() computes Gaussian derivatives
%
%   This code makes use of the 'Color Descriptor' implementation by Joost van de Weijer. 
%   (http://cat.uab.es/~joost/software.html)
%
%============================================================================
function [H]= gDer(f,sigma, iorder,jorder)
    % compute Gaussian derivatives

    break_off_sigma = 3.;
    filtersize = floor(break_off_sigma*sigma+0.5);

    f=fill_border(f,filtersize);

    x=-filtersize:1:filtersize;

    Gauss=1/(sqrt(2 * pi) * sigma)* exp((x.^2)/(-2 * sigma * sigma) );

    switch(iorder)
    case 0
        Gx= Gauss/sum(Gauss);
    case 1
        Gx  =  -(x/sigma^2).*Gauss;
        Gx  =  Gx./(sum(sum(x.*Gx)));
    case 2
        Gx = (x.^2/sigma^4-1/sigma^2).*Gauss;
        Gx = Gx-sum(Gx)/size(x,2);
        Gx = Gx/sum(0.5*x.*x.*Gx);
    end
    H = filter2(Gx,f);

    switch(jorder)
    case 0
        Gy= Gauss/sum(Gauss);
    case 1
        Gy  =  -(x/sigma^2).*Gauss;
        Gy  =  Gy./(sum(sum(x.*Gy)));
    case 2
        Gy = (x.^2/sigma^4-1/sigma^2).*Gauss;
        Gy = Gy-sum(Gy)/size(x,2);
        Gy = Gy/sum(0.5*x.*x.*Gy);
    end
    H = filter2(Gy',H);

    H=H(filtersize+1:size(H,1)-filtersize,filtersize+1:size(H,2)-filtersize);
end

function out=fill_border(in,bw)

    hh=size(in,1);
    ww=size(in,2);
    dd=size(in,3);

    if(dd==1)
        out=zeros(hh+bw*2,ww+bw*2);

        out(1:bw,1:bw)=ones(bw,bw).*in(1,1);
        out(bw+hh+1:2*bw+hh,1:bw)=ones(bw,bw).*in(hh,1);
        out(1:bw,bw+1+ww:2*bw+ww)=ones(bw,bw).*in(1,ww);
        out(bw+hh+1:2*bw+hh,bw+1+ww:2*bw+ww)=ones(bw,bw).*in(hh,ww);
        out( bw+1:bw+hh,bw+1:bw+ww )= in;
        out(1:bw,bw+1:bw+ww)=ones(bw,1)*in(1,:);
        out(bw+hh+1:2*bw+hh,bw+1:bw+ww)=ones(bw,1)*in(hh,:);
        out(bw+1:bw+hh,1:bw)=in(:,1)*ones(1,bw);
        out(bw+1:bw+hh,bw+ww+1:2*bw+ww)=in(:,ww)*ones(1,bw);
    else
        out=zeros(hh+bw*2,ww+bw*2,dd);
        for ii = 1:dd
            out(1:bw,1:bw,ii)=ones(bw,bw).*in(1,1,ii);
            out(bw+hh+1:2*bw+hh,1:bw,ii)=ones(bw,bw).*in(hh,1,ii);
            out(1:bw,bw+1+ww:2*bw+ww,ii)=ones(bw,bw).*in(1,ww,ii);
            out(bw+hh+1:2*bw+hh,bw+1+ww:2*bw+ww,ii)=ones(bw,bw).*in(hh,ww,ii);
            out( bw+1:bw+hh,bw+1:bw+ww,ii )= in(:,:,ii);
            out(1:bw,bw+1:bw+ww,ii)=ones(bw,1)*in(1,:,ii);
            out(bw+hh+1:2*bw+hh,bw+1:bw+ww,ii)=ones(bw,1)*in(hh,:,ii);
            out(bw+1:bw+hh,1:bw,ii)=in(:,1,ii)*ones(1,bw);
            out(bw+1:bw+hh,bw+ww+1:2*bw+ww,ii)=in(:,ww,ii)*ones(1,bw);
        end
    end
end
