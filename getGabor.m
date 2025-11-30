
function GaborReal = getGabor(ke_h,ke_w)

    par.raT        =    0.9;            % Gabor kernel's energy preseving ratio
    par.Kmax       =    pi/2;           % Gabor kernel's para, default(pi/2)
    par.f          =    sqrt(2);        % Gabor kernel's para, default(sqrt(2))

    par.sigma      =    pi;             % Gabor kernel's para, default(pi or 1.5pi)
    
    par.ke_h=ke_h;
    par.ke_w=ke_w;

    [ GaborReal, GaborImg ]  =   MakeAllGaborKernal( par.ke_h, par.ke_w ,par.Kmax, par.f, par.sigma);

    [hh,ww,cc]=size(GaborReal);
    GaborReal = reshape(GaborReal,[hh*ww,cc])';
end







function [ GaborReal, GaborImg ] = MakeAllGaborKernal( GaborH, GaborW, Kmax, f, sigma)
    %%% GaborReal, [GaborH,GaborW,40] 40个Gabor模板实部
    %%% GaborImg,  [GaborH,GaborW,40] 40个Gabor模板虚部
    %%% 缺省输入前2个参数,后面参数 Kmax=2.5*pi/2, f=sqrt(2), sigma=1.5*pi;
    %%% GaborH, GaborW, Gabor模板大小
    %%% U,方向因子{0,1,2,3,4,5,6,7}
    %%% V,大小因子{0,1,2,3,4}
    %%% Kmax,f,sigma 其他参数

    %%
    GaborReal = zeros( GaborH, GaborW, 40 );
    GaborImg = zeros( GaborH, GaborW, 40 );

    for v = 0 : 4
        for u = 0 : 7
            [ GaborReal(:,:,v*8+u+1), GaborImg(:,:,v*8+u+1) ] = MakeGaborKernal( GaborH, GaborW, u, v, Kmax, f, sigma);
        end
    end
end





function [GaborReal, GaborImg] = MakeGaborKernal(GaborH, GaborW, U, V, Kmax,f,sigma)
    %%% function [GaborReal, GaborImg] = MakeGaborKernal[GaborH, GaborW, U, V]
    %%% 用以生成 Gabor 核
    %%% GaborReal: 核实部 GaborImg: 虚部
    %%% GaborH,GaborW: Gabor窗口 高宽.
    %%% U,V: 方向 大小
    %%%            ||Ku,v||^2
    %%% G(Z) = ---------------- exp(-||Ku,v||^2 * Z^2)/(2*sigma*sigma)(exp(i*Ku,v*Z)-exp(-sigma*sigma/2))
    %%%          sigma*sigma

    HarfH = fix(GaborH/2);
    HarfW = fix(GaborW/2);

    Qu = pi*U/8;

    sqsigma = sigma*sigma;
    %%%% Kv = 2.5*pi*(2^(-(V+2)/2));
    Kv = Kmax/(f^V);

    postmean = exp(-sqsigma/2);

    GaborReal = zeros(length(-HarfH : HarfH),length(-HarfW : HarfW));
    GaborImg = zeros(length(-HarfH : HarfH),length(-HarfW : HarfW));
    for j = -HarfH : HarfH
        for i =  -HarfW : HarfW

            tmp1 = exp(-(Kv*Kv*(j*j+i*i)/(2*sqsigma)));

            tmp2 = cos(Kv*cos(Qu)*i+Kv*sin(Qu)*j) - postmean;
            %%%%%%%tmp3 = sin(Kv*cos(Qu)*i+Kv*sin(Qu)*j) - exp(-sqsigma/2);
            tmp3 = sin(Kv*cos(Qu)*i+Kv*sin(Qu)*j);

            GaborReal(j+HarfH+1, i+HarfW+1) = Kv*Kv*tmp1*tmp2/sqsigma;
            GaborImg(j+HarfH+1, i+HarfW+1) = Kv*Kv*tmp1*tmp3/sqsigma;
        end
    end
end