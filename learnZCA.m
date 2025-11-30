function [P,C] = learnZCA(InImg, PatchSize) 

    ImgZ = length(InImg);

    %% get patches
    numper = ceil(200000/ImgZ);
    patches = zeros(numper*ImgZ,PatchSize^2);
    for i = 1:ImgZ  
        im = im2col_mean_removal(InImg{i},[PatchSize PatchSize]); 
        indm=randperm(size(im,2));
        thepatches=im(:,indm(1:numper));
        patches((i-1)*numper+1:i*numper,:)=thepatches';
        clear thepatches im;
    end
    
    %% get P and C
    %%%%%% normalize for contrast
    patches = bsxfun(@rdivide, bsxfun(@minus, patches, mean(patches,2)), sqrt(var(patches,[],2)+10));

    %%%%%% ZCA whitening (with low-pass)
    CovMatr = cov(patches);
    P = mean(patches);
    [V,D] = eig(CovMatr);
    C = V * diag(sqrt(1./(diag(D) + 0.1))) * V';