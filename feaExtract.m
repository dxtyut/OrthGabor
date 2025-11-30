function f = feaExtract(ImgCells,V,P,C,PatchSize,paraBlkSize)

    NumImg = length(ImgCells);
    ImgIdx = (1:NumImg)';
    NumFilters = 8;
    
    [ImgCells, ImgIdx] = PCANet_output(ImgCells, ImgIdx, PatchSize, NumFilters, V, P, C);  

    [f] = PCANet_HashingHist(paraBlkSize,ImgIdx,ImgCells);
end


function [OutImg,OutImgIdx] = PCANet_output(InImg, InImgIdx, PatchSize, NumFilters, V, P, C)
    ImgZ = length(InImg);
    mag = (PatchSize-1)/2;
    OutImg = cell(NumFilters*ImgZ,1); 
    
    cnt = 0;
    for i = 1:ImgZ
        [ImgX, ImgY, NumChls] = size(InImg{i});
        img = zeros(ImgX+PatchSize-1,ImgY+PatchSize-1, NumChls);
        img((mag+1):end-mag,(mag+1):end-mag,:) = InImg{i};     
        im = im2col_mean_removal(img,[PatchSize PatchSize]); 
        %%%%% collect all the patches of the ith image in a matrix, and perform patch mean removal

        %%%%%% normalize for contrast, and ZCA
        patches=im';
        patches = bsxfun(@rdivide, bsxfun(@minus, patches, mean(patches,2)), sqrt(var(patches,[],2)+10));    
        patches = bsxfun(@minus, patches, P) * C;
        im=patches';
        clear patches;

        for j = 1:NumFilters
            cnt = cnt + 1;
            OutImg{cnt} = reshape(V(:,j)'*im,ImgX,ImgY);  % convolution output
        end
        InImg{i} = [];
    end
    OutImgIdx = kron(InImgIdx,ones(NumFilters,1)); 
end


function [f] = PCANet_HashingHist(paraBlkSize,ImgIdx,OutImg)

    NumImg = max(ImgIdx);
    f = cell(NumImg,1);
    map_weights = 2.^((8-1):-1:0); 
    %%%% weights for binary to decimal conversion

    for Idx = 1:NumImg

        Idx_span = find(ImgIdx == Idx);

        T = 0;
        for j = 1:8
            T = T + map_weights(j)*Heaviside(OutImg{Idx_span(j)}); 
            %%% weighted combination; hashing codes to decimal number conversion

            OutImg{Idx_span(j)} = [];
        end

        BlkOverLapRatio = 0;
        stride = round((1-BlkOverLapRatio)*paraBlkSize); 
        blkwise_fea = sparse(histc(im2col_general(T,paraBlkSize,stride),(0:2^8-1)')); 
        %%%%% calculate histogram for each local block in "T"

        feas = bsxfun(@times, blkwise_fea,2^8./sum(blkwise_fea)); 
        f{Idx} = vec([feas]');
    end

    f = [f{:}];
end


function X = Heaviside(X) %%% binary quantization
    X = sign(X);
    X(X<=0) = 0;
end

function x = vec(X) %%% vectorization
    x = X(:);
end