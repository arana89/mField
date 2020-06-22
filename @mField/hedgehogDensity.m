function obj = hedgehogDensity(obj)
%% load properties from mField class
    if ~obj.isNormalized
        obj = normalize(obj);
    end
    nDim = obj.nDim;
    mx = obj.mx;
    my = obj.my;
    mz = obj.mz;
    %%
    [duX, duY, duZ] = gradient(mx);
    [dvX, dvY, dvZ] = gradient(my);
    [dwX, dwY, dwZ] = gradient(mz);

    dMdX = cat(1,duX(:)',dvX(:)',dwX(:)');

    dMdY = cat(1,duY(:)',dvY(:)',dwY(:)');

    dMdZ = cat(1,duZ(:)',dvZ(:)',dwZ(:)');

    H_density = (3 / (4 * pi)) .* dot(dMdX, cross(dMdY,dMdZ));

    obj.tHedgehog = reshape(H_density,nDim);

    if nargout == 0
        clear obj
    end
end
