function obj = machNumber(obj)
%contributed by Ezio Iacocca
%% load properties from mField class
    if ~obj.isNormalized
        obj = normalize(obj);
    end
    showPlot = obj.showPlot;
    mx = obj.mx;
    my = obj.my;
    mz = obj.mz;
    if isempty(obj.Aex) || isempty(obj.Ms)
        out = input('Please specify [Aex(J/m) Ms(A/m)]:\n');
        Aex = out(1);
        Ms = out(2);
        obj.Aex = Aex;
        obj.Ms = Ms;
    end
    Ms = obj.Ms;
    Aex = obj.Aex;
    mu0 = obj.mu0;
    if isempty(obj.cellSize)
        cS = input('Please specify cell size [X(m) Y(m) Z(m)]:\n');
        obj.cellSize = cS;
    end
    cellSize = obj.cellSize;
%%
    lex=sqrt(2*Aex/(mu0*Ms^2));
    Ux = cell(3,1);
    Uy = cell(3,1);
    Uz = cell(3,1);
    tMach = cell(3,1);
% Loop through dimensions
    for I = 1:3
        switch I
            case 1
                Re = my;
                Im = mz;
                N = mx;
            case 2
                Re = mz;
                Im = mx;
                N = my;
            case 3
                Re = mx;
                Im = my;
                N = mz;
        end
% Compute the fluid flow
        [Rex,Rey,Rez] = gradient(Re);
        [Imx,Imy,Imz] = gradient(Im);

        Ux{I} = (Im.*Rex - Re.*Imx)./ (1 - N.^2);
        Uy{I} = (Im.*Rey - Re.*Imy) ./ (1 - N.^2);
        Uz{I} = (Im.*Rez - Re.*Imz) ./ (1 - N.^2);

% Scale discretization
        Ux{I} = Ux{I} * lex / cellSize(1);
        Uy{I} = Uy{I} * lex / cellSize(2);
        Uz{I} = Uz{I} * lex / cellSize(3);

% Compute Mach numbers
        U = sqrt(Ux{I}.^2 + Uy{I}.^2 + Uz{I}.^2);
        tMach{I} = U.*sqrt((1+3*N.^2)./(1-N.^2));
% compute fluid velocity
        U2{I}(:,:,:) = sqrt(Ux{I}.^2 + Uy{I}.^2 + Uz{I}.^2); %per-spin basis
%         tMach2{I} = U2.*sqrt((1+3*N.^2)./(1-N.^2));
        tFluidVel{I} = U2{I};
    end
    obj.tMach = tMach;
    obj.tFluidVel = tFluidVel;
    if nargout == 0
        clear obj
    end
end