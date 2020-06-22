function obj = hopfNumber(obj)
%% load properties from mField class
    if ~obj.isNormalized
        obj = normalize(obj);
    end
    mx = obj.mx;
    my = obj.my;
    mz = obj.mz;
    nDim = obj.nDim;
    showPlot = obj.showPlot;
    n_center = round((nDim(1)+1)/2);
%% compute B from m

    A =       [0 0 0;
               0 0 1;
               0 -1 0];

    B =        [0 0 -1;
                0 0 0;
                1 0 0];

    C =        [0 1 0;
                -1 0 0;
                 0 0 0];

    eps_3d = cat(3,A,B,C); 

    [Dmx_Dx, Dmx_Dy, Dmx_Dz] = gradient(mx);
    [Dmy_Dx, Dmy_Dy, Dmy_Dz] = gradient(my);
    [Dmz_Dx, Dmz_Dy, Dmz_Dz] = gradient(mz);

    Dm_Dx = cat(1,Dmx_Dx(:)',Dmy_Dx(:)', Dmz_Dx(:)');
    Dm_Dy = cat(1,Dmx_Dy(:)',Dmy_Dy(:)', Dmz_Dy(:)');
    Dm_Dz = cat(1,Dmx_Dz(:)',Dmy_Dz(:)', Dmz_Dz(:)');

    m = cat(1, mx(:)', my(:)', mz(:)');
    Dm = cat(3,Dm_Dx,Dm_Dy,Dm_Dz);
    B = zeros(size(m));

    for ii = 1:3
        B_buffer = zeros(1,length(m));
        for jj = 1:3
            for kk = 1:3
                B_buffer = B_buffer + eps_3d(ii,jj,kk) .* ...
                    dot(m,cross(Dm(:,:,jj),Dm(:,:,kk)));
            end
        end
        B(ii,:) = B_buffer;
    end

    B_scale = 8*pi;
    B_x = reshape(B(1,:),nDim) ./ B_scale ;
    B_y = reshape(B(2,:),nDim) ./ B_scale ;
    B_z = reshape(B(3,:),nDim) ./ B_scale ;
%% plot slices and toroid    
    if showPlot
        figure;
        subplot(1,3,1);
        quiver(squeeze(B_x(:,:,n_center)),squeeze(B_y(:,:,n_center)));axis image; title('B_{XY}');
        subplot(1,3,2);
        quiver(squeeze(B_y(:,n_center,:)),squeeze(B_z(:,n_center,:)));axis image; title('B_{YZ}');
        m_r = sqrt(mx.^2 + my.^2 + mz.^2);
        m_phi = atan2(my,mx);
        m_theta = acos(mz ./ m_r);
        % toroid
        subplot(1,3,3);
        hold on;
        fvc = isosurface(m_theta,pi/2,abs(m_phi));
%         fvc = isosurface(m_theta,pi/2,m_phi);
        p = patch(fvc,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.6);
        isonormals(m_theta,p);
        camlight
        lighting none;
        hold off;
        axis image;
        view(45,45);
    end
%%
    Bx_FFT = My_FFTN(B_x);
    Bx_FFT = reshape(Bx_FFT,1,numel(Bx_FFT));
    By_FFT = My_FFTN(B_y);
    By_FFT = reshape(By_FFT,1,numel(By_FFT));
    Bz_FFT = My_FFTN(B_z);
    Bz_FFT = reshape(Bz_FFT,1,numel(Bz_FFT));
    B_k_FFT2 = vertcat(Bx_FFT,By_FFT,Bz_FFT);
    nn = size(B_x,1);
    n_center = round((nn+1)/2);
    n_max = n_center-1;
    kvec = 2*pi.* ((1:nn) - n_center) ./ n_max;
    [kx,ky,kz] = meshgrid(kvec,kvec,kvec);
    k = cat(1,kx(:)',ky(:)',kz(:)');
    k2 = dot(k,k);
    lll = dot(B_k_FFT2,cross(k,B_k_FFT2,1),1) ./ k2 ;
    inds = find(sqrt(k2) <= (2*pi) & k2 ~= 0);
    obj.tHopfion = sum(lll(inds)) .* 1i ./ numel(inds);
    if nargout == 0
        clear obj
    end
end