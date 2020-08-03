classdef mField < handle
    properties (Constant)
        mu0 = 4e-7*pi;
    end
    properties 
        % fields
        mx
        my
        mz
        Bx
        By
        Bz
        
        % geometry
        nDim
        cellSize %physical voxel dimensions (3x1 array)
        
        % physical
        Ms %saturation magnetization
        gamma
        Aex %exch const
        
        
        % topological metrics
        tHedgehog
        tHopfion
        tVortex
        tMach
        tMach2
        
        % boolean flags
        showPlot = 0;
        isNormalized = 0;
    end
    methods
        %long functions defined in separate .m files
        obj = hopfNumber(obj) 
        obj = hedgehogDensity(obj)
        obj = machNumber(obj);
        %short functions
        function obj = normalize(obj)
            mNorm = sqrt(obj.mx.^2 + obj.my.^2 + obj.mz.^2);
            mxb = zeros(size(obj.mx));
            myb = zeros(size(obj.my));
            mzb = zeros(size(obj.mz));
            inds = find(mNorm ~= 0);
            mxb(inds) = obj.mx(inds) ./ mNorm(inds);
            mxb(isnan(mxb)) = 0;
            myb(inds) = obj.my(inds) ./ mNorm(inds);
            myb(isnan(myb)) = 0;
            mzb(inds) = obj.mz(inds) ./ mNorm(inds);
            mzb(isnan(mzb)) = 0;
            obj.mx = mxb;
            obj.my = myb;
            obj.mz = mzb;
            obj.isNormalized = 1;
        end
        
        function obj = computeB(obj)
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

            [Dmx_Dx, Dmx_Dy, Dmx_Dz] = gradient(obj.mx);
            [Dmy_Dx, Dmy_Dy, Dmy_Dz] = gradient(obj.my);
            [Dmz_Dx, Dmz_Dy, Dmz_Dz] = gradient(obj.mz);

            Dm_Dx = cat(1,Dmx_Dx(:)',Dmy_Dx(:)', Dmz_Dx(:)');
            Dm_Dy = cat(1,Dmx_Dy(:)',Dmy_Dy(:)', Dmz_Dy(:)');
            Dm_Dz = cat(1,Dmx_Dz(:)',Dmy_Dz(:)', Dmz_Dz(:)');

            m = cat(1, obj.mx(:)', obj.my(:)', obj.mz(:)');
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
            obj.Bx = reshape(B(1,:),obj.nDim) ./ B_scale ;
            obj.By = reshape(B(2,:),obj.nDim) ./ B_scale ;
            obj.Bz = reshape(B(3,:),obj.nDim) ./ B_scale ;
            if nargout == 0
                clear obj
            end
        end
        
        function obj = interpField(obj,interpFac) 
            %interpolate to improve numerical accuracy in certain calculations (i.e. hopfNumber)
            if ~obj.isNormalized
                obj.normalize;
            end
            nn = obj.nDim(1);
            n_center = round((nn+1)/2);
            rvec1 = ((1:nn) - n_center) ./ (n_center-1);
            nn2 = round(interpFac*nn);
            n_center2 = round((nn2 + 1)/2);
            rvec2 = ((1:nn2) - n_center2) ./ (n_center2-1);
            
            [x1,y1,z1] = meshgrid(rvec1);
            [x2,y2,z2] = meshgrid(rvec2);
            
            mx_i = interp3(x1,y1,z1,obj.mx,x2,y2,z2,'linear',0);
            my_i = interp3(x1,y1,z1,obj.my,x2,y2,z2,'linear',0);
            mz_i = interp3(x1,y1,z1,obj.mz,x2,y2,z2,'linear',0);

            obj.mx = mx_i;
            obj.my = my_i;
            obj.mz = mz_i;
            obj.nDim = size(mx_i);
%             obj.isNormalized = 0;
            obj.tHedgehog = [];
            obj.tHopfion = [];
            obj.tVortex = [];
            if nargout == 0
                clear obj
            end
        end
        %constructor
        function obj = mField(varargin) % mField(mx,my,mz)
            if nargin > 1
                obj.mx = varargin{1};
                obj.my = varargin{2};
                obj.mz = varargin{3};
                obj.nDim = size(obj.mx);
            end
        end
        function obj = loadOVF(obj,fname,varargin) %load cubic geometry from ovf file
            fileID = fopen(fname, 'r');
            while 1
                tline = fgetl(fileID);
                if ~isempty(str2num(tline)), break, end
            end
            sizem = [3 Inf];
            formatSpec = '%f';
            m = fscanf(fileID,formatSpec, sizem);
            m = m';
            firstline = str2num(tline);
            m = cat(1,firstline,m);
            fclose(fileID);
            if nargin > 2
                ndim = varargin{1};
            else
                nn = nthroot(length(m),3);
                ndim = [nn nn nn];
            end
            obj.mx = reshape(m(:,1),ndim);
            obj.my = reshape(m(:,2),ndim);
            obj.mz = reshape(m(:,3),ndim);

            obj.nDim = ndim;
            if nargout == 0
                clear obj
            end
        end
        function obj = streamB(obj,varargin) %plot B field streamlines w/ randomly generated seed points
            n = obj.nDim;
            [sx,sy,sz] = meshgrid(1:n(1),1:n(2),1:n(3));
            
            if (isempty(obj.Bx) || isempty(obj.By) || isempty(obj.Bz))
                obj.computeB;
            end
            
            bx = obj.Bx;
            by = obj.By;
            bz = obj.Bz;
            
            if nargin > 1
                Nseed = varargin{1};
            else
                Nseed = round(0.005 * prod(n));
            end
            
            xstart = n(1)*rand(Nseed,1);
            ystart = n(2)*rand(Nseed,1);
            zstart = n(3)*rand(Nseed,1);
            
            figure;
            streamline(stream3(sx,sy,sz,bx,by,bz,xstart,ystart,zstart));
            title(sprintf('B field with %d seed points',Nseed));
            if nargout == 0
                clear obj
            end
        end       
    end
end
