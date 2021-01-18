function obj = vortexDensity(obj)
% computes 3 components of vortex density
    nDim = obj.nDim;
    vD.x = zeros(nDim);
    vD.y = zeros(nDim);
    vD.z = zeros(nDim);
    
    if ~obj.isNormalized
        obj = normalize(obj);
    end
    
    mx = obj.mx;
    my = obj.my;
    mz = obj.mz;
    
    [duX, duY, duZ] = gradient(mx);
    [dvX, dvY, dvZ] = gradient(my);
    [dwX, dwY, dwZ] = gradient(mz);

    dMdX = cat(1,duX(:)',dvX(:)', dwX(:)');
    dMdY = cat(1,duY(:)',dvY(:)', dwY(:)');
    dMdZ = cat(1,duZ(:)',dvZ(:)', dwZ(:)');

    vdX = cross(dMdY,dMdZ);
    vD.x = reshape(vdX(1,:), nDim);

    vdY = cross(dMdZ,dMdX);
    vD.y = reshape(vdY(2,:), nDim);
    
    vdZ = cross(dMdX,dMdY);
    vD.z = reshape(vdZ(3,:), nDim);
    
    
%     for v_component = 1:3 % dimensions
%         switch v_component
%             case 1 %y
%                 m1 = 'mz';
%                 m2 = 'mx';
%                 d3 = nDim(1);
%                 v_expr = 'vD.y';
%                 slice_expr = '(ii,:,:)';
%             case 2 %x
% %                 m1 = 'my';
% %                 m2 = 'mz';
% %                 d3 = nDim(2);
% %                 v_expr = 'vD.x';
% %                 slice_expr = '(:,ii,:)';
%                  for ii = 1:nDim(2)
%                     [duX, duY] = gradient(squeeze(my(:,ii,:)));
%                     [dvX, dvY] = gradient(squeeze(mz(:,ii,:)));
%                     dZ = zeros(size(duX));
%                     dMdY = cat(1,duX(:)',dvX(:)', dZ(:)');
%                     dMdZ = cat(1,duY(:)',dvY(:)', dZ(:)');
%                     V_density = cross(dMdY,dMdZ);
%                     vD.x(:,ii,:) = reshape((1/pi).*V_density(3,:),size(my(:,ii,:)));
%                 end
%             case 3 %z
% %                 m1 = 'mx';
% %                 m2 = 'my';
% %                 d3 = nDim(3);
% %                 v_expr = 'vD.z';
% %                 slice_expr = '(:,:,ii)';
%                 for ii = 1:nDim(3)
%                     [duX, duY] = gradient(mx(:,:,ii));
%                     [dvX, dvY] = gradient(my(:,:,ii));
%                     dZ = zeros(size(duX));
%                     dMdX = cat(1,duX(:)',dvX(:)', dZ(:)');
%                     dMdY = cat(1,duY(:)',dvY(:)', dZ(:)');
%                     V_density = cross(dMdX,dMdY);
%                     vD.z(:,:,ii) = reshape((1/pi).*V_density(3,:),size(mx(:,:,ii)));
%                 end
%         end
% %               for ii = 1:d3
% %                 eval(sprintf('[duA, duB] = gradient(%s%s);',m1,slice_expr));
% %                 eval(sprintf('[dvA, dvB] = gradient(%s%s);',m2,slice_expr));
% %                 dC = zeros(size(duA));
% %                 dMdA = cat(1,duA(:)',dvA(:)', dC(:)');
% %                 dMdB = cat(1,duB(:)',dvB(:)', dC(:)');
% %                 V_slice = cross(dMdA,dMdB);
% % %                 vD{v_component}(:,:,ii) = reshape((1/pi).*V_slice(v_component,:),size(d1,d2));
% % %                 eval(sprintf('vD{v_component}%s = reshape((1/pi).*V_slice(v_component,:),size(dC))',slice_expr));
% %                 eval(sprintf('%s%s = reshape((1/pi).*V_slice(v_component,:),size(dC));',v_expr,slice_expr));
% %               end
%     end
        
    obj.tVortexDensity = vD;
    
    if nargout == 0
        clear obj
    end
end