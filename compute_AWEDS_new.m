function aweds = compute_AWEDD(shape, param)

addpath(genpath('utils'));
addpath(genpath('tools'));
addpath(genpath('curvature/'));
% refidx=4168; %鍙傝?冪偣鐨勭储寮?
% load('./data_mesh/data_coe.mat');% 11,15,19,23,27,31,35
%-----------------------------options-----------------------%
% type = 'conformal';  %default
options.symmetrize = 1;
options.normalize = 0;
options.n_eigen = param.n_eigen;
options.curv_smooth = param.curv_smooth;
options.n_eigen = param.n_eigen;
options.alpha = param.alpha;
options.angle = 0;   %init
%----------------------------type-----------------------------------%
type.wavelet='mexican_hat';  % mexican_hat
type.Nscales=param.Nscales;%31 default
type.collection=2;

DESC=cell(1,size(param.angles,2));

[descriptor,weight]=compute_AWEDS_single(shape,param,options,type);
DESC{1}=descriptor;
%% replace it with Aniostropic
for an=2:numel(param.angles)
%     param.k_scale=param.k_scale+an-1;
    options.angle=pi*param.angles(an)/180;
    %采用取各项同性的小波基进行加权
      descriptor=compute_AWEDS_single_2(shape,param,options,type,weight);  
    
    DESC{an}=descriptor;

end

DESC= cell2mat(DESC);
aweds=DESC;
% aweds=normalize(DESC,'L2',2);

end




function [descriptor,weight]=compute_AWEDD_single(shape,param,options,type)
if ~isfield(shape, 'VERT')
    shape.VERT = [shape.X,shape.Y,shape.Z];
end
% FACE = shape.TRIV;
Nscales=type.Nscales;
K = param.n_eigen; %take first 300 eigfunctions
weds_scale=param.weds_scale; %128
k_scale=param.k_scale;
[V, D, L, A] = calc_ALB(shape.VERT, shape.TRIV, options);
lmax=D(end);
lmax=lmax*1.01;
fprintf('Measuring largest eigenvalue, lmax = %g\n',lmax);
D=diag(D);%300*300
%%
%fprintf('Designing transform in spectral domain\n');
designtype=type.wavelet;%default
collection_type =type.collection;%default
%designtype='meyer';
[g,~]=sgwt_filter_design(lmax,Nscales,'designtype',designtype);
%[g,~]=my_mexican_hat_design(lmax,Nscales,coe_z);
%----------------initialize------------------------%
Wnk=zeros(size(L,1),3, Nscales+1);% 6890x3x32  rij?
% Tnk=zeros(size(L,1),size(L,1)); %6890*6890
Win=zeros(length(k_scale),size(L,1),size(L,1)); % wij??? 3*6890*6890
%Win=zeros(size(L,1),size(L,1));
clk=zeros(K,Nscales+1);%K=300,the count of eigvalues!
des=zeros(size(L,1),Nscales+1,3);%6890*32*3

for k=1:numel(g)
    clk(:,k)=g{k}(diag(D));  %
    clk(:,k)=normalize(clk(:,k),'L1',1);
    %clk(:,k)=normalize(clk(:,k),'L1',1);  %娴嬭瘯涓?涓嬶紒2021骞?3鏈?15
    Wnk(:,:,k)=V*diag(clk(:,k))*V'*A*shape.VERT;
    index = find(k_scale==k,1);
    if ~isempty(index)
        Tnk=V*diag(clk(:,k))*V'*A;
        if collection_type == 1
            WinT = abs(Tnk);
            Win(index,:,:) = ( (WinT-repmat(min(WinT),size(WinT,1),1))./repmat(max(WinT)-min(WinT),size(WinT,1),1) );
        elseif collection_type == 2
            WinT = Tnk.^2;
            Win(index,:,:) = (WinT./repmat(sqrt(sum(WinT.^2,1)),size(WinT,1),1));
        end
    end
end

weight.clk=clk;
weight.V=V;
weight.D=D;
weight.A=A;

for j=1:3
    wl=V'*A*squeeze(Wnk(:,j,:)).*clk;  %Remove singleton dimensions
    wl=diag(sum(wl,2));
    des(:,:,j)=(clk'*wl*(D.^2)*V'*A)'.*squeeze(Wnk(:,j,:));% 
end


if collection_type == 1
    descriptor=sum(des,3);
elseif collection_type == 2    %run this step!
    descriptor=sqrt(abs(sum(des,3))); %6890*32
end


desc = [];
for i=1:length(k_scale)  %k_scale=24,16,8
    T = squeeze(Win(i,:,:)); %6890*6890  
    desc = [desc, (T'*descriptor)];  %6890*32,
    %desc =[desc, descriptor];
end
%         T = squeeze(Win(i,:,:)); %6890*6890  鍒犻櫎缁村害涓?1鐨勭煩闃碉紝3*1*2  to 3*2
%         T=Win;
%         desc = [desc, (T'*descriptor)];  %6890*32,骞挎挱锛熸嫾鍑戞崲鎴?6890*96  ,6890*128
descriptor = desc;
Fn = size(descriptor,2);
%鏀瑰彉閲囨牱鏂规硶
descriptor = descriptor(:,floor(linspace(1,Fn,weds_scale)));% downsample  鍙互鏀硅繘鍚楋紵杩欎釜閲囨牱
end







