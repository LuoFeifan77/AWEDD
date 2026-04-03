function aweds = compute_AWEDS(shape, param)

addpath(genpath('utils'));
addpath(genpath('tools'));
addpath(genpath('curvature/'));
% refidx=4168; %参考点的索引
% load('./data_mesh/data_coe.mat');% 11,15,19,23,27,31,35
%-----------------------------options-----------------------%
% type = 'conformal';  %default
options.symmetrize = 1;
options.normalize = 0;
options.n_eigen = param.n_eigen;
options.curv_smooth = param.curv_smooth;
options.n_eigen = param.n_eigen;
options.alpha = param.alpha;
options.angle = 0;   %why angle =0!
%----------------------------type-----------------------------------%
type.wavelet='mexican_hat';  % mexican_hat
type.Nscales=param.Nscales;%31 default
type.collection=2;

DESC=cell(1,size(param.angles,2));

%% replace it with Aniostropic
for an=1:numel(param.angles)
    
    options.angle=pi*param.angles(an)/180;
    DESC{an}=compute_AWEDS_single(shape,param,options,type);
    
end

DESC= cell2mat(DESC);
aweds=DESC;

% aweds=normalize(DESC,'L2',2);

end




function descriptor=compute_AWEDS_single(shape,param,options,type)
if ~isfield(shape, 'VERT')
    shape.VERT = [shape.X,shape.Y,shape.Z];
end
% FACE = shape.TRIV; %直接使用.mat文件
Nscales=type.Nscales;
K = param.n_eigen; %take first 300 eigfunctions
weds_scale=param.weds_scale; %128
k_scale=param.k_scale;
[V, D, L, A] = calc_ALB(shape.VERT, shape.TRIV, options);%part 1    D is not diagonal！
lmax=D(end);
lmax=lmax*1.01;%  why?
fprintf('Measuring largest eigenvalue, lmax = %g\n',lmax);

D=diag(D);%300*300

% if weds_scale <= 96 % 96
%     sample = 3 + 2;
% else
%     sample = ceil(weds_scale/32 ) + 2; % ceil means round up to an integer!  ceil(weds_scale/16) + 2;
% end
% if sample <= Nscales   %16  防止跑出最大维度1024=32*32  Ndim=xdim*samples
%     k_scale = floor(linspace((Nscales+1) ,1,sample));% floor means round down to an integer! floor(linspace(16,1,sample))
%     
%     k_scale = k_scale(2:end-1); % m  St_s   25    19    13     7
% else
%     k_scale = floor(linspace((Nscales+1) ,1,sample-2)); %24,16,8  %floor(linspace(16,1,sample-2));
% end

%%

%fprintf('Designing transform in spectral domain\n');
designtype=type.wavelet;%default
collection_type =type.collection;%default
%designtype='meyer';
[g,~]=sgwt_filter_design(lmax,Nscales,'designtype',designtype);
%[g,~]=my_mexican_hat_design(lmax,Nscales,coe_z);
designtype='meyer';
Ns=9;
[f,~]=sgwt_filter_design(lmax,Ns,'designtype',designtype);
%----------------initialize------------------------%
Wnk=zeros(size(L,1),3, Nscales+1);% 6890x3x32  rij?
Tnk=zeros(size(L,1),size(L,1)); %6890*6890
Win=zeros(length(k_scale),size(L,1),size(L,1)); % wij??? 3*6890*6890
%Win=zeros(size(L,1),size(L,1));
clk=zeros(K,Nscales+1);%K=300,the count of eigvalues!
dlk=zeros(K,Ns+1);
des=zeros(size(L,1),Nscales+1,3);%6890*32*3


for k=1:numel(g)
    clk(:,k)=g{k}(diag(D));  %
    %clk(:,k)=normalize(clk(:,k),'L1',1);  %测试一下！2021年3月15
    Wnk(:,:,k)=V*diag(clk(:,k))*V'*A*shape.VERT;

end


for k=1:numel(f)
    index = find(k_scale==k,1);
    if ~isempty(index)
       dlk(:,k)=f{k}(diag(D));
%        dlk(:,k)=normalize(clk(:,k),'L1',1);
       Tnk=V*diag(dlk(:,k))*V'*A;
        if collection_type == 1
            WinT = abs(Tnk);
            Win(index,:,:) = ( (WinT-repmat(min(WinT),size(WinT,1),1))./repmat(max(WinT)-min(WinT),size(WinT,1),1) );% what is this?
            %                     Win = ( (WinT-repmat(min(WinT),size(WinT,1),1))./repmat(max(WinT)-min(WinT),size(WinT,1),1) );
        elseif collection_type == 2   %run this step!
            WinT = Tnk.^2;
            Win(index,:,:) = (WinT./repmat(sqrt(sum(WinT.^2,1)),size(WinT,1),1));
            %                     Win = (WinT./repmat(sqrt(sum(WinT.^2,1)),size(WinT,1),1));
        end
    end
end    
    
for j=1:3
    wl=V'*A*squeeze(Wnk(:,j,:)).*clk;  %Remove singleton dimensions
    wl=diag(sum(wl,2));
    des(:,:,j)=(clk'*wl*(D.^2)*V'*A)'.*squeeze(Wnk(:,j,:));% WENS的能量
end

if collection_type == 1
    descriptor=sum(des,3);
elseif collection_type == 2    %run this step!
    descriptor=sqrt(abs(sum(des,3))); %6890*32
end

desc = [];
for i=1:length(k_scale)  %k_scale=24,16,8
    T = squeeze(Win(i,:,:)); %6890*6890  删除维度为1的矩阵，3*1*2  to 3*2
    desc = [desc, (T'*descriptor)];  %6890*32,广播？拼凑换成6890*96  ,6890*128
    %desc =[desc, descriptor];
end

%         T = squeeze(Win(i,:,:)); %6890*6890  删除维度为1的矩阵，3*1*2  to 3*2
%         T=Win;
%         desc = [desc, (T'*descriptor)];  %6890*32,广播？拼凑换成6890*96  ,6890*128

descriptor = desc;
Fn = size(descriptor,2);
%改变采样方法
descriptor = descriptor(:,floor(linspace(1,Fn,weds_scale)));% downsample  可以改进吗？这个采样
end
% descriptor=normalize(desc,'L2',2);

