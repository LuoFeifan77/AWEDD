 function descriptor = waveletEnergyDecompositionSignature(shape_dir, shape_name, ref, shape_ref, if_output, output_dir, collection_type,param)

addpath(genpath('utils'));
addpath(genpath('tools'));
addpath(genpath('curvature/'));
%直接使用.mat文件
shape_fullname = [shape_dir, '/', shape_name];% shape_dir=./data_mesh/  shape_name=tr_reg_000.off
... shape_fullname=./data_mesh/tr_reg_000.off 
if shape_name(end-2:end) == 'off'
    shape=read_shape(shape_fullname);
elseif shape_name(end-2:end) == 'obj'
    shape=read_shape_obj(shape_fullname);
elseif shape_name(end-2:end) == 'ply'
    shape=read_shape_ply(shape_fullname);
end


VERT = shape.VERT; FACE = shape.TRIV;%shape is struct! why T?
% [m,n]=size(V)
% [~, shape.n] = size(V);  % the conut of column 
options.symmetrize = 1;
options.normalize = 0;
options.n_eigen = param.n_eigen;

% type = 'conformal';  %default 
K = param.n_eigen; %take first 300 eigfunctions  
weds_scale=param.weds_scale; %128
%----------------------options---------------%
options.curv_smooth = param.curv_smooth;
options.n_eigen = param.n_eigen;
options.alpha = param.alpha;
options.angle = 0;   %why angle =0!
% refidx=4168; %参考点的索引
load('./data_mesh/data_coe.mat');% 11,15,19,23,27,31,35

for i=3:size(data_coe,2)-4  %6*7  选取Nscales=19
    Nscales=data_coe(1,i);
    coe_z=data_coe(2:end,i);
    DESC=cell(1,size(param.angles,2));
    sample=floor(linspace((Nscales+1),1,10));
    sample=sample(2:end);
    %% replace it with Aniostropic
    for an=1:numel(param.angles)
        options.angle=pi*param.angles(an)/180;
        [V, D, L, A] = calc_ALB(VERT, FACE, options);%part 1    D is not diagonal！
        D=diag(D);%300*300
        
%         if weds_scale <= 3*(Nscales+1)  % 96
%             sample = 3 + 2;
%         else
%             sample = ceil(weds_scale/(Nscales+1) ) + 2; % ceil means round up to an integer!  ceil(weds_scale/16) + 2;
%         end
%         if sample <= (Nscales+1)   %16  防止跑出最大维度1024=32*32  Ndim=xdim*samples
%             k_scale = floor(linspace((Nscales+1) ,1,sample));% floor means round down to an integer! floor(linspace(16,1,sample))
%             
%             k_scale = k_scale(2:end-1); % m  St_s   25    19    13     7
%         else
%             k_scale = floor(linspace((Nscales+1) ,1,sample-2)); %24,16,8  %floor(linspace(16,1,sample-2));
%         end
%         %换掉采样的情况，使用非线性采样的方法
%         N=cos(pi*(0:(Nscales+1))/(2*(Nscales+1)));
        
        k_scale=sample(an);
        
        if ref
            %[L2,A2] = compute_mesh_laplacian_plusA_half(shape_ref.VERT',shape_ref.TRIV',type,options);
            [~, evals,~, ~] = calc_ALB(shape_ref.VERT,shape_ref.TRIV, options);%part 2
            %     opts=struct('tol',5e-3,'p',10,'disp',0);
            lmax=max( evals);
        else
            lmax=max(diag(D));
        end
        % Nscales=32;% K=31 in the paper if replace it with other number,A,B,C,D,E need to change!
        %%
        lmax=lmax*1.01;%  why?
        %         fprintf('Measuring largest eigenvalue, lmax = %g\n',lmax);
        %         fprintf('Designing transform in spectral domain\n');
        % designtype='mexican_hat';%default
        designtype='meyer';
        [g,~]=sgwt_filter_design(lmax,Nscales,'designtype',designtype);
        %         [g,~]=my_mexican_hat_design(lmax,Nscales,coe_z);
        
        % [g,~]=sgwt_filter_design(lmax,Nscales,'designtype',designtype);
        % g - wavelet and scaling function kernel function handles,and g is struct
        % t - values of scale parameters used for wavelet kernels  why didn't use it
        
        % initialize
        Wnk=zeros(size(L,1),3, Nscales+1);% 6890x3x32  rij?
        Tnk=zeros(size(L,1),size(L,1)); %6890*6890
%         Win=zeros(length(k_scale),size(L,1),size(L,1)); % wij??? 3*6890*6890
        Win=zeros(size(L,1),size(L,1));
        clk=zeros(K,Nscales+1);%K=300,the count of eigvalues!
        des=zeros(size(L,1),Nscales+1,3);%6890*32*3
        
        for k=1:numel(g)  %numel return the quality of the element of the g
            clk(:,k)=g{k}(diag(D));  % g{1} is scaling function kernel, g{2} ... g{Nscales+1} are wavelet kernels
%             clk(:,k)=normalize(clk(:,k),'L1',1);  %测试一下！2021年3月15
            Wnk(:,:,k)=V*diag(clk(:,k))*V'*A*shape.VERT;%shape.VERT means vertices(vertex)
            Tnk=V*diag(clk(:,k))*V'*A;% what is this ?
            index = find(k_scale==k,1);  %return indix!   %wede_scale=96,3,2,1;wede_scale=128,4,3,2,1
            if ~isempty(index)
                if collection_type == 1
                    WinT = abs(Tnk);
%                     Win(index,:,:) = ( (WinT-repmat(min(WinT),size(WinT,1),1))./repmat(max(WinT)-min(WinT),size(WinT,1),1) );% what is this?
                    Win = ( (WinT-repmat(min(WinT),size(WinT,1),1))./repmat(max(WinT)-min(WinT),size(WinT,1),1) );
                elseif collection_type == 2   %run this step!
                    WinT = Tnk.^2;
%                     Win(index,:,:) = (WinT./repmat(sqrt(sum(WinT.^2,1)),size(WinT,1),1));
                    Win = (WinT./repmat(sqrt(sum(WinT.^2,1)),size(WinT,1),1));
                end
            end
        end
        
        for i=1:3
            wl=V'*A*squeeze(Wnk(:,i,:)).*clk;  %Remove singleton dimensions
            wl=diag(sum(wl,2));
            des(:,:,i)=(clk'*wl*(D.^2)*V'*A)'.*squeeze(Wnk(:,i,:));% WENS的能量
        end
        
        if collection_type == 1
            descriptor=sum(des,3);
        elseif collection_type == 2    %run this step!
            descriptor=sqrt(abs(sum(des,3))); %6890*32
        end
%           不加权值        
%         desc = [];
% %         for i=1:length(k_scale)  %k_scale=24,16,8
% %             T = squeeze(Win(i,:,:)); %6890*6890  删除维度为1的矩阵，3*1*2  to 3*2
% %             desc = [desc, (T'*descriptor)];  %6890*32,广播？拼凑换成6890*96  ,6890*128
% %             %             desc =[desc, descriptor];
% %         end
% %          T = squeeze(Win(i,:,:)); %6890*6890  删除维度为1的矩阵，3*1*2  to 3*2
%          T=Win;
%          desc = [desc, (T'*descriptor)];  %6890*32,广播？拼凑换成6890*96  ,6890*128
%         
%         descriptor = desc;
%         Fn = size(descriptor,2);
%         %改变采样方法
%         if Fn>=weds_scale
%             descriptor = descriptor(:,floor(linspace(1,Fn,weds_scale)));% downsample  可以改进吗？这个采样
%         end
        DESC{an}=descriptor;
        %         DESC{an}=normalize(descriptor,'L2',2);
        % temp=DESC{an};
        % save([output_dir, '/', shape_name(1:end-3),num2str(an-1, '%01d'),'.', 'mat'],'temp','-mat');
        % DESC= cell2mat(DESC);
    end
    % mesh visualization
    % DESC= cell2mat(DESC); %  Convert a cell to an array
    % dist=pdist2(DESC,DESC(refidx,:)); %
    % mplot_mesh(shape.VERT,shape.TRIV,dist);
    % colorbar;
    
    if if_output  % if_output = true;
        %     Winn.V = V;
        %     Winn.clk = clk;
        %     Winn.A = diag(full(A));
        %     Winn.D = diag(D);
        %     save([output_dir, '/', shape_name(1:end-3), 'mat'],'-struct','Winn');
        %     for an=1:numel(param.angles)
        %     save([output_dir, '/', shape_name(1:end-6),num2str(an, '%03d'), 'mat'],'DESC{an}','-mat');
        %     end
        DESC= cell2mat(DESC);
        DESC=normalize(DESC,'L2',2);
         save([output_dir, '/', shape_name(1:end-4),'_','descriptor','.', 'mat'],'DESC','-mat');
        
        %     save([output_dir, '/', shape_name(1:end-4),'_',num2str(param.alpha),'_','descriptor','.', 'mat'],'DESC','-mat');%测试alpha
        %     save([output_dir, '/', shape_name(1:end-4),'_',num2str(param.alpha),'_',num2str(size(param.angles,2)),'_','descriptor','.', 'mat'],'DESC','-mat');%保持维度不变，测试angles
%         save([output_dir, '/', shape_name(1:end-4),'_','Ns','_',num2str(Nscales),'_','descriptor','.', 'mat'],'DESC','-mat');%保持维度不变，测试angles
        %     save([output_dir, '/', shape_name(1:end-6),num2str(an, '%03d'), 'mat'],'-struct','Winn');
        %     dlmwrite([output_dir, '/', shape_name(1:end-6),num2str(an, '%03d'), 'txt'], descriptor,'delimiter',' ', 'precision','%4.6e');
        %       dlmwrite([output_dir, '/', shape_name(1:end-6),num2str(0, '%03d'), 'txt'], DESC,'delimiter',' ', 'precision','%4.6e');
    end
end