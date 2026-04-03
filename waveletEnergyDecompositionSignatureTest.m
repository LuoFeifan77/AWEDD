function descriptor = waveletEnergyDecompositionSignatureTest(shape_dir, shape_name, ref, shape_ref, if_output, output_dir, collection_type,param)

addpath(genpath('utils'));
addpath(genpath('tools'));
addpath(genpath('curvature/'));
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
% type = 'conformal';  %default 
K = param.n_eigen; %take first 300 eigfunctions  
weds_scale=param.weds_scale; %128
%----------------------options---------------%
options.symmetrize = 1;
options.normalize = 0;
options.n_eigen = param.n_eigen;
options.curv_smooth = param.curv_smooth;
options.n_eigen = param.n_eigen;
options.alpha = param.alpha;
options.angle = 0;   %why angle =0!
% refidx=4168; %参考点的索引
load('./data_mesh/data_coe.mat');% 11,15,19,23,27,31,35
%% aweds part
% for i=3:size(data_coe,2)-4  %6*7  选取Nscales=19
%     Nscales=data_coe(1,i);
%     coe_z=data_coe(2:end,i);
    DESC=cell(1,size(param.angles,2));
%     sample=floor(linspace((Nscales+1),1,10));
%     sample=sample(2:end);
    Nscales=9;
    %% replace it with Aniostropic
    for an=1:numel(param.angles)
        options.angle=pi*param.angles(an)/180;
        [V, D, L, A] = calc_ALB(VERT, FACE, options);%part 1    D is not diagonal！
        D=diag(D);%300*300
        %         k_scale=sample(an);
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
        % initialize
        Wnk=zeros(size(L,1),3, Nscales+1);% 6890x3x32  rij?
        %         Tnk=zeros(size(L,1),size(L,1)); %6890*6890
        %         % Win=zeros(length(k_scale),size(L,1),size(L,1)); % wij??? 3*6890*6890
        %         Win=zeros(size(L,1),size(L,1));
        clk=zeros(K,Nscales+1);%K=300,the count of eigvalues!
        des=zeros(size(L,1),Nscales+1,3);%6890*32*3
        
        for k=1:numel(g)  %numel return the quality of the element of the g
            clk(:,k)=g{k}(diag(D));  % g{1} is scaling function kernel, g{2} ... g{Nscales+1} are wavelet kernels
            %             clk(:,k)=normalize(clk(:,k),'L1',1);  %测试一下！2021年3月15
            Wnk(:,:,k)=V*diag(clk(:,k))*V'*A*shape.VERT;%shape.VERT means vertices(vertex)
%             Tnk=V*diag(clk(:,k))*V'*A;% what is this ?
            %             index = find(k_scale==k,1);  %return indix!   %wede_scale=96,3,2,1;wede_scale=128,4,3,2,1
            %             if ~isempty(index)
            %                 if collection_type == 1
            %                     WinT = abs(Tnk);
            %                     %                     Win(index,:,:) = ( (WinT-repmat(min(WinT),size(WinT,1),1))./repmat(max(WinT)-min(WinT),size(WinT,1),1) );% what is this?
            %                     Win = ( (WinT-repmat(min(WinT),size(WinT,1),1))./repmat(max(WinT)-min(WinT),size(WinT,1),1) );
            %                 elseif collection_type == 2   %run this step!
            %                     WinT = Tnk.^2;
            %                     %                     Win(index,:,:) = (WinT./repmat(sqrt(sum(WinT.^2,1)),size(WinT,1),1));
            %                     Win = (WinT./repmat(sqrt(sum(WinT.^2,1)),size(WinT,1),1));
            %                 end
            %             end
            %         end
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
            
            DESC{an}=descriptor;
            
        
    end
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
    end




