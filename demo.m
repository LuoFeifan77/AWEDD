clc; clf; clear;
% to get Aniostropic WEDN
addpath(genpath('utils/'));%add path ,utils box!
addpath(genpath('tools/'));
addpath(genpath('curvature/'));

%%
% 测试
mesh_dir =  './TOSCA_SET/DESCS/AWEDS';
% num_model =1; % two shapes!

if_output = true;
% output_dir = './data_wavelet/test_all';%change  path of save
output_dir = './data_wavelet/test_meyer_wavelet_desc_Ns_10';
collection_type = 2;  % 1; diffrence between there!
% -----------------editable parameters------------
% the number of eigenfunctions and scales for WEDS
n_angles=8;%[10,12,14,16,1]; % default %n_angles=[4,6,8];
num_eigs=300; % default 300
weds_scale=192./n_angles; % weds_scale = 64;% difine 128? %keep it dont change!  24
alpha=12;%default 12% alpha=10; % value in [5,20] is ok, default 10; 6，8，10，12，14，16，18
curv_smooth=10; % default 10 is ok for all shapes
%----------------------parameters---------------------------% 
param.curv_smooth=curv_smooth;
param.n_eigen=num_eigs;
param.alpha=alpha;
angles=linspace(0,180,n_angles+1);
param.angles=angles(1:end-1);
param.weds_scale=weds_scale;% i wanna laugh!

%pass ,There are partial matching,I can try later!
ref = false;  % generating desc for partial matching   % descriptors
if ref   
    ref_name = 'cuts_victoria_shape_0.off'; % for example
    ref_fullname = [mesh_dir, '/', ref_name];
    if ref_name(end-2:end) == 'off'
        shape_ref=read_shape(ref_fullname);
    elseif ref_name(end-2:end) == 'obj'
        shape_ref=read_shape_obj(ref_fullname);
    elseif ref_name(end-2:end) == 'ply'
        shape_ref=read_shape_ply(ref_fullname);
    end
else
    shape_ref = 0;
end

    
for k = 0:num_model
    % parfor k = 0:num_model
    s_name = ['tr_reg_', num2str(k, '%03d'), '.off'];%num2str(0, '%03d')=000 num2str(1, '%03d')=001
    WEDS = waveletEnergyDecompositionSignatureTest(mesh_dir, s_name, ref, shape_ref, if_output, output_dir, collection_type,param);
    fprintf('这是第%d个模型,还剩下%d个模型未完成\n',k+1,num_model-k);
end
elapsed_time = toc(time_start);
fprintf(1,'%3.2fs\n',elapsed_time); % 打印运行时间

