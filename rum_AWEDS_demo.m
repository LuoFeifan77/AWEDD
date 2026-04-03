clc; clf; clear;
% to get Aniostropic WEDN
addpath(genpath('utils/'));%add path ,utils box!
addpath(genpath('tools/'));
addpath(genpath('curvature/'));
addpath(genpath('data_mesh/')); 
addpath('gspbox/');
%%
collection_type = 2;  % 1; diffrence between there!
% -----------------editable parameters------------
% the number of eigenfunctions and scales for WEDS
n_angles=8; % default 8
num_eigs=300; % default 300
weds_scale=30; % weds_scale = 64;% difine 128?   [12,24,48,64,128] 
% alpha=4;%[2,4,8,10,20,30,50,100]; % [1,2,4,8,10,20,30,50,100];%default   alpha=4
alpha=4;
curv_smooth=10; % default 10 is ok for all shapes
Nscales=31;
k_sample=6; 
% ad=1:3;
%----------------------parameters---------------------------% 
param.curv_smooth=curv_smooth;
param.n_eigen=num_eigs;
param.alpha=alpha;
angles=linspace(0,180,n_angles+1);
param.angles=angles(1:end-1);
param.weds_scale=weds_scale;% i wanna laugh!
param.Nscales=31;
Ns=9;

% -------------------DATA PATH & SAVE PATH------------------------------%
% DATA PATH
SHAPE_DIR_path{1}           =fullfile('../data_mesh/');  %
DESC_DIR_path{1}           =fullfile('../AWEDD/');  %

for t=1:1
    %--------------改这个区�?---------------------%
    SHAPE_DIR=SHAPE_DIR_path{t};
    DESC_DIR=DESC_DIR_path{t};
%     param.alpha=alpha(t);
%     DESC_DIR=fullfile('G:\test_param\alpha\non_norm_part',['alpha_',num2str(alpha(t))]);
    %---------------改变下采样的情况--------------------%
    if weds_scale <= 96 % 96
        sample = k_sample + 2;
    else
        sample = ceil(weds_scale/32 ) + 2; % ceil means round up to an integer!  ceil(weds_scale/16) + 2;
    end
    if sample <= Ns   %16  防止跑出�?大维�?1024=32*32  Ndim=xdim*samples
        
        k_scale = floor(linspace((Ns+1) ,1,sample));% 原始采样情况
        k_scale = k_scale(2:end-1); % m  St_s   25    19    13     7
        
    else
        k_scale = floor(linspace((Ns+1) ,1,sample-2)); %24,16,8  %floor(linspace(16,1,sample-2));
    end
    
    param.k_scale=k_scale;%选择加权小波基的个数
    
    
    %  DESC_DIR=fullfile(DESC_DIR,['k_sample_',num2str(k_sample(t))]);
    
    if ~exist(DESC_DIR, 'dir')
        mkdir(DESC_DIR);
    end
    
    % dispay infos
    fprintf('[i] compute AWEDS_non_norm descriptors:\n');
    
    %-----------------处理不同类型的数据文�?-----------------------%
    datatype='off';
    switch datatype
        case 'mat'
            SHAPES = dir(fullfile(SHAPE_DIR, '*.mat'));  %数据 .mat
        case 'off'
            SHAPES = dir(fullfile(SHAPE_DIR, '*.off'));  %数据 .off
        otherwise
            error('Unknown design type');
    end
    
    SHAPES = natsort({SHAPES.name}');
    
    time_start = tic;
    
    for s = 1:length(SHAPES) %运行后面90-99个模�?
        
        shapename = SHAPES{s};
        fprintf(1, '  %-30s \t', shapename);
        % Load shape
        switch datatype
            case 'mat'
                load(fullfile(SHAPE_DIR, shapename),'shape');%读取 .mat文件
            case 'off'
                shape=read_shape(fullfile(SHAPE_DIR, shapename)); %读取 .off文件
            otherwise
                error('Unknown design type');
        end
        
        
        if exist(fullfile(DESC_DIR, shapename), 'file')
            fprintf(1, '\t\t file exist...\n ')
            continue;
        end
        
        time_start = tic;
        desc = compute_AWEDS(shape,param);
        elapsed_time = toc(time_start);
        
        fprintf('%3.2fs\n',elapsed_time);
        
        
        if shapename(end-2:end)=='off'
            shapename=[shapename(1:end-3),'mat'];
        end
        
        save(fullfile(DESC_DIR,shapename),'desc');%
        fprintf('这是�?%d个模�?,还剩�?%d个模型未完成\n',s,length(SHAPES)-s);
        
    end
end
