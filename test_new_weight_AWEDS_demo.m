clc; clf; clear;
% to get Aniostropic WEDN
addpath(genpath('utils/'));%add path ,utils box!
addpath(genpath('tools/'));
addpath(genpath('curvature/'));
%%
collection_type = 2;  % 1; diffrence between there!
% -----------------editable parameters------------
%鏁版嵁娴嬭瘯 浣跨敤鏈?浼樺弬鏁帮紒
% the number of eigenfunctions and scales for WEDS
n_angles=8; % default 8
num_eigs=300; % default 300
weds_scale=30; % weds_scale = 64;% difine 128?   [12,24,48,64,128]  缁存寔鎻忚堪瀛愮殑缁村害涓?192缁?
% alpha=4;%[2,4,8,10,20,30,50,100]; % [1,2,4,8,10,20,30,50,100];%default   alpha=4
alpha=4;
curv_smooth=10; % default 10 is ok for all shapes
Nscales=31;
k_sample=6;  %鍏ㄩ噰鏍?
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

% --------------------娴嬭瘯TOSCA鏁版嵁闆?------------------------------%
% 鏂囦欢鐩綍
% DATA_DIR =  './TOSCA_SET/DESCS/AWEDS';
% DESC_DIR = fillfile('E:\MatlabWorkPlace\TwoType_WEDS-master\AWEDS-master\TOSCA_SET\DESCS\AWEDS'); %鍒涘缓鎻忚堪瀛愮洰褰?
% DATA_ROOT_DIR       =fullfile('E:\MatlabWorkPlace\TwoType_WEDS-master\AWEDS-master\TOSCA_SET');
% SHAPE_DIR           =fullfile(DATA_ROOT_DIR,'TOSCA_DATA','shapes');  %
%  娴嬭瘯remesh FAUST鏁版嵁闆?
SHAPE_DIR_path{1}           =fullfile('I:\dataset\Faust_remesh\Faust_5k_off');  %
DESC_DIR_path{1}           =fullfile('G:\Test_Data\test_remesh_FAUST\AWEDS_non_norm_desc968');  %

%  娴嬭瘯remesh SCAPE鏁版嵁闆?
SHAPE_DIR_path{2}           =fullfile('I:\dataset\Scape_remesh\SCAPE_offset');  %
DESC_DIR_path{2}           =fullfile('G:\Test_Data\test_remesh_SCAPE\AWEDS_deta_meyer_weight_desc');  %

%  娴嬭瘯remesh TOSCA鏁版嵁闆?
SHAPE_DIR_path{3}           =fullfile('I:\dataset\TOSCA_remesh\TOSCA_offset');  %
DESC_DIR_path{3}           =fullfile('G:\Test_Data\test_remesh_TOSCA\AWEDS_deta_meyer_weight_desc');  %
%--------------------娴嬭瘯鍙傛暟alpha--------------------------
% SHAPE_DIR=fullfile('E:\MatlabWorkPlace\TwoType_WEDS-master\DATA_AND_TEST\Data-Faust\offset'); 
% SHAPE_DIR_path{4}           =fullfile('G:\data_new\SHERC19\off_2');  %
% DESC_DIR_path{4}           =fullfile('G:\Test_Data\test_SHREC\AWEDS_new_meyer_weight_desc');  %

%--------------------娴嬭瘯鍙傛暟alpha--------------------------
SHAPE_DIR_path{4}           =fullfile('G:\remesh_data(faust,scape,tosca,shrec19)\TOSCA_nonIsometric\vtx_5k');  %
DESC_DIR_path{4}           =fullfile('G:\Test_Data\test_remesh_TOSCA_non_iso\AWEDS_deta_meyer_weight_desc');  %


for t=2:4
    %--------------鏀硅繖涓尯鍩?---------------------%
    SHAPE_DIR=SHAPE_DIR_path{t};
    DESC_DIR=DESC_DIR_path{t};
%     param.alpha=alpha(t);
%     DESC_DIR=fullfile('G:\test_param\alpha\non_norm_part',['alpha_',num2str(alpha(t))]);
    %---------------鏀瑰彉涓嬮噰鏍风殑鎯呭喌--------------------%
    if weds_scale <= 96 % 96
        sample = k_sample + 2;
    else
        sample = ceil(weds_scale/32 ) + 2; % ceil means round up to an integer!  ceil(weds_scale/16) + 2;
    end
    if sample <= Ns   %16  闃叉璺戝嚭鏈?澶х淮搴?1024=32*32  Ndim=xdim*samples
        
        k_scale = floor(linspace((Ns+1) ,1,sample));% 鍘熷閲囨牱鎯呭喌
        k_scale = k_scale(2:end-1); % m  St_s   25    19    13     7
        
    else
        k_scale = floor(linspace((Ns+1) ,1,sample-2)); %24,16,8  %floor(linspace(16,1,sample-2));
    end
    
    param.k_scale=k_scale;%閫夋嫨鍔犳潈灏忔尝鍩虹殑涓暟
    
    
    %  DESC_DIR=fullfile(DESC_DIR,['k_sample_',num2str(k_sample(t))]);
    
    if ~exist(DESC_DIR, 'dir')
        mkdir(DESC_DIR);
    end
    
    % dispay infos
    fprintf('[i] compute AWEDS_non_norm descriptors:\n');
    
    %-----------------澶勭悊涓嶅悓绫诲瀷鐨勬暟鎹枃浠?-----------------------%
    datatype='off';
    switch datatype
        case 'mat'
            SHAPES = dir(fullfile(SHAPE_DIR, '*.mat'));  %鏁版嵁 .mat
        case 'off'
            SHAPES = dir(fullfile(SHAPE_DIR, '*.off'));  %鏁版嵁 .off
        otherwise
            error('Unknown design type');
    end
    
    SHAPES = natsort({SHAPES.name}');
    
    time_start = tic;
    
    for s = 1:length(SHAPES) %杩愯鍚庨潰90-99涓ā鍨?
        
        shapename = SHAPES{s};
        fprintf(1, '  %-30s \t', shapename);
        % Load shape
        switch datatype
            case 'mat'
                load(fullfile(SHAPE_DIR, shapename),'shape');%璇诲彇 .mat鏂囦欢
            case 'off'
                shape=read_shape(fullfile(SHAPE_DIR, shapename)); %璇诲彇 .off鏂囦欢
            otherwise
                error('Unknown design type');
        end
        
        
        if exist(fullfile(DESC_DIR, shapename), 'file')
            fprintf(1, '\t\t file exist...\n ')
            continue;
        end
        
        time_start = tic;
        desc = compute_new_weight_AWEDS(shape,param);
        elapsed_time = toc(time_start);
        
        fprintf('%3.2fs\n',elapsed_time);
        
        
        if shapename(end-2:end)=='off'
            shapename=[shapename(1:end-3),'mat'];
        end
        
        save(fullfile(DESC_DIR,shapename),'desc');%瀛樺偍鏂囦欢
        fprintf('杩欐槸绗?%d涓ā鍨?,杩樺墿涓?%d涓ā鍨嬫湭瀹屾垚\n',s,length(SHAPES)-s);
        
    end
    %     fprintf('杩欐槸绗?%d涓枃浠?,杩樺墿涓?%d涓枃浠舵湭瀹屾垚\n',t,length(alpha)-t);
end