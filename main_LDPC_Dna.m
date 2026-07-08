
% LDPC码在DNA纳米孔信道下的matlab仿真
% 论文: [LDPC_Codes_for_Portable_DNA_Storage]

clear all
close all
clc

%% 添加工作路径
addpath('Encoder')
addpath('Decoder')
addpath('Dna')
addpath('date')
addpath('diag')
addpath('Low_Density_Parity_Check_-LDPC-_Codes_-_MATLAB_Simulation-main')

% %% H矩阵生成

 A = IEEE80216e(768, '1/2');
 
 A =double(A);


[DiagA,GT,H,rankA,N0,K0,M0]=Diag(A) ;%时间长，提前准备GT， load
[N1,N2]=size(GT');%N1 信息位 N2 码长




lambda=0.002:0.005:0.02;
BER = zeros(1, length(lambda));

avgTime = zeros(1, length(lambda));

% 日志记录到 mydiary.txt
clock;


for lambda_i = 1:1:length(lambda)
    disp(['lambda=' num2str(lambda(lambda_i)) ' is simulating...']);
    % 设定停止条件
    if lambda(lambda_i) == 0.00015
        maxBlocks = 3000;
    else
        maxBlocks = 50000;
    end
    times = zeros(1, maxBlocks);
    % 设定译码算法最大迭代次数
    iterMax=10;
    ErrorBits_SP1 = 0; 
    ErrorBits_SP2 = 0;

    blocks_SP=0;

    % 编码 
    for i = 1:1:maxBlocks

        
        % 生成矩阵编码（s --> x）
        s1 = randi([0, 1], 1, N1);
        x1 = mod(s1 * GT', 2);

        s2 = randi([0, 1], 1, N1);
        x2 = mod(s2 * GT', 2);
 
        

        if sum(mod(H*(x1'), 2)) > 0
            %sprintf('the '+ num2str(i) + ' th encoding is not right');
            disp('encoding is not right');
            continue;
        end
        if sum(mod(H*(x2'), 2)) > 0
            %sprintf('the '+ num2str(i) + ' th encoding is not right');
            disp('encoding is not right');
            continue;
        end

        % ACGT调制
        dna_seq=binary_to_acgt(x1,x2);

        % 纳米孔信道

        y = dna_channel(dna_seq,lambda(lambda_i));

        % 译码端接收
        [llr_1, llr_2] = calculate_dna_llr_count(dna_seq,y,lambda(lambda_i));  %文献《An_Asymmetric-Error-Aware_LDPC_Decoding_Algorithm_for_DNA_Storage》
        %[llr_1, llr_2] = calculate_dna_llr(y,lambda(lambda_i));               %文献《Error_Rate-Based_Log-Likelihood_Ratio_Processing_for_Low-Density_Parity-Check_Codes_in_DNA_Storage》

        %译码
         tic;
        %[v_SP1,v_SP2] = LDPCDecoder_dna_fangfa1( H, y,llr_1, llr_2,iterMax );
        [v_SP1,v_SP2] = LDPCDecoder_dna( H, llr_1, llr_2, iterMax );
        times(i) = toc;
        result1 = v_SP1(1:N1);
        result2 = v_SP2(1:N1);
        errorbits_SP1 = sum(s1 ~= result1);
        errorbits_SP2 = sum(s2 ~= result2);
        ErrorBits_SP1 = ErrorBits_SP1 + errorbits_SP1;
        ErrorBits_SP2 = ErrorBits_SP2 + errorbits_SP2;

        blocks_SP = blocks_SP + 1;
       
        fprintf("\n LAMBDA-%d round-%d Errors - %d,%d ,  BER - %d",lambda(lambda_i) ,i,errorbits_SP1,errorbits_SP2, (ErrorBits_SP1+ErrorBits_SP2)/(N1 * blocks_SP * 2));
    end

    BER(1, lambda_i) = (ErrorBits_SP1+ErrorBits_SP2)/(N1 * blocks_SP * 2);
    avgTime(lambda_i) = mean(times);

end


% 绘制  BER


figure('numbertitle','off','name','BER of  Decode algorithms 1')
semilogy(lambda, BER(1, :), 'K-^', 'LineWidth', 1.0, 'MarkerSize', 6); hold on; % SP 三角marker 黑线


xlabel('lambda'); ylabel('BER');

legend('BER',  'BER - ldpc')
grid on;



diary off
