% LDPC码在DNA纳米孔信道下的matlab仿真
% 实现功能: 在论文[LDPC_Codes_for_Portable_DNA_Storage]的基础上添加rs级联结构，并且在译码过程中添加反馈结构


clear all
close all
clc

%% 预定义变量
  k            = 51;  
 %k=44;
  n             = 63; 

%k            = 103;       
%n             = 127; 

M= n+1;
rscode=n*log2(n+1);
rsmescode=k*log2(n+1);


%% 添加工作路径
addpath('Encoder')
addpath('Decoder')
addpath('Dna')
addpath('date')
addpath('diag')
addpath('Low_Density_Parity_Check_-LDPC-_Codes_-_MATLAB_Simulation-main')

%% 定义rs编解码器
rsEncoder = comm.RSEncoder(n,k,'BitInput',true);
rsDecoder = comm.RSDecoder(n,k,'BitInput',true);



%% H矩阵生成
%load("LDPC(1000,500)_(6.3).mat") 
%load("HI_QC.mat") 

%A =double(HI_QC);

%A =double(H);
%A = IEEE80216e(2304, '1/2');
%A = IEEE80216e(1440, '2/3B');
A = IEEE80216e(768, '1/2');
A =double(A);

RR=rank(A);
[DiagA,GT,H,rankA,N0,K0,M0]=Diag(A);
[N1,N2]=size(GT');%N1 信息位 N2 码长

%% 

lambda=0.03:0.005:0.03;
BER = zeros(4, length(lambda));
avgTime = zeros(1, length(lambda));
avgIter = zeros(1, length(lambda));
% 日志记录到 mydiary.txt
diary 'mylog.txt'
clock;

%rs级联ldpc
for lambda_i = 1:1:length(lambda)
    disp(['lambda=' num2str(lambda(lambda_i)) ' is simulating...']);
    % 设定停止条件
    if lambda(lambda_i) <= 0.03
        maxBlocks = 1500;
    else
        maxBlocks = 150;
    end

    % 设定译码算法最大迭代次数
    iterMax = 30;
    ErrorBits_SP1 = 0; 
    ErrorBits_SP2 = 0;
    blocks_SP=0;
    times = zeros(1, maxBlocks);
    totalIterSum = 0;
    % 编码 - 纳米孔信道
    for i = 1:1:maxBlocks

        recordStr = [];

        % 产生信号
        s1 = randi([0 1],k*log2(n+1),1);
        s2 = randi([0 1],k*log2(n+1),1);
        %rs编码
        Rs_Enc1  = rsEncoder(s1);
        Rs_Enc2  = rsEncoder(s2);
        combinedRs = [Rs_Enc1', Rs_Enc2'];
        %补数据
        message1    = [Rs_Enc1',(rand(1,N1-rscode)>=0.5)];
        message2    = [Rs_Enc2',(rand(1,N1-rscode)>=0.5)];

        %ldpc编码
        x1 = mod(message1 * GT', 2);
        x2 = mod(message2 * GT', 2);

        if sum(mod(H*(x1'), 2)) > 0
            sprintf('the '+ num2str(i) + ' th encoding is not right');
            continue;
        end
        if sum(mod(H*(x2'), 2)) > 0
            sprintf('the '+ num2str(i) + ' th encoding is not right');
            continue;
        end

        % ACGT调制
        dna_seq=binary_to_acgt(x1,x2);

        % 纳米孔信道
        
        y = dna_channel(dna_seq,lambda(lambda_i));

        % 译码端接收
        [llr_1, llr_2] = calculate_dna_llr(y,lambda(lambda_i));

        % %译码

      

        %提出的回馈级联
        tic;
        [Rs_Dec1,Rs_Dec2,results] = LDPCDecoder_dna_iter( H, y,llr_1, llr_2, iterMax,lambda(lambda_i),n,k);
        times(i) = toc;
        itertotal = results(1) + results(2);
        totalIterSum = totalIterSum + itertotal;

        errorbits_SP1 = sum(s1 ~= Rs_Dec1);
        errorbits_SP2 = sum(s2 ~= Rs_Dec2);
        ErrorBits_SP1 = ErrorBits_SP1 + errorbits_SP1;
        ErrorBits_SP2 = ErrorBits_SP2 + errorbits_SP2;

        blocks_SP = blocks_SP + 1;
        fprintf("\n LAMBDA-%d round-%d Errors - %d,%d ,  BER - %d",lambda(lambda_i) ,i,errorbits_SP1,errorbits_SP2, (ErrorBits_SP1+ErrorBits_SP2)/(rsmescode * blocks_SP * 2));
    end
    avgIter(lambda_i) = totalIterSum / blocks_SP;
    avgTime(lambda_i) = mean(times);
    BER(1, lambda_i) = (ErrorBits_SP1+ErrorBits_SP2)/(rsmescode * blocks_SP * 2);

end





% 绘制  BER
xlswrite('./BERofFourAlgorithm.xlsx', BER);


figure('numbertitle','off','name','BER of  Decode algorithms 1')
semilogy(lambda, BER(1, :), 'K-^', 'LineWidth', 1.0, 'MarkerSize', 6); hold on; % SP 三角marker 黑线
xlabel('lambda'); ylabel('BER');
legend('BER - rs级联ldpc',  'BER - ldpc')
grid on;
diary off
