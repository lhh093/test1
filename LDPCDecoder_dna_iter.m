function [ o1,o2,results ] = LDPCDecoder_dna_iter( H, y,LLR_y1, LLR_y2, iterMax,lambda,N,K)
% 方法1译码的matlab仿真
% 复现论文: [LDPC_Codes_for_Portable_DNA_Storage]
% 编码者: 刘海航 华南师范大学        
% 时间: [25.3.24]

%   H为校验矩阵，用以进行判决；LLR_y为接收到的信号初始置信度；iterMax为最大迭代次数；返回v为解码后的信息序列估计值

addpath('Dna')
a=0.1;
rscode=N*log2(N+1);
rsmescode=K*log2(N+1);
%注释
rsEncoder = comm.RSEncoder(N,K,'BitInput',true);
rsdecoder = comm.RSDecoder(N, K, ...
           'BitInput', true, ...
           'NumCorrectedErrorsOutputPort', true);




%初始化第一码字 u 和 v 矩阵
ch1=zeros(size(LLR_y1));
Uji1 = zeros(size(H));
Uji1_iter = zeros(size(H));
Vij1 = zeros(size(H'));
[VerificationNodes1, VariableNodes1] = size(H);
N1=VariableNodes1-VerificationNodes1;

%初始化第二码字 u 和 v 矩阵
ch2=zeros(size(LLR_y2));
Uji2 = zeros(size(H));
Uji2_iter = zeros(size(H));
Vij2 = zeros(size(H'));
[VerificationNodes2, VariableNodes2] = size(H);
results = zeros(1, 2); 

%初始化Uji
for i = 1:1:VariableNodes1
            idx1 = find(H(:, i) == 1);
            for k = 1:1:length(idx1)
                 Vij1(i, idx1(k)) =  LLR_y1(i);
                 Vij2(i, idx1(k)) =  LLR_y2(i);
            end   
end 
T=1;
while T <= 2
    
    [y1, y2] = dna_to_binary(y);
    [x1, x2] = dna_to_binary(y);
    for iter = 1:1:iterMax
    
        
        %停止迭代判决
        if all(mod(H*(x1'), 2) == 0) && all(mod(H*(x2'), 2) == 0)
            fprintf("\n 迭代提前结束 %d",iter);
            break;
        end
        
        % Message passing between two codes
        for i = 1:VariableNodes1
                idx1 = find(H(:, i) == 1);
                %if (x1(i) == 1 && x2(i) == 0) || (x1(i) == 0 && x2(i) == 1)
                if (y1(i) == 1 && y2(i) == 0) || (y1(i) == 0 && y2(i) == 1)
                    % 当组合为 (1, 0) 或 (0, 1) 时，更新交换信息ch
                        %tmp1=a*sum(Uji2(idx1, i));
                        ch1(i)=LLR_y1(i)-a*sum(Uji2(idx1, i));
                        %tmp2=a*sum(Uji1(idx1, i));
                        ch2(i)=LLR_y2(i)-a*sum(Uji1(idx1, i));
                else
                        ch1(i)=LLR_y1(i);
                        ch2(i)=LLR_y2(i);
                end
        end
        
    
        %Check node update
        for j = 1:1:VerificationNodes1
            idx1 = find(H(j, :) == 1);
            for k = 1:1:length(idx1)
                multipleVal1 = 1;
                multipleVal2 = 1;
                for l = 1:1:length(idx1)
                     if k == l
                          continue;
                     end
                     tmp=Vij1(idx1(l), j);
                     tmp = min(tmp, 20);
                     tmp = max(tmp, -20);
                     multipleVal1 = multipleVal1 * tanh(tmp/2);
                     tmp=Vij2(idx1(l), j);
                     tmp = min(tmp, 20);
                     tmp = max(tmp, -20);
                     multipleVal2 = multipleVal2 * tanh(tmp/2);
                end
                Uji1(j, idx1(k)) = 2 * atanh(multipleVal1);
                Uji2(j, idx1(k)) = 2 * atanh(multipleVal2);
                
    
            end
        end
        
        %Variable node update
        for i = 1:1:VariableNodes1
             idx1 = find(H(:, i) == 1);
             for k = 1:1:length(idx1)
                 Vij1(i, idx1(k)) =  ch1(i) + sum(Uji1(idx1, i)) - Uji1(idx1(k), i);
                 Vij2(i, idx1(k)) =  ch2(i) + sum(Uji2(idx1, i)) - Uji2(idx1(k), i);
             end
        end
        %判决
        for i = 1:1:VariableNodes1
            idx1 = find(H(:, i) == 1);
            addVal = sum(Uji1(idx1, i)) + ch1(i);
            if(addVal < 0)
                x1(i) = 1;
            else
                x1(i) = 0;
            end
        end
        %第二判决
        for i = 1:1:VariableNodes1
            idx1 = find(H(:, i) == 1);
            addVal = sum(Uji2(idx1, i)) + ch2(i);
            if(addVal < 0)
                x2(i) = 1;
            else
                x2(i) = 0;
            end
        end
    end

    results(T) =iter;
    %rs硬译码部分
    

    % 检查输入数组长度是否相同
    if length(x1) ~= length(x2)
        error('输入数组x1和x2长度必须相同');
    end
    
    
    v1 = x1(1:N1);
    v2 = x2(1:N1);
    tmp=N1-rscode;
    v1(end-tmp+1:end) = [];
    v2(end-tmp+1:end) = [];
    [Rs_Dec1, errnum1]=rsdecoder(v1');
    [Rs_Dec2, errnum2]=rsdecoder(v2');
       


   
  
    if errnum1 == -1 && errnum2 == -1
        disp('\n rs译码错误');
        break;  % 满足条件则跳出循环
    end
    
    T=T+1;
    %重新赋值LLR
    Rs_Enc1  = rsEncoder(Rs_Dec1);
    Rs_Enc2  = rsEncoder(Rs_Dec2);
    dna_seq=binary_to_acgt(Rs_Enc1,Rs_Enc2);
    y(1:length(dna_seq)) = dna_seq;
    
    [Rs_llr1, Rs_llr2] = calculate_dna_llr(dna_seq,lambda);
    % 加权
    weighted_Rs_llr1 = 11 * Rs_llr1;
    weighted_Rs_llr2 = 11 * Rs_llr2;
    % 赋值到 LLR_y1 和 LLR_y2
    LLR_y1(1:length(Rs_llr1)) = weighted_Rs_llr1;
    LLR_y2(1:length(Rs_llr2)) = weighted_Rs_llr2;
    
    

    Uji1 = zeros(size(H));
    Uji2 = zeros(size(H));
    Vij2 = zeros(size(H'));
    for i = 1:1:VariableNodes1
            idx1 = find(H(:, i) == 1);
            for k = 1:1:length(idx1)
                 Vij1(i, idx1(k)) =  LLR_y1(i);
                 Vij2(i, idx1(k)) =  LLR_y2(i);
            end 

    end

end
o1= Rs_Dec1;
o2= Rs_Dec2;