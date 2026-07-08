function [ Rs_Dec ] = LDPCDecoder_SP_iter( H, LLR_y, iterMax ,N,K,sigma)
% 复现论文: 在论文[RS-LDPC 级联码的新型编译码技术]
% 编码者: 刘海航 华南师范大学        
% 时间: [25.6.15]

U0i = LLR_y;
Uji = zeros(size(H));
Vij = zeros(size(H'));
[VerificationNodes, VariableNodes] = size(H);
x = zeros(size(LLR_y));

rsencoder = comm.RSEncoder(N,K,'BitInput',true);
rsdecoder = comm.RSDecoder(N, K, ...
           'BitInput', true, ...
           'NumCorrectedErrorsOutputPort', true);
T=1;



while T <= 2
    U0i = LLR_y;
    Uji = zeros(size(H));
    Vij = zeros(size(H'));
    for iter = 1:1:iterMax
       %disp(['the ' num2str(iter) '-th iteration of SP'])
        % 求解Vij矩阵
        for i = 1:1:VariableNodes
            idx = find(H(:, i) == 1);
            for k = 1:1:length(idx)
                 Vij(i, idx(k)) =  U0i(i) + sum(Uji(idx, i)) - Uji(idx(k), i);
            end
        end
        
    %     % 求解Uji矩阵
        for j = 1:1:VerificationNodes
            idx1 = find(H(j, :) == 1);
            for k = 1:1:length(idx1)
                multipleVal = 1;
                
                for l = 1:1:length(idx1)
                     if k == l
                          continue;
                     end
                     tmp=Vij(idx1(l), j);
                     tmp = min(tmp, 20);
                     tmp = max(tmp, -20);
                     multipleVal = multipleVal * tanh(tmp/2);
                     
                end
                Uji(j, idx1(k)) = 2 * atanh(multipleVal);
            end
        end


        %判决
        for i = 1:1:length(x)
            idx = find(H(:, i) == 1);
            addVal = sum(Uji(idx, i)) + U0i(i);
            if(addVal < 0)
                x(i) = 1;
            else
                x(i) = 0;
            end
        end
        
        %如果校验关系满足 break;
        %否则继续迭代
        if mod(H*(x'), 2) == 0
            fprintf("\n 迭代提前结束 %d",iter);
            break;
        end
    end
    
    v = x(1:384);
    %删除补充比特
    v(end-6+1:end) = [];
    %rs译码
    [Rs_Dec,errnum1]=rsdecoder(v');
    if errnum1 == -1
        disp('rs译码错误');
        break;  % 满足条件则跳出循环
    end
    Rs_Eec=rsencoder(Rs_Dec);
    LLR = 10*2*(1-2*Rs_Eec)/(sigma^2);
    LLR_y(1:length(LLR)) = LLR;
    T=  T+1;

end

%v = x;