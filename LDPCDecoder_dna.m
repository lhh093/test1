function [ v1,v2 ] = LDPCDecoder_dna( H, LLR_y1, LLR_y2, iterMax )
%LDPCDecoder_SP LDPC使用和积算法解码
%   H为校验矩阵，用以进行判决；LLR_y为接收到的信号初始置信度；iterMax为最大迭代次数；返回v为解码后的信息序列估计值
a=0.1;

%初始化第一码字 u 和 v 矩阵
U0i1 = LLR_y1;
Uji1 = zeros(size(H));
Vij1 = zeros(size(H'));
[VerificationNodes1, VariableNodes1] = size(H);
x1 = zeros(size(LLR_y1));

%初始化第二码字 u 和 v 矩阵
U0i2 = LLR_y2;
Uji2 = zeros(size(H));
Vij2 = zeros(size(H'));
[VerificationNodes2, VariableNodes2] = size(H);
x2 = zeros(size(LLR_y2));


for iter = 1:1:iterMax

    % % 求解Vij矩阵，存储初始或上一次迭代矩阵llr信息
    % for i = 1:1:VariableNodes1
    %     idx1 = find(H(:, i) == 1);
    %     for k = 1:1:length(idx1)
    %          Vij1(i, idx1(k)) =  U0i1(i) + sum(Uji1(idx1, i)) - Uji1(idx1(k), i);
    %          Vij2(i, idx1(k)) =  U0i2(i) + sum(Uji2(idx1, i)) - Uji2(idx1(k), i);
    %     end
    % end


    % 求解Vij矩阵，dna修改
    if iter == 1
        for i = 1:1:VariableNodes1
            idx1 = find(H(:, i) == 1);
            for k = 1:1:length(idx1)
                 Vij1(i, idx1(k)) =  U0i1(i) ;
                 Vij2(i, idx1(k)) =  U0i2(i) ;
            end
        end
    elseif iter > 1 && iter <= iterMax
        for i = 1:length(x1)
            idx1 = find(H(:, i) == 1);
            if (x1(i) == 1 && x2(i) == 0) || (x1(i) == 0 && x2(i) == 1)
                % 当组合为 (1, 0) 或 (0, 1) 时，执行算法1
                for k = 1:1:length(idx1)
                    Vij1(i, idx1(k)) =  U0i1(i) + sum(Uji1(idx1, i)) - Uji1(idx1(k), i) - a*sum(Uji2(idx1, i));
                    Vij2(i, idx1(k)) =  U0i2(i) + sum(Uji2(idx1, i)) - Uji2(idx1(k), i) - a*sum(Uji1(idx1, i));
                end
            else
                % 否则，执行算法2
                for k = 1:1:length(idx1)
                    Vij1(i, idx1(k)) =  U0i1(i) + sum(Uji1(idx1, i)) - Uji1(idx1(k), i);
                    Vij2(i, idx1(k)) =  U0i2(i) + sum(Uji2(idx1, i)) - Uji2(idx1(k), i);
                end
            end
        end
    end

    
    % 求解Uji矩阵，行信息更新
    for j = 1:1:VerificationNodes1
        idx1 = find(H(j, :) == 1);
        for k = 1:1:length(idx1)
            multipleVal1 = 2*atanh(prod(tanh(Vij1(idx1, j)/2))/tanh(Vij1(idx1(k), j)/2));
            multipleVal2 = 2*atanh(prod(tanh(Vij2(idx1, j)/2))/tanh(Vij2(idx1(k), j)/2));

            %第一码字
            if multipleVal1 == inf || multipleVal1 == -inf
                if multipleVal1 == inf
                    Uji1(j, idx1(k)) = 10;
                    %disp(['>> Uji is inf when j = ' num2str(j) ', i = ' num2str(idx1(k))]);
                else
                    Uji1(j, idx1(k)) = -10;
                   % disp(['>> Uji is -inf when j = ' num2str(j) ', i = ' num2str(idx1(k))]);
                end
%                 prodOfSign = prod( sign(Vij(idx, j)) ) / sign(Vij(idx(k), j));
%                 if k == 1
%                     minOfVal = min(abs(Vij(idx(2:end), j)));
%                 elseif k == length(idx)
%                     minOfVal = min(abs(Vij(idx(1:k-1), j)));
%                 else
%                     minOfVal = min(min( abs(Vij(idx(1:k-1), j)) ), min( abs(Vij(idx(k+1:end), j)) ) );
%                 end
%                 Uji(j, idx(k)) = prodOfSign * minOfVal;
%                 
            else
                Uji1(j, idx1(k)) = multipleVal1;
            end
            
            %第二码字
            if multipleVal2 == inf || multipleVal2 == -inf
                if multipleVal2 == inf
                    Uji2(j, idx1(k)) = 10;
                   % disp(['>> Uji is inf when j = ' num2str(j) ', i = ' num2str(idx1(k))]);
                else
                    Uji2(j, idx1(k)) = -10;
                    %disp(['>> Uji is -inf when j = ' num2str(j) ', i = ' num2str(idx1(k))]);
                end
%                 prodOfSign = prod( sign(Vij(idx, j)) ) / sign(Vij(idx(k), j));
%                 if k == 1
%                     minOfVal = min(abs(Vij(idx(2:end), j)));
%                 elseif k == length(idx)
%                     minOfVal = min(abs(Vij(idx(1:k-1), j)));
%                 else
%                     minOfVal = min(min( abs(Vij(idx(1:k-1), j)) ), min( abs(Vij(idx(k+1:end), j)) ) );
%                 end
%                 Uji(j, idx(k)) = prodOfSign * minOfVal;
%                 
            else
                Uji2(j, idx1(k)) = multipleVal2;
            end

        end
    end
    
    %第一判决
    for i = 1:1:length(x1)
        idx1 = find(H(:, i) == 1);
        addVal = sum(Uji1(idx1, i)) + U0i1(i);
        if(addVal < 0)
            x1(i) = 1;
        else
            x1(i) = 0;
        end
    end

    %第二判决
    for i = 1:1:length(x2)
        idx1 = find(H(:, i) == 1);
        addVal = sum(Uji2(idx1, i)) + U0i2(i);
        if(addVal < 0)
            x2(i) = 1;
        else
            x2(i) = 0;
        end
    end
    
    %如果校验关系满足 break;
    %否则继续迭代
    if all(mod(H*(x1'), 2) == 0) && all(mod(H*(x2'), 2) == 0)
        break;
    end
end

 % v1 = x1(1009:end);
 % v2 = x2(1009:end);
 v1 = x1;
 v2 = x2;
 % v1 = x1(1:384);
 % v2 = x2(1:384);
 % v1 = x1(1:1008);
 % v2 = x2(1:1008);