function [ v ] = LDPCDecoder_SP( H, LLR_y, iterMax )
%LDPCDecoder_SP LDPC使用和积算法解码
%   H为校验矩阵，用以进行判决；LLR_y为接收到的信号初始置信度；iterMax为最大迭代次数；返回v为解码后的信息序列估计值

%初始化 u 和 v 矩阵
U0i = LLR_y;
Uji = zeros(size(H));
Vij = zeros(size(H'));
[VerificationNodes, VariableNodes] = size(H);
x = zeros(size(LLR_y));

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
%     for j = 1:1:VerificationNodes
%         idx = find(H(j, :) == 1);
%         for k = 1:1:length(idx)
%             multipleVal = 2*atanh(prod(tanh(Vij(idx, j)/2))/tanh(Vij(idx(k), j)/2));
% 
%             if multipleVal == inf || multipleVal == -inf
%                 if multipleVal == inf
%                     Uji(j, idx(k)) = 10;
%                     disp(['>> Uji is inf when j = ' num2str(j) ', i = ' num2str(idx(k))]);
%                 else
%                     Uji(j, idx(k)) = -10;
%                     disp(['>> Uji is -inf when j = ' num2str(j) ', i = ' num2str(idx(k))]);
%                 end
% %                 prodOfSign = prod( sign(Vij(idx, j)) ) / sign(Vij(idx(k), j));
% %                 if k == 1
% %                     minOfVal = min(abs(Vij(idx(2:end), j)));
% %                 elseif k == length(idx)
% %                     minOfVal = min(abs(Vij(idx(1:k-1), j)));
% %                 else
% %                     minOfVal = min(min( abs(Vij(idx(1:k-1), j)) ), min( abs(Vij(idx(k+1:end), j)) ) );
% %                 end
% %                 Uji(j, idx(k)) = prodOfSign * minOfVal;
% %                 
%             else
%                 Uji(j, idx(k)) = multipleVal;
%             end
% 
%         end
%     end

     %new
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
                
                
                % %%%源程序
                %             multipleVal1 = 2*atanh(prod(tanh(Vij1(idx1, j)/2))/tanh(Vij1(idx1(k), j)/2));
                %             multipleVal2 = 2*atanh(prod(tanh(Vij2(idx1, j)/2))/tanh(Vij2(idx1(k), j)/2));
                % 
                % %第一码字
                % if multipleVal1 == inf || multipleVal1 == -inf
                %     if multipleVal1 == inf
                %         Uji1(j, idx1(k)) = 10;
                %         %disp(['>> Uji is inf when j = ' num2str(j) ', i = ' num2str(idx1(k))]);
                %     else
                %         Uji1(j, idx1(k)) = -10;
                %         %disp(['>> Uji is -inf when j = ' num2str(j) ', i = ' num2str(idx1(k))]);
                %     end              
                % else
                %     Uji1(j, idx1(k)) = multipleVal1;
                % end
                % 
                % %第二码字
                % if multipleVal2 == inf || multipleVal2 == -inf
                %     if multipleVal2 == inf
                %         Uji2(j, idx1(k)) = 10;
                %         %disp(['>> Uji is inf when j = ' num2str(j) ', i = ' num2str(idx1(k))]);
                %     else
                %         Uji2(j, idx1(k)) = -10;
                %         %disp(['>> Uji is -inf when j = ' num2str(j) ', i = ' num2str(idx1(k))]);
                %     end
                % else
                %     Uji2(j, idx1(k)) = multipleVal2;
                % end
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
        break;
    end
end

v = x(1:384);
%v = x;