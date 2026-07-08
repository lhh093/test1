
function [someSNR, ave_BER50] = BER(m,n)
rows=m;
cols=n;
rate=(cols-rows)/cols;
tic
cycle=10;
amp=1;
H=genH(rows,cols);
someSNR=[0.9:0.3:3];
ave_BER50=zeros(1,length(someSNR));
count=1;

for S_num=1:length(someSNR) 
    total_num50=0;    
    SNR=someSNR(S_num);
    EbNo=10.^(SNR/10);
    sigma=1/sqrt(2*rate*EbNo);
    for i=1:cycle
        s=round(rand(1, cols-rows));      
        %use s and H to encode 
        [code,P,rearranged_cols]=ldpc_encode(s,H);    
        tx_code=bpsk(code,amp);                                %编码后的码字     
        re_waveform=tx_code+randn(1,cols)*sigma;               %加噪声后的码字      
        
        %LDPC Decoding 译码后输出的码字
        [vhat50]=bp_decoder50(re_waveform,H,rate,EbNo);        %采用BP译码
        de_mes50= extract_mesg(vhat50,rearranged_cols);
        errors50=find(s~=de_mes50);
        total_num50=total_num50+length(errors50);              %译码后错误总数        
       
    end
    count=count*1.5;
    ave_BER50(S_num)=total_num50/(count*cycle*(cols-rows));
end
toc