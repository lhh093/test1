rows=512;
cols=1024;
rate=(cols-rows)/cols;
H=genH(rows,cols);
s=round(rand(1, cols-rows));      
 %use s and H to encode 
[code,P,rearranged_cols]=ldpc_encode(s,H); 